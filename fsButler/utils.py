import re
import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.analysis.utils as utils
import lsst.afw.display.ds9 as ds9

"""
Utility functions to process data elements delivered by fsButler
"""

#TODO: Make these configurable options

_fixedFields = ["id", "coord"]

_fixedPatterns = []

_suffixableFields = ["parent",
                     "deblend.nchild",
                     "classification.extendedness",
                     "flags.pixel.bad",
                     #"flux.kron.radius",
                     "flags.pixel.edge",
                     #"flags.pixel.interpolated.any",
                     #"flags.pixel.interpolated.center",
                     #"flags.pixel.saturated.any",
                     "flags.pixel.saturated.center",
                     "flux.psf",
                     "flux.psf.err",
                     "cmodel.flux",
                     "cmodel.flux.err"]

_suffixablePatterns = ["flux.zeromag*"]
                       #"flux.psf*",
                       #"cmodel*"]

_suffixRegex = re.compile(r'(_[grizy])$')

_zeroMagField = afwTable.Field["F"]("flux.zeromag",
                                    "The flux corresponding to zero magnitude")
_zeroMagErrField = afwTable.Field["F"]("flux.zeromag.err",
                                       "The flux error corresponding to zero magnitude")
_stellarField = afwTable.Field["Flag"]("stellar",
                                       "If true, the object is known to be a star if false it's known not to be a star.")

def _getFilterSuffix(filterSuffix):
    if filterSuffix == 'HSC-G':
        return '_g'
    elif filterSuffix == 'HSC-R':
        return '_r'
    elif filterSuffix == 'HSC-I':
        return '_i'
    elif filterSuffix == 'HSC-Z':
        return '_z'
    elif filterSuffix == 'HSC-Y':
        return '_y'
    else:
        return filterSuffix

def _suffixOrder(suffix):
    if suffix == '_g':
        return 1
    if suffix == '_r':
        return 2
    if suffix == '_i':
        return 3
    if suffix == '_z':
        return 4
    if suffix == '_y':
        return 5

def getCatSuffixes(cat):
    suffixes = []
    for schemaItem in cat.getSchema():
        fieldName = schemaItem.getField().getName()
        match = _suffixRegex.search(fieldName)
        if match:
            suffix = match.group(1)
            if suffix not in suffixes:
                suffixes.append(suffix) 
    suffixes.sort(key=_suffixOrder)
    return suffixes
    
def createSchemaMapper(cat, cat2=None, filterSuffix=None, withZeroMagFlux=False, withStellar=False):

    if cat2 and filterSuffix:
        raise ValueError("Can't use filterSuffix for two catalogs")

    suffixes = getCatSuffixes(cat)
    if len(suffixes) > 0 and filterSuffix:
        raise ValueError("Can't add a suffix to a catalog that already has one")

    schema = cat.getSchema()
    scm = afwTable.SchemaMapper(schema)

    # First fixed fields and patterns
    for f in _fixedFields: 
        scm.addMapping(schema.find(f).getKey())
    for p in _fixedPatterns:
        for f in schema.extract(p):
            scm.addMapping(schema.find(f).getKey())

    # Now suffixable fields and patterns
    if filterSuffix:
        suffix = _getFilterSuffix(filterSuffix)
    for f in _suffixableFields:
        if filterSuffix:
            field = schema.find(f).getField()
            newField = field.copyRenamed(f + suffix)
            scm.addMapping(schema.find(f).getKey(), newField)
        else:
            if len(suffixes) == 0:
                scm.addMapping(schema.find(f).getKey())
            else:
                for s in suffixes:
                    scm.addMapping(schema.find(f+s).getKey())
    for p in _suffixablePatterns:
        for f in schema.extract(p):
            if filterSuffix:
                field = schema.find(f).getField()
                newField = field.copyRenamed(f + suffix)
                scm.addMapping(schema.find(f).getKey(), newField)
            else:
                # The extract command gets all the suffixes for me
                scm.addMapping(schema.find(f).getKey())

    if cat2:
        suffixes2 = getCatSuffixes(cat2)
        schema2 = cat2.getSchema()
        for f in _suffixableFields:
            for s in suffixes2:
                field = schema2.find(f+s).getField()
                scm.addOutputField(field)
        for p in _suffixablePatterns:
            for f in schema2.extract(p):
                # The extract command gets the suffixes for me
                field = schema2.find(f).getField()
                scm.addOutputField(field)

    if withZeroMagFlux:
        if filterSuffix:
            scm.addOutputField(_zeroMagField.copyRenamed("flux.zeromag"+suffix))
            scm.addOutputField(_zeroMagErrField.copyRenamed("flux.zeromag.err"+suffix))
        else:
            if len(suffixes) == 0:
                scm.addOutputField(_zeroMagField)
            else:
                for s in suffixes:
                    scm.addOutputField(_zeroMagField.copyRenamed("flux.zeromag"+s))
                    scm.addOutputField(_zeroMagErrField.copyRenamed("flux.zeromag.err"+s))

    if withStellar:
        scm.addOutputField(_stellarField)

    return scm

def goodSources(cat):
    # Get the list of sources with bad flags
    bad = reduce(lambda x, y: np.logical_or(x, cat.get(y)),
                 ["flags.pixel.edge",
                  "flags.pixel.bad",
                  "flags.pixel.saturated.center"],
                  False)
    good = np.logical_not(bad)
    # Get rid of objects that have children, i.e. the deblender thinks it's a set of objects
    good = np.logical_and(good, cat.get("deblend.nchild") == 0)
    return good

def strictMatch(cat1, cat2, matchRadius=1*afwGeom.arcseconds):
    """
    Match two catalogs using a one to one relation where each match is the closest
    object
    """
    
    #import ipdb; ipdb.set_trace()

    matched = afwTable.matchRaDec(cat1, cat2, matchRadius, True)

    bestMatches = {}
    for m1, m2, d in matched:
        id = m2.getId()
        if id not in bestMatches:
            bestMatches[id] = (m1, m2, d)
        else:
            if d < bestMatches[id][1]:
                bestMatches[id] = (m1, m2, d)

    scm = createSchemaMapper(cat1, cat2)
    schema = scm.getOutputSchema()
    cat = afwTable.SimpleCatalog(schema)
    cat2Fields = []
    schema2 = cat2.getSchema()
    suffixes = getCatSuffixes(cat2)
    for suffix in suffixes:
        cat2Fields.extend(schema2.extract("*" + suffix).keys())
    for id in bestMatches:
        m1, m2, d = bestMatches[id]
        record = cat.addNew()
        record.assign(m1, scm)
        for f in cat2Fields:
            record.set(f, m2.get(f))
    return cat

def buildXY(hscCat, sgTable, matchRadius=1*afwGeom.arcseconds):

    print "Matching with HST catalog"
    matchedSG = afwTable.matchRaDec(sgTable, hscCat, matchRadius, False)
    print "Found {0} matches with HST objects".format(len(matchedSG))
    
    # Build truth table
    stellar = {}
    for m1, m2, d in matchedSG:
        id = m2.getId()
        isStar = (m1.get("mu.class") == 2)
        if id not in stellar:
            stellar[id] = [isStar, d]
        else:
            if d < stellar[id][1]:
                stellar[id] = [isStar, d] # Only keep closest for now
    print "Of which I picked {0}".format(len(stellar)) 

    scm = createSchemaMapper(hscCat, withStellar=True)
    schema = scm.getOutputSchema()
    cat = afwTable.SourceCatalog(schema)

    for m in hscCat:
        id = m.getId()
        if id in stellar:
            record = cat.addNew()
            record.assign(m, scm)
            record.set('stellar', stellar[id][0])

    return cat

def displayObject(objId, fsButler, prefix=''):

    dataType = prefix + '_calexp_sub'
    info = utils.makeMapperInfo(fsButler.butler)
    dataId = info.splitCoaddId(objId)
    dataId.pop('objId')
    src = fsButler.butler.get('deepCoadd_src', **dataId)
    src = src[objId == src.get("id")][0]
    bbox = src.getFootprint().getBBox()
    im = fsButler.butler.get(dataType, bbox=bbox, imageOrigin="PARENT", **dataId)
    ds9.mtv(im, frame=0)
