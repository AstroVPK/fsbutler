import re
import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

"""
Utility functions to process data elements delivered by fsButler
"""

_fixedFields = ["id", "coord"]

_fixedPatterns = []

_suffixableFields = ["parent",
                     "deblend.nchild",
                     "classification.extendedness",
                     "flags.pixel.bad",
                     "flux.kron.radius",
                     "flags.pixel.edge",
                     "flags.pixel.interpolated.any",
                     "flags.pixel.interpolated.center",
                     "flags.pixel.saturated.any",
                     "flags.pixel.saturated.center"]

_suffixablePatterns = ["flux.psf*", "cmodel*"]


_suffixRegex = re.compile(r'(_[grizy])$')

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

def getCatSuffixes(cat):
    suffixes = []
    for schemaItem in cat.getSchema():
        fieldName = schemaItem.getField().getName()
        match = _suffixRegex.search(fieldName)
        if match:
            suffix = match.group(1)
            if suffix not in suffixes:
                suffixes.append(suffix) 
    return suffixes
    
def createSchemaMapper(cat, cat2=None, filterSuffix=None, withZeroMagFlux=False):

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
        scm.addOutputField(afwTable.Field["F"]("flux.zeromag",
                                               "The flux corresponding to zero magnitude"))

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
    cat = afwTable.SourceCatalog(schema)
    cat2Fields = []
    schema2 = cat2.getSchema()
    suffixes = getCatSuffixes(cat2)
    for suffix in suffixes:
        cat2Fields.extend(schema2.extract("*" + suffix + "$").keys())
    for id in bestMatches:
        m1, m2, d = bestMatches[id]
        record = cat.addNew()
        record.assign(m1, scm)
        for f in cat2Fields:
            record.set(f, m2.get(f))

    return cat

def buildXY(hscCat, sgTable):
    #print "Loading {0}...".format(hscSources)
    #cat = afwTable.SourceCatalog.readFits(hscSources)
    #print "Loading {0}...".format(selectSG)
    #sgTable = afwTable.SimpleCatalog.readFits(selectSG)
    # Alexie's catalog is in degrees and LSST catalogs are in radians
    #sgTable["coord.ra"][:]  = np.radians(sgTable["coord.ra"])
    #sgTable["coord.dec"][:] = np.radians(sgTable["coord.dec"])
    
    #print "Matching catalogs..."
    matchRadius = 1*afwGeom.arcseconds
    matchedSG = afwTable.matchRaDec(sgTable, hscCat, matchRadius, False)
    
    #print "Building truth table..."
    stellar = {}
    for m1, m2, d in matchedSG:
        id = m2.getId()
        isStar = (m1.get("mu.class") == 2)
        if id not in stellar:
            stellar[id] = [isStar, d]
        else:
            if d < stellar[id][1]:
                stellar[id] = [isStar, d] # Only keep closest for now
    
    #print "Building arrays to train classifier..."
    X = []; Y = []
    for m in hscCat:
        id = m.getId()
        if id in stellar:
            ra = m.get("coord.ra")
            dec = m.get("coord.dec")
            cmodelFlux = m.get("cmodel.flux")
            psfFlux = m.get("flux.psf")
            fluxMag0 = m.get("flux.zeromag")
            if cmodelFlux <= 0.0 or psfFlux <= 0.0 or fluxMag0 <= 0.0:
                continue
            cmodelMag = -2.5*np.log10(cmodelFlux/fluxMag0)
            psfMag = -2.5*np.log10(psfFlux/fluxMag0)
            extendedness = -2.5*np.log10(psfFlux/cmodelFlux)
            # Keep only the objects that have no nans
            if not np.isnan(cmodelMag) and not np.isnan(extendedness):
                X.append([id, ra, dec, psfMag, cmodelMag, extendedness])
                Y.append(stellar[id][0])
    X = np.array(X)
    Y = np.array(Y, dtype=int)
    return X, Y
