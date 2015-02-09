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
                     "flags.pixel.interpolated.any",
                     "flags.pixel.interpolated.center",
                     "flags.pixel.saturated.any",
                     "flags.pixel.saturated.center"]

_suffixablePatterns = ["flux.zeromag*",
                       "flux.psf*",
                       "cmodel*",
                       "centroid*",
                       "seeing*",
                       "exptime*"]

_suffixRegex = re.compile(r'(_[grizy])$')

_zeroMagField = afwTable.Field["F"]("flux.zeromag",
                                    "The flux corresponding to zero magnitude.")
_zeroMagErrField = afwTable.Field["F"]("flux.zeromag.err",
                                       "The flux error corresponding to zero magnitude.")
_stellarField = afwTable.Field["Flag"]("stellar",
                                       "If true, the object is known to be a star if false it's known not to be a star.")
_seeingField = afwTable.Field["F"]("seeing",
                                    "The PSF FWHM.")
_exptimeField = afwTable.Field["F"]("exptime",
                                    "Exposure time.")

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
    
def createSchemaMapper(cat, cat2=None, filterSuffix=None, withZeroMagFlux=False,
                       withStellar=False, withSeeing=False, withExptime=False):

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
                scm.addOutputField(_zeroMagErrField)
            else:
                for s in suffixes:
                    scm.addOutputField(_zeroMagField.copyRenamed("flux.zeromag"+s))
                    scm.addOutputField(_zeroMagErrField.copyRenamed("flux.zeromag.err"+s))

    if withSeeing:
        if filterSuffix:
            scm.addOutputField(_seeingField.copyRenamed("seeing"+suffix))
        else:
            if len(suffixes) == 0:
                scm.addOutputField(_seeingField)
            else:
                for s in suffixes:
                    scm.addOutputField(_seeingField.copyRenamed("seeing"+s))

    if withExptime:
        if filterSuffix:
            scm.addOutputField(_exptimeField.copyRenamed("exptime"+suffix))
        else:
            if len(suffixes) == 0:
                scm.addOutputField(_exptimeField)
            else:
                for s in suffixes:
                    scm.addOutputField(_exptimeField.copyRenamed("exptime"+s))

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
            if d < bestMatches[id][2]:
                bestMatches[id] = (m1, m2, d)

    scm = createSchemaMapper(cat1, cat2)
    schema = scm.getOutputSchema()
    cat = afwTable.SimpleCatalog(schema)
    cat.reserve(len(bestMatches))
    cat2Fields = []; cat2Keys = []; catKeys = []
    schema2 = cat2.getSchema()
    suffixes = getCatSuffixes(cat2)
    for suffix in suffixes:
        cat2Fields.extend(schema2.extract("*" + suffix).keys())
    for f in cat2Fields:
        cat2Keys.append(schema2.find(f).key)
        catKeys.append(schema.find(f).key)
    for id in bestMatches:
        m1, m2, d = bestMatches[id]
        record = cat.addNew()
        record.assign(m1, scm)
        for i in range(len(cat2Keys)):
            record.set(catKeys[i], m2.get(cat2Keys[i]))
    return cat

def buildXY(hscCat, sgTable, matchRadius=1*afwGeom.arcseconds):

    print "Matching with HST catalog"
    matchedSG = afwTable.matchRaDec(sgTable, hscCat, matchRadius, False)
    print "Found {0} matches with HST objects".format(len(matchedSG))
    
    # Build truth table
    stellar = {}
    classKey = sgTable.getSchema().find('mu.class').key
    for m1, m2, d in matchedSG:
        id = m2.getId()
        isStar = (m1.get(classKey) == 2)
        if id not in stellar:
            stellar[id] = [isStar, d, m2]
        else:
            if d < stellar[id][1]:
                stellar[id] = [isStar, d, m2] # Only keep closest for now
    print "Of which I picked {0}".format(len(stellar)) 

    scm = createSchemaMapper(hscCat, withStellar=True)
    schema = scm.getOutputSchema()
    cat = afwTable.SourceCatalog(schema)
    cat.reserve(len(stellar))
    stellarKey = schema.find('stellar').key

    for id in stellar:
        isStar, d, m2 = stellar[id]
        record = cat.addNew()
        record.assign(m2, scm)
        record.set(stellarKey, isStar)

    return cat

def displayObject(objId, fsButler, prefix='', frame=None):
    #TODO: Enable single exposure objects
    #dataType = prefix + '_calexp_sub'
    dataType = prefix + '_sub'
    info = utils.makeMapperInfo(fsButler.butler)
    dataId = info.splitCoaddId(objId)
    dataId.pop('objId')
    src = fsButler.butler.get('deepCoadd_src', **dataId)
    src = src[objId == src.get("id")][0]
    bbox = src.getFootprint().getBBox()
    im = fsButler.butler.get(dataType, bbox=bbox, imageOrigin="PARENT", **dataId)
    ds9.mtv(im, frame=frame)
    return im

def showCoaddInputs(objId, fsButler, coaddType="deepCoadd"):
    """Show the inputs for the specified object Id, optionally at the specified position
    @param fsButler    Butler to provide inputs
    @param coaddType Type of coadd to examine
    """
    info = utils.makeMapperInfo(fsButler.butler)
    dataId = info.splitCoaddId(objId)
    dataId.pop('objId')
    src = fsButler.butler.get('deepCoadd_src', **dataId)
    src = src[objId == src.get("id")][0]
    pos = src.getCentroid()
    coadd = fsButler.butler.get(coaddType, **dataId)
    visitInputs = coadd.getInfo().getCoaddInputs().visits
    ccdInputs = coadd.getInfo().getCoaddInputs().ccds
    posSky = coadd.getWcs().pixelToSky(pos)

    psf = coadd.getPsf()
    sigmaCoadd = psf.computeShape(pos).getDeterminantRadius()

    print "%6s %3s %7s %5s %5s" % ("visit", "ccd", "exptime", "FWHM", "weight")

    totalExpTime = 0.0
    expTimeVisits = set()

    for i in range(len(ccdInputs)):
        input = ccdInputs[i]
        ccd = input.get("ccd")
        v = input.get("visit")
        bbox = input.getBBox()
        # It's quicker to not read all the pixels, so just read 1
        calexp = fsButler.butler.get("calexp_sub", bbox=afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(1, 1)),
                            visit=int(v), ccd=ccd)
        calib = calexp.getCalib()
        psf = calexp.getPsf()
        pos = calexp.getWcs().skyToPixel(posSky)
        sigma = psf.computeShape(pos).getDeterminantRadius()
        exptime = calib.getExptime()
        weight = "%5.2f" % (input.get("weight"))
        if v not in expTimeVisits:
            totalExpTime += exptime
            expTimeVisits.add(v)
        print  "%6s %3s %7.0f %5.2f %5s" % (v, ccd, exptime, sigma, weight)
    print "Total Exposure time {0}".format(totalExpTime)
    print "Coadd FWHM {0}".format(sigmaCoadd)
