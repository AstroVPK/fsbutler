import re
import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.coord as afwCoord
import lsst.analysis.utils as utils
import lsst.afw.display.ds9 as ds9

from lsst.pex.exceptions import LsstCppException
from lsst.meas.extensions.multiShapelet import FitPsfModel

"""
Utility functions to process data elements delivered by fsButler
"""

#TODO: Make these configurable options

_fixedFields = ["id", "coord"]

_fixedPatterns = []

_suffixableFields = ["parent",
                     "classification.extendedness",
                     "flags.pixel.bad",
                     #"flux.kron.radius",
                     "flags.pixel.edge",
                     "flags.pixel.interpolated.any",
                     "flags.pixel.interpolated.center",
                     "flags.pixel.saturated.any",
                     "flags.pixel.saturated.center"]

_suffixablePatterns = ["flux.zeromag*",
                       "multishapelet.psf*",
                       "shape*",
                       "flux.psf*",
                       "flux.kron*",
                       "cmodel*",
                       "centroid*",
                       "seeing*",
                       "exptime*",
                       "multId*",
                       "dGauss.radInner*",
                       "dGauss.radOuter*",
                       "dGauss.ampRat*",
                       "dGauss.qInner*",
                       "dGauss.qOuter*",
                       "dGauss.thetaInner*",
                       "dGauss.thetaOuter*"]

_suffixRegex = re.compile(r'(_[grizy])$')
_bandRegex = re.compile(r'(\.[grizy])$')

_zeroMagField = afwTable.Field["F"]("flux.zeromag",
                                    "The flux corresponding to zero magnitude.")
_zeroMagErrField = afwTable.Field["F"]("flux.zeromag.err",
                                       "The flux error corresponding to zero magnitude.")
_stellarField = afwTable.Field["Flag"]("stellar",
                                       "If true, the object is known to be a star if false it's known not to be a star.")
_magAutoField = afwTable.Field["F"]("mag.auto",
                                    "The magnitude computed by SExtractor in the HST catalog.")
_seeingField = afwTable.Field["F"]("seeing",
                                    "The PSF FWHM.")
_exptimeField = afwTable.Field["F"]("exptime",
                                    "Exposure time.")
_idField = afwTable.Field["L"]("multId",
                               "Multiple id, this field is in place to keep track of the ids of the matches in other catalogs/bands")
_dGaussRadInner = afwTable.Field["F"]("dGauss.radInner",
                                      "Determinant radius of the inner Gaussian in the double Gaussian fit to the PSF.")
_dGaussRadOuter = afwTable.Field["F"]("dGauss.radOuter",
                                      "Determinant radius of the outer Gaussian in the double Gaussian fit to the PSF.")
_dGaussAmpRat = afwTable.Field["F"]("dGauss.ampRat",
                                    "Peak amplitudes ratio of the two Gaussians in the double Gaussian fit to the PSF.")
_dGaussQInner = afwTable.Field["F"]("dGauss.qInner",
                                    "Ellipticity of the inner Gaussian in the double Gaussian fit to the PSF.")
_dGaussQOuter = afwTable.Field["F"]("dGauss.qOuter",
                                    "Ellipticity of the outer Gaussian in the double Gaussian fit to the PSF.")
_dGaussThetaInner = afwTable.Field["F"]("dGauss.thetaInner",
                                        "Inclination angle of the inner Gaussian in the double Gaussian fit to the PSF.")
_dGaussThetaOuter = afwTable.Field["F"]("dGauss.thetaOuter",
                                        "Inclination angle of the outer Gaussians in the double Gaussian fit to the PSF.")


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

def _bandOrder(suffix):
    if suffix == 'g':
        return 1
    if suffix == 'r':
        return 2
    if suffix == 'i':
        return 3
    if suffix == 'z':
        return 4
    if suffix == 'y':
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
    
def getCatBands(cat):
    bands = []
    for schemaItem in cat.getSchema():
        fieldName = schemaItem.getField().getName()
        match = _bandRegex.search(fieldName)
        if match:
            band = match.group(1)[-1]
            if band not in bands:
                bands.append(band) 
    bands.sort(key=_bandOrder)
    return bands

def createSchemaMapper(cat, cat2=None, filterSuffix=None, withZeroMagFlux=False,
                       withStellar=False, withSeeing=False, withExptime=False, withDGaussPsf=False):

    if cat2 is not None and filterSuffix:
        raise ValueError("Can't use filterSuffix for two catalogs")

    suffixes = getCatSuffixes(cat)
    if len(suffixes) > 0 and filterSuffix is not None:
        raise ValueError("Can't add a suffix to a catalog that already has suffixes")

    schema = cat.getSchema()
    scm = afwTable.SchemaMapper(schema)

    # First fixed fields and patterns
    for f in _fixedFields: 
        scm.addMapping(schema.find(f).getKey())
    for p in _fixedPatterns:
        for f in schema.extract(p):
            scm.addMapping(schema.find(f).getKey())

    # Now suffixable fields and patterns
    if filterSuffix is not None:
        suffix = _getFilterSuffix(filterSuffix)
        scm.addOutputField(_idField.copyRenamed("multId"+suffix))
    for f in _suffixableFields:
        if filterSuffix is not None:
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

    if cat2 is not None:
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
        scm.addOutputField(_magAutoField)

    if withDGaussPsf:
        if filterSuffix:
            scm.addOutputField(_dGaussRadInner.copyRenamed("dGauss.radInner"+suffix))
            scm.addOutputField(_dGaussRadOuter.copyRenamed("dGauss.radOuter"+suffix))
            scm.addOutputField(_dGaussAmpRat.copyRenamed("dGauss.ampRat"+suffix))
            scm.addOutputField(_dGaussQInner.copyRenamed("dGauss.qInner"+suffix))
            scm.addOutputField(_dGaussQOuter.copyRenamed("dGauss.qOuter"+suffix))
            scm.addOutputField(_dGaussThetaInner.copyRenamed("dGauss.thetaInner"+suffix))
            scm.addOutputField(_dGaussThetaOuter.copyRenamed("dGauss.thetaOuter"+suffix))
        else:
            if len(suffixes) == 0:
                scm.addOutputField(_dGaussRadInner)
                scm.addOutputField(_dGaussRadOuter)
                scm.addOutputField(_dGaussAmpRat)
                scm.addOutputField(_dGaussQInner)
                scm.addOutputField(_dGaussQOuter)
                scm.addOutputField(_dGaussThetaInner)
                scm.addOutputField(_dGaussThetaOuter)
            else:
                for s in suffixes:
                    scm.addOutputField(_dGaussRadInner.copyRenamed("dGauss.radInner"+s))
                    scm.addOutputField(_dGaussRadOuter.copyRenamed("dGauss.radOuter"+s))
                    scm.addOutputField(_dGaussAmpRat.copyRenamed("dGauss.ampRat"+s))
                    scm.addOutputField(_dGaussQInner.copyRenamed("dGauss.qInner"+s))
                    scm.addOutputField(_dGaussQOuter.copyRenamed("dGauss.qOuter"+s))
                    scm.addOutputField(_dGaussThetaInner.copyRenamed("dGauss.thetaInner"+s))
                    scm.addOutputField(_dGaussThetaOuter.copyRenamed("dGauss.thetaOuter"+s))

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

def strictMatch(cat1, cat2, matchRadius=1*afwGeom.arcseconds, includeMismatches=True,
                multiMeas=False):
    """
    Match two catalogs using a one to one relation where each match is the closest
    object
    """
    
    mc = afwTable.MatchControl()
    mc.includeMismatches = includeMismatches
    mc.findOnlyClosest = True

    #matched = afwTable.matchRaDec(cat1, cat2, matchRadius, True)
    matched = afwTable.matchRaDec(cat1, cat2, matchRadius, mc)

    bestMatches = {}
    noMatch = []
    for m1, m2, d in matched:
        if m2 is None:
            noMatch.append(m1)
        else:
            if not multiMeas:
                id = m2.getId()
                if id not in bestMatches:
                    bestMatches[id] = (m1, m2, d)
                else:
                    if d < bestMatches[id][2]:
                        bestMatches[id] = (m1, m2, d)
            else:
                id = m1.getId()
                bestMatches[id] = (m1, m2, d)

    if includeMismatches:
        print "{0} objects from {1} in the first catalog had no match in the second catalog.".format(len(noMatch), len(cat1))
        print "{0} objects from the first catalog with a match in the second catalog were not the closest match.".format(len(matched) - len(noMatch) - len(bestMatches))

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

def matchMultiBand(butler, dataType, filters=['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y'],
                   multiMeas=False, quick=False, **kargs):
    cats = []
    for f in filters:
        if quick:
            ids = butler.fetchIds(dataType, filter=f)
            cat = butler.fetchDataset(dataType, filterSuffix=f, **ids[0])
        else:
            cat = butler.fetchDataset(dataType, filterSuffix=f, filter=f, **kargs)
        cats.append(cat)

    matched = cats[0]

    for i in range(1, len(filters)):
        matched = strictMatch(matched, cats[i], multiMeas=multiMeas)
    
    return matched

def matchCats(cat1, cat2, matchRadius=1*afwGeom.arcseconds, includeMismatches=False,
              multiMeas=False, suffix='.2'):
    """
    Match to catalogs and return a catalog with the fields of the two catalogs
    """

    mc = afwTable.MatchControl()
    mc.includeMismatches = includeMismatches
    mc.findOnlyClosest = True

    matched = afwTable.matchRaDec(cat1, cat2, matchRadius, mc)

    haveCentroid = {}
    for m1, m2, d in matched:
        haveCentroid[m1.getId()] = (m1, m2, d)

    bestMatches = {}
    if includeMismatches:
        noMatch = []
    for m1, m2, d in matched:
        if m2 is None:
            noMatch.append(m1)
        else:
            if not multiMeas:
                id = m2.getId()
                if id not in bestMatches:
                    bestMatches[id] = (m1, m2, d)
                else:
                    if d < bestMatches[id][2]:
                        bestMatches[id] = (m1, m2, d)
            else:
                id = m1.getId()
                bestMatches[id] = (m1, m2, d)

    if includeMismatches:
        print "{0} objects from {1} in the first catalog had no match in the second catalog.".format(len(noMatch), len(cat1))
        print "{0} objects from the first catalog with a match in the second catalog were not the closest match.".format(len(matched) - len(noMatch) - len(bestMatches))

    if includeMismatches and not multiMeas:
        nMatches = len(cat1)
        print "I found {0} matches".format(len(bestMatches))
    else:
        nMatches = len(bestMatches)
        print "I found {0} matches".format(nMatches)

    schema1 = cat1.getSchema(); schema2 = cat2.getSchema()
    names1 = cat1.schema.getNames(); names2 = cat2.schema.getNames()

    schema = afwTable.SimpleTable.makeMinimalSchema()

    catKeys = []; cat1Keys = []; cat2Keys = []
    for name in names1:
        cat1Keys.append(schema1.find(name).getKey())
        if name not in ['id', 'coord']:
            catKeys.append(schema.addField(schema1.find(name).getField()))
        else:
            catKeys.append(schema.find(name).getKey())
    for name in names2:
        cat2Keys.append(schema2.find(name).getKey())
        if name not in schema1.getNames():
            catKeys.append(schema.addField(schema2.find(name).getField()))
        else:
            catKeys.append(schema.addField(schema2.find(name).getField().copyRenamed(name+suffix)))

    cat = afwTable.SimpleCatalog(schema)
    cat.reserve(nMatches)

    if includeMismatches and not multiMeas:
        for m1 in cat1:
            id1 = m1.getId()
            record = cat.addNew()
            for i in range(len(cat1Keys)):
                record.set(catKeys[i], m1.get(cat1Keys[i]))
            if id1 in haveCentroid:
                m2 = haveCentroid[id1][1]
                if m2 is not None:
                    id2 = m2.getId()
                    if id2 in bestMatches:
                        if bestMatches[id2][0] == m1:
                            for i in range(len(cat1Keys), len(catKeys)):
                                record.set(catKeys[i], m2.get(cat2Keys[i-len(cat1Keys)]))
    else:
        for id in bestMatches:
            m1, m2, d = bestMatches[id]
            record = cat.addNew()
            for i in range(len(cat1Keys)):
                record.set(catKeys[i], m1.get(cat1Keys[i]))
            for i in range(len(cat1Keys), len(catKeys)):
                record.set(catKeys[i], m2.get(cat2Keys[i-len(cat1Keys)]))

    return cat

def buildPermissiveXY(butler, dataType, filters=['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y'],
                      multiMeas=False, quick=False,
                      selectSG="/tigress/garmilla/data/cosmos_sg_all.fits",
                      matchRadius=1*afwGeom.arcseconds, **kargs):

    patch = None

    sgTable = afwTable.SimpleCatalog.readFits(selectSG)
    sgTable["coord.ra"][:]  = np.radians(sgTable["coord.ra"])
    sgTable["coord.dec"][:] = np.radians(sgTable["coord.dec"])
    if quick:
        indexes = np.random.choice(len(sgTable), 10000, replace=False)
        quickCat = afwTable.SimpleCatalog(sgTable.getSchema())
        quickCat.reserve(10000)
        for i in indexes:
            record = quickCat.addNew()
            record.assign(sgTable[i])
        sgTable = quickCat
    matches = []; idKeys = []

    for f in filters:
        if quick:
            if patch is None:
                ids = butler.fetchIds(dataType, filter=f)
                patch = ids[0]['patch']
            cat = butler.fetchDataset(dataType, filterSuffix=f, filter=f, patch=patch, **kargs)
        else:
            cat = butler.fetchDataset(dataType, filterSuffix=f, filter=f, **kargs)

        suffix = _getFilterSuffix(f)
        match = matchCats(sgTable, cat, matchRadius=matchRadius, multiMeas=False,
                          includeMismatches=True, suffix=suffix)
        matches.append(match)
        idKeys.append(match.schema.find('id').getKey())

    nMatches = len(sgTable)

    outputSchema = afwTable.SimpleTable.makeMinimalSchema()
    outputSchema.addField(sgTable.schema.find('mag.auto').getField())
    outputSchema.addField(sgTable.schema.find('mu.class').getField())

    schema = afwTable.SimpleTable.makeMinimalSchema()

    for match in matches:
        schema = match.getSchema()
        names = schema.getNames()
        for name in names:
            if name not in ['id', 'coord', 'mag.auto', 'mu.class']:
                outputSchema.addField(schema.find(name).getField())

    outputCat = afwTable.SimpleCatalog(outputSchema)
    outputCat.reserve(nMatches)
    for i in range(nMatches):
        outputCat.addNew()

    outputCat.get('id')[:] = sgTable.get('id')
    outputCat.get('coord.ra')[:] = sgTable.get('coord.ra')
    outputCat.get('coord.dec')[:] = sgTable.get('coord.dec')
    for match in matches:
        schema = match.getSchema()
        names = schema.getNames()
        for name in names:
            if name not in ['id', 'coord']:
                try:
                    outputCat.get(name)[:] = match.get(name)
                except LsstCppException:
                    keyOut = outputCat.schema.find(name).getKey()
                    keyMatch = match.schema.find(name).getKey()
                    for i in range(nMatches):
                        outputCat[i].set(keyOut, match[i].get(keyMatch))
                except ValueError:
                    keyOut = outputCat.schema.find(name).getKey()
                    keyMatch = match.schema.find(name).getKey()
                    for i in range(nMatches):
                        outputCat[i].set(keyOut, match[i].get(keyMatch))

    return outputCat

def buildXY(hscCat, sgTable, matchRadius=1*afwGeom.arcseconds, includeMismatches=True,
            multiMeas=False):

    mc = afwTable.MatchControl()
    mc.includeMismatches = includeMismatches
    mc.findOnlyClosest = True

    print "Matching with HST catalog"
    matchedSG = afwTable.matchRaDec(hscCat, sgTable, matchRadius, mc)
    print "Found {0} matches with HST objects".format(len(matchedSG))
    
    # Build truth table
    stellar = {}
    classKey = sgTable.getSchema().find('mu.class').key
    magAutoKey = sgTable.getSchema().find('mag.auto').key
    noMatch = []
    for m1, m2, d in matchedSG:
        if m2 is None:
            noMatch.append(m1.getId())
        else:
            if not multiMeas:
                id = m2.getId()
                isStar = (m2.get(classKey) == 2)
                magAuto = m2.get(magAutoKey)
                if id not in stellar:
                    stellar[id] = [isStar, magAuto, d, m1]
                else:
                    if d < stellar[id][2]:
                        stellar[id] = [isStar, magAuto, d, m1] # Only keep closest for now
            else:
                id = m1.getId()
                isStar = (m2.get(classKey) == 2)
                magAuto = m2.get(magAutoKey)
                stellar[id] = [isStar, magAuto, d, m1]

    if includeMismatches:
        print "{0} objects from {1} in the HSC catalog had no match in the HST catalog.".format(len(noMatch), len(hscCat))
        print "{0} objects from the HSC catalog with a match in the HST catalog were not the closest match.".format(len(matchedSG) - len(noMatch) - len(stellar))

    print "Of which I picked {0}".format(len(stellar)) 

    scm = createSchemaMapper(hscCat, withStellar=True)
    schema = scm.getOutputSchema()
    cat = afwTable.SourceCatalog(schema)
    cat.reserve(len(stellar))
    stellarKey = schema.find('stellar').key
    magAutoKey = schema.find('mag.auto').key

    for id in stellar:
        isStar, magAuto, d, m2 = stellar[id]
        record = cat.addNew()
        record.assign(m2, scm)
        record.set(stellarKey, isStar)
        record.set(magAutoKey, magAuto)

    if includeMismatches:
        return cat, noMatch

    return cat

def getNoMatchCat(butler, dataType, filters=['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y'],
                  quick=False, selectSG="/tigress/garmilla/data/cosmos_sg_all.fits",
                  matchRadius=1*afwGeom.arcseconds, mode='hsc', **kargs):

    mc = afwTable.MatchControl()
    mc.includeMismatches = True
    mc.findOnlyClosest = True

    sgTable = afwTable.SimpleCatalog.readFits(selectSG)
    sgTable["coord.ra"][:]  = np.radians(sgTable["coord.ra"])
    sgTable["coord.dec"][:] = np.radians(sgTable["coord.dec"])

    outputCats = []
    for f in filters:
        cat = butler.fetchDataset(dataType, filterSuffix=f, filter=f, **kargs)
        if mode == 'hsc':
            schema = cat.getSchema()
        elif mode == 'hst':
            schema = sgTable.getSchema()
        outputCat = afwTable.SimpleCatalog(schema)
        scm = afwTable.SchemaMapper(schema)
        suffix = _getFilterSuffix(f)
        if mode == 'hsc':
            matched = afwTable.matchRaDec(cat, sgTable, matchRadius, mc)
        elif mode == 'hst':
            matched = afwTable.matchRaDec(sgTable, cat, matchRadius, mc)
        for m1, m2, d in matched:
            if m2 is None:
                record = outputCat.addNew()
                record.assign(m1, scm)
        outputCats.append(outputCat)

    result = outputCats[0]
    for i in range(1, len(outputCats)):
        result.extend(outputCats[i], deep=False)

    return result

def buildCatFromIds(objIds, fsButler, dataType='deepCoadd'):
    info = utils.makeMapperInfo(fsButler.butler)
    cat = None
    for objId in objIds:
        if 'Coadd' in dataType or 'coadd' in dataType:
            dataId = info.splitCoaddId(objId)
        else:
            dataId = info.splitExposureId(objId)
        dataId.pop('objId')
        src = fsButler.butler.get(dataType+'_src', immediate=True, **dataId)
        record = src[objId == src.get("id")][0]
        if cat is None:
            schema = record.getSchema()
            cat = afwTable.SourceCatalog(schema)
            cat.reserve(len(objIds))
        cat.append(record)
    return cat

def getRecord(objId, fsButler, dataType='deepCoadd'):
    info = utils.makeMapperInfo(fsButler.butler)
    if 'Coadd' in dataType or 'coadd' in dataType:
        dataId = info.splitCoaddId(objId)
    else:
        dataId = info.splitExposureId(objId)
    dataId.pop('objId')
    src = fsButler.butler.get(dataType+'_src', immediate=True, **dataId)
    record = src[objId == src.get("id")][0]
    return record

def getParent(objId, fsButler, dataType='deepCoadd'):
    record = getRecord(objId, fsButler, dataType='deepCoadd')
    info = utils.makeMapperInfo(fsButler.butler)
    parentId = record.getParent()
    if parentId == 0:
        print "This object has no parent"
        return None
    if 'Coadd' in dataType or 'coadd' in dataType:
        dataId = info.splitCoaddId(parentId)
    else:
        dataId = info.splitExposureId(parentId)
    dataId.pop('objId')
    src = fsButler.butler.get(dataType+'_src', immediate=True, **dataId)
    parent = src[objId == src.get("id")][0]
    return parent

def getMultId(cat):
   bands = getCatBands(cat)
   multIds = np.zeros((len(cat),), dtype=[(b, 'int64') for b in bands])
   for b in bands:
       multIds[b] = cat.get('multId.'+b)
   return multIds

def extractDoubleGaussianModel(ctrl, record):
    model = FitPsfModel(ctrl, record)
    shapelets = model.asMultiShapelet(record.getCentroid())
    shapelet1 = shapelets.getElements()[0]
    shapelet2 = shapelets.getElements()[1]
    ampRat = shapelet2.evaluate()(record.getCentroid())/shapelet1.evaluate()(record.getCentroid())
    ellipse1 = afwGeom.ellipses.SeparableConformalShearDeterminantRadius(shapelet1.getEllipse().getCore())
    ellipse2 = afwGeom.ellipses.SeparableConformalShearDeterminantRadius(shapelet2.getEllipse().getCore())
    radInner = ellipse1.getRadius().getValue()
    radOuter = ellipse2.getRadius().getValue()
    qInner = ellipse1.getEllipticity().getAxisRatio()
    qOuter = ellipse2.getEllipticity().getAxisRatio()
    thetaInner = ellipse1.getEllipticity().getTheta()
    thetaOuter = ellipse2.getEllipticity().getTheta()
    return radInner, radOuter, ampRat, qInner, qOuter, thetaInner, thetaOuter

def displayObject(objId, fsButler, dataType='deepCoadd', suffix='', nPixel=15, frame=None):
    #TODO: Enable single exposure objects
    info = utils.makeMapperInfo(fsButler.butler)
    if 'Coadd' in dataType or 'coadd' in dataType:
        dataId = info.splitCoaddId(objId)
    else:
        dataId = info.splitExposureId(objId)
    dataId.pop('objId')
    src = fsButler.butler.get(dataType+'_src', immediate=True, **dataId)
    src = src[objId == src.get("id")][0]
    coord = src.getCoord()
    de = fsButler.butler.get(dataType, **dataId)
    pixel = de.getWcs().skyToPixel(coord)
    pixel = afwGeom.Point2I(pixel)
    bbox = afwGeom.Box2I(pixel, pixel)
    bbox.grow(nPixel)
    im = afwImage.ExposureF(de, bbox, afwImage.PARENT)
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

def getCoaddCutOut(fsButler, ra, dec, nPixel=15, filter='HSC-R'):
    """
    Return a cutout centered at `ra` `dec` (equatiorial coordinates in radians) with a bounding
    box size of nPixel x nPixel in filter `filter`. This function also returns a catalog of all
    the records with centroids contained in the bounding box.
    """

    catSubImage = []

    skyMap = fsButler.butler.get("deepCoadd_skyMap", immediate=True)
    coord = afwCoord.IcrsCoord(afwGeom.Point2D(ra, dec), afwGeom.radians)

    for tract, patch in skyMap.findClosestTractPatchList([coord]):
        coadd = fsButler.butler.get("deepCoadd", tract=tract.getId(),
                                    patch="%d,%d" % patch[0].getIndex(),
                                    filter=filter, immediate=True)
        pixel = afwGeom.Point2I(coadd.getWcs().skyToPixel(coord))
        bbox = afwGeom.Box2I(pixel, pixel)
        bbox.grow(nPixel)
        bbox.clip(coadd.getBBox(afwImage.PARENT))
        if bbox.isEmpty():
            continue
        origin = afwGeom.Extent2D(bbox.getMin())
        subImage = afwImage.ExposureF(coadd, bbox, afwImage.PARENT)
        cat = fsButler.butler.get("deepCoadd_src", tract=tract.getId(),
                                  patch="%d,%d" % patch[0].getIndex(),
                                  filter=filter, immediate=True)
        for record in cat:
            centroid = afwGeom.PointI(record.getCentroid())
            if bbox.contains(centroid):
                catSubImage.append(record)

    return subImage, origin, catSubImage

def displayCutout(subImage, origin=None, catSubImage=None, frame=0):
    """
    Display a subImage. If catSubImage is not None then the centroids of the records in the catalog
    are also displayed, this catalog is meant to be the catalog of records in the subImage's
    bounding box.
    """

    ds9.mtv(subImage, frame=frame)

    if origin is not None and catSubImage is not None:
        for s in catSubImage:
            centroid = s.getCentroid() - origin
            if s.get('deblend.nchild') == 0:
                ds9.dot('+', centroid[0], centroid[1], ctype=ds9.GREEN, frame=frame)
            else:
                ds9.dot('+', centroid[0], centroid[1], ctype=ds9.RED, frame=frame)
            if s.getParent() == 0:
                ds9.dot('o', centroid[0], centroid[1], ctype=ds9.CYAN, frame=frame)
