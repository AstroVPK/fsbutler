import os
import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

from . import utils

"""
Butler class that works around the LSST butler by looking at the file system.
"""

def _concatenateCats(cats):
    """
    Concatenate a list of catalogs into a single catalog
    """

    cat = cats[0]

    for i in range(1, len(cats)):
        cat.extend(cats[i], deep=False)

    return cat

class fsButler(object):

    _filters = ['HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'HSC-Y']

    def __init__(self, dataRoot):
        
        self.dataRoot = dataRoot

        import lsst.daf.persistence as dafPersist
        self.butler = dafPersist.Butler(dataRoot)

    @staticmethod
    def singleExpIds(dataRoot, filter=None, visit=None, ccd=None):
        """
        Returns a list of data Ids for single exposures in a given dataRoot directory
    
        Arguments
        dataRoot: The directory where the rerun data is stored
    
        Keywords
        filter: If None return the ids with all the filters, if present return only the ids
                with a given filter
        visit: If None return the ids with all the visits, if present return only the ids
                with a given visit 
        ccd: If None return the ids with all the ccds, if present return only the ids
                with a given ccd
        """
    
        if filter and visit and ccd:
            return [{'filter' : filter, 'visit' : visit, 'ccd' : ccd}]

        import re
    
        subDirs = os.listdir(dataRoot)
        dataDirs = []
        for d in subDirs:
            if d.isdigit():
                dataDirs.append(d)
    
        dataIds = []
        for d in dataDirs:
            cDir = os.path.join(dataRoot, d)
            if filter:
                filterPath = os.path.join(cDir, filter)
                if os.path.isdir(filterPath):
                    dirFilters = [filter]
                else:
                    continue
            else:
                dirFilters = os.listdir(cDir)
            for df in dirFilters:
                outputPath = os.path.join(os.path.join(cDir, df), 'output')
                outputs = os.listdir(outputPath)
                for o in outputs:
                    if re.match(r"^SRC-", o):
                        match = re.match(r"^SRC-([0-9]{7})-([0-9]{3}).fits", o)    

                        if match == None:
                            print "WARNING: Failed to read visit and ccd numbers from\
                                  {0}".format(os.path.join(outputPath, o))
                            continue

                        visitNumber = int(match.group(1))
                        ccdNumber = int(match.group(2))
                        if visit == None and ccd == None:
                            dataIds.append({'filter' : df, 'visit' : visitNumber, 'ccd' : ccdNumber})
                        elif visit != None and ccd == None:
                            if visitNumber == visit:
                                dataIds.append({'filter' : df, 'visit' : visitNumber, 'ccd' : ccdNumber})
                        elif visit == None and ccd != None:
                            if ccdNumber == ccd:
                                dataIds.append({'filter' : df, 'visit' : visitNumber, 'ccd' : ccdNumber})
                        elif visit != None and ccd != None:
                            if visitNumber == visit and ccdNumber == ccd:
                                dataIds.append({'filter' : df, 'visit' : visitNumber, 'ccd' : ccdNumber})

        return dataIds
    
    @staticmethod
    def deepCoaddIds(dataRoot, filter=None, tract=None, patch=None):
        """
        Returns a list of data Ids for deep coadds in a given dataRoot directory
    
        Arguments
        dataRoot: The directory where the rerun data is stored
    
        Keywords
        filter: If None return the ids with all the filters, if present return only the ids
                with a given filter
        tract: If None return the ids with all the tracts, if present return only the ids
                with a given tract
        patch: If None return the ids with all the patches, if present return only the ids
                with a given patch
        """
    
        if filter and tract and patch:
            return [{'filter' : filter, 'tract' : tract, 'patch' : patch}]

        deepCoaddPath = os.path.join(dataRoot,'deepCoadd-results')
    
        if filter:
            filters = [filter]
        else:
            filters = os.listdir(deepCoaddPath)
    
        dataIds = []
        for f in filters:
            filterPath = os.path.join(deepCoaddPath, f)
    
            if tract:
                tracts = [tract]
            else: 
                tracts = os.listdir(filterPath)
    
            for t in tracts:
                tractPath = os.path.join(filterPath, t)
    
                if patch:
                    patches = [patch]
                else:
                    patches = os.listdir(tractPath)
    
                for p in patches:
                    dataIds.append({'filter' : f, 'tract' : int(t), 'patch' : p})
    
        return dataIds
    
    @staticmethod
    def getIds(dataRoot, dataType, **dataId):
        """
        Returns a list of data Ids for data of type `dataType` in a given dataRoot directory
    
        Arguments
        dataRoot: The directory where the rerun data is stored
        dataType: The type of data element we want to fetch
    
        Keywords
        dataId: Dictionary of keywords that specify a dataId, if None all the data elements of type
                `dataType` will be returned. If the id is complete it will only return a single data
                element.
        """

        if dataType == 'src' or dataType == 'calexp_md':
            dataIds = fsButler.singleExpIds(dataRoot, **dataId)
        elif dataType == 'deepCoadd' or dataType == 'deepCoadd_src' or\
             dataType == 'deepCoadd_calexp_md' or dataType == 'deepCoadd_meas':
            dataIds = fsButler.deepCoaddIds(dataRoot, **dataId)
        else:
            raise ValueError("Data type {0} is not implemented".format(dataType))

        return dataIds

    def fetchIds(self, dataType, **dataId):
        """
        Returns a list of data Ids for data of type `dataType` in a given dataRoot directory
    
        Arguments
        dataType: The type of data element we want to fetch
    
        Keywords
        dataId: Dictionary of keywords that specify a dataId, if None all the data elements of type
                `dataType` will be returned. If the id is complete it will only return a single data
                element.
        """

        return self.getIds(self.dataRoot, dataType, **dataId)

    @staticmethod
    def _getCalexpType(dataType):
        if dataType == 'src':
            return 'calexp'
        elif dataType == 'deepCoadd_src' or dataType == 'deepCoadd_meas':
            return 'deepCoadd_calexp_md'
        elif dataType == 'calexp_md' or dataType == 'deepCoadd_calexp_md':
            return dataType
        elif dataType == 'calexp' or dataType == 'deepCoadd':
            return dataType
        else:
            raise ValueError("Unkown dataType")

    def _getCoaddExpTime(self, coadd):
        """
        Get the total exposure time associated with a coadd patch by adding up
        the exposure times of the individual visits used to build the coadd.
        """
        ccdInputs = coadd.getInfo().getCoaddInputs().ccds

        totalExpTime = 0.0
        expTimeVisits = set()

        for i in range(len(ccdInputs)):
            input = ccdInputs[i]
            ccd = input.get("ccd")
            v = input.get("visit")
            bbox = input.getBBox()
            # It's quicker to not read all the pixels, so just read 1
            calexp = self.butler.get("calexp_sub", bbox=afwGeom.BoxI(afwGeom.PointI(0, 0), afwGeom.ExtentI(1, 1)),
                                visit=int(v), ccd=ccd)
            calib = calexp.getCalib()
            exptime = calib.getExptime()
            if v not in expTimeVisits:
                totalExpTime += exptime
                expTimeVisits.add(v)
        return totalExpTime

    def _getCalibData(self, dataType, **id):
        """
        Get the zero magnitude flux for the given data type and data id
        It also returns the psf associated with the data id in case it's needed
        later (e.g. when withSeeing=True)
        """
        calexpType = self._getCalexpType(dataType)
        if self.butler.datasetExists(calexpType, **id):
            calexp = self.butler.get(calexpType, **id)
            try:
                calib = afwImage.Calib(calexp)
            except NotImplementedError:
                calib = calexp.getCalib()
        else:
            calexpType = dataType[:-4]
            if calexpType[-1] == '_':
                calexpType = calexpType[:-1]
            calexp = self.butler.get(calexpType, **id)
            calib = calexp.getCalib()
        psf = calexp.getPsf()
        fluxMag0, fluxMag0Err = calib.getFluxMag0()
        if 'deepCoadd' in dataType:
            exptime = self._getCoaddExpTime(calexp)
        else:
            exptime = calib.getExptime()
        return fluxMag0, fluxMag0Err, psf, exptime

    def fetchDataset(self, dataType='src', flags=None, immediate=True, withZeroMagFlux=True,
                     filterSuffix=None, scm=None, withSeeing=True, seeingAtPos=False, 
                     withExptime=True, withDGaussPsf=False, **dataId):
        """
        Returns the union of all the data elements of type `dataType` that match the id `dataId`
    
        Keywords
        dataType: The type of data element we want to fetch
        flags: Flags for the data id
        immediate: If True the butler makes sure it returns the actual data and not a proxy of it
        withZeroMagFlux: If true, an extra pair of columns will be added to the catalog for the zero
                         magnitude flux and its error estimate.
        filterSuffix: If present, a suffix corresponding to the filter will be appended to suffixable
                      fields. For example if dataId has `filter='HSC-I'`, you can set
                      `filterSiffix=i` to append a `.i` to all the column names.
        scm: If None, generate a new schema mapper
        dataId: Dictionary of keywords that specify a dataId, if None all the data elements of type
                `dataType` will be returned. If the id is complete it will only return a single data
                element.
        withSeeing: If true add a column for the FWHM of the PSF at each source's location.
        seeingAtPos: If true compute the seeing at each object's position, if False simply use
                     the seeing at the average position of the patch (in the case of coadds) or
                     ccd (in the case of single exposures).
        withExptime: If true add a column for the exposure time.
        withDGaussPSF: If true carry out a double Gaussian fit to the PSF and store the parameters.
        """
    
        dataIds = self.getIds(self.dataRoot, dataType, **dataId)
    
        if filterSuffix is not None:
            suffix = utils._getFilterSuffix(filterSuffix)

        dataset = []
        if withDGaussPsf:
            cfg = self.butler.get("multiband_config", immediate=True)
            ctrl = cfg.measureCoaddSources.measurement.algorithms["multishapelet.psf"].makeControl()
        for id in dataIds:
            if self.butler.datasetExists(dataType, **id):
                if flags is not None:
                    dataElement = self.butler.get(dataType, flags=flags, immediate=immediate, **id)
                else:
                    dataElement = self.butler.get(dataType, immediate=immediate, **id)
                if withZeroMagFlux or withSeeing or withExptime:
                    fluxMag0, fluxMag0Err, psf, exptime = self._getCalibData(dataType, **id)
                    # Compute seeing at average position
                    seeingAvgPos = psf.computeShape().getDeterminantRadius()
                if isinstance(dataElement, afwTable.SourceCatalog):
                    if scm == None:
                        scm = utils.createSchemaMapper(dataElement, filterSuffix=filterSuffix,
                                                       withZeroMagFlux=withZeroMagFlux,
                                                       withSeeing=withSeeing,
                                                       withExptime=withExptime,
                                                       withDGaussPsf=withDGaussPsf)
                    outputSchema = scm.getOutputSchema()
                    outputCat = afwTable.SimpleCatalog(outputSchema)
                    good = utils.goodSources(dataElement)
                    outputCat.reserve(np.sum(good))
                    if filterSuffix is not None:
                        idKey = outputSchema.find('multId'+suffix).key
                    if withZeroMagFlux:
                        if filterSuffix is not None:
                            zeroKey = outputSchema.find('flux.zeromag'+suffix).key
                            zeroErrKey = outputSchema.find('flux.zeromag.err'+suffix).key
                        else:
                            zeroKey = outputSchema.find('flux.zeromag').key
                            zeroErrKey = outputSchema.find('flux.zeromag.err').key
                    if withSeeing:
                        if filterSuffix is not None:
                            seeingKey = outputSchema.find('seeing'+suffix).key
                        else:
                            seeingKey = outputSchema.find('seeing').key
                    if withExptime:
                        if filterSuffix is not None:
                            exptimeKey = outputSchema.find('exptime'+suffix).key
                        else:
                            exptimeKey = outputSchema.find('exptime').key
                    if withDGaussPsf:
                        if filterSuffix is not None:
                            dGaussRadInnerKey = outputSchema.find('dGauss.radInner'+suffix).key
                            dGaussRadOuterKey = outputSchema.find('dGauss.radOuter'+suffix).key
                            dGaussAmpRatKey = outputSchema.find('dGauss.ampRat'+suffix).key
                            dGaussQInnerKey = outputSchema.find('dGauss.qInner'+suffix).key
                            dGaussQOuterKey = outputSchema.find('dGauss.qOuter'+suffix).key
                            dGaussThetaInnerKey = outputSchema.find('dGauss.thetaInner'+suffix).key
                            dGaussThetaOuterKey = outputSchema.find('dGauss.thetaOuter'+suffix).key
                        else:
                            dGaussRadInnerKey = outputSchema.find('dGauss.radInner').key
                            dGaussRadOuterKey = outputSchema.find('dGauss.radOuter').key
                            dGaussAmpRatKey = outputSchema.find('dGauss.ampRat').key
                            dGaussQInnerKey = outputSchema.find('dGauss.qInner').key
                            dGaussQOuterKey = outputSchema.find('dGauss.qOuter').key
                            dGaussThetaInnerKey = outputSchema.find('dGauss.thetaInner').key
                            dGaussThetaOuterKey = outputSchema.find('dGauss.thetaOuter').key
                    for i, record in enumerate(dataElement):
                        if good[i]:
                            outputRecord = outputCat.addNew()
                            outputRecord.assign(record, scm)
                            if filterSuffix is not None:
                                outputRecord.set(idKey, record.getId())
                            if withZeroMagFlux:
                                outputRecord.set(zeroKey, fluxMag0)
                                outputRecord.set(zeroErrKey, fluxMag0Err)
                            if withSeeing:
                                if seeingAtPos:
                                    pos = record.getCentroid()
                                    try:
                                        seeing = psf.computeShape(pos).getDeterminantRadius()
                                    except:
                                        seeing = psf.computeShape().getDeterminantRadius()
                                else:
                                    seeing = seeingAvgPos
                                outputRecord.set(seeingKey, seeing)
                            if withExptime:
                                outputRecord.set(exptimeKey, exptime)
                            if withDGaussPsf:
                                radInner, radOuter, ampRat, qInner, qOuter, thetaInner, thetaOuter\
                                = utils.extractDoubleGaussianModel(ctrl, record)
                                outputRecord.set(dGaussRadInnerKey, radInner)
                                outputRecord.set(dGaussRadOuterKey, radOuter)
                                outputRecord.set(dGaussAmpRatKey, ampRat)
                                outputRecord.set(dGaussQInnerKey, qInner)
                                outputRecord.set(dGaussQOuterKey, qOuter)
                                outputRecord.set(dGaussThetaInnerKey, thetaInner)
                                outputRecord.set(dGaussThetaOuterKey, thetaOuter)
                    dataset.append(outputCat)
                else:
                    dataset.append(dataElement)
            else:
                print "WARNING: The data id {0} does not exist for data type {1}".format(id,dataType)
    
        if len(dataset) == 0:
            return None

        if len(dataset) == 1:
            return dataset[0]
    
        return _concatenateCats(dataset)
