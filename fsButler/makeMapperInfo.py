from __future__ import print_function
import re
import numpy as np

import lsst.afw.image as afwImage

# Inherited from lsst.analysis.utils
# The lsst.analysis package may have become bit-rotted. Rather than fix the package, we have copied the bits
# of lsst.analysis.utils that are used by sgs-fsbutler into this file.


def makeMapperInfo(butler):
    """Return an object with extra per-mapper information (e.g. which fields fully specify an exposure)"""

    if not butler:
        return None

    mapper = butler.mapper

    class MapperInfo(object):
        @staticmethod
        def getColorterm(filterName):
            return None

        def getId(self, src, field="objId"):  # can't be static as it calls derived function splitId
            idDict = self.splitId(src.getId())

            return idDict[field] if field else idDict

        @staticmethod
        def canonicalFiltername(filterName):
            return filterName

        @staticmethod
        def idMask(dataId):
            return 0x0

    class LsstSimMapperInfo(MapperInfo):
        def __init__(self, Mapper):
            LsstSimMapperInfo.Mapper = Mapper

        @staticmethod
        def getFields(dataType):
            fields = ["visit", "filter", "raft", "sensor", ]
            if dataType == "raw":
                fields += ["snap", "channel", ]

            return fields

        @staticmethod
        def getTrimmableData():
            """Return a list of data products that needs to be trimmed"""
            return ("raw", "flat", "bias", "dark",)

        @staticmethod
        def dataIdToTitle(dataIds, rerunName=None):
            try:
                dataIds[0]
            except TypeError:
                dataIds = [dataIds]

            filters = set()
            sensors = set()
            rafts = set()
            visits = set()
            for dataId in dataIds:
                dataId = dataId.copy()
                for k, v in dataId.items():
                    if isinstance(v, np.int32):
                        dataId[k] = int(v)

                if dataId.get("sensor") is None:
                    dataId["sensor"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except Exception as err:
                        filters.add("?")
                    sensors.add("all")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except Exception as err:
                        filters.add("?")

                    try:
                        sensors.add(dataId["sensor"])
                    except TypeError:
                        for c in dataId["sensor"]:
                            sensors.add(c)

                if dataId.get("raft") is None:
                    did = dataId.copy()
                    did["raft"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **did)).getName())
                    except Exception as err:
                        filters.add("?")
                    rafts.add("all")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except Exception as err:
                        filters.add("?")

                    try:
                        rafts.add(dataId["raft"])
                    except TypeError:
                        for c in dataId["raft"]:
                            rafts.add(c)

                try:
                    visits.add(dataId["visit"])
                except TypeError:
                    for v in dataId["visit"]:
                        visits.add(v)

            sensors = sorted(list(sensors))
            rafts = sorted(list(rafts))
            visits = sorted(list(visits))
            filters = sorted(list(filters))

            if len(visits) > 1 and len(filters) > 1:
                print >> sys.stderr, \
                      "I don't know how to make a title out of multiple visits and filters: %s %s" % \
                      (visits, filters)
                visits = visits[0:1]

            title = "%s R%s S%s [%s]" % (getNameOfSet(visits),
                                         getNameOfSRSet(rafts, 5, ['0,0', '4,0', '0,4', '4,4']),
                                         getNameOfSRSet(sensors, 3), ", ".join(filters))
            if rerunName:
                title += " %s" % rerunName

            return title

        @staticmethod
        def exposureToStr(exposure):
            ccdId = cameraGeom.cast_Ccd(exposure.getDetector()).getId().getName()
            visit = exposure.getMetadata().get("OBSID")

            return "%s %s" % (visit, ccdId)

        assembleCcd = staticmethod(assembleCcdLsst)

        @staticmethod
        def getInButler(dataRoot, registry, butler=None):
            inputRoot = os.path.join(os.path.split(dataRoot)[0], "input")
            return dafPersist.ButlerFactory(mapper=LsstSimMapperInfo.Mapper(root=inputRoot,
                                                                            registry=registry)).create()

        @staticmethod
        def splitId(oid, asDict=True):
            """Split an ObjectId into visit, raft, sensor, and objId"""
            objId = int((oid & 0xffff) - 1)     # Should be the same value as was set by apps code
            oid >>= 16
            raftSensorId = oid & 0x1ff
            oid >>= 9
            visit = int(oid)

            raftId, sensorId = int(raftSensorId//10), int(raftSensorId%10)
            raft = "%d,%d" % (raftId//5, raftId%5)
            sensor = "%d,%d" % (sensorId//3, sensorId%3)

            if asDict:
                return dict(visit=visit, raft=raft, sensor=sensor, objId=objId)
            else:
                return visit, raft, sensor, objId

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def splitSdssCcdExposureId(oid, hasFilter=True, asDict=True):
        """Split an ObjectId into run, camcol, [filter], objId

    If hasFilter is True, the ObjectId encodes a filtername
        """
        nbits = 26                      # number of bits reserved for objId
        oid = long(oid)

        omask = 0xffffffffffffffff << nbits
        objId = int(oid & ~omask)     # Should be the same value as was set by apps code
        oid >>= nbits

        field = int(oid % 10000)
        oid //= 10000
        camcol = int(oid % 10)
        oid //= 10
        filter = int(oid % 10)
        oid //= 10
        run = int(oid)

        if hasFilter:
            filterNames = [k for k, v in butler.mapper.filterIdMap.items() if v == filter]
            try:
                filter = filterNames[0]
                assert len(filterNames) == 1
            except IndexError:
                raise RuntimeError("Invalid filter index %d" % filter)

            if asDict:
                return dict(run=run, camcol=camcol, filter=filter, field=field, objId=objId)
            else:
                return run, camcol, filter, field, objId
        else:
            if asDict:
                return dict(run=run, camcol=camcol, field=field, objId=objId)
            else:
                return run, camcol, field, objId

    def splitSdssCoaddId(oid, hasFilter=True, asDict=True):
        """Split an ObjectId into tract, patch, [filter], objId

    If hasFilter is True, the ObjectId encodes a filtername
        """
        nbits = 34                  # number of bits used by patch etc. part of ID
        if hasFilter:
            nbits += 3                  # add 3 bits for filters
        nbits = 64 - nbits          # length
        oid = long(oid)

        omask = 0xffffffffffffffff << nbits
        objId = int(oid & ~omask)     # Should be the same value as was set by apps code
        oid >>= nbits
        if hasFilter:
            filter = int(oid & 0x7)
            oid >>= 3
        patchY = int(oid & 0x1fff)
        oid >>= 13
        patchX = int(oid & 0x1fff)
        oid >>= 13
        tract = int(oid)

        patch = "%d,%d" % (patchX, patchY)

        if hasFilter:
            filterNames = [k for k, v in butler.mapper.filterIdMap.items() if v == filter]
            try:
                filter = filterNames[0]
                assert len(filterNames) == 1
            except IndexError:
                raise RuntimeError("Invalid filter index %d" % filter)

            if asDict:
                return dict(tract=tract, patch=patch, filter=filter, objId=objId)
            else:
                return tract, patch, filter, objId
        else:
            if asDict:
                return dict(tract=tract, patch=patch, objId=objId)
            else:
                return tract, patch, objId

    class SdssMapperInfo(MapperInfo):
        def __init__(self, Mapper):
            SdssMapperInfo.Mapper = Mapper

        @staticmethod
        def getFields(dataType):
            if _prefix_ in ("", "forced",) or dataType in ("coaddTempExp",):
                fields = ["run", "filter", "camcol"]

                if dataType not in ("flat",):
                    fields.append("field")
            elif _prefix_ in ("deepCoadd", "deepCoaddForced", "goodSeeingCoadd",):
                fields = ["patch", "tract", "filter"]
            else:
                raise RuntimeError("I don't know which fields I need to read %s data" % _prefix_)
                pass

            return fields

        @staticmethod
        def getTrimmableData():
            """Return a list of data products that needs to be trimmed"""
            return ("raw",)

        @staticmethod
        def photometricTransform(desiredBand, primaryMag, secondaryMag):
            """Return the primary/secondary magnitude transformed into the desiredBand"""
            return SdssMapperInfo._Colorterm.transformMags(desiredBand, primaryMag, secondaryMag)

        @staticmethod
        def dataIdToTitle(dataIds, rerunName=None):
            try:
                dataIds[0]
            except TypeError:
                dataIds = [dataIds]

            runs = set()
            filters = set()
            camcols = set()
            fields = set()
            for dataId in dataIds:
                dataId = dataId.copy()
                for k, v in dataId.items():
                    if isinstance(v, np.int32):
                        dataId[k] = int(v)

                if dataId.get("camcol") is None:
                    dataId["camcol"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except Exception as err:
                        filters.add(dataId.get("filter", "?"))
                    if _prefix_ in ("", "forced", "deepCoaddForced"):
                        camcols.add("(all)")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except Exception as err:
                        filters.add("?")

                    try:
                        camcols.add(dataId["camcol"])
                    except TypeError:
                        for c in dataId["camcol"]:
                            camcols.add(c)

                for k in ["run", "patch", "tract"]:
                    try:
                        runs.add(dataId[k])
                    except KeyError:
                        pass
                    except TypeError:
                        for f in dataId[k]:
                            fields.add(f)

                for k in ["field", ]:
                    try:
                        runs.add(dataId[k])
                    except KeyError:
                        pass
                    except TypeError:
                        for f in dataId["field"]:
                            fields.add(f)

            runs = sorted(list(runs))
            fields = sorted(list(fields))
            camcols = sorted(list(camcols))
            filters = sorted(list(filters))

            if len(runs) > 1 and len(filters) > 1:
                print >> sys.stderr, \
                      "I don't know how to make a title out of multiple runs and filters: %s %s" % \
                      (runs, filters)
                runs = runs[0:1]

            nameOfFilters = "".join(filters)
            if len(filters) > 1:
                nameOfFilters = "[%s]" % nameOfFilters
            title = "%s %s%s %s" % (getNameOfSet(runs), nameOfFilters, getNameOfSet(camcols),
                                    getNameOfSet(fields))
            if rerunName:
                title += " %s" % rerunName

            return title

        @staticmethod
        def exposureToStr(exposure):
            try:
                ccdId = cameraGeom.cast_Ccd(exposure.getDetector()).getId().getSerial()
                visit = re.sub(r"^SUPA", "", exposure.getMetadata().get("FRAMEID"))
            except AttributeError:
                return "??"

            return "%s %s" % (visit, ccdId)

        @staticmethod
        def getInButler(dataRoot, registry, butler=None):
            return butler

        @staticmethod
        def splitId(oid, hasFilter=True, asDict=True):
            """Split an ObjectId into run, camcol, [filter], field, objId or tract, patch, [filter], objId

        If hasFilter is True, the ObjectId encodes a filtername
            """

            if _prefix_ in ("goodSeeingCoadd",):
                return splitSdssCoaddId(oid, hasFilter=hasFilter, asDict=asDict)
            else:
                return splitSdssCcdExposureId(oid, hasFilter=hasFilter, asDict=asDict)

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    from lsst.meas.photocal.colorterms import Colorterm
    from lsst.obs.suprimecam.colorterms import colortermsData

    class SubaruMapperInfo(MapperInfo):
        def __init__(self, Mapper):
            SubaruMapperInfo.Mapper = Mapper

            if False:
                SubaruMapperInfo._Colorterm = Colorterm
                SubaruMapperInfo.getColorterm = lambda x, y: Colorterm.getColorterm(y)
                SubaruMapperInfo._Colorterm.setColorterms(colortermsData, "Hamamatsu")

        @staticmethod
        def getFields(dataType):
            if _prefix_ in ("",):
                fields = ["visit", "ccd"]
                if dataType not in ("flat",):
                    fields.append("filter")
            elif _prefix_ in ("deepCoadd", "deepCoaddForced", "forced",):
                fields = ["tract", "patch", "filter"]
            elif _prefix_ in ("stack",):
                fields = ["stack", "patch", "filter"]
            else:
                raise RuntimeError("I don't know which fields I need to read %s data" % _prefix_)

            return fields

        @staticmethod
        def getTrimmableData():
            """Return a list of data products that needs to be trimmed"""
            return ("raw",)

        @staticmethod
        def photometricTransform(desiredBand, primaryMag, secondaryMag):
            """Return the primary/secondary magnitude transformed into the desiredBand"""
            return SubaruMapperInfo._Colorterm.transformMags(desiredBand, primaryMag, secondaryMag)

        @staticmethod
        def dataIdToTitle(dataIds, rerunName=None):
            try:
                dataIds[0]
            except TypeError:
                dataIds = [dataIds]

            if _prefix_ == "stack":
                title = []
                for did in dataIds:
                    title.append("stack %(stack)d patch %(patch)d filter %(filter)s" % did)

                return "[%s]" % "], [".join(title)

            filters = set()
            ccds = set()
            visits = set()
            for dataId in dataIds:
                dataId = dataId.copy()
                for k, v in dataId.items():
                    if isinstance(v, np.int32):
                        dataId[k] = int(v)

                if dataId.get("ccd") is None:
                    dataId["ccd"] = 0
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except Exception as err:
                        filters.add(dataId.get("filter", "?"))
                    ccds.add("(all)")
                else:
                    try:
                        filters.add(afwImage.Filter(butler.get(dtName("calexp", True), **dataId)).getName())
                    except Exception as err:
                        filters.add("?")

                    try:
                        ccds.add(dataId["ccd"])
                    except TypeError:
                        for c in dataId["ccd"]:
                            ccds.add(c)
                try:
                    visits.add(dataId["visit"])
                except TypeError:
                    for v in dataId["visit"]:
                        visits.add(v)
                except KeyError:
                    pass

            ccds = sorted(list(ccds))
            filters = sorted(list(filters))
            visits = sorted(list(visits))

            if len(visits) > 1 and len(filters) > 1:
                print >> sys.stderr, \
                      "I don't know how to make a title out of multiple visits and filters: %s %s" % \
                      (visits, filters)
                visits = visits[0:1]

            title = "%s CCD%s [%s]" % (getNameOfSet(visits), getNameOfSet(ccds), ", ".join(filters))
            if rerunName:
                title += " %s" % rerunName

            return title

        @staticmethod
        def exposureToStr(exposure):
            try:
                ccdId = cameraGeom.cast_Ccd(exposure.getDetector()).getId().getSerial()
                visit = re.sub(r"^SUPA", "", exposure.getMetadata().get("FRAMEID"))
            except AttributeError:
                return "??"

            return "%s %s" % (visit, ccdId)

        assembleCcd = staticmethod(assembleCcdSubaru)

        @staticmethod
        def getInButler(dataRoot, registry, butler=None):
            return butler

        @staticmethod
        def idMask(dataId):
            return 0x0

        @staticmethod
        def splitId(oid, asDict=True):
            """Split an ObjectId into visit, ccd, and objId.
            See obs/subaru/python/lsst/obs/suprimecam/suprimecamMapper.py"""
            oid = np.array(oid, dtype='int64')
            objId = np.bitwise_and(oid, 0xffff)  # Should be the same value as was set by apps code
            oid = np.right_shift(oid, 22).astype('int32')

            if _prefix_ == "stack":
                print("Warning: not vectorized")
                nfilter = len(butler.mapper.filters)
                nPatches = 1000000L

                ifilter = oid % nfilter
                oid //= nfilter

                patch = oid % nPatches
                oid //= nPatches

                stack = int(oid)

                filter = [k for k, v in butler.mapper.filterIdMap.items() if v == ifilter][0]

                if asDict:
                    return dict(stack=stack, patch=patch, filter=filter, objId=objId)
                else:
                    return stack, patch, filter, objId

            else:
                oid = np.right_shift(oid, 10).astype('int32')
                ccd = oid % 10
                oid //= 10
                visit = oid

                if visit.size == 1:     # sqlite doesn't like numpy types
                    visit = int(visit)
                    ccd = int(ccd)
                    objId = int(objId)

                if asDict:
                    return dict(visit=visit, ccd=ccd, objId=objId)
                else:
                    return visit, ccd, objId

        @staticmethod
        def canonicalFiltername(filterName):
            mat = re.search(r"W-J-(.)", filterName)
            if mat:
                return mat.group(1)

            mat = re.search(r"W-S-(.)\+", filterName)
            if mat:
                return mat.group(1).lower()

            return filterName

    class SubaruMapperInfoMit(SubaruMapperInfo):
        def __init__(self, Mapper):
            SubaruMapperInfo.__init__(self, None)
            SubaruMapperInfoMit.Mapper = Mapper
            SubaruMapperInfo._Colorterm.setColorterms(colortermsData, "MIT")

    from lsst.obs.hsc.colorterms import colortermsData

    class HscMapperInfo(SubaruMapperInfo):
        def __init__(self, Mapper):
            SubaruMapperInfo.__init__(self, None)
            HscMapperInfo.Mapper = Mapper

            if False:
                HscMapperInfo._Colorterm.setColorterms(colortermsData, "Hamamatsu")

        @staticmethod
        def exposureToStr(exposure):
            try:
                ccdId = cameraGeom.cast_Ccd(exposure.getDetector()).getId().getSerial()
                visit = re.sub(r"^HSC", "", exposure.getMetadata().get("FRAMEID"))
            except AttributeError:
                return "??"

            return "%s %s" % (visit, ccdId)

        @staticmethod
        def splitId(oid, asDict=True):
            """Split an ObjectId into (visit, ccd, objId) or (tract, patch, [filter], objId)

        If hasFilter is True, the ObjectId encodes a filtername
            """

            if _prefix_ in ("deepCoadd", "deepCoaddForced", "forced",):
                return HscMapperInfo.splitCoaddId(oid, asDict=asDict, hasFilter=True)
            elif _prefix_ in ("chisqCoadd",):
                return HscMapperInfo.splitCoaddId(oid, asDict=asDict, hasFilter=False)
            elif _prefix_ in ("",):
                return HscMapperInfo.splitExposureId(oid, asDict=asDict)
            else:
                raise RuntimeError("Please teach HscMapperInfo how to process splitId on a %s" % _prefix_)

        @staticmethod
        def splitExposureId(oid, asDict=True):
            """Split an ObjectId (maybe an numpy array) into visit, ccd, and objId.
            See obs/subaru/python/lsst/obs/hscSim/hscMapper.py"""
            oid = np.array(oid, dtype='int64')
            objId = np.bitwise_and(oid, 2**32 - 1)  # Should be the same value as was set by apps code
            oid = np.right_shift(oid, 32).astype('int32')

            ccd = (oid % 200).astype('int32')
            oid //= 200
            visit = oid.astype('int32')

            if visit.size == 1:     # sqlite doesn't like numpy types
                visit = int(visit)
                ccd = int(ccd)
                objId = int(objId)

            if asDict:
                return dict(visit=visit, ccd=ccd, objId=objId)
            else:
                return visit, ccd, objId

        @staticmethod
        def splitCoaddId(oid, asDict=True, hasFilter=True):
            """Split an ObjectId (maybe an numpy array) into tract, patch, [filter], and objId.
            See obs/subaru/python/lsst/obs/hscSim/hscMapper.py"""
            mapper = HscMapperInfo.Mapper

            oid = np.array(oid, dtype='int64')
            objId = np.bitwise_and(oid, 2**mapper._nbit_id - 1)
            oid >>= mapper._nbit_id

            if hasFilter:
                filterId = np.bitwise_and(oid, 2**mapper._nbit_filter - 1).astype('int32')
                oid >>= mapper._nbit_filter

                filterName = np.empty(oid.size, "a6")

                if filterId.size == 1:
                    filterId = [int(filterId)]  # as you can't iterate over a length-1 np array

                for fid in set(filterId):
                    name = afwImage.Filter(int(fid)).getName()

                    filesystemName = "HSC-%s" % name.upper()  # name mapper needs
                    try:
                        afwImage.Filter(filesystemName)
                        name = filesystemName
                    except Exception as err:
                        pass

                    filterName[filterId == fid] = name
            else:
                filterName = None

            patchY = np.bitwise_and(oid, 2**mapper._nbit_patch - 1).astype('int32')
            oid >>= mapper._nbit_patch
            patchX = np.bitwise_and(oid, 2**mapper._nbit_patch - 1).astype('int32')
            oid >>= mapper._nbit_patch
            add = np.core.defchararray.add  # why isn't this easier to find?
            patch = add(add(patchX.astype(str), ","), patchY.astype(str))
            patch.shape = filterName.shape  # why do I have to do this?

            tract = oid.astype('int32')

            if oid.size == 1:     # sqlite doesn't like numpy types
                filterName = str(filterName[0])
                tract = int(tract)
                patch = str(patch[0])
                objId = int(objId)

            if asDict:
                return {"filter": filterName, "tract": tract, "patch": patch, "objId": objId}
            else:
                return filterName, tract, patch, objId

    # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    if isinstance(mapper, LsstSimMapper):
        return LsstSimMapperInfo(LsstSimMapper)
    elif isinstance(butler.mapper, SdssMapper):
        return SdssMapperInfo(SdssMapper)
    elif isinstance(butler.mapper, SuprimecamMapper):
        return SubaruMapperInfo(SuprimecamMapper)
    elif isinstance(butler.mapper, SuprimecamMapperMit):
        return SubaruMapperInfoMit(SuprimecamMapperMit)
    elif isinstance(butler.mapper, HscMapper):
        return HscMapperInfo(HscMapper)
    else:
        raise RuntimeError("Impossible mapper")
