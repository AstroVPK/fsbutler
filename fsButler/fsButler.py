import os

"""
Butler class that works around the LSST butler by looking at the file system.
"""

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
    
    
    def fetchDataset(self, dataType='src', flags=None, immediate=True, **dataId):
        """
        Returns a list with all the data elemnts of type `dataType` that match the id `dataId`
    
        Keywords
        dataType: The type of data element we want to fetch
        flags: Flags for the data id
        immediate: If True the butler makes sure it returns the actual data and not a proxy of it
        dataId: Dictionary of keywords that specify a dataId, if None all the data elements of type
                `dataType will be returned. If the id is complete it will only return a single data
                element.
        """
    
        if dataType == 'src' or dataType == 'calexp_md':
            dataIds = self.singleExpIds(self.dataRoot, **dataId)
        elif dataType == 'deepCoadd_src' or dataType == 'deepCoadd_calexp_md':
            dataIds = self.deepCoaddIds(self.dataRoot, **dataId)
        else:
            raise ValueError("Data type {0} is not implemented".format(dataType))
    
        dataset = []
        for id in dataIds:
            if self.butler.datasetExists(dataType, **id):
                dataElement = self.butler.get(dataType, flags=flags, immediate=immediate, **id)
                dataset.append(dataElement)
            else:
                print "WARNING: The data id {0} does not exist".format(id)
    
        if len(dataset) == 0:
            return None

        if len(dataset) == 1:
            return dataset[0]
    
        return dataset
