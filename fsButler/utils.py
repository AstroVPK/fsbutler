import numpy as np

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

"""
Utility functions to process data elements delivered by fsButler
"""

def createSchemaMapper(cat):
    # Get the table associated with the data
    table = cat.getTable()
    # Get the schema of the source catalog
    schema = table.getSchema()
    # Generate a schema mapper to build the schema of our output catalog
    scm = afwTable.SchemaMapper(schema)
    # Columns that we are interested in
    columns = ["id", "coord", "parent", "deblend.nchild", "classification.extendedness",
               "flux.kron.radius",
               "flags.pixel.bad",
               "flags.pixel.edge",
               "flags.pixel.interpolated.any",
               "flags.pixel.interpolated.center",
               "flags.pixel.saturated.any",
               "flags.pixel.saturated.center"]
    for f in schema.extract("flux.psf*"):
            columns.append(f)
    # Check if the catalog has composite model measurements
    hasCmodel = True if len(schema.extract("cmodel*")) > 0 else False
    if hasCmodel:
        for f in schema.extract("cmodel*"):
            columns.append(f)
    else:
        print "WARNING: This catalog has no cmodel measurements"
    # Add the mappings for the columns we are interested in
    for f in columns: 
        scm.addMapping(schema.find(f).getKey())
    # Add mappings for the zeroMag flux
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

#def buildXY(hscSources="/tigress/garmilla/data/sgClassCOSMOS.fits", selectSG="/tigress/garmilla/data/cosmos_sg_all.fits"):
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
