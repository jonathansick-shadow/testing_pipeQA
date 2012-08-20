import sys, os, re

import QaDataUtils as qaDataUtils


###################################################
# Factory for QaData
###################################################
def makeQaData(label, rerun=None, retrievalType=None, camera=None, **kwargs):
    """Factory to make a QaData object for either Butler data, or Database data.

    @param label         identifier for the data - either a directory in TESTBED_PATH or a database name
    @param rerun         data rerun to retrieve
    @param retrievalType 'butler', 'db', or None (will search first for butler, then database)
    @param camera        Specify which camera is to be used
    """
    
    if retrievalType is None:
        
        # if there's only one possibility, use that
        
        # see if there's a testbed directory called 'label'
        # if TESTBED_PATH isn't set, skip this and assume it's a db
        validButler = False
        if os.environ.has_key('TESTBED_PATH'):
            testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label, raiseOnFailure=False)
            if (not testbedDir is None) and (not testdataDir is None):
                validButler = True

        # see if we can connect to a database with name 'label'
        # NOTE: must update if/when non-lsst databases get used
        validDb = True

        try:
            if not camera is None and re.search("^(suprimecam|hsc)", camera):
                from HscDatabaseQuery import DbInterface, DatabaseIdentity
                dbInterface = DbInterface(DatabaseIdentity(label))
            else:
                from DatabaseQuery import LsstSimDbInterface, DatabaseIdentity
                dbInterface = LsstSimDbInterface(DatabaseIdentity(label))
        except Exception, e:
            validDb = False

        if validButler and not validDb:
            retrievalType = 'butler'
        if validDb and not validButler:
            retrievalType = 'db'
        if validDb and validButler:
            raise Exception("The label "+label+" is present as both a testbed directory and a database."\
                            "Please specify retrievalType='butler', or retrievalType='db'.")
        if not validDb and not validButler:
            raise Exception("Unable to find "+label+" as a testbed directory or a database.")


    # handle specially requested camera via retrievalType
    #if re.search("(lsstsim|suprimecam|cfht|hsc|sdss)", retrievalType):
    #    return makeButlerQaData(label, rerun, camera=retrievalType, **kwargs)


    print "RetrievalType=", retrievalType
    print "camera=", camera
    
    if re.search("^[Bb]utler$", retrievalType):
        from ButlerQaData  import makeButlerQaData
        return makeButlerQaData(label, rerun, camera=camera, **kwargs)
    
    if re.search("^([Dd][Bb]|[Dd]atabase)$", retrievalType):


        import CameraInfo as qaCamInfo
        cameraInfos = {
    #       "cfht": qaCamInfo.CfhtCameraInfo(), # XXX CFHT camera geometry is currently broken following #1767
            "hsc" : qaCamInfo.HscCameraInfo(),
            "suprimecam": qaCamInfo.SuprimecamCameraInfo(),
            "suprimecam-old": qaCamInfo.SuprimecamCameraInfo(True),
            "sdss" : qaCamInfo.SdssCameraInfo(),
            "coadd"          : qaCamInfo.CoaddCameraInfo(),
            "lsstsim": qaCamInfo.LsstSimCameraInfo(),
            }


        cameraToUse = None
        if not camera is None:
            cameraToUse = cameraInfos[camera]
        else:
            cameraToUse = cameraInfos['lsstsim']
            camera = 'lsstsim'
            
        if re.search("^(hsc|suprimecam|suprimecam-old)$", camera):
            from HscDbQaData      import HscDbQaData
            return HscDbQaData(label, rerun, cameraToUse)
        else:
            from DbQaData      import DbQaData
            return DbQaData(label, rerun, cameraToUse)




