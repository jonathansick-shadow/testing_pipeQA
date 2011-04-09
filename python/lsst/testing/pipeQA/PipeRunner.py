import lsst.pex.policy                   as pexPolicy
import lsst.daf.persistence              as dafPersist
import lsst.daf.base                     as dafBase
import lsst.afw.detection                as afwDet

class PipeRunner(object):
    """Manage a full run of data through Pipette."""
    
    def __init__(self):
        self.testdataList = []
        pass


    def addTestData(self, testdata):
        """Add a dataset by name.  It must exist within a directory on TESTBED_PATH."""
        self.testdataList.append(testdata)

        
    def run(self, **kwargs):
        """Run all the data we have through Pipette."""
        for testdata in self.testdataList:
            testdata.run(kwargs)

    def getUncaughtExceptionDict(self):
        """Get a dictionary of all exceptions coming out of Pipette."""
        # merge exceptionDicts from all test data sets into one dictionary
        uncaughtExceptionDict = {}
        for testdata in self.testdataList:
            uncaughtExceptionDict = dict(uncaughtExceptionDict.items() +
                                         testdata.getUncaughtExceptionDict().items())
        return uncaughtExceptionDict

    
    def getLogFiles(self):
        """Get output log files written by Pipette."""
        fileList = []
        for testdata in self.testdataList:
            fileList += testdata.getLogFiles()
        return fileList

    
    def getEupsSetupFiles(self):
        """Get files listing the EUPS setups used when Pipette ran."""
        fileList = []
        for testdata in self.testdataList:
            fileList += testdata.getEupsSetupFiles()
        return fileList
    
    

    def getSourceSet(self, dataId):
        """Get sources for requested data as one sourceSet. """
        
        ss = []
        for testdata in self.testdataList:
            ss += testdata.getSourceSet(dataId)
        return ss

    
