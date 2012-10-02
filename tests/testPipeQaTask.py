import os, shutil
import re
import unittest
import eups
import lsst.utils.tests as tests
from lsst.testing.pipeQA.analysis.PipeQaTask import PipeQaTask

class PipeQaTestCases(unittest.TestCase):
    """For testing purposes we will disable all tests except for ZptFit"""
    def setUp(self):
        self.qaTask       = PipeQaTask()
        self.testDatabase = "rplante_PT1_2_u_pt12prod_im3000"  # LSST PT 1.2 series 3000
        self.testVisit1   = "885335881"                        # r-band
        self.testVisit2   = "885335891"                        # r-band
        self.testFilt     = "r"
        self.testRaft     = "2,2"
        self.testCcd      = "1,1"

        self.testDatabase = "krughoff_S12_lsstsim_u_krughoff_2012_0706_183555"  # lsstSim S21 schema
        self.testVisit1   = "899551571"                        # z-band
        self.testVisit2   = "899553091"                        # r-band
        self.testFilt1    = "z"
        self.testFilt2    = "r"
        self.testRaft     = "2,2"
        self.testCcd      = "1,1"

        self.wwwRerun     = os.path.join(eups.productDir("testing_pipeQA"), "tests", "www_rerun")
        os.environ["WWW_RERUN"] = self.wwwRerun
        if os.path.isdir(self.wwwRerun):
            shutil.rmtree(self.wwwRerun)
            

    def disableTasks(self):
        disArgs = ["--config"]
        disArgs.append("doAstromQa=False")
        disArgs.append("doCompleteQa=False")
        disArgs.append("doEmptySectorQa=False")
        disArgs.append("doPhotCompareQa=False")
        disArgs.append("doPsfShapeQa=False")
        disArgs.append("doVignettingQa=False")
        disArgs.append("doVisitQa=False")
        disArgs.append("doPerformanceQa=False")
        #disArgs.append("doZptFitQa=False")
        return disArgs

    def validateFiles(self, visit, filt, raft, ccd, checkFpa = True):
        if not os.path.isdir(self.wwwRerun):
            self.fail()

        testDir = os.path.join(self.wwwRerun, "test_%s-%s_testPipeQaTask.ZeropointFitQaTask" % (visit, filt))

        if not os.path.isdir(testDir):
            self.fail()
        
        perCcd = "zeropointFit-R:%s_S:%s--%s%s.png" % (raft, ccd, re.sub(",", "", raft), re.sub(",", "", ccd))
        if not os.path.isfile(os.path.join(testDir, perCcd)):
            self.fail()

        if checkFpa:
            if not os.path.isfile(os.path.join(testDir, "zeropoint.png")):
                self.fail()
        
    def tearDown(self):
        del self.qaTask

        try:
            shutil.rmtree(self.wwwRerun) # just in case
        except:
            pass
    
    def testBasic(self):
        os.mkdir(self.wwwRerun)

        args = ["-e", "-v", self.testVisit1, "-r", self.testRaft, "-c", self.testCcd, self.testDatabase]
        for disArg in self.disableTasks():
            args.append(disArg)

        self.qaTask.parseAndRun(args)
        self.validateFiles(self.testVisit1, self.testFilt1, self.testRaft, self.testCcd)

        shutil.rmtree(self.wwwRerun)

    def testAll(self):
        os.mkdir(self.wwwRerun)

        args = ["-e", "-v", self.testVisit1, "-r", self.testRaft, "-c", self.testCcd, self.testDatabase]
        self.qaTask.parseAndRun(args)
        self.validateFiles(self.testVisit1, self.testFilt1, self.testRaft, self.testCcd)

        shutil.rmtree(self.wwwRerun)


    def testRegexp(self):
        os.mkdir(self.wwwRerun)
        
        args = ["-e", "-v", self.testVisit1, "-r", self.testRaft, "-c", re.sub(",1", ".*", self.testCcd), self.testDatabase]
        for disArg in self.disableTasks():
            args.append(disArg)

        self.qaTask.parseAndRun(args)
        for ccd in (0, 1, 2):
            self.validateFiles(self.testVisit1, self.testFilt1, self.testRaft, re.sub(",1", ",%d" % (ccd), self.testCcd))

        shutil.rmtree(self.wwwRerun)

    def testMemoryOpt(self):
        #  -b ccd  will run 1 ccd at a time and free memory after each
        #  -k is needed to write the values to cache so they can be retrieved
        #     when the final summary figure is made (i.e., since we freed them with -b ccd)
        #  -d delays making the summary figure until all ccds have run
        #     otherwise we'd waste cycles making the figure and overwriting it repeatedly
        #  -f will fork the process before plotting. When all the data are loaded for
        #     the final summary figures, the memory footprint grows and the os can't get it back from the PID
        #     so the plot() method is run as a separate PID that dies and returns control to pipeQa.py 
        os.mkdir(self.wwwRerun)

        args = ["-b", "ccd", "-k", "-d", "-f", "-e", "-v", self.testVisit1, "-r", self.testRaft, "-c", re.sub(",1", ".*", self.testCcd), self.testDatabase]
        for disArg in self.disableTasks():
            args.append(disArg)

        self.qaTask.parseAndRun(args)
        for ccd in (0, 1, 2):
            self.validateFiles(self.testVisit1, self.testFilt1, self.testRaft, re.sub(",1", ",%d" % (ccd), self.testCcd))

        shutil.rmtree(self.wwwRerun)

    def testMulti(self):
        # -g 5:n says 'group all visits matching '888.*' in groups of 5, and run the n'th one
        #        so the first example runs the first 5 visits, the second one runs the next 5 visits
        # --noWwwCache is essential if multiple pipeQas will be writing to the same place.
        #              more than 2 or 3, and there will be sqlite conflicts 
        os.mkdir(self.wwwRerun)

        args = ["--noWwwCache", "-f", "-d", "-b", "ccd", "-k", "-e", "-g", "1:0", "-v", "899.*", "-r", self.testRaft, "-c", self.testCcd, self.testDatabase]
        for disArg in self.disableTasks():
            args.append(disArg)
        self.qaTask.parseAndRun(args)
        self.validateFiles(self.testVisit1, self.testFilt1, self.testRaft, self.testCcd)

        args = ["--noWwwCache", "-f", "-d", "-b", "ccd", "-k", "-e", "-g", "1:1", "-v", "899.*", "-r", self.testRaft, "-c", self.testCcd, self.testDatabase]
        for disArg in self.disableTasks():
            args.append(disArg)
        self.qaTask.parseAndRun(args)
        self.validateFiles(self.testVisit2, self.testFilt2, self.testRaft, self.testCcd)

        shutil.rmtree(self.wwwRerun)
        

    def testExcept(self):
        os.mkdir(self.wwwRerun)
        args = ["-e", "-b", "invalid", "-v", self.testVisit1, "-r", self.testRaft, "-c", self.testCcd, self.testDatabase]
        try:
            self.qaTask.parseAndRun(args)
        except:
            pass
        else:
            self.fail()
        shutil.rmtree(self.wwwRerun)
#####
        
def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(PipeQaTestCases)
    suites += unittest.makeSuite(tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(doExit=False):
    """Run the tests"""
    tests.run(suite(), doExit)

if __name__ == "__main__":
    run(True)
