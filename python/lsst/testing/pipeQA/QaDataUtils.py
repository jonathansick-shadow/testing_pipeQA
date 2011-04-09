import os, re

def findDataInTestbed(label, raiseOnFailure=True):
    """Scan TESTBED_PATH directories to find a testbed dataset with a given name"""
    
    # If label==None -> use TESTBOT_DIR (_DIR, not _PATH ... only one dataset can be run)
    if re.search("^testBot", label):
        testbotDir = os.getenv("TESTBOT_DIR")
        testbedDir, testdataDir  = os.path.split(testbotDir)
        return testbedDir, testbotDir

    
    # otherwise, get a specific test-data set from one of the testbed directories
    testbedPath = os.getenv("TESTBED_PATH")
    if testbedPath is not None:
        testbedDirs = testbedPath.split(":")
    else:
        raise Exception("Must specify environment variable TESTBED_PATH.")
    
    #############################
    # find the label in the testbed path
    testbedDir = None
    testdataDir = None
    for tbDir in testbedDirs:
        path = os.path.join(tbDir, label)
        if os.path.exists(path):
            testbedDir = tbDir
            testdataDir = path
                
    if raiseOnFailure and (testbedDir is None):
        msg = "Testbed %s was not found in any of the TESTBED_PATH directories:\n" % (label)
        msg += "\n".join(testbedDirs) + "\n"
        raise Exception(msg)
    
    return testbedDir, testdataDir


