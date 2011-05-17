import sys, os, re

from DatabaseQuery import LsstSimDbInterface, DatabaseIdentity

import QaDataUtils as qaDataUtils

from ButlerQaData  import makeButlerQaData
from DbQaData      import makeDbQaData

###################################################
# Factory for QaData
###################################################
def makeQaData(label, rerun=None, retrievalType=None, **kwargs):
    """ """
    
    if retrievalType is None:
        
        # if there's only one possibility, use that
        
        # see if there's a testbed directory called 'label'
        validButler = False
        testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label, raiseOnFailure=False)
        if (not testbedDir is None) and (not testdataDir is None):
            validButler = True

        # see if we can connect to a database with name 'label'
        # NOTE: must update if/when non-lsst databases get used
        validDb = True
        try:
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

            
    if re.search("[Bb]utler", retrievalType):
        return makeButlerQaData(label, rerun, **kwargs)

    if re.search("^([Dd][Bb]|[Dd]atabase)$", retrievalType):
        return makeDbQaData(label, rerun, **kwargs)

