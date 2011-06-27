import os
import lsst.testing.pipeQA.figures as qaFig
import numpy
import cPickle as pickle

import lsst.testing.pipeQA.TestCode as testCode

class QaAnalysis(object):
    """Baseclass for analysis classes."""

    def __init__(self, testLabel=None, useCache=False):
        """
        @param testLabel   A name for this kind of analysis test.
        """
        
        self.testSets = {}
        self.testLabel = testLabel

        # if we're not going to use the cached values
        # we'll have to clean the output directory on our first call
        self.useCache = useCache
        self.clean    = not useCache
        
    
    def getTestSet(self, data, dataId, label=None):
        """Get a TestSet object in the correct group.

        @param data    a QaData object
        @param dataId  a dataId dictionary
        @param label   a label for particular TestSet
        """
        
        group = dataId['visit']
        filter = data.getFilterBySensor(dataId)
        # all sensors have the same filter, so just grab one
        key = filter.keys()[0]
        filterName = filter[key].getName()

        if not label is None:
            label = self.__class__.__name__ + "."+ label
        else:
            label = self.__class__.__name__

        tsId = group + "-" + filterName
        if not self.testSets.has_key(tsId):
            self.testSets[tsId] = testCode.TestSet(label, group=tsId, clean=self.clean)
            self.testSets[tsId].addMetadata('dataset', data.getDataName())
            self.testSets[tsId].addMetadata('visit-filter', tsId)

        return self.testSets[tsId]


    def __str__(self):
        testLabel = ""
        if not self.testLabel is None:
            testLabel = "."+self.testLabel
        return self.__class__.__name__ + testLabel
    

    ##########################################
    # pure virtual methods
    #########################################
    def free(self):
        """Method to free attributes to minimize memory consumption."""
        pass
    
    def test(self):
        """Method to perform tests. """
        return []

    def plot(self):
        """Method to make figures."""
        return []

