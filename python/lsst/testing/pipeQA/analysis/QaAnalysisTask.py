import os
import numpy
import cPickle as pickle

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures as qaFig

class QaAnalysisConfig(pexConfig.Config):
    shapeAlgorithm = pexConfig.ChoiceField(
        dtype = str,
        doc = "Shape Algorithm to load.",
        default = "HSM_REGAUSS",
        allowed = {
            "HSM_REGAUSS": "Hirata-Seljac-Mandelbaum Regaussianization",
            "HSM_BJ" : "Hirata-Seljac-Mandelbaum Bernstein&Jarvis",
            "HSM_LINEAR" : "Hirata-Seljac-Mandelbaum Linear",
            "HSM_KSB" : "Hirata-Seljac-Mandelbaum KSB",
            "HSM_SHAPELET" : "Hirata-Seljac-Mandelbaum SHAPELET"
        }
    )

    doZptFitQa = pexConfig.Field(dtype = bool, doc = "Photometric Zeropoint: qaAnalysis.ZeropointFitQa", default = True)
    doEmptySectorQa = pexConfig.Field(dtype = bool, doc = "Empty Sectors: qaAnalysis.EmptySectorQaAnalysis", default = True)
    doAstromQa = pexConfig.Field(dtype = bool, doc = "Astrometric Error: qaAnalysis.AstrometricErrorQaAnalysis", default = True)
    doPhotCompareQa = pexConfig.Field(dtype = bool, doc = "Photometric Error: qaAnalysis.PhotCompareQaAnalysis", default = True)
    doPsfShapeQa = pexConfig.Field(dtype = bool, doc = "Psf Shape: qaAnalysis.PsfShapeQaAnalysis", default = True)
    doCompleteQa = pexConfig.Field(dtype = bool, doc = "Photometric Depth: qaAnalysis.CompletenessQa", default = True)
    doVignettingQa = pexConfig.Field(dtype = bool, doc = "Vignetting testing: qaAnalysis.VignettingQa", default = True)
    doVisitQa = pexConfig.Field(dtype = bool, doc = "Visit to visit: qaAnalysis.VisitToVisitPhotQaAnalysis and qaAnalysis.VisitToVisitAstromQaAnalysis", default = False)

class QaAnalysisTask(pipeBase.Task):
    """Baseclass for analysis classes."""
    ConfigClass  = QaAnalysisConfig
    _DefaultName = "qaAnalysis"

    def __init__(self, testLabel=None, useCache=False, wwwCache=True, delaySummary=False, *args, **kwargs):
        """
        @param testLabel   A name for this kind of analysis test.
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.testSets = {}
        self.testLabel = testLabel

        # if we're not going to use the cached values
        # we'll have to clean the output directory on our first call
        self.useCache = useCache
        self.clean    = not useCache
        self.wwwCache = wwwCache
        self.delaySummary = delaySummary
        
    
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
	filterName = '?'
	if not filter[key] is None:
	    filterName = filter[key].getName()

        if not label is None:
            label = self.__class__.__name__ + "."+ label
        else:
            label = self.__class__.__name__

        tsId = group + "-" + filterName
        if not self.testSets.has_key(tsId):
            self.testSets[tsId] = testCode.TestSet(label, group=tsId, clean=self.clean,
                                                   wwwCache=self.wwwCache)
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
    def run(self, doTest = True, doPlot = True, *args, **kwargs):
        """Base method to perform pipeQA testing or plotting"""
        if doTest:
            self.test(*args, **kwargs)
        if doPlot:
            self.plot(*args, **kwargs)

    def free(self):
        """Method to free attributes to minimize memory consumption."""
        pass
    
    def test(self):
        """Method to perform tests. """
        return []

    def plot(self):
        """Method to make figures."""
        return []

