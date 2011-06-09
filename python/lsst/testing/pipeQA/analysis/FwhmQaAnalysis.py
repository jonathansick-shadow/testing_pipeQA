import sys, os, re
import lsst.meas.algorithms        as measAlg
import lsst.testing.pipeQA.figures as qaFig
import numpy

import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.font_manager as fm
from matplotlib.collections import LineCollection


class FwhmQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
        qaAna.QaAnalysis.__init__(self)

    def free(self):

        # please free all large data structures here.
        del self.calexpDict
        
    def test(self, data, dataId):
        
        # get data
        self.calexpDict        = data.getCalexpBySensor(dataId)

        for k,v in self.calexpDict.items():
            print k, v['fwhm']


    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)


    
