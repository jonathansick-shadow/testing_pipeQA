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


class MultiEpochAstrometryQa(qaAna.QaAnalysis):

    def __init__(self):
        qaAna.QaAnalysis.__init__(self)

    def free(self):
        pass
        
    def test(self, data, dataId):
        
        # get data
        self.clusterList   = data.getSourceClusters(dataId)

        raStd = []
        decStd = []

        for ss in self.clusterList:
            ra, dec = [], []
            for s in ss:
                ra.append(s.getRa())
                dec.append(s.getDec())
            ra, dec = numpy.array(ra), numpy.array(dec)

            raStd.append(206265.0*numpy.cos(dec.mean())*ra.std())
            decStd.append(206265.0*dec.std())

        self.raStd = numpy.array(raStd)
        self.decStd = numpy.array(decStd)
        
        self.limits = [0, 0.2]

        testSet = self.getTestSet(data, dataId)

        visits = data.getVisits(dataId)
        testSet.addMetadata("visits", "\n".join(visits))
        
        # add tests for acceptible numpy of empty sectors
        testRa = testCode.Test("median RA stdev", numpy.median(self.raStd), self.limits,
                               "median RA standard deviation [arcsec]", areaLabel='all')
        testDec = testCode.Test("median Dec stdev", numpy.median(self.decStd), self.limits,
                                "median Declination standard deviation [arcsec]", areaLabel='all')
        for test in [testRa, testDec]:
            testSet.addTest(test)
            

    def plot(self, data, dataId):

        testSet = self.getTestSet(data, dataId)

        fig = qaFig.QaFigure()
        ax = fig.fig.add_subplot(111)

        ax.scatter(self.raStd, self.decStd, s=1.0, color='k')
        ax.set_xlabel("$\Delta$ R.A.", size='x-small')
        ax.set_ylabel("$\Delta$ Decl.", size='x-small')
        ax.set_xlim([0, 0.5])
        ax.set_ylim([0, 0.5])
        for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_size('xx-small')


        testSet.addFigure(fig, "multiEpochAstrom.png", "Std.Devs. of RA/Dec for matched objects.",
                          areaLabel='all')
