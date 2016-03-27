import sys
import os
import re
import numpy
import time

import lsst.meas.algorithms as measAlg
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .QaAnalysisTask import QaAnalysisTask
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil
import QaPlotUtils as qaPlotUtil

import lsst.testing.pipeQA.source as pqaSource

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.font_manager as fm


class PhotCompareQaConfig(pexConfig.Config):

    cameras = pexConfig.ListField(dtype = str, doc = "Cameras to run PhotCompareQaTask",
                                  default = ("lsstSim", "hscSim", "suprimecam", "cfht", "sdss", "coadd"))
    magCut = pexConfig.Field(dtype = float, doc = "Faintest magnitude for establishing photometric RMS",
                             default = 20.0)
    deltaMin = pexConfig.Field(dtype = float, doc = "Minimum allowed delta", default = -0.02)
    deltaMax = pexConfig.Field(dtype = float, doc = "Maximum allowed delta", default = 0.02)
    rmsMax = pexConfig.Field(dtype = float, doc = "Maximum allowed photometric RMS on bright end",
                             default = 0.02)
    derrMax = pexConfig.Field(dtype = float, doc = "Maximum allowed error bar underestimate on bright end",
                              default = 0.01)
    slopeMinSigma = pexConfig.Field(dtype = float,
                                    doc = "Minimum (positive valued) std.devs. of slope below slope=0", default = 3.5)
    slopeMaxSigma = pexConfig.Field(dtype = float,
                                    doc = "Maximum std.dev. of slope above slope=0", default = 3.5)

    compareTypes = pexConfig.ListField(dtype = str,
                                       doc = "Photometric Error: qaAnalysis.PhotCompareQaAnalysis",
                                       default = ("psf cat", "psf ap", "psf mod",
                                                  "ap cat", "psf inst", "inst cat", "mod cat", "mod inst"))

# allowed = {
#    "psf cat"  : "Compare Psf magnitudes to catalog magnitudes",
#    "psf ap"   : "Compare Psf and aperture magnitudes",
#    "psf mod"  : "Compare Psf and model magnitudes",
#    "ap cat"   : "Compare Psf and model magnitudes",
#    "psf inst" : "Compare PSF and instrument magnitudes",
#    "inst cat" : "Compare Inst (Gaussian) and catalog magnitudes",
#    "mod cat"  : "Compare model and catalog magnitudes",
#    "mod inst" : "Separate stars/gxys for model and inst (Gaussian) magnitudes"
# }
#

    starGalaxyToggle = pexConfig.ListField(dtype = str, doc = "Make separate figures for stars and galaxies.",
                                           default = ("mod cat", "inst cat", "ap cat", "psf cat"))

# allowed = {
#    "psf cat"  : "Separate stars/gxys for Psf magnitudes to catalog magnitudes",
#    "psf ap"   : "Separate stars/gxys for Psf and aperture magnitudes",
#    "psf mod"  : "Separate stars/gxys for Psf and model magnitudes",
#    "ap cat"   : "Separate stars/gxys for Psf and model magnitudes",
#    "psf inst" : "Separate stars/gxys for PSF and instrument magnitudes",
#    "inst cat" : "Separate stars/gxys for Inst (Gaussian) and catalog magnitudes",
#    "mod cat"  : "Separate stars/gxys for model and catalog magnitudes",
#    "mod inst" : "Separate stars/gxys for model and inst (Gaussian) magnitudes"
# }
#


class PhotCompareQaTask(QaAnalysisTask):

    ConfigClass = PhotCompareQaConfig
    _DefaultName = "photCompareQa"

    def __init__(self, magType1, magType2, starGalaxyToggle, **kwargs):
        testLabel = magType1+"-"+magType2
        QaAnalysisTask.__init__(self, testLabel, **kwargs)

        self.magCut = self.config.magCut
        self.deltaLimits = [self.config.deltaMin, self.config.deltaMax]
        self.rmsLimits = [0.0, self.config.rmsMax]
        self.derrLimits = [0.0, self.config.derrMax]
        self.slopeLimits = [-self.config.slopeMinSigma, self.config.slopeMaxSigma]
        self.starGalaxyToggle = starGalaxyToggle  # not from config!

        self.sCatDummy = pqaSource.Catalog()
        self.srefCatDummy = pqaSource.RefCatalog()

        def magType(mType):
            if re.search("(psf|PSF)", mType):
                return "psf"
            elif re.search("^ap", mType):
                return "ap"
            elif re.search("^mod", mType):
                return "mod"
            elif re.search("^cat", mType):
                return "cat"
            elif re.search("^inst", mType):
                return "inst"

        self.magType1 = magType(magType1)
        self.magType2 = magType(magType2)

        self.description = """
         For each CCD, the difference between magnitude1 and magnitude2 are
         plotted as a function of magnitude1, and also versus x,y pixel coordinate 
         (color and point size represent delta-mag and mag1).  We make comparisons between
         aperture magnitudes (ap), reference catalog magnitudes (cat), psf 
         magnitudes (psf), multifit model magnitudes (mod), and
         gaussian model magnitudes (inst).  The magnitudes used for a given test are shown
         in the test name, eg.: psf-cat.  The width of the bright end of
         this distribution (stars are plotted in red) reflects the systematic floor
         in these measurement comparisons.  For psf magnitudes, this is
         typically 1-2% for PT1.2 measurements.  The summary figure showing all
         data from all CCDs in the FPA shows stars (red) and galaxies (blue), and 
         on the bottom the average photometric error bar of the stars and the empirical
         RMS in bins of 0.5 mags.  The bottom panel also shows the additional error
         that must be added in quadrature to the photometric error bars to match the empirical
         scatter (derr).  The FPA figures on the right show the the mean magnitude 
         offset, slope in this offset as a function of magnitude, and width (standard deviation) 
         of this distribution for the bright stars (mag < 20).  The Derr plot shows the value
         of derr for all CCDs.  The per-CCD figures (everything, galaxies, stars) show the scatter
         on a per-CCD basis, as well as the distribution of the residuals across the focal plane.
         The per-CCD derr figure compares the error bars v. magnitude with the empirical RMS.
        """

    def _getFlux(self, mType, s, sref):

        # if the source isn't valid, return NaN
        if not hasattr(s, 'getId') or not hasattr(sref, 'getId'):
            return numpy.NaN

        if mType == "psf":
            return s.getD(self.sCatDummy.PsfFluxKey)
        elif mType == "ap":
            return s.getD(self.sCatDummy.ApFluxKey)
        elif mType == "mod":
            return s.getD(self.sCatDummy.ModelFluxKey)
        elif mType == "cat":
            return sref.getD(self.srefCatDummy.PsfFluxKey)
        elif mType == "inst":
            return s.getD(self.sCatDummy.InstFluxKey)

    def _getFluxErr(self, mType, s, sref):

        # if the source isn't valid, return NaN
        if not hasattr(s, 'getId') or not hasattr(sref, 'getId'):
            return numpy.NaN

        if mType == "psf":
            return s.getD(self.sCatDummy.PsfFluxErrKey)
        elif mType == "ap":
            return s.getD(self.sCatDummy.ApFluxErrKey)
        elif mType == "mod":
            return s.getD(self.sCatDummy.ModelFluxErrKey)
        elif mType == "cat":
            return 0.0
        elif mType == "inst":
            return s.getD(self.sCatDummy.InstFluxErrKey)

    def free(self):
        del self.x
        del self.y
        del self.mag
        del self.diff
        del self.filter
        del self.detector
        if self.matchListDictSrc is not None:
            del self.matchListDictSrc
        if self.ssDict is not None:
            del self.ssDict
        del self.means
        del self.medians
        del self.stds
        del self.trend
        del self.star

    def test(self, data, dataId):

        # get data
        self.detector = data.getDetectorBySensor(dataId)
        self.filter = data.getFilterBySensor(dataId)

        self.derr = raftCcdData.RaftCcdVector(self.detector)
        self.diff = raftCcdData.RaftCcdVector(self.detector)
        self.mag = raftCcdData.RaftCcdVector(self.detector)
        self.x = raftCcdData.RaftCcdVector(self.detector)
        self.y = raftCcdData.RaftCcdVector(self.detector)
        self.star = raftCcdData.RaftCcdVector(self.detector)

        filter = None

        self.matchListDictSrc = None
        self.ssDict = None

        # if we're asked to compare catalog fluxes ... we need a matchlist
        if self.magType1 == "cat" or self.magType2 == "cat":
            self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
            for key in self.matchListDictSrc.keys():
                raft = self.detector[key].getParent().getId().getName()
                ccd = self.detector[key].getId().getName()
                filter = self.filter[key].getName()

                matchList = self.matchListDictSrc[key]['matched']
                # qaAnaUtil.isStar(matchList)

                for m in matchList:
                    sref, s, dist = m

                    f1 = self._getFlux(self.magType1, s, sref)
                    f2 = self._getFlux(self.magType2, s, sref)
                    df1 = self._getFluxErr(self.magType1, s, sref)
                    df2 = self._getFluxErr(self.magType2, s, sref)

                    #badFlags = pqaSource.INTERP_CENTER | pqaSource.SATUR_CENTER | pqaSource.EDGE
                    intcen = s.getD(self.sCatDummy.FlagPixInterpCenKey)
                    satcen = s.getD(self.sCatDummy.FlagPixSaturCenKey)
                    edge = s.getD(self.sCatDummy.FlagPixEdgeKey)

                    if data.cameraInfo.name == 'coadd':
                        flagit = (satcen or edge)  # coadds have excessive area covered by InterpCen flags
                    else:
                        flagit = (intcen or satcen or edge)

                    if (f1 > 0.0 and f2 > 0.0 and not flagit):
                        m1 = -2.5*numpy.log10(f1)
                        m2 = -2.5*numpy.log10(f2)
                        dm1 = 2.5 / numpy.log(10.0) * df1 / f1
                        dm2 = 2.5 / numpy.log(10.0) * df2 / f2

                        star = 0 if s.getD(self.sCatDummy.ExtendednessKey) else 1

                        if numpy.isfinite(m1) and numpy.isfinite(m2):
                            self.derr.append(raft, ccd, numpy.sqrt(dm1**2 + dm2**2))
                            self.diff.append(raft, ccd, m1 - m2)
                            self.mag.append(raft, ccd, m1)
                            self.x.append(raft, ccd, s.getD(self.sCatDummy.XAstromKey))
                            self.y.append(raft, ccd, s.getD(self.sCatDummy.YAstromKey))
                            self.star.append(raft, ccd, star)

        # if we're not asked for catalog fluxes, we can just use a sourceSet
        else:
            self.ssDict = data.getSourceSetBySensor(dataId)
            for key, ss in self.ssDict.items():
                raft = self.detector[key].getParent().getId().getName()
                ccd = self.detector[key].getId().getName()

                filter = self.filter[key].getName()

                # qaAnaUtil.isStar(ss)  # sets the 'STAR' flag
                for s in ss:
                    f1 = self._getFlux(self.magType1, s, s)
                    f2 = self._getFlux(self.magType2, s, s)
                    df1 = self._getFluxErr(self.magType1, s, s)
                    df2 = self._getFluxErr(self.magType2, s, s)
                    intcen = s.getD(self.sCatDummy.FlagPixInterpCenKey)
                    satcen = s.getD(self.sCatDummy.FlagPixSaturCenKey)
                    edge = s.getD(self.sCatDummy.FlagPixEdgeKey)

                    if data.cameraInfo.name == 'coadd':
                        flagit = (satcen or edge)  # coadds have excessive area covered by InterpCen flags
                    else:
                        flagit = (intcen or satcen or edge)

                    if ((f1 > 0.0 and f2 > 0.0) and not flagit):

                        m1 = -2.5*numpy.log10(f1)  # self.calib[key].getMagnitude(f1)
                        m2 = -2.5*numpy.log10(f2)  # self.calib[key].getMagnitude(f2)
                        dm1 = 2.5 / numpy.log(10.0) * df1 / f1
                        dm2 = 2.5 / numpy.log(10.0) * df2 / f2

                        star = 0 if s.getD(self.sCatDummy.ExtendednessKey) else 1

                        # if star:
                        #    print "isStar: ", star
                        if numpy.isfinite(m1) and numpy.isfinite(m2):
                            self.derr.append(raft, ccd, numpy.sqrt(dm1**2 + dm2**2))
                            self.diff.append(raft, ccd, m1 - m2)
                            self.mag.append(raft, ccd, m1)
                            self.x.append(raft, ccd, s.getD(self.sCatDummy.XAstromKey))
                            self.y.append(raft, ccd, s.getD(self.sCatDummy.YAstromKey))
                            self.star.append(raft, ccd, star)

        testSet = self.getTestSet(data, dataId, label=self.magType1+"-"+self.magType2)

        testSet.addMetadata('magType1', self.magType1)
        testSet.addMetadata('magType2', self.magType2)
        testSet.addMetadata({"Description": self.description})

        self.means = raftCcdData.RaftCcdData(self.detector)
        self.medians = raftCcdData.RaftCcdData(self.detector)
        self.stds = raftCcdData.RaftCcdData(self.detector)
        self.derrs = raftCcdData.RaftCcdData(self.detector)
        self.trend = raftCcdData.RaftCcdData(self.detector, initValue=[0.0, 0.0])

        self.dmagMax = 0.4
        allMags = numpy.array([])
        allDiffs = numpy.array([])

        for raft, ccd in self.mag.raftCcdKeys():
            dmag0 = self.diff.get(raft, ccd)
            mag0 = self.mag.get(raft, ccd)
            derr0 = self.derr.get(raft, ccd)
            star = self.star.get(raft, ccd)

            wGxy = numpy.where((mag0 > 10) & (mag0 < self.magCut) &
                               (star == 0) & (numpy.abs(dmag0) < 1.0))[0]
            w = numpy.where((mag0 > 10) & (mag0 < self.magCut) & (star > 0))[0]

            mag = mag0[w]
            dmag = dmag0[w]
            derr = derr0[w]

            allMags = numpy.append(allMags, mag)
            allDiffs = numpy.append(allDiffs, dmag)

            # already using NaN for 'no-data' for this ccd
            #  (because we can't test for 'None' in a numpy masked_array)
            # unfortunately, these failures will have to do
            mean = 99.0
            median = 99.0
            std = 99.0
            derrmed = 99.0
            n = 0
            lineFit = [[99.0, 0.0, 0.0, 0.0]]*3
            lineCoeffs = [[99.0, 0.0]]*3
            if len(dmag) > 0:
                stat = afwMath.makeStatistics(dmag, afwMath.NPOINT | afwMath.MEANCLIP |
                                              afwMath.STDEVCLIP | afwMath.MEDIAN)
                mean = stat.getValue(afwMath.MEANCLIP)
                median = stat.getValue(afwMath.MEDIAN)
                std = stat.getValue(afwMath.STDEVCLIP)
                n = stat.getValue(afwMath.NPOINT)

                derrmed = afwMath.makeStatistics(derr, afwMath.MEDIAN).getValue(afwMath.MEDIAN)

                # get trendlines for stars/galaxies
                # for alldata, use trendline for stars
                if len(dmag) > 1:
                    lineFit[0] = qaAnaUtil.robustPolyFit(mag, dmag, 1)
                    lineCoeffs[0] = lineFit[0][0], lineFit[0][2]
                if len(wGxy) > 1:
                    lineFit[1] = qaAnaUtil.robustPolyFit(mag0[wGxy], dmag0[wGxy], 1)
                    lineCoeffs[1] = lineFit[1][0], lineFit[1][2]
                lineFit[2] = lineFit[0]
                lineCoeffs[2] = lineCoeffs[0]

            tag = self.magType1+"_vs_"+self.magType2
            dtag = self.magType1+"-"+self.magType2
            self.means.set(raft, ccd, mean)
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "mean "+tag  # +" " + areaLabel
            comment = "mean "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmag), n)
            testSet.addTest(testCode.Test(label, mean, self.deltaLimits, comment, areaLabel=areaLabel))

            self.medians.set(raft, ccd, median)
            label = "median "+tag  # +" "+areaLabel
            comment = "median "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmag), n)
            testSet.addTest(testCode.Test(label, median, self.deltaLimits, comment, areaLabel=areaLabel))

            self.stds.set(raft, ccd, std)
            label = "stdev "+tag  # +" " + areaLabel
            comment = "stdev of "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmag), n)
            testSet.addTest(testCode.Test(label, std, self.rmsLimits, comment, areaLabel=areaLabel))

            self.derrs.set(raft, ccd, derrmed)
            label = "derr "+tag
            comment = "add phot err in quad for "+dtag + \
                " (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmag), n)
            testSet.addTest(testCode.Test(label, derrmed, self.derrLimits, comment, areaLabel=areaLabel))

            self.trend.set(raft, ccd, lineFit)
            label = "slope "+tag  # +" " + areaLabel
            slopeLimits = self.slopeLimits[0]*lineFit[0][1], self.slopeLimits[1]*lineFit[0][1]
            comment = "slope of "+dtag+" (mag lt %.1f, nstar/clip=%d/%d) limits=(%.1f,%.1f)sigma" % \
                      (self.magCut, len(dmag), n, self.slopeLimits[0], self.slopeLimits[1])
            testSet.addTest(testCode.Test(label, lineCoeffs[0][0], slopeLimits, comment,
                                          areaLabel=areaLabel))

        # do a test of all CCDs for the slope ... suffering small number problems
        #  on indiv ccds and could miss a problem

        lineFit = [99.0, 0.0, 0.0, 0.0]
        lineCoeffs = [99.0, 0.0]
        if len(allDiffs) > 1:
            lineFit = qaAnaUtil.robustPolyFit(allMags, allDiffs, 1)
            lineCoeffs = lineFit[0], lineFit[2]
        label = "slope"
        slopeLimits = self.slopeLimits[0]*lineFit[1], self.slopeLimits[1]*lineFit[1]
        comment = "slope for all ccds (mag lt %.1f, nstar=%d) limits=(%.1f,%.1f sigma)" % (
            self.magCut, len(allDiffs), self.slopeLimits[0], self.slopeLimits[1])
        testSet.addTest(testCode.Test(label, lineCoeffs[0], slopeLimits, comment, areaLabel="all"))

    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId, label=self.magType1+"-"+self.magType2)
        testSet.setUseCache(self.useCache)

        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        xlim = [14.0, 25.0]
        ylimStep = 0.4
        ylim = [-ylimStep, ylimStep]
        aspRatio = (xlim[1]-xlim[0])/(ylim[1]-ylim[0])

        tag1 = "m$_{\mathrm{"+self.magType1.upper()+"}}$"
        tag = "m$_{\mathrm{"+self.magType1.upper()+"}}$ - m$_{\mathrm{"+self.magType2.upper()+"}}$"
        dtag = self.magType1+"-"+self.magType2
        wtag = self.magType1+"minus"+self.magType2

        if self.showFpa:
            # fpa figure
            meanFilebase = "mean" + wtag
            stdFilebase = "std"+wtag
            derrFilebase = "derr"+wtag
            slopeFilebase = "slope"+wtag
            meanData, meanMap = testSet.unpickle(meanFilebase, default=[None, None])
            stdData, stdMap = testSet.unpickle(stdFilebase, default=[None, None])
            derrData, derrMap = testSet.unpickle(derrFilebase, default=[None, None])
            slopeData, slopeMap = testSet.unpickle(slopeFilebase, default=[None, None])

            meanFig = qaFig.FpaQaFigure(data.cameraInfo, data=meanData, map=meanMap)
            stdFig = qaFig.FpaQaFigure(data.cameraInfo, data=stdData, map=stdMap)
            derrFig = qaFig.FpaQaFigure(data.cameraInfo, data=derrData, map=derrMap)
            slopeFig = qaFig.VectorFpaQaFigure(data.cameraInfo, data=slopeData, map=slopeMap)

            for raft, ccd in self.means.raftCcdKeys():

                meanFig.data[raft][ccd] = self.means.get(raft, ccd)
                stdFig.data[raft][ccd] = self.stds.get(raft, ccd)
                derrFig.data[raft][ccd] = self.derrs.get(raft, ccd)
                slope = self.trend.get(raft, ccd)[0]

                if slope is not None and not slope[1] == 0:
                    # aspRatio will make the vector have the same angle as the line in the figure
                    slopeSigma = slope[0]/slope[1]
                    slopeFig.data[raft][ccd] = [numpy.arctan2(aspRatio*slope[0], 1.0), None, slopeSigma]
                else:
                    slopeSigma = None
                    slopeFig.data[raft][ccd] = [None, None, None]

                if self.means.get(raft, ccd) is not None:
                    meanFig.map[raft][ccd] = "mean=%.4f" % (self.means.get(raft, ccd))
                    stdFig.map[raft][ccd] = "std=%.4f" % (self.stds.get(raft, ccd))
                    derrFig.map[raft][ccd] = "derr=%.4f" % (self.derrs.get(raft, ccd))
                    fmt0, fmt1, fmtS = "%.4f", "%.4f", "%.1f"
                    if slope[0] is None:
                        fmt0 = "%s"
                    if slope[1] is None:
                        fmt1 = "%s"
                    if slopeSigma is None:
                        fmtS = "%s"
                    fmt = "slope="+fmt0+"+/-"+fmt1+"("+fmtS+"sig)"
                    slopeFig.map[raft][ccd] = fmt % (slope[0], slope[1], slopeSigma)

            blue, red = '#0000ff', '#ff0000'

            testSet.pickle(meanFilebase, [meanFig.data, meanFig.map])
            testSet.pickle(stdFilebase, [stdFig.data, stdFig.map])
            testSet.pickle(derrFilebase, [derrFig.data, derrFig.map])
            testSet.pickle(slopeFilebase, [slopeFig.data, slopeFig.map])

            if not self.delaySummary or isFinalDataId:
                self.log.log(self.log.INFO, "plotting FPAs")
                meanFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[-0.03, 0.03],
                                   title="Mean "+tag, cmapOver=red, cmapUnder=blue, failLimits=self.deltaLimits)
                testSet.addFigure(meanFig, "f01"+meanFilebase+".png",
                                  "mean "+dtag+" mag   (brighter than %.1f)" % (self.magCut), navMap=True)
                del meanFig

                stdFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 0.03],
                                  title="Stdev "+tag, cmapOver=red, failLimits=self.rmsLimits)
                testSet.addFigure(stdFig, "f02"+stdFilebase+".png",
                                  "stdev "+dtag+" mag  (brighter than %.1f)" % (self.magCut), navMap=True)
                del stdFig

                derrFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 0.01],
                                   title="Derr "+tag, cmapOver=red, failLimits=self.derrLimits)
                testSet.addFigure(derrFig, "f03"+derrFilebase+".png",
                                  "derr "+dtag+" mag (brighter than %.1f)" % (self.magCut), navMap=True)
                del derrFig

                cScale = 2.0
                slopeFig.makeFigure(cmap="RdBu_r",
                                    vlimits=[cScale*self.slopeLimits[0], cScale*self.slopeLimits[1]],
                                    title="Slope "+tag, failLimits=self.slopeLimits)
                testSet.addFigure(slopeFig, "f04"+slopeFilebase+".png",
                                  "slope "+dtag+" mag (brighter than %.1f)" % (self.magCut), navMap=True)
                del slopeFig
            else:
                del meanFig
                del stdFig
                del derrFig
                del slopeFig

        #############################################
        #
        #xmin, xmax = self.mag.summarize('min', default=0.0), self.mag.summarize('max', default=25.0)
        #ymin, ymax = self.diff.summarize('min', default=-1.0), self.diff.summarize('max', default=1.0)
        #xrang = xmax-xmin
        #xmin, xmax = int(xmin-0.05*xrang), int(xmax+0.05*xrang)+1
        #yrang = ymax-ymin
        #ymin, ymax = ymin-0.05*yrang, ymax+0.05*yrang
        xlim2 = xlim  # [xmin, xmax]
        ylim2 = [-2.0, 2.0]  # [ymin, ymax]

        figsize = (6.5, 3.75)

        figbase = "diff_" + dtag

        shelfData = {}

        for raft, ccd in self.mag.raftCcdKeys():
            mag0 = self.mag.get(raft, ccd)
            diff0 = self.diff.get(raft, ccd)
            star0 = self.star.get(raft, ccd)
            derr0 = self.derr.get(raft, ccd)

            self.log.log(self.log.INFO, "plotting %s" % (ccd))

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            statBlurb = "Points used for statistics/trendline shown in red."

            dataDict = {
                'mag0': mag0,
                'diff0': diff0,
                'star0': star0,
                'derr0': derr0,
                'areaLabel': areaLabel,
                'raft': raft,
                'ccd': ccd,
                'figsize': figsize,
                'xlim': xlim,
                'ylim': ylim,
                'xlim2': xlim2,
                'ylim2': ylim2,
                'ylimStep': ylimStep,
                'tag1': tag1,
                'tag': tag,

                'x': self.y.get(raft, ccd),
                'y': self.x.get(raft, ccd),
                'trend': self.trend.get(raft, ccd),
                'magCut': self.magCut,
                'summary': False,
            }

            masterToggle = None
            if self.starGalaxyToggle:

                masterToggle = '0_stars'

                import PhotCompareQaAnalysisPlot as plotModule
                label = areaLabel
                pngFile = figbase+".png"

                ##############################################
                caption = dtag + " vs. " + self.magType1 + "(stars). "+statBlurb
                toggle = '0_stars'
                if self.lazyPlot.lower() in ['sensor', 'all']:
                    testSet.addLazyFigure(dataDict, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="stars standard",
                                          toggle=toggle,
                                          masterToggle=masterToggle)
                else:
                    testSet.cacheLazyData(dataDict, pngFile, areaLabel=label, toggle=toggle,
                                          masterToggle=masterToggle)
                    dataDict['mode'] = 'stars'
                    dataDict['figType'] = 'standard'
                    fig = plotModule.plot(dataDict)
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label, toggle=toggle)
                    del fig

                ##############################################
                caption = dtag + " vs. " + self.magType1 + "(galaxies). "+statBlurb
                toggle = '1_galaxies'
                if self.lazyPlot.lower() in ['sensor', 'all']:
                    testSet.addLazyFigure(dataDict, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="galaxies standard",
                                          toggle=toggle,
                                          masterToggle=masterToggle)
                else:
                    testSet.cacheLazyData(dataDict, pngFile, areaLabel=label, toggle=toggle,
                                          masterToggle=masterToggle)
                    dataDict['mode'] = 'galaxies'
                    dataDict['figType'] = 'standard'
                    fig = plotModule.plot(dataDict)
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label, toggle=toggle)
                    del fig

                ##############################################
                caption = dtag + " vs. " + self.magType1 + "(everything). "+statBlurb
                toggle = '2_everything'
                if self.lazyPlot.lower() in ['sensor', 'all']:
                    testSet.addLazyFigure(dataDict, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="all standard",
                                          toggle=toggle,
                                          masterToggle=masterToggle)
                else:
                    testSet.cacheLazyData(dataDict, pngFile, areaLabel=label, toggle=toggle,
                                          masterToggle=masterToggle)
                    dataDict['mode'] = 'all'
                    dataDict['figType'] = 'standard'
                    fig = plotModule.plot(dataDict)
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label, toggle=toggle)
                    del fig

                ##############################################
                caption = dtag + " vs. " + self.magType1 + "(star derr). "+statBlurb
                toggle = '3_derr'
                if self.lazyPlot.lower() in ['sensor', 'all']:
                    testSet.addLazyFigure(dataDict, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="stars derr",
                                          toggle=toggle,
                                          masterToggle=masterToggle)
                else:
                    testSet.cacheLazyData(dataDict, pngFile, areaLabel=label, toggle=toggle,
                                          masterToggle=masterToggle)
                    dataDict['mode'] = 'stars'
                    dataDict['figType'] = 'derr'
                    fig = plotModule.plot(dataDict)
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label, toggle=toggle)
                    del fig

            else:

                import PhotCompareQaAnalysisPlot as plotModule
                label = areaLabel
                caption = dtag + " vs. " + self.magType1 + ". "+statBlurb
                pngFile = figbase+".png"

                if self.lazyPlot.lower() in ['sensor', 'all']:
                    testSet.addLazyFigure(dataDict, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="all standard")
                else:
                    testSet.cacheLazyData(dataDict, pngFile, areaLabel=label)
                    dataDict['mode'] = 'all'
                    dataDict['figType'] = 'standard'
                    fig = plotModule.plot(dataDict)
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                    del fig

        if not self.delaySummary or isFinalDataId:
            self.log.log(self.log.INFO, "plotting Summary figure")

            import PhotCompareQaAnalysisPlot as plotModule
            label = 'all'
            caption = dtag+" vs. "+self.magType1
            pngFile = figbase+".png"

            if self.lazyPlot in ['all']:
                testSet.addLazyFigure({}, pngFile, caption,
                                      plotModule, areaLabel=label, plotargs="all summary",
                                      masterToggle=masterToggle)
            else:
                dataDict, isSummary = qaPlotUtil.unshelveGlob(figbase+"-all.png", testSet=testSet)
                dataDict['mode'] = 'all'
                dataDict['figType'] = 'summary'
                fig = plotModule.plot(dataDict)
                testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                del fig




