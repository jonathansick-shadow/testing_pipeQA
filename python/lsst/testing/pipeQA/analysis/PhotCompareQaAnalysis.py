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

class PhotCompareQaAnalysis(qaAna.QaAnalysis):

    def __init__(self, magType1, magType2, magCut,
                 deltaMin, deltaMax, rmsMax, slopeMinSigma, slopeMaxSigma,
                 **kwargs):
        testLabel = magType1+"-"+magType2
        qaAna.QaAnalysis.__init__(self, testLabel, **kwargs)

        self.magCut = magCut
        self.deltaLimits = [deltaMin, deltaMax]
        self.rmsLimits = [0.0, rmsMax]
        self.slopeLimits = [-slopeMinSigma, slopeMaxSigma]
        
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
         plotted as a function of magnitude1.  We make comparisons between
         aperture magnitudes and reference catalog magnitudes (ap-cat), psf and
         aperture magnitudes (psf-ap), psf and reference catalog magnitudes
         (psf-cat), psf and multifit model magnitudes (psf-inst), and psf and
         gaussian model magnitudes (psf-mod).  The width of the bright end of
         this distribution (stars plotted in red) reflects the systematic floor
         in these measurement comparisons.  For psf magnitudes, this is
         typically 1-2% for PT1.2 measurements.  The summary FPA figures show
         the mean magnitude offset, slope in this offset as a function of
         magnitude, and width (standard deviation) of this distribution for the
         bright stars.
        """

    def _getFlux(self, mType, s, sref):

        # if the source isn't valid, return NaN
        if not hasattr(s, 'getId') or not hasattr(sref, 'getId'):
            return numpy.NaN
            
        
        if mType=="psf":
            return s.getPsfFlux()
        elif mType=="ap":
            return s.getApFlux()
        elif mType=="mod":
            return s.getModelFlux()
        elif mType=="cat":
            return sref.getPsfFlux()
        elif mType=="inst":
            return s.getInstFlux()

    def free(self):
        del self.x
        del self.y
        del self.mag
        del self.diff
        del self.filter
        del self.detector
        if not self.matchListDictSrc is None:
            del self.matchListDictSrc
        if not self.ssDict is None:
            del self.ssDict
        del self.means
        del self.medians
        del self.stds
        del self.trend
        del self.star
        

    def test(self, data, dataId):
        
        # get data
        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        
        self.diff = raftCcdData.RaftCcdVector(self.detector)
        self.mag  = raftCcdData.RaftCcdVector(self.detector)
        self.x    = raftCcdData.RaftCcdVector(self.detector)
        self.y    = raftCcdData.RaftCcdVector(self.detector)
        self.star = raftCcdData.RaftCcdVector(self.detector)

        filter = None
        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER

        self.matchListDictSrc = None
        self.ssDict = None

        # if we're asked to compare catalog fluxes ... we need a matchlist
        if  self.magType1=="cat" or self.magType2=="cat":
            self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
            for key in self.matchListDictSrc.keys():
                raft = self.detector[key].getParent().getId().getName()
                ccd  = self.detector[key].getId().getName()
                filter = self.filter[key].getName()

                matchList = self.matchListDictSrc[key]['matched']
                #qaAnaUtil.isStar(matchList)

                for m in matchList:
                    sref, s, dist = m
                    
                    f1 = self._getFlux(self.magType1, s, sref)
                    f2 = self._getFlux(self.magType2, s, sref)
                    flags = s.getFlagForDetection()

                    if (f1 > 0.0 and f2 > 0.0  and not flags & badFlags):
                        m1 = -2.5*numpy.log10(f1)
                        m2 = -2.5*numpy.log10(f2)

                        star = flags & measAlg.Flags.STAR
                        
                        if numpy.isfinite(m1) and numpy.isfinite(m2):
                            self.diff.append(raft, ccd, m1 - m2)
                            self.mag.append(raft, ccd, m1)
                            self.x.append(raft, ccd, s.getXAstrom())
                            self.y.append(raft, ccd, s.getYAstrom())
                            self.star.append(raft, ccd, star)

        # if we're not asked for catalog fluxes, we can just use a sourceSet
        else:
            self.ssDict        = data.getSourceSetBySensor(dataId)
            for key, ss in self.ssDict.items():
                raft = self.detector[key].getParent().getId().getName()
                ccd  = self.detector[key].getId().getName()

                filter = self.filter[key].getName()

                qaAnaUtil.isStar(ss)  # sets the 'STAR' flag
                for s in ss:
                    f1 = self._getFlux(self.magType1, s, s)
                    f2 = self._getFlux(self.magType2, s, s)
                    flags = s.getFlagForDetection()
                    
                    if ((f1 > 0.0 and f2 > 0.0) and not flags & badFlags):
                        m1 = -2.5*numpy.log10(f1) #self.calib[key].getMagnitude(f1)
                        m2 = -2.5*numpy.log10(f2) #self.calib[key].getMagnitude(f2)

                        star = flags & measAlg.Flags.STAR
                        #if star:
                        #    print "isStar: ", star
                        
                        self.diff.append(raft, ccd, m1 - m2)
                        self.mag.append(raft, ccd, m1)
                        self.x.append(raft, ccd, s.getXAstrom())
                        self.y.append(raft, ccd, s.getYAstrom())
                        self.star.append(raft, ccd, star)
                    
        testSet = self.getTestSet(data, dataId, label=self.magType1+"-"+self.magType2)
        testSet.addMetadata('magType1', self.magType1)
        testSet.addMetadata('magType2', self.magType2)

        self.means = raftCcdData.RaftCcdData(self.detector)
        self.medians = raftCcdData.RaftCcdData(self.detector)
        self.stds  = raftCcdData.RaftCcdData(self.detector)
        self.trend = raftCcdData.RaftCcdData(self.detector, initValue=[0.0, 0.0])
        
        self.dmagMax = 0.4
        allMags = numpy.array([])
        allDiffs = numpy.array([])

        for raft,  ccd in self.mag.raftCcdKeys():
            dmag = self.diff.get(raft, ccd)
            mag = self.mag.get(raft, ccd)
            star = self.star.get(raft, ccd)
            w = numpy.where((mag > 10) & (mag < self.magCut) & (star > 0))[0] # & (abs(dmag) < self.dmagMax) )[0]

            mag = mag[w]
            dmag = dmag[w]

            allMags = numpy.append(allMags, mag)
            allDiffs = numpy.append(allDiffs, dmag)
             
            # already using NaN for 'no-data' for this ccd
            #  (because we can't test for 'None' in a numpy masked_array)
            # unfortunately, these failures will have to do
            mean = 99.0
            median = 99.0
            std = 99.0
            n = 0
            lineFit = [99.0, 0.0, 0.0, 0.0]
            lineCoeffs = [99.0, 0.0]
            if len(dmag) > 0:
                stat = afwMath.makeStatistics(dmag, afwMath.NPOINT | afwMath.MEANCLIP |
                                              afwMath.STDEVCLIP | afwMath.MEDIAN)
                mean = stat.getValue(afwMath.MEANCLIP)
                median = stat.getValue(afwMath.MEDIAN)
                std = stat.getValue(afwMath.STDEVCLIP)
                n = stat.getValue(afwMath.NPOINT)

                if len(dmag) > 1:
                    lineFit = qaAnaUtil.robustPolyFit(mag, dmag, 1)
                    lineCoeffs = lineFit[0], lineFit[2]
                    

            tag = self.magType1+"_vs_"+self.magType2
            dtag = self.magType1+"-"+self.magType2
            self.means.set(raft, ccd, mean)
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "mean "+tag # +" " + areaLabel
            comment = "mean "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmag),n)
            testSet.addTest( testCode.Test(label, mean, self.deltaLimits, comment, areaLabel=areaLabel))

            self.medians.set(raft, ccd, median)
            label = "median "+tag #+" "+areaLabel
            comment = "median "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmag), n)
            testSet.addTest( testCode.Test(label, median, self.deltaLimits, comment, areaLabel=areaLabel))

            self.stds.set(raft, ccd, std)
            label = "stdev "+tag #+" " + areaLabel
            comment = "stdev of "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmag), n)
            testSet.addTest( testCode.Test(label, std, self.rmsLimits, comment, areaLabel=areaLabel))

            self.trend.set(raft, ccd, lineFit)
            label = "slope "+tag #+" " + areaLabel
            slopeLimits = self.slopeLimits[0]*lineFit[1], self.slopeLimits[1]*lineFit[1]
            comment = "slope of "+dtag+" (mag lt %.1f, nstar/clip=%d/%d) limits=(%.1f,%.1f)sigma" % (self.magCut, len(dmag), n, self.slopeLimits[0], self.slopeLimits[1])
            testSet.addTest( testCode.Test(label, lineCoeffs[0], slopeLimits, comment, areaLabel=areaLabel))


        # do a test of all CCDs for the slope ... suffering small number problems
        #  on indiv ccds and could miss a problem
        if len(allDiffs) > 1:
            lineFit = qaAnaUtil.robustPolyFit(allMags, allDiffs, 1)
            lineCoeffs = lineFit[0], lineFit[2]
        label = "slope"
        slopeLimits = self.slopeLimits[0]*lineFit[1], self.slopeLimits[1]*lineFit[1]
        comment = "slope for all ccds (mag lt %.1f, nstar=%d) limits=(%.1f,%.1f sigma)" % (self.magCut, len(allDiffs), self.slopeLimits[0], self.slopeLimits[1])
        testSet.addTest( testCode.Test(label, lineCoeffs[0], slopeLimits, comment, areaLabel="all"))


    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId, label=self.magType1+"-"+self.magType2)
        testSet.setUseCache(self.useCache)

        xlim = [14.0, 25.0]
        ylimStep = 0.4
        ylim = [-ylimStep, ylimStep]
	aspRatio = (xlim[1]-xlim[0])/(ylim[1]-ylim[0])

        tag1 = "m$_{\mathrm{"+self.magType1.upper()+"}}$"
        tag = "m$_{\mathrm{"+self.magType1.upper()+"}}$ - m$_{\mathrm{"+self.magType2.upper()+"}}$"
        dtag = self.magType1+"-"+self.magType2
        wtag = self.magType1+"minus"+self.magType2

        # fpa figure
        meanFilebase = "mean" + wtag
        stdFilebase  = "std"+wtag
        slopeFilebase  = "slope"+wtag
        meanData, meanMap   = testSet.unpickle(meanFilebase, default=[None, None])
        stdData, stdMap     = testSet.unpickle(stdFilebase, default=[None, None])
        slopeData, slopeMap = testSet.unpickle(slopeFilebase, default=[None, None])

        meanFig  = qaFig.FpaQaFigure(data.cameraInfo, data=meanData, map=meanMap)
        stdFig   = qaFig.FpaQaFigure(data.cameraInfo, data=stdData, map=stdMap)
        slopeFig = qaFig.VectorFpaQaFigure(data.cameraInfo, data=slopeData, map=slopeMap)

        for raft, ccd in self.means.raftCcdKeys():

            meanFig.data[raft][ccd] = self.means.get(raft, ccd)
            stdFig.data[raft][ccd] = self.stds.get(raft, ccd)
            slope = self.trend.get(raft, ccd)

            if not slope is None and not slope[1] == 0:
                # aspRatio will make the vector have the same angle as the line in the figure
                slopeSigma = slope[0]/slope[1]
                slopeFig.data[raft][ccd] = [numpy.arctan2(aspRatio*slope[0],1.0), None, slopeSigma]
            else:
                slopeSigma = None
                slopeFig.data[raft][ccd] = [None, None, None]

            if not self.means.get(raft, ccd) is None:
                meanFig.map[raft][ccd] = "mean=%.4f" % (self.means.get(raft, ccd))
                stdFig.map[raft][ccd] = "std=%.4f" % (self.stds.get(raft, ccd))
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

        meanFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[-0.03, 0.03],
                           title="Mean "+tag, cmapOver=red, cmapUnder=blue, failLimits=self.deltaLimits)
        testSet.addFigure(meanFig, meanFilebase+".png",
                          "mean "+dtag+" mag   (brighter than %.1f)" % (self.magCut), navMap=True)
        testSet.pickle(meanFilebase, [meanFig.data, meanFig.map])
        
        stdFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 0.03],
                          title="Stdev "+tag, cmapOver=red, failLimits=self.rmsLimits)
        testSet.addFigure(stdFig, stdFilebase+".png",
                          "stdev "+dtag+" mag  (brighter than %.1f)" % (self.magCut), navMap=True)
        testSet.pickle(stdFilebase, [stdFig.data, stdFig.map])
        

        cScale = 2.0
        slopeFig.makeFigure(cmap="RdBu_r", vlimits=[cScale*self.slopeLimits[0], cScale*self.slopeLimits[1]],
                            title="Slope "+tag, failLimits=self.slopeLimits)
        testSet.addFigure(slopeFig, slopeFilebase+".png",
                          "slope "+dtag+" mag (brighter than %.1f)" % (self.magCut), navMap=True)
        testSet.pickle(slopeFilebase, [slopeFig.data, slopeFig.map])


        #############################################
        #
        figsize = (6.5, 3.75)

        conv = colors.ColorConverter()
        red = conv.to_rgba('r')
        black = conv.to_rgba('k')
        size = 1.0
        
        xmin, xmax = self.mag.summarize('min', default=0.0), self.mag.summarize('max', default=25.0)
        ymin, ymax = self.diff.summarize('min', default=-1.0), self.diff.summarize('max', default=1.0)
        xrang = xmax-xmin
        xmin, xmax = int(xmin-0.05*xrang), int(xmax+0.05*xrang)+1
        yrang = ymax-ymin
        ymin, ymax = ymin-0.05*yrang, ymax+0.05*yrang
        xlim2 = xlim        #[xmin, xmax]
        ylim2 = [-2.0, 2.0] #[ymin, ymax]

        # grab any cached values
        figbase = "diff_" + dtag
        allMags, allDiffs, allStars, allColor, allLabels = \
                 testSet.unpickle(figbase, default=[{},{},{},{},{}])

        colorId = {}
        i = 0
        for k in data.cameraInfo.detectors.keys():
            colorId[k] = i
            i += 1

        ccdKeysThisRun = list(zip(*self.mag.raftCcdKeys())[1])
        ccdKeysCache   = allMags.keys()
        allCcds = set(ccdKeysCache + ccdKeysThisRun)
        colorIdList = []
        for ccd in allCcds:
            colorIdList.append(colorId[ccd])
        minColorId = numpy.min(colorIdList)
        maxColorId = numpy.max(colorIdList)
        #nKeys = len(colorId.keys())
        norm = colors.Normalize(vmin=minColorId, vmax=maxColorId)
        sm = cm.ScalarMappable(norm, cmap=cm.jet)

        for raft, ccd in self.mag.raftCcdKeys():
            mag  = self.mag.get(raft, ccd)
            diff = self.diff.get(raft, ccd)
            x    = self.x.get(raft, ccd)
            y    = self.y.get(raft, ccd)
            star = self.star.get(raft, ccd)
            lineFit = self.trend.get(raft, ccd)
            trendCoeffs = lineFit[0], lineFit[2]
            trendCoeffsLo = lineFit[0]+lineFit[1], lineFit[2]-lineFit[3]
            trendCoeffsHi = lineFit[0]-lineFit[1], lineFit[2]+lineFit[3]
            
            if len(mag) == 0:
                mag = numpy.array([xmax])
                diff = numpy.array([0.0])
                x    = numpy.array([0.0])
                y    = numpy.array([0.0])
                star = numpy.array([0])

            whereCut = numpy.where((mag < self.magCut) & (star > 0))[0] # & (abs(diff) < self.dmagMax))[0]
            print "plotting ", ccd

            #################
            # data for one ccd
            fig = qaFig.QaFigure(size=figsize)
            fig.fig.subplots_adjust(left=0.125, bottom=0.125)
            ax_1 = fig.fig.add_subplot(121)
            ax_2 = fig.fig.add_subplot(122)
            clr = [black] * len(mag)
            clr = numpy.array(clr)
            if len(whereCut) > 0:
                clr[whereCut] = [red] * len(whereCut)

            xTrend = numpy.array(xlim)
            ax_1.plot(xTrend, numpy.array([0.0, 0.0]), "-k", lw=1.0)
            
            for ax in [ax_1, ax_2]:
                ax.scatter(mag, diff, size, color=clr, label=ccd)
                ax.set_xlabel(tag1)

            lineVals = numpy.lib.polyval(trendCoeffs, xTrend)
            lineValsLo = numpy.lib.polyval(trendCoeffsLo, xTrend)
            lineValsHi = numpy.lib.polyval(trendCoeffsHi, xTrend)
            

            if abs(trendCoeffs[0] - 99.0) > 1.0e-6:
                lmin, lmax = lineVals.min(), lineVals.max()
                if lmin < ylim[0]:
                    ylim[0] = -(int(abs(lmin)/ylimStep) + 1)*ylimStep
                if lmax > ylim[1]:
                    ylim[1] = (int(lmax/ylimStep) + 1)*ylimStep

                ax_1.plot(xTrend, lineVals, "-r", lw=1.0)
                ax_1.plot(xTrend, lineValsLo, "--r", lw=1.0)
                ax_1.plot(xTrend, lineValsHi, "--r", lw=1.0)
                ax_2.plot(xTrend, lineVals, "-r", lw=1.0)
                ax_2.plot(xTrend, lineValsLo, "--r", lw=1.0)
                ax_2.plot(xTrend, lineValsHi, "--r", lw=1.0)
            
            ax_2.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]],
                      [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], '-k')
            ax_1.set_ylabel(tag)
            ax_1.set_xlim(xlim)
            ax_2.set_xlim(xlim2)
            ax_1.set_ylim(ylim)
            ax_2.set_ylim(ylim2)

            # move the y axis on right panel
            ax_2.yaxis.set_label_position('right')
            ax_2.yaxis.set_ticks_position('right')

            dmag = 0.1
            ddiff1 = 0.02
            ddiff2 = ddiff1*(ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]) # rescale for larger y range
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            for j in xrange(len(mag)):
                info = "nolink:x:%.2f_y:%.2f" % (x[j], y[j])
                area = (mag[j]-dmag, diff[j]-ddiff1, mag[j]+dmag, diff[j]+ddiff1)
                fig.addMapArea(areaLabel, area, info, axes=ax_1)
                area = (mag[j]-dmag, diff[j]-ddiff2, mag[j]+dmag, diff[j]+ddiff2)
                fig.addMapArea(areaLabel, area, info, axes=ax_2)
                

            testSet.addFigure(fig, figbase+".png",
                              dtag+" vs. "+self.magType1 + ". Point used for statistics shown in red.",
                              areaLabel=areaLabel)


            # append values to arrays for a plot showing all data
            allMags[ccd] = mag  
            allDiffs[ccd] = diff
            allStars[ccd] = star
            allColor[ccd] = colorId[ccd]*numpy.ones(len(mag))
            allLabels[ccd] = [areaLabel] * len(mag)


        # stash values
        testSet.pickle(figbase, [allMags, allDiffs, allStars, allColor, allLabels])

        withDelete = True
        allMags   = qaAnaUtil.dictToList(allMags, withDelete=withDelete)
        allDiffs  = qaAnaUtil.dictToList(allDiffs, withDelete=withDelete)
        allStars  = qaAnaUtil.dictToList(allStars, withDelete=withDelete)
        allColor  = qaAnaUtil.dictToList(allColor, withDelete=withDelete)
        allLabels = qaAnaUtil.dictToList(allLabels, withDelete=withDelete)

        allColor = sm.to_rgba(allColor)
        
        # dmag vs mag
        fig0 = qaFig.QaFigure(size=figsize)
        fig0.fig.subplots_adjust(left=0.125, bottom=0.125)
        ax0_1 = fig0.fig.add_subplot(121)
        ax0_2 = fig0.fig.add_subplot(122)

        w = numpy.where( (allMags < self.magCut) & (allStars > 0))[0] # & (abs(allDiffs) < self.dmagMax))[0]

        if len(w) > 0:
            lineFit = qaAnaUtil.robustPolyFit(allMags[w], allDiffs[w], 1)
            trendCoeffs = lineFit[0], lineFit[2]
            trendCoeffsLo = lineFit[0]+lineFit[1], lineFit[2]-lineFit[3]
            trendCoeffsHi = lineFit[0]-lineFit[1], lineFit[2]+lineFit[3]
        else:
            trendCoeffs = [0.0, 0.0]
            trendCoeffsLo = [0.0, 0.0]
            trendCoeffsHi = [0.0, 0.0]

        ####################
        # data for all ccds
        if len(allMags) == 0:
            allMags = numpy.array([xmax])
            allDiffs = numpy.array([0.0])
            allColor = [black]
            allLabels = ["no_valid_data"]
            trendCoeffs = [0.0, 0.0]
            trendCoeffsLo = [0.0, 0.0]
            trendCoeffsHi = [0.0, 0.0]

        allColor = numpy.array(allColor)
        for ax in [ax0_1, ax0_2]:
            ax.plot(xlim2, [0.0, 0.0], "-k", lw=1.0)  # show an x-axis at y=0
            ax.scatter(allMags, allDiffs, size, color=allColor)
            # 99 is the 'no-data' values
            if abs(trendCoeffs[0] - 99.0) > 1.0e-6:
                ax.plot(xlim2, numpy.polyval(trendCoeffs, xlim2), "-k", lw=1.0)
                ax.plot(xlim2, numpy.polyval(trendCoeffsLo, xlim2), "--k", lw=1.0)
                ax.plot(xlim2, numpy.polyval(trendCoeffsHi, xlim2), "--k", lw=1.0)
        ax0_1.set_xlim(xlim)
        ax0_2.set_xlim(xlim2)
        ax0_1.set_ylim(ylim)
        ax0_2.set_ylim(ylim2)

        dmag = 0.1
        ddiff1 = 0.01
        ddiff2 = ddiff1*(ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]) # rescale for larger y range
        for j in xrange(len(allMags)):
            area = (allMags[j]-dmag, allDiffs[j]-ddiff1, allMags[j]+dmag, allDiffs[j]+ddiff1)
            fig0.addMapArea(allLabels[j], area, "%.3f_%.3f"% (allMags[j], allDiffs[j]), axes=ax0_1)
            area = (allMags[j]-dmag, allDiffs[j]-ddiff2, allMags[j]+dmag, allDiffs[j]+ddiff2)
            fig0.addMapArea(allLabels[j], area, "%.3f_%.3f"% (allMags[j], allDiffs[j]), axes=ax0_2)



        # move the yaxis ticks/labels to the other side
        ax0_2.yaxis.set_label_position('right')
        ax0_2.yaxis.set_ticks_position('right')

        ax0_2.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]],
                   [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], '-k')
        ax0_2.set_xlim(xlim2)
        ax0_2.set_ylim(ylim2)

        for ax in [ax0_1, ax0_2]:
            ax.set_xlabel(tag1)
            ax.set_ylabel(tag)

        testSet.addFigure(fig0, figbase+".png", dtag+" vs. "+self.magType1, areaLabel="all")

        
