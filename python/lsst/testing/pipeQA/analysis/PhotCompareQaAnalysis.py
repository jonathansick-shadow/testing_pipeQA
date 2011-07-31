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
                 deltaMin, deltaMax, rmsMax, slopeMinSigma, slopeMaxSigma, starGalaxyToggle,
                 **kwargs):
        testLabel = magType1+"-"+magType2
        qaAna.QaAnalysis.__init__(self, testLabel, **kwargs)

        self.magCut = magCut
        self.deltaLimits = [deltaMin, deltaMax]
        self.rmsLimits = [0.0, rmsMax]
        self.slopeLimits = [-slopeMinSigma, slopeMaxSigma]
        self.starGalaxyToggle = starGalaxyToggle
        
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
         this distribution (stars plotted in red) reflects the systematic floor
         in these measurement comparisons.  For psf magnitudes, this is
         typically 1-2% for PT1.2 measurements.  The summary figure showing all
         data from all CCDs in the FPA shows stars (red) and galaxies (green), but 
         does not include the x,y panel.  The summary FPA figures show
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
        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER | measAlg.Flags.EDGE

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
        testSet.addMetadata({"Description": self.description})

        self.means = raftCcdData.RaftCcdData(self.detector)
        self.medians = raftCcdData.RaftCcdData(self.detector)
        self.stds  = raftCcdData.RaftCcdData(self.detector)
        self.trend = raftCcdData.RaftCcdData(self.detector, initValue=[0.0, 0.0])
        
        self.dmagMax = 0.4
        allMags = numpy.array([])
        allDiffs = numpy.array([])

        for raft,  ccd in self.mag.raftCcdKeys():
            dmag0 = self.diff.get(raft, ccd)
            mag0 = self.mag.get(raft, ccd)
            star = self.star.get(raft, ccd)

            wGxy = numpy.where((mag0 > 10) & (mag0 < self.magCut) &
                               (star == 0) & (numpy.abs(dmag0) < 1.0))[0]
            w = numpy.where((mag0 > 10) & (mag0 < self.magCut) & (star > 0))[0]

            mag = mag0[w]
            dmag = dmag0[w]

            allMags = numpy.append(allMags, mag)
            allDiffs = numpy.append(allDiffs, dmag)
             
            # already using NaN for 'no-data' for this ccd
            #  (because we can't test for 'None' in a numpy masked_array)
            # unfortunately, these failures will have to do
            mean = 99.0
            median = 99.0
            std = 99.0
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
            slopeLimits = self.slopeLimits[0]*lineFit[0][1], self.slopeLimits[1]*lineFit[0][1]
            comment = "slope of "+dtag+" (mag lt %.1f, nstar/clip=%d/%d) limits=(%.1f,%.1f)sigma" % \
                      (self.magCut, len(dmag), n, self.slopeLimits[0], self.slopeLimits[1])
            testSet.addTest( testCode.Test(label, lineCoeffs[0][0], slopeLimits, comment,
                                           areaLabel=areaLabel))


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

        isFinalDataId = False
        if len(data.brokenDataIdList) > 0 and data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        xlim = [14.0, 25.0]
        ylimStep = 0.4
        ylim = [-ylimStep, ylimStep]
	aspRatio = (xlim[1]-xlim[0])/(ylim[1]-ylim[0])

        tag1 = "m$_{\mathrm{"+self.magType1.upper()+"}}$"
        tag  = "m$_{\mathrm{"+self.magType1.upper()+"}}$ - m$_{\mathrm{"+self.magType2.upper()+"}}$"
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
            slope = self.trend.get(raft, ccd)[0]

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


        testSet.pickle(meanFilebase, [meanFig.data, meanFig.map])
        testSet.pickle(stdFilebase, [stdFig.data, stdFig.map])
        testSet.pickle(slopeFilebase, [slopeFig.data, slopeFig.map])

        if not self.delaySummary or isFinalDataId:
            print "plotting FPAs"
            meanFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[-0.03, 0.03],
                               title="Mean "+tag, cmapOver=red, cmapUnder=blue, failLimits=self.deltaLimits)
            testSet.addFigure(meanFig, meanFilebase+".png",
                              "mean "+dtag+" mag   (brighter than %.1f)" % (self.magCut), navMap=True)
            del meanFig
            
            stdFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 0.03],
                              title="Stdev "+tag, cmapOver=red, failLimits=self.rmsLimits)
            testSet.addFigure(stdFig, stdFilebase+".png",
                              "stdev "+dtag+" mag  (brighter than %.1f)" % (self.magCut), navMap=True)
            del stdFig
            
            cScale = 2.0
            slopeFig.makeFigure(cmap="RdBu_r",
                                vlimits=[cScale*self.slopeLimits[0], cScale*self.slopeLimits[1]],
                                title="Slope "+tag, failLimits=self.slopeLimits)
            testSet.addFigure(slopeFig, slopeFilebase+".png",
                              "slope "+dtag+" mag (brighter than %.1f)" % (self.magCut), navMap=True)
            del slopeFig
        else:
            del meanFig
            del stdFig
            del slopeFig


        #############################################
        #
        #xmin, xmax = self.mag.summarize('min', default=0.0), self.mag.summarize('max', default=25.0)
        #ymin, ymax = self.diff.summarize('min', default=-1.0), self.diff.summarize('max', default=1.0)
        #xrang = xmax-xmin
        #xmin, xmax = int(xmin-0.05*xrang), int(xmax+0.05*xrang)+1
        #yrang = ymax-ymin
        #ymin, ymax = ymin-0.05*yrang, ymax+0.05*yrang
        xlim2 = xlim        #[xmin, xmax]
        ylim2 = [-2.0, 2.0] #[ymin, ymax]

        figsize = (6.5, 3.75)

        figbase = "diff_" + dtag

        shelfData = {}

        for raft, ccd in self.mag.raftCcdKeys():
            mag0  = self.mag.get(raft, ccd)
            diff0 = self.diff.get(raft, ccd)
            star0 = self.star.get(raft, ccd)
            
            print "plotting ", ccd

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            statBlurb = "Points used for statistics/trendline shown in red."

            if self.starGalaxyToggle:
                fig = self.regularFigure([mag0, diff0, star0,
                                          areaLabel, raft, ccd, figsize,
                                          xlim, ylim, xlim2, ylim2, ylimStep,
                                          tag1, tag, 'stars'])
                testSet.addFigure(fig, figbase+".png",
                                  dtag+" vs. "+self.magType1 + "(stars). "+statBlurb,
                                  areaLabel=areaLabel, toggle='stars')
                del fig

                fig = self.regularFigure([mag0, diff0, star0,
                                          areaLabel, raft, ccd, figsize,
                                          xlim, ylim, xlim2, ylim2, ylimStep,
                                          tag1, tag, 'galaxies'])
                testSet.addFigure(fig, figbase+".png",
                                  dtag+" vs. "+self.magType1 + "(galaxies). "+statBlurb,
                                  areaLabel=areaLabel, toggle='galaxies')
                del fig
                fig = self.regularFigure([mag0, diff0, star0,
                                          areaLabel, raft, ccd, figsize,
                                          xlim, ylim, xlim2, ylim2, ylimStep,
                                          tag1, tag, 'all'])
                testSet.addFigure(fig, figbase+".png",
                                  dtag+" vs. "+self.magType1 + ". "+statBlurb,
                                  areaLabel=areaLabel, toggle="everything")
                del fig
                
            else:
                fig = self.regularFigure([mag0, diff0, star0,
                                          areaLabel, raft, ccd, figsize,
                                          xlim, ylim, xlim2, ylim2, ylimStep, 
                                          tag1, tag, 'all'])
                testSet.addFigure(fig, figbase+".png",
                                  dtag+" vs. "+self.magType1 + ". "+statBlurb,
                                  areaLabel=areaLabel)
                del fig
                
            
            # append values to arrays for a plot showing all data
            shelfData[ccd] = [mag0, diff0, star0, [areaLabel]*len(mag0)]
            
        # stash values
        if self.useCache:
            testSet.shelve(figbase, shelfData)

        if not self.delaySummary or isFinalDataId:
            print "plotting Summary figure"

            self.summaryFigure([
                data, figbase, testSet, shelfData, figsize,
                xlim, ylim, xlim2, ylim2,
                tag1, tag, dtag,
                ])






    def regularFigure(self, args):

        mag0, diff0, star0, areaLabel, raft, ccd, figsize, xlim, ylim, xlim2, ylim2, ylimStep, \
              tag1, tag, mode = args

        conv = colors.ColorConverter()
        size = 2.0
        
        red = conv.to_rgba('r')
        black = conv.to_rgba('k')

        x0    = self.x.get(raft, ccd)
        y0    = self.y.get(raft, ccd)
        modeIdx = {'all': 0, 'galaxies': 1, 'stars': 2 }
        lineFit = self.trend.get(raft, ccd)[modeIdx[mode]]
            
        trendCoeffs = lineFit[0], lineFit[2]
        trendCoeffsLo = lineFit[0]+lineFit[1], lineFit[2]-lineFit[3]
        trendCoeffsHi = lineFit[0]-lineFit[1], lineFit[2]+lineFit[3]
	#print trendCoeffs        
        if len(mag0) == 0:
            mag0 = numpy.array([xlim[1]])
            diff0 = numpy.array([0.0])
            x0    = numpy.array([0.0])
            y0    = numpy.array([0.0])
            star0 = numpy.array([0])

        #################
        # data for one ccd
        if mode == 'fourPanel':
            figsize = (6.5, 5.0)
        
        fig = qaFig.QaFigure(size=figsize)
        fig.fig.subplots_adjust(left=0.09, right=0.93, bottom=0.125)

        if mode == 'fourPanel':
            ax_1s = fig.fig.add_subplot(221)
            ax_2s = fig.fig.add_subplot(222)
            ax_1g = fig.fig.add_subplot(223)
            ax_2g = fig.fig.add_subplot(224)
            axSets        = [ [ax_1s, ax_2s], [ax_1g, ax_2g] ]
            starGalLabels = ["stars", "galaxies"]
            whereStarGals = [numpy.where(star0 == 0)[0], numpy.where(star0 > 0)[0] ]
        else:
            ax_1 = fig.fig.add_subplot(131)
            ax_2 = fig.fig.add_subplot(132)
            ax_3 = fig.fig.add_subplot(133)
            axSets        = [ [ax_1, ax_2, ax_3] ]
            starGalLabels = [mode]
            if mode == 'stars':
                whereStarGals = [numpy.where(star0 > 0)[0] ]
            if mode == 'galaxies':
                whereStarGals = [numpy.where(star0 == 0)[0] ]
            if mode == 'all':
                starGalLabels = ["all data"]
                whereStarGals = [numpy.where(star0 > -1)[0] ]
        
        for iSet in range(len(axSets)):
            ax_1, ax_2, ax_3   = axSets[iSet]
            starGalLabel = starGalLabels[iSet]
            whereStarGal = whereStarGals[iSet]

            mag  = mag0[whereStarGal]
            diff = diff0[whereStarGal]
            x    = x0[whereStarGal]
            y    = y0[whereStarGal]
            star = star0[whereStarGal]

            whereCut = numpy.where((mag < self.magCut))[0]
            whereOther = numpy.where((mag > self.magCut))[0]

            xTrend = numpy.array(xlim)
            ax_1.plot(xTrend, numpy.array([0.0, 0.0]), "-k", lw=1.0)

            ax_1.text(1.02*xlim[0], 0.87*ylim[1], starGalLabel, size='x-small', horizontalalignment='left')
            for ax in [ax_1, ax_2]:
                ax.plot(mag[whereOther], diff[whereOther], "k.", ms=size, label=ccd)
                ax.plot(mag[whereCut], diff[whereCut], "r.", ms=size, label=ccd)
                ax.set_xlabel(tag1, size='small')


            norm = colors.Normalize(vmin=-0.05, vmax=0.05, clip=True)
            cdict = {'red': ((0.0, 0.0, 0.0),
                             (0.5, 0.0, 0.0),
                             (1.0, 1.0, 1.0)),
                     'green': ((0.0, 0.0, 0.0),
                               (0.5, 0.0, 0.0),
                               (1.0, 0.0, 0.0)),
                     'blue': ((0.0, 1.0, 1.0),
                              (0.5, 0.0, 0.0),
                              (1.0, 0.0, 0.0))}
            my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)
            
            midSize = 2.0
            maxSize = 6*midSize
            minSize = midSize/2
            magLo = xlim[0] if xlim[0] > mag0.min() else mag0.min()
            magHi = self.magCut + 2.0 #xlim[1] if xlim[1] < mag0.max() else mag0.max()
            sizes = minSize + (maxSize - minSize)*(magHi - mag)/(magHi - magLo)
            sizes = numpy.clip(sizes, minSize, maxSize)
            xyplot = ax_3.scatter(x, y, s=sizes, c=diff, marker='o',
                                  cmap=my_cmap, norm=norm, edgecolors='none')
            ax_3.set_xlim([x0.min(), x0.max()])
            ax_3.set_ylim([y0.min(), y0.max()])
            ax_3.set_xlabel("x", size='x-small')
            #ax_3.set_ylabel("y", labelpad=20, y=1.0, size='x-small', rotation=0.0)
            cb = fig.fig.colorbar(xyplot)
            cb.ax.set_xlabel("$\delta$ mag", size="x-small")
            cb.ax.xaxis.set_label_position('top')
            for tick in cb.ax.get_yticklabels():
                tick.set_size("x-small")
            for t in ax_3.get_xticklabels() + ax_3.get_yticklabels():
                t.set_rotation(45.0)
                t.set_size('xx-small')
            #for t in ax_3.get_yticklabels():
            #    t.set_size('xx-small')
                
            
            lineVals = numpy.lib.polyval(trendCoeffs, xTrend)
            lineValsLo = numpy.lib.polyval(trendCoeffsLo, xTrend)
            lineValsHi = numpy.lib.polyval(trendCoeffsHi, xTrend)


            if abs(trendCoeffs[0] - 99.0) > 1.0e-6:
                lmin, lmax = lineVals.min(), lineVals.max()
                if False:
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
            ax_1.set_ylabel(tag, size="small")
            ax_1.set_xlim(xlim)
            ax_2.set_xlim(xlim2)
            ax_1.set_ylim(ylim)
            ax_2.set_ylim(ylim2)

            # move the y axis on right panel
            #ax_3.yaxis.set_label_position('right')
            #ax_3.yaxis.set_ticks_position('right')

            dmag = 0.1
            ddiff1 = 0.02
            ddiff2 = ddiff1*(ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]) # rescale for larger y range
            
            for j in xrange(len(mag)):
                info = "nolink:x:%.2f_y:%.2f" % (x[j], y[j])
                area = (mag[j]-dmag, diff[j]-ddiff1, mag[j]+dmag, diff[j]+ddiff1)
                fig.addMapArea(areaLabel, area, info, axes=ax_1)
                area = (mag[j]-dmag, diff[j]-ddiff2, mag[j]+dmag, diff[j]+ddiff2)
                fig.addMapArea(areaLabel, area, info, axes=ax_2)

            for ax in [ax_1, ax_2]:
                for tick in ax.get_xticklabels() + ax.get_yticklabels():
                    tick.set_size("x-small")

        return fig



                    
    def summaryFigure(self, summaryArgs):

        data, figbase, testSet, shelfData, figsize, xlim, ylim, xlim2, ylim2, tag1, tag, dtag = summaryArgs

        size = 2.0
        
        # unstash the values
        if self.useCache:
            shelfData = testSet.unshelve(figbase)

        #colorId = {}
        #i = 0
        #for k in sorted(data.cameraInfo.sensors.keys()):
        #    colorId[k] = i
        #    i += 1

        allMags = numpy.array([])
        allDiffs = numpy.array([])
        allStars = numpy.array([])
        #allColor = numpy.array([])
        allLabels = numpy.array([])
        allCcds = []
        for k,v in shelfData.items():
            allCcds.append(k)
            mags, diffs, stars, arealabels = v
            allMags = numpy.append(allMags, mags)
            allDiffs = numpy.append(allDiffs, diffs)
            allStars = numpy.append(allStars, stars)
            #allColor = numpy.append(allColor, colorId[k]*numpy.ones(len(mags)))
            allLabels = numpy.append(allLabels, arealabels)

        # figure out the color map
        #colorIdList = []
        #for ccd in allCcds:
        #    clr = colorId[ccd]
        #    colorIdList.append(clr)
        #minColorId = numpy.min(colorIdList)
        #maxColorId = numpy.max(colorIdList)
        #norm = colors.Normalize(vmin=minColorId, vmax=maxColorId)
        #sm = cm.ScalarMappable(norm, cmap=cm.jet)
        #allColor = sm.to_rgba(allColor)

        # dmag vs mag
        fig0 = qaFig.QaFigure(size=figsize)
        fig0.fig.subplots_adjust(left=0.125, bottom=0.125)
        ax0_1 = fig0.fig.add_subplot(121)
        ax0_2 = fig0.fig.add_subplot(122)

        w = numpy.where( (allMags < self.magCut) & (allStars > 0))[0]
        wStar = numpy.where(allStars > 0)[0]
        wGxy = numpy.where(allStars == 0)[0]

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
            allMags = numpy.array([xlim[1]])
            allDiffs = numpy.array([0.0])
            #allColor = [black]
            allLabels = ["no_valid_data"]
            trendCoeffs = [0.0, 0.0]
            trendCoeffsLo = [0.0, 0.0]
            trendCoeffsHi = [0.0, 0.0]

        #allColor = numpy.array(allColor)
        for ax in [ax0_1, ax0_2]:
            ax.plot(xlim2, [0.0, 0.0], "-k", lw=1.0)  # show an x-axis at y=0
            ax.plot(allMags[wGxy], allDiffs[wGxy], "g.", ms=size, label="galaxies") #0.75, alpha=0.5)
            ax.plot(allMags[wStar], allDiffs[wStar], "r.", ms=size, label="stars") #, alpha=0.5)
            # 99 is the 'no-data' values
            if abs(trendCoeffs[0] - 99.0) > 1.0e-6:
                ax.plot(xlim2, numpy.polyval(trendCoeffs, xlim2), "-k", lw=1.0)
                ax.plot(xlim2, numpy.polyval(trendCoeffsLo, xlim2), "--k", lw=1.0)
                ax.plot(xlim2, numpy.polyval(trendCoeffsHi, xlim2), "--k", lw=1.0)

        # show outliers railed at the ylims
        wStarOutlier = numpy.where((numpy.abs(allDiffs) > ylim2[1]) & (allStars > 0))[0]
        wGxyOutlier = numpy.where((numpy.abs(allDiffs) > ylim2[1]) & (allStars == 0))[0]
        clip = 0.99
        ax0_2.plot(allMags[wStarOutlier], numpy.clip(allDiffs[wStarOutlier], clip*ylim2[0], clip*ylim2[1]),
                   'r.', ms=size)
        ax0_2.plot(allMags[wGxyOutlier], numpy.clip(allDiffs[wGxyOutlier], clip*ylim2[0], clip*ylim2[1]),
                   'g.', ms=size)
        ax0_2.legend(prop=fm.FontProperties(size="x-small"), loc="upper left")
        
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


        del allMags
        del allDiffs
        #del allColor
        del allLabels

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

        del fig0

