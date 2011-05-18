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

    def __init__(self, magType1='psf', magType2='aperture', cut=20.0):
        testLabel = magType1+"-"+magType2
        qaAna.QaAnalysis.__init__(self, testLabel)

        self.cut = cut
        
        def magType(mType):
            if re.search("(psf|PSF)", mType):
                return "psf"
            elif re.search("^ap", mType):
                return "ap"
            elif re.search("^mod", mType):
                return "mod"
            elif re.search("^cat", mType):
                return "cat"

        self.magType1 = magType(magType1)
        self.magType2 = magType(magType2)


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

    def free(self):
        del self.x
        del self.y
        del self.mag
        del self.diff
        del self.filter
        del self.detector
        if not self.matchListDict is None:
            del self.matchListDict
        if not self.ssDict is None:
            del self.ssDict
        del self.means
        del self.medians
        del self.stds
        

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

        self.matchListDict = None
        self.ssDict = None

        # if we're asked to compare catalog fluxes ... we need a matchlist
        if  self.magType1=="cat" or self.magType2=="cat":
            self.matchListDict = data.getMatchListBySensor(dataId)
            for key, matchList in self.matchListDict.items():
                raft = self.detector[key].getParent().getId().getName()
                ccd  = self.detector[key].getId().getName()
                filter = self.filter[key].getName()

                qaAnaUtil.isStar(matchList)

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
        self.deltaLimits = [-0.02, 0.02]
        self.rmsLimits = [0.0, 0.02]
        for raft,  ccd in self.mag.raftCcdKeys():
            dmag = self.diff.get(raft, ccd)
            mag = self.mag.get(raft, ccd)
            star = self.star.get(raft, ccd)
            w = numpy.where((mag > 10) & (mag < self.cut) & (star > 0))
            dmag = dmag[w]

            if len(dmag) > 0:
                stat = afwMath.makeStatistics(dmag, afwMath.NPOINT | afwMath.MEANCLIP |
                                              afwMath.STDEVCLIP | afwMath.MEDIAN)
                mean = stat.getValue(afwMath.MEANCLIP)
                median = stat.getValue(afwMath.MEDIAN)
                std = stat.getValue(afwMath.STDEVCLIP)
                n = stat.getValue(afwMath.NPOINT)

            else:
                # already using NaN for 'no-data' for this ccd
                #  (because we can't test for 'None' in a numpy masked_array)
                # unfortunately, these failures will have to do
                mean = 99.0
                median = 99.0
                std = 99.0
                n = 0

            tag = self.magType1+"_vs_"+self.magType2
            dtag = self.magType1+"-"+self.magType2
            self.means.set(raft, ccd, mean)
            areaLabel = re.sub("\s+", "_", ccd)
            label = "mean "+tag +" " + areaLabel
            comment = "mean "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.cut, len(dmag),n)
            testSet.addTest( testCode.Test(label, mean, self.deltaLimits, comment), areaLabel=areaLabel )

            self.medians.set(raft, ccd, median)
            label = "median "+tag+" "+areaLabel
            comment = "median "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.cut, len(dmag), n)
            testSet.addTest( testCode.Test(label, median, self.deltaLimits, comment), areaLabel=areaLabel )

            self.stds.set(raft, ccd, std)
            label = "stdev "+tag+" " + areaLabel
            comment = "stdev of "+dtag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.cut, len(dmag), n)
            testSet.addTest( testCode.Test(label, std, self.rmsLimits, comment), areaLabel=areaLabel )
                


    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId, label=self.magType1+"-"+self.magType2)

        # fpa figure
        meanFig = qaFig.FpaQaFigure(data.cameraInfo.camera)
        stdFig = qaFig.FpaQaFigure(data.cameraInfo.camera)
        for raft, ccdDict in meanFig.data.items():
            for ccd, value in ccdDict.items():
                meanFig.data[raft][ccd] = self.means.get(raft, ccd)
                stdFig.data[raft][ccd] = self.stds.get(raft, ccd)
                if not self.means.get(raft, ccd) is None:
                    meanFig.map[raft][ccd] = "mean=%.4f" % (self.means.get(raft, ccd))
                    stdFig.map[raft][ccd] = "std=%.4f" % (self.stds.get(raft, ccd))

        tag1 = "m$_{\mathrm{"+self.magType1.upper()+"}}$"
        tag = "m$_{\mathrm{"+self.magType1.upper()+"}}$ - m$_{\mathrm{"+self.magType2.upper()+"}}$"
        dtag = self.magType1+"-"+self.magType2
        wtag = self.magType1+"minus"+self.magType2
        deepPink = '#ff1493'
        darkViolet = '#9400d3'
        blue = '#0000ff'
        red = '#ff0000'
        meanFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[-0.03, 0.03],
                           title="Mean "+tag, cmapOver=red, cmapUnder=blue, failLimits=self.deltaLimits)
        testSet.addFigure(meanFig, "mean"+wtag+".png", "mean "+dtag+" mag   (brighter than %.1f)" % (self.cut),
                          navMap=True)
        stdFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 0.03],
                          title="Stdev "+tag, cmapOver=red, failLimits=self.rmsLimits)
        testSet.addFigure(stdFig, "std"+wtag+".png", "stdev "+dtag+" mag  (brighter than %.1f)" % (self.cut),
                          navMap=True)
        



        #############################################
        #
        figsize = (6.5, 3.75)
        nKeys = len(self.mag.raftCcdKeys())
        norm = colors.Normalize(vmin=0, vmax=nKeys)
        sm = cm.ScalarMappable(norm, cmap=cm.jet)

        xlim = [14.0, 25.0]
        ylim = [-0.4, 0.4]

        conv = colors.ColorConverter()
        red = conv.to_rgba('r')
        black = conv.to_rgba('k')
        size = 1.0
        
        i = 0
        xmin, xmax = self.mag.summarize('min', default=0.0), self.mag.summarize('max', default=25.0)
        ymin, ymax = self.diff.summarize('min', default=-1.0), self.diff.summarize('max', default=1.0)
        xrang = xmax-xmin
        xmin, xmax = xmin-0.05*xrang, xmax+0.05*xrang
        yrang = ymax-ymin
        ymin, ymax = ymin-0.05*yrang, ymax+0.05*yrang
        xlim2 = [xmin, xmax]
        ylim2 = [ymin, ymax]

        allMags = numpy.array([])
        allDiffs = numpy.array([])
        allColor = [] #numpy.array([])
        allLabels = []
        for raft, ccd in self.mag.raftCcdKeys():
            mag  = self.mag.get(raft, ccd)
            diff = self.diff.get(raft, ccd)
            x    = self.x.get(raft, ccd)
            y    = self.y.get(raft, ccd)
            star = self.star.get(raft, ccd)
            
            if len(mag) == 0:
                mag = numpy.array([xmax])
                diff = numpy.array([0.0])
                x    = numpy.array([0.0])
                y    = numpy.array([0.0])
                star = numpy.array([0])
                
            whereCut = numpy.where((mag < self.cut) & (star > 0))

            print "plotting ", ccd

            #################
            # data for one ccd
            fig = qaFig.QaFigure(size=figsize)
            fig.fig.subplots_adjust(left=0.125, bottom=0.125)
            ax_1 = fig.fig.add_subplot(121)
            ax_2 = fig.fig.add_subplot(122)
            clr = [black] * len(mag)
            clr = numpy.array(clr)
            clr[whereCut] = [red] * len(whereCut)

            for ax in [ax_1, ax_2]:
                ax.scatter(mag, diff, size, color=clr, label=ccd)
                ax.set_xlabel(tag1)
                
            ax_2.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]],
                      [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], '-k')
            ax_1.set_ylabel(tag)
            ax_1.set_xlim(xlim)
            ax_2.set_xlim(xlim2)
            ax_1.set_ylim(ylim)
            ax_2.set_ylim(ylim2)

            # move the y axis on right panel
            ax_2dummy = ax_2.twinx()
            ax_2dummy.set_ylim(ax_2.get_ylim())
            ax_2.set_yticks([])
            ax_2dummy.set_ylabel(tag)

            dmag = 0.1
            ddiff1 = 0.02
            ddiff2 = ddiff1*(ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]) # rescale for larger y range
            areaLabel = re.sub("\s+", "_", ccd)
            for j in xrange(len(mag)):
                info = "nolink:x:%.2f_y:%.2f" % (x[j], y[j])
                area = (mag[j]-dmag, diff[j]-ddiff1, mag[j]+dmag, diff[j]+ddiff1)
                fig.addMapArea(areaLabel, area, info, axes=ax_1)
                area = (mag[j]-dmag, diff[j]-ddiff2, mag[j]+dmag, diff[j]+ddiff2)
                fig.addMapArea(areaLabel, area, info, axes=ax_2)
                

            testSet.addFigure(fig, "diff_"+dtag+".png",
                              dtag+" vs. "+self.magType1 + ". Point used for statistics shown in red.",
                              areaLabel=areaLabel)


            # append values to arrays for a plot showing all data
            allMags = numpy.append(allMags, mag)
            allDiffs = numpy.append(allDiffs, diff)
            color = [sm.to_rgba(i)] * len(mag)
            allColor += color
            allLabels += [areaLabel] * len(mag)
            i += 1


        # dmag vs mag
        fig0 = qaFig.QaFigure(size=figsize)
        fig0.fig.subplots_adjust(left=0.125, bottom=0.125)
        ax0_1 = fig0.fig.add_subplot(121)
        ax0_2 = fig0.fig.add_subplot(122)

        
        ####################
        # data for all ccds
        if len(allMags) == 0:
            allMags = numpy.array([xmax])
            allDiffs = numpy.array([0.0])
            allColor = [black]
            allLabels = ["no_valid_data"]
            
        allColor = numpy.array(allColor)
        for ax in [ax0_1, ax0_2]:
            ax.scatter(allMags, allDiffs, size, color=allColor)
        ax0_1.set_xlim(xlim)
        ax0_2.set_xlim(xlim2)
        ax0_1.set_ylim(ylim)
        ax0_2.set_ylim(ylim2)

        dmag = 0.1
        ddiff1 = 0.02
        ddiff2 = ddiff1*(ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]) # rescale for larger y range
        for j in xrange(len(allMags)):
            area = (allMags[j]-dmag, allDiffs[j]-ddiff1, allMags[j]+dmag, allDiffs[j]+ddiff1)
            fig0.addMapArea(allLabels[j], area, "%.3f_%.3f"% (allMags[j], allDiffs[j]), axes=ax0_1)
            area = (allMags[j]-dmag, allDiffs[j]-ddiff2, allMags[j]+dmag, allDiffs[j]+ddiff2)
            fig0.addMapArea(allLabels[j], area, "%.3f_%.3f"% (allMags[j], allDiffs[j]), axes=ax0_2)



        # move the yaxis ticks/labels to the other side
        ax0_2dummy = ax0_2.twinx()
        ax0_2dummy.set_ylim(ax0_2.get_ylim())
        ax0_2.set_yticks([])

        ax0_2.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]],
                   [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], '-k')
        ax0_2.set_xlim(xlim2)
        ax0_2.set_ylim(ylim2)

        for ax in [ax0_1, ax0_2dummy]:
            ax.set_xlabel(tag1)
            ax.set_ylabel(tag)

            
        testSet.addFigure(fig0, "diff_"+dtag+".png", dtag+" vs. "+self.magType1, areaLabel="all")

            
