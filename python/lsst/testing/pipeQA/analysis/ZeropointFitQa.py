import re
import numpy as num

import lsst.meas.algorithms         as measAlg
import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import QaAnalysisUtils as qaAnaUtil
import RaftCcdData as raftCcdData

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse

class ZeropointFitQa(qaAna.QaAnalysis):
    def __init__(self, medOffsetMin, medOffsetMax, figsize=(5.0,5.0), **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.figsize = figsize
        self.limits = [medOffsetMin, medOffsetMax]

        self.description = """
         For each CCD, the central panel shows the instrumental magnitude of
         matched stars and galaxies, plotted as a function of the catalog
         magnitude of the reference objects they were matched to.  The fitted
         zeropoint is shown as the dashed line.  The bottom panel shows a
         histogram of the "matched" and "orphan" source as a function of
         instrumental magnitude, while the left panel shows a histogram of the
         "matched" and "unmatched" entries in the reference catalog.  The top
         panel shows the scatter of the matched stars and galaxies around the
         zeropoint fit (similar to pipeQa.PhotCompareQaAnalysis.psf-cat) with
         the median offset of star from the zeropoint indicated with the dotted
         line.  The summary FPA figures show the median offset of stars from
         the zeropoint across the focal plane, as well as the fitted zeropoint.
        """

    def free(self):
        del self.detector
        del self.filter
        del self.calib
        del self.matchListDictSrc

        del self.orphan
        del self.matchedStar
        del self.undetectedStar
        del self.matchedGalaxy
        del self.undetectedGalaxy
        
        del self.zeroPoint
        del self.medOffset
        
    def test(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)
        self.fluxType         = fluxType
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        self.calib            = data.getCalibBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')

        # Ignore blends in this analysis
        self.orphan           = raftCcdData.RaftCcdVector(self.detector)
        self.matchedStar      = raftCcdData.RaftCcdVector(self.detector)
        self.undetectedStar   = raftCcdData.RaftCcdVector(self.detector)
        self.matchedGalaxy    = raftCcdData.RaftCcdVector(self.detector)
        self.undetectedGalaxy = raftCcdData.RaftCcdVector(self.detector)

        # Results
        self.zeroPoint        = raftCcdData.RaftCcdData(self.detector)
        self.medOffset        = raftCcdData.RaftCcdData(self.detector)

        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER

        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            fmag0 = self.calib[key].getFluxMag0()[0]
            if fmag0 <= 0.0:
                self.zeroPoint.set(raftId, ccdId, NaN)
                continue
            zpt = -2.5*num.log10(fmag0)
            self.zeroPoint.set(raftId, ccdId, zpt)

            if self.matchListDictSrc.has_key(key):
                # Matched
                mdict    = self.matchListDictSrc[key]['matched']
                stars    = []
                galaxies = []

                for m in mdict:
                    sref, s, dist = m
                    if fluxType == "psf":
                        fref  = sref.getPsfFlux()
                        f     = s.getPsfFlux()
                        ferr  = s.getPsfFluxErr()
                    else:
                        fref  = sref.getPsfFlux()
                        f     = s.getApFlux()
                        ferr  = s.getApFluxErr()

                    # un-calibrate the magnitudes
                    f *= fmag0
                    
                    flags = s.getFlagForDetection()
                    if (fref > 0.0 and f > 0.0 and not flags & badFlags):
                        mrefmag  = -2.5*num.log10(fref)
                        mimgmag  = -2.5*num.log10(f)
                        mimgmerr =  2.5 / num.log(10.0) * ferr / f
    
                        star = flags & measAlg.Flags.STAR
                        
                        if num.isfinite(mrefmag) and num.isfinite(mimgmag):
                            if star > 0:
                                stars.append((mrefmag, mimgmag, mimgmerr))
                            else:
                                galaxies.append((mrefmag, mimgmag, mimgmerr))
                self.matchedStar.set(raftId, ccdId, {"Refmag": num.array([x[0] for x in stars]),
                                                     "Imgmag": num.array([x[1] for x in stars]),
                                                     "Imgerr": num.array([x[2] for x in stars])})
                self.matchedGalaxy.set(raftId, ccdId, {"Refmag": num.array([x[0] for x in galaxies]),
                                                       "Imgmag": num.array([x[1] for x in galaxies]),
                                                       "Imgerr": num.array([x[2] for x in galaxies])})
            
                # Non-detections
                undetectedStars = []
                undetectedGalaxies = []
                for nondet in self.matchListDictSrc[key]['undetected']:
                    mag = nondet.getMag(filterName)
                    if nondet.getIsStar():
                        undetectedStars.append(mag)
                    else:
                        undetectedGalaxies.append(mag)
                self.undetectedStar.set(raftId, ccdId, num.array(undetectedStars))
                self.undetectedGalaxy.set(raftId, ccdId, num.array(undetectedGalaxies))

                # Orphans
                orphans = []
                for orphan in self.matchListDictSrc[key]['orphan']:
                    if self.fluxType == "psf":
                        f = orphan.getPsfFlux()
                    else:
                        f = orphan.getApFlux()
                    if f > 0.0:
                        # un-calibrate the magnitudes
                        f *= fmag0
                        
                        orphans.append(-2.5 * num.log10(f))
                        
                self.orphan.set(raftId, ccdId, num.array(orphans))

                # Metrics
                offset      = num.array(self.matchedStar.get(raftId, ccdId)["Imgmag"]) # make a copy
                offset     -= self.zeroPoint.get(raftId, ccdId)
                offset     -= self.matchedStar.get(raftId, ccdId)["Refmag"]
                med         = num.median(offset) 
                self.medOffset.set(raftId, ccdId, med)
                
                areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
                label = "median offset from zeropoint"
                comment = "Median offset of calibrated stellar magnitude to zeropoint fit"
                test = testCode.Test(label, med, self.limits, comment, areaLabel=areaLabel)
                testSet.addTest(test)

            
    def plot(self, data, dataId, showUndefined=False):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)

        # fpa figure
        zpts = []
        zptBase = "zeropoing"
        zptData, zptMap = testSet.unpickle(zptBase, default=[None, None])
        zptFig = qaFig.FpaQaFigure(data.cameraInfo, data=zptData, map=zptMap)

        offsetBase = "medZeropointOffset"
        offsetData, offsetMap = testSet.unpickle(offsetBase, default=[None, None])
        offsetFig = qaFig.FpaQaFigure(data.cameraInfo, data=offsetData, map=offsetMap)

        for raft, ccdDict in zptFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.zeroPoint.get(raft, ccd) is None:
                    zpt = self.zeroPoint.get(raft, ccd)
                    zpts.append(zpt)
                    zptFig.data[raft][ccd] = zpt
                    zptFig.map[raft][ccd] = 'zpt=%.2f' % (zpt)

                    offset = self.medOffset.get(raft, ccd)
                    offsetFig.data[raft][ccd] = offset
                    offsetFig.map[raft][ccd] = 'offset=%.2f' % (offset)
                    
        blue = '#0000ff'
        red = '#ff0000'
        
        zptFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=[num.min(zpts), num.max(zpts)],
                          title="Zeropoint", cmapOver=red, 
                          cmapUnder=blue)
        testSet.addFigure(zptFig, zptBase+".png", "Photometric zeropoint", 
                          navMap=True)
        testSet.pickle(zptBase, [zptFig.data, zptFig.map])
        
        offsetFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=self.limits,
                             title="Med offset from Zpt Fit", cmapOver=red, failLimits=self.limits,
                             cmapUnder=blue)
        testSet.addFigure(offsetFig, offsetBase + ".png", "Median offset from photometric zeropoint", 
                          navMap=True)
        testSet.pickle(offsetBase, [offsetFig.data, offsetFig.map])
        

        # Each CCD
        for raft, ccd in self.zeroPoint.raftCcdKeys():
            if self.zeroPoint.get(raft, ccd) == 0.0:
                continue
            
            print "Plotting", ccd

            # Just to get the histogram results
            fig = qaFig.QaFigure(size=self.figsize)

            legLines  = []
            legLabels = []
            
            axis = fig.fig.add_axes([0.225, 0.225, 0.675, 0.550])
    
            # Plot all matched galaxies
            mrefGmag  = self.matchedGalaxy.get(raft, ccd)["Refmag"]
            mimgGmag  = self.matchedGalaxy.get(raft, ccd)["Imgmag"]
            mimgGmerr = self.matchedGalaxy.get(raft, ccd)["Imgerr"]
            mimgGplot = axis.plot(mimgGmag, mrefGmag, '.', color='g', mfc='g', mec='g',
                                  alpha=0.5, zorder=10, label = 'Matched Galaxies')

            ################    --------------------------------  ##################
            # can't fig.savefig() with this Ellipse stuff ??? ... i love matplotlib
            if False:
                for i in range(len(mrefGmag)):
                    a = Ellipse(xy=num.array([mimgGmag[i], mrefGmag[i]]),
                                width=mimgGmerr[i]/2., height=mimgGmerr[i]/2.,
                                alpha=0.5, fill=True, ec='g', fc='g', zorder=10)
                    axis.add_artist(a)
            else:
                axis.plot(mimgGmag, mrefGmag, 'g.', zorder=10, alpha=0.5)
            #########################################################################

                
            legLines.append(mimgGplot[0])
            legLabels.append("Matched Galaxies")
                        
            # Plot all matched stars
            mrefSmag  = self.matchedStar.get(raft, ccd)["Refmag"]
            mimgSmag  = self.matchedStar.get(raft, ccd)["Imgmag"]
            mimgSmerr = self.matchedStar.get(raft, ccd)["Imgerr"]
            mimgSplot = axis.plot(mimgSmag, mrefSmag, '.', color='b', mfc='b', mec='b',
                                  alpha=0.5, zorder=12, label = 'Matched Stars')

            ################    --------------------------------  ##################
            # can't fig.savefig() with this Ellipse stuff ??? ... i love matplotlib
            if False:
                for i in range(len(mrefSmag)):
                    a = Ellipse(xy=num.array([mimgSmag[i], mrefSmag[i]]),
                                width=mimgSmerr[i]/2., height=mimgSmerr[i]/2.,
                                alpha=0.5, fill=True, ec='b', fc='b', zorder=12)
                    axis.add_artist(a)
            else:
                axis.plot(mimgSmag, mrefSmag, "b.", zorder=12, alpha=0.5)
            ########################################################################
                
            legLines.append(mimgSplot[0])
            legLabels.append("Matched Stars")
    
            if len(mrefGmag) == 0 and len(mrefSmag) == 0:
                xmin, xmax, ymin, ymax = -15, -8, 16, 28
            else:
                xmin, xmax, ymin, ymax = axis.axis()
    
            # Plot zpt
            xzpt = num.array((xmin, xmax))
            pzpt = axis.plot(xzpt, xzpt - self.zeroPoint.get(raft, ccd), 'k--', label = 'Zeropoint')
            legLines.append(pzpt)
            legLabels.append("Zeropoint")
    
            # Unmatched objects
            urefmag     = num.concatenate((self.undetectedStar.get(raft, ccd), self.undetectedGalaxy.get(raft, ccd)))
            uimgmag     = self.orphan.get(raft, ccd)

            # Unmatched & matched reference objects
            ax2  = fig.fig.add_axes([0.1,   0.225, 0.125, 0.550], sharey=axis)
            if len(urefmag) > 0:
                nu, bu, pu = ax2.hist(urefmag, bins=num.arange(ymin, ymax, 0.25),
                                      orientation='horizontal', facecolor = 'r',
                                      log = True, alpha = 0.5, zorder = 1)
                legLines.append(pu[0])
                legLabels.append("Unmatched Sources")
                
            if len(mrefGmag) > 0 and len(mrefSmag) > 0:
                ax2.hist(num.concatenate((mrefGmag,mrefSmag)), bins=num.arange(ymin, ymax, 0.25),
                         orientation='horizontal', log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
            elif len(mrefGmag):
                ax2.hist(mrefGmag, bins=num.arange(ymin, ymax, 0.25),
                         orientation='horizontal', log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
            elif len(mrefSmag) > 0:
                ax2.hist(mrefSmag, bins=num.arange(ymin, ymax, 0.25),
                         orientation='horizontal', log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
            ax2.set_xlabel('N', fontsize = 10)
            ax2.set_ylabel('Reference catalog mag', fontsize = 10)

    
            # Unmatched & matched stellar objects
            ax3  = fig.fig.add_axes([0.225, 0.1,   0.675, 0.125], sharex=axis)
            ax3.get_yaxis().set_ticks_position('right')
            ax3.get_yaxis().set_label_position('right')
            if len(mimgGmag) > 0 and len(mimgSmag) > 0:
                nm, bm, pm = ax3.hist(num.concatenate((mimgGmag,mimgSmag)), bins=num.arange(xmin, xmax, 0.25),
                                      log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
                legLines.append(pm[0])
                legLabels.append("Matched Sources")
            elif len(mimgGmag) > 0:
                nm, bm, pm = ax3.hist(mimgGmag, bins=num.arange(xmin, xmax, 0.25),
                                      log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
                legLines.append(pm[0])
                legLabels.append("Matched Sources")
            elif len(mimgSmag) > 0:
                nm, bm, pm = ax3.hist(mimgSmag, bins=num.arange(xmin, xmax, 0.25),
                                      log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
                legLines.append(pm[0])
                legLabels.append("Matched Sources")

            if len(uimgmag) > 0:
                ax3.hist(uimgmag, bins=num.arange(xmin, xmax, 0.25),
                         log = True, facecolor = 'r', alpha = 0.5, zorder = 1)
            ax3.set_xlabel('Image instrumental %s mag' % (self.fluxType), fontsize = 10)
            ax3.set_ylabel('N', rotation = 180, fontsize = 10)
    
            # Mag - Zpt
            ax4  = fig.fig.add_axes([0.225, 0.775, 0.675, 0.125], sharex=axis)
            if len(mimgSmag):
                mimgSeb = ax4.errorbar(mimgSmag, mimgSmag - self.zeroPoint.get(raft, ccd) - mrefSmag,
                                       yerr = mimgSmerr,
                                       fmt = 'bo', ms = 2, alpha = 0.25, capsize = 0, elinewidth = 0.5)
                mimgSeb[2][0].set_alpha(0.25) # alpha for error bars
    
            if len(mimgGmag):
                mimgGeb = ax4.errorbar(mimgGmag, mimgGmag - self.zeroPoint.get(raft, ccd) - mrefGmag,
                                       yerr = mimgGmerr,
                                       fmt = 'go', ms = 2, alpha = 0.25, capsize = 0, elinewidth = 0.5)
                mimgGeb[2][0].set_alpha(0.25) # alpha for error bars
    
            ax4.get_yaxis().set_ticks_position('right')
            ax4.get_yaxis().set_label_position('right')
            ax4.set_ylabel('Cal-Ref', fontsize = 10, rotation = 270)
            ax4.axhline(y = 0, c='k', linestyle='--', alpha = 0.25)
    
            # Cleaning up figure
            qaFigUtils.qaSetp(axis.get_xticklabels()+axis.get_yticklabels(), visible=False)
            qaFigUtils.qaSetp(ax2.get_xticklabels()+ax2.get_yticklabels(), fontsize = 8)
            qaFigUtils.qaSetp(ax3.get_xticklabels()+ax3.get_yticklabels(), fontsize = 8)
            qaFigUtils.qaSetp(ax4.get_xticklabels(), visible=False)
            qaFigUtils.qaSetp(ax4.get_yticklabels(), fontsize = 8)
    
            fig.fig.legend(legLines, legLabels,
                           numpoints=1, prop=FontProperties(size='x-small'), loc = 'center right')
            label = data.cameraInfo.getDetectorName(raft, ccd)
            fig.fig.suptitle('%s' % (label), fontsize = 12)
    
            numerator   = (mimgSmag - self.zeroPoint.get(raft, ccd) - mrefSmag)
            denominator = mimgSmerr
            med         = num.median(numerator) 
            ax4.axhline(y = med, c='k', linestyle=':', alpha = 0.5)
    
            # Final axis limits
            ax2.set_xlim(0.75, 999)
            ax3.set_ylim(0.75, 999)
            ax4.set_ylim(-0.24, 0.24)

            #axis.axis((xmax, xmin, ymax, ymin))
            axis.set_ylim(26, 14)
            axis.set_xlim(26 + self.zeroPoint.get(raft, ccd), 14 + self.zeroPoint.get(raft, ccd))
    
            label = data.cameraInfo.getDetectorName(raft, ccd)
            testSet.addFigure(fig, "zeropointFit.png", "Zeropoint fit "+label, areaLabel=label)
    
