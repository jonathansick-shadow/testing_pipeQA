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
    def __init__(self, medOffsetMin, medOffsetMax, figsize=(5.0,5.0)):
        qaAna.QaAnalysis.__init__(self)
        self.figsize = figsize
        self.limits = [medOffsetMin, medOffsetMax]


    def free(self):
        del self.detector
        del self.filter
        del self.calib
        del self.matchListDict
        del self.ssDict
        del self.sroDict
        
        del self.zeroPoint
        del self.medOffset
        
        del self.matchedGal
        del self.matchedStar
        del self.unmatchedRef
        del self.unmatchedImg


    def test(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)
        self.fluxType      = fluxType
        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        self.calib         = data.getCalibBySensor(dataId)
        self.matchListDict = data.getMatchListBySensor(dataId)
        self.ssDict        = data.getSourceSetBySensor(dataId)
        self.sroDict       = data.getRefObjectSetBySensor(dataId)

        
        self.zeroPoint     = raftCcdData.RaftCcdData(self.detector)
        self.medOffset     = raftCcdData.RaftCcdData(self.detector)

        self.matchedGal    = raftCcdData.RaftCcdData(self.detector)
        self.matchedStar   = raftCcdData.RaftCcdData(self.detector)
        self.unmatchedRef  = raftCcdData.RaftCcdData(self.detector)
        self.unmatchedImg  = raftCcdData.RaftCcdData(self.detector)

        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER

        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()
            
            zpt = -2.5*num.log10(self.calib[key].getFluxMag0())[0]
            self.zeroPoint.set(raftId, ccdId, zpt)

            mrefSmag, mimgSmag, mimgSmerr = [], [], []
            mrefGmag, mimgGmag, mimgGmerr = [], [], []

            # Unmatched detections (found by never inserted ... false positives)
            sids = {}
            refIds = {}

            if self.matchListDict.has_key(key):
                matchList = self.matchListDict[key]
                for m in matchList:
                    sref, s, dist = m

                    sids[s.getId()] = 1
                    refIds[sref.getId()] = 1
    
                    #oid = sref.getId()
    
                    if fluxType == "psf":
                        fref  = sref.getPsfFlux()
                        f     = s.getPsfFlux()
                        ferr  = s.getPsfFluxErr()
                    else:
                        fref  = sref.getPsfFlux()
                        f     = s.getApFlux()
                        ferr  = s.getApFluxErr()
                        
                    flags = s.getFlagForDetection()
    
                    if (fref > 0.0 and f > 0.0  and not flags & badFlags):
                        mrefmag  = -2.5*num.log10(fref)
                        mimgmag  = -2.5*num.log10(f)
                        mimgmerr =  2.5 / num.log(10.0) * ferr / f
    
                        star = flags & measAlg.Flags.STAR
    
                        if num.isfinite(mrefmag) and num.isfinite(mimgmag):
                            if star > 0:
                                mrefSmag.append(mrefmag)
                                mimgSmag.append(mimgmag)
                                mimgSmerr.append(mimgmerr)
                            else:
                                mrefGmag.append(mrefmag)
                                mimgGmag.append(mimgmag)
                                mimgGmerr.append(mimgmerr)

            mrefSmag = num.array(mrefSmag)
            mimgSmag = num.array(mimgSmag)
            mimgSmerr = num.array(mimgSmerr)
            self.matchedStar.set(raftId, ccdId, {"Refmag": mrefSmag,
                                                 "Imgmag": mimgSmag,
                                                 "Imgerr": mimgSmerr})
            mrefGmag = num.array(mrefGmag)
            mimgGmag = num.array(mimgGmag)
            mimgGmerr = num.array(mimgGmerr)
            self.matchedGal.set(raftId, ccdId, {"Refmag": mrefGmag,
                                                "Imgmag": mimgGmag,
                                                "Imgerr": mimgGmerr})

            unmatchedImg = []
            if self.ssDict.has_key(key):
                for s in self.ssDict[key]:
                    if not sids.has_key(s.getId()):
                        if self.fluxType == 'psf':
                            f = s.getPsfFlux()
                        else:
                            f = s.getApFlux()
                        unmatchedImg.append(-2.5*num.log10(f))
            uimgmag = num.array(unmatchedImg)
            self.unmatchedImg.set(raftId, ccdId, uimgmag)

            # Unmatched reference objects (inserted, but not found ... false negatives)
            unmatchedRef = []
            if self.sroDict.has_key(key):
                for sro in self.sroDict[key]:
                    if not refIds.has_key(sro.refObjectId):
                        mag = sro.getMag(filterName)
                        unmatchedRef.append(mag)
            urefmag = num.array(unmatchedRef)
            self.unmatchedRef.set(raftId, ccdId, urefmag)
                        

            # Metrics
            numerator   = mimgSmag - mrefSmag
            med         = num.median(numerator) 
            self.medOffset.set(raftId, ccdId, med)
            
            areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
            label = "median offset from zeropoint"
            comment = "Median offset of calibrated stellar magnitude to zeropoint fit"
            test = testCode.Test(label, med, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)

            
    def plot(self, data, dataId, showUndefined=False):
        testSet = self.getTestSet(data, dataId)

        # fpa figure
        zpts = []
        zptFig = qaFig.FpaQaFigure(data.cameraInfo)
        offsetFig = qaFig.FpaQaFigure(data.cameraInfo)
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
        testSet.addFigure(zptFig, "zeropoint.png", "Photometric zeropoint", 
                          navMap=True)
        
        offsetFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=self.limits,
                             title="Med offset from Zpt Fit", cmapOver=red, failLimits=self.limits,
                             cmapUnder=blue)
        testSet.addFigure(offsetFig, "medZeropointOffset.png", "Median offset from photometric zeropoint", 
                          navMap=True)

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
            mrefGmag  = self.matchedGal.get(raft, ccd)["Refmag"]
            mimgGmag  = self.matchedGal.get(raft, ccd)["Imgmag"]
            mimgGmerr = self.matchedGal.get(raft, ccd)["Imgerr"]
            mimgGplot = axis.plot(mimgGmag, mrefGmag, '.', color='g', mfc='g', mec='g',
                                  alpha=0.5, zorder=10, label = 'Matched Galaxies')

            ################    --------------------------------  ##################
            # can't fig.savefig() with this Ellipse stuff ??? ... i love matplotlib
            if False:
                for i in range(len(mrefGmag)):
                    print mimgGmag[i], mrefGmag[i]
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
            urefmag     = self.unmatchedRef.get(raft, ccd)
            uimgmag     = self.unmatchedImg.get(raft, ccd)

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
                mimgSeb = ax4.errorbar(mimgSmag, mimgSmag - mrefSmag,
                                       yerr = mimgSmerr,
                                       fmt = 'bo', ms = 2, alpha = 0.25, capsize = 0, elinewidth = 0.5)
                mimgSeb[2][0].set_alpha(0.25) # alpha for error bars
    
            if len(mimgGmag):
                mimgGeb = ax4.errorbar(mimgGmag, mimgGmag - mrefGmag,
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
    
            numerator   = (mimgSmag - mrefSmag)
            denominator = mimgSmerr
            med         = num.median(numerator) 
            ax4.axhline(y = med, c='k', linestyle=':', alpha = 0.5)
    
            # Final axis limits
            ax2.set_xlim(0.75, 999)
            ax3.set_ylim(0.75, 999)
            ax4.set_ylim(-0.24, 0.24)
            axis.axis((xmax, xmin, ymax, ymin))
    
            label = data.cameraInfo.getDetectorName(raft, ccd)
            testSet.addFigure(fig, "zeropointFit.png", "Zeropoint fit "+label, areaLabel=label)
    
