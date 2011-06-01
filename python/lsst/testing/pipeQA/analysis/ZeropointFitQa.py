import re
import numpy as num

import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse

class ZeropointFitQa(qaAna.QaAnalysis):
    def __init__(self, medOffsetMin, medOffsetMax, figsize=(5.0,5.0)):
        qaAna.QaAnalysis.__init__(self)
        self.figsize = figsize
        self.limits = [medOffsetMin, medOffsetMax]

    def test(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)
        self.fluxType      = fluxType
        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        
        self.zeroPoint     = raftCcdData.RaftCcdData(self.detector)
        self.medOffset     = raftCcdData.RaftCcdData(self.detector)

        self.matchedGal    = raftCcdData.RaftCcdData(self.detector)
        self.matchedStar   = raftCcdData.RaftCcdData(self.detector)
        self.unmatchedRef  = raftCcdData.RaftCcdData(self.detector)
        self.unmatchedImg  = raftCcdData.RaftCcdData(self.detector)
            
        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()
            
            print "Running", ccdId

            # Get zeropoint
            zptsql      = 'select scienceCcdExposureId,fluxMag0 from Science_Ccd_Exposure'
            zptsql     += ' where (visit = %s)' % (dataId['visit'])
            zptsql     += ' and (raftName = "%s")' % (re.sub("R:", "", raftId))
            zptsql     += ' and (ccdName = "%s")' % (ccdId[-3:])
            zptresults  = data.dbInterface.execute(zptsql)
            if len(zptresults) != 1:
                # throw exception or something
                return
            sceId, fluxMag0 = zptresults[0]
            zpt = -2.5 * num.log10(fluxMag0)
            self.zeroPoint.set(raftId, ccdId, zpt)
            
            # Select all reference stars within this field
            catsql      = 'select sro.refObjectId, sro.isStar'
            catsql     += ' from SimRefObject as sro, Science_Ccd_Exposure as sce'
            catsql     += ' where (sce.scienceCcdExposureId = %d)' % (sceId)
            catsql     += ' and qserv_ptInSphPoly(sro.ra, sro.decl,'
            catsql     += ' concat_ws(" ", sce.llcRa, sce.llcDecl, sce.lrcRa, sce.lrcDecl, '
            catsql     += ' sce.urcRa, sce.urcDecl, sce.ulcRa, sce.ulcDecl))'
            catresults  = data.dbInterface.execute(catsql)

            # Sort by star,gal
            refAll = {'s': [], 'g': []}
            for result in catresults:
                oid, isStar = result
                if isStar: refAll['s'].append(oid)
                else: refAll['g'].append(oid)
    
            # Select all matched galaxies
            if len(refAll['g']) == 0:
                self.matchedGal.set(raftId, ccdId, {"Refmag": num.array(()),
                                                    "Imgmag": num.array(()),
                                                    "Imgerr": num.array(())})
            else:
                mrefGsql  = 'select sro.%sMag,' % (filterName)
                mrefGsql += ' s.%sFlux, s.%sFluxSigma' % (fluxType, fluxType)
                mrefGsql += ' from SimRefObject as sro, RefObjMatch as rom, Source as s'
                mrefGsql += ' where (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
                mrefGsql += ' and (s.scienceCcdExposureId = %d)' % (sceId)
                mrefGsql += ' and (s.objectID is not NULL)'
                mrefGsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['g'])))
                mrefGresults = data.dbInterface.execute(mrefGsql)
                mrefGmag  = num.array([x[0] for x in mrefGresults])
                mimgGflu  = num.array([x[1] for x in mrefGresults])
                mimgGferr = num.array([x[2] for x in mrefGresults])
                mimgGmag  = -2.5 * num.log10(mimgGflu)
                mimgGmerr =  2.5 / num.log(10.0) * mimgGferr / mimgGflu
                self.matchedGal.set(raftId, ccdId, {"Refmag": mrefGmag,
                                                    "Imgmag": mimgGmag,
                                                    "Imgerr": mimgGmerr})
            
            # Select all matched stars
            mrefSsql  = 'select sro.%sMag,' % (filterName)
            mrefSsql += ' s.%sFlux, s.%sFluxSigma' % (fluxType, fluxType)
            mrefSsql += ' from SimRefObject as sro, RefObjMatch as rom, Source as s'
            mrefSsql += ' where (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
            mrefSsql += ' and (s.scienceCcdExposureId = %d)' % (sceId)
            mrefSsql += ' and (s.objectID is not NULL)'
            mrefSsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['s'])))
            mrefSresults  = data.dbInterface.execute(mrefSsql)
            mrefSmag  = num.array([x[0] for x in mrefSresults])
            mimgSflu  = num.array([x[1] for x in mrefSresults])
            mimgSferr = num.array([x[2] for x in mrefSresults])
            mimgSmag  = -2.5 * num.log10(mimgSflu)
            mimgSmerr =  2.5 / num.log(10.0) * mimgSferr / mimgSflu
            self.matchedStar.set(raftId, ccdId, {"Refmag": mrefSmag,
                                                 "Imgmag": mimgSmag,
                                                 "Imgerr": mimgSmerr})
    
            ####
    
            # Umatched reference catalog objects
            urefsql  = 'select sro.%sMag from SimRefObject as sro, RefObjMatch as rom' % (filterName)
            urefsql += ' where (sro.refObjectId = rom.refObjectId)'
            urefsql += ' and (rom.objectId is NULL)'
            urefsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['g'] + refAll['s'])))
            urefresults = data.dbInterface.execute(urefsql)
            urefmag     = num.array([x[0] for x in urefresults])
            self.unmatchedRef.set(raftId, ccdId, urefmag)
    
            # Unmatched detections
            uimgsql  = 'select %sFlux from Source ' % (fluxType)
            uimgsql += ' where (scienceCcdExposureId = %d)' % (sceId)
            uimgsql += ' and (objectId is NULL)'
            uimgresults  = data.dbInterface.execute(uimgsql)
            uimgmag      = -2.5 * num.log10( num.array([x[0] for x in uimgresults]) )
            self.unmatchedImg.set(raftId, ccdId, uimgmag)

            # Metrics
            numerator   = (mimgSmag - zpt) - mrefSmag
            denominator = mimgSmerr
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
            for i in range(len(mrefGmag)):
                a = Ellipse(xy=num.array([mimgGmag[i], mrefGmag[i]]),
                            width=mimgGmerr[i]/2., height=mimgGmerr[i]/2.,
                            alpha=0.5, fill=True, ec='g', fc='g', zorder=10)
                axis.add_artist(a)            
            legLines.append(mimgGplot[0])
            legLabels.append("Matched Galaxies")
            
            # Plot all matched stars
            mrefSmag  = self.matchedStar.get(raft, ccd)["Refmag"]
            mimgSmag  = self.matchedStar.get(raft, ccd)["Imgmag"]
            mimgSmerr = self.matchedStar.get(raft, ccd)["Imgerr"]
            mimgSplot = axis.plot(mimgSmag, mrefSmag, '.', color='b', mfc='b', mec='b',
                                  alpha=0.5, zorder=12, label = 'Matched Stars')
            for i in range(len(mrefSmag)):
                a = Ellipse(xy=num.array([mimgSmag[i], mrefSmag[i]]),
                            width=mimgSmerr[i]/2., height=mimgSmerr[i]/2.,
                            alpha=0.5, fill=True, ec='b', fc='b', zorder=12)
                axis.add_artist(a)            
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
                                      orientation='horizontal', facecolor = 'r', log = True, alpha = 0.5, zorder = 1)
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
                
            ax3.hist(uimgmag, bins=num.arange(xmin, xmax, 0.25),
                     log = True, facecolor = 'r', alpha = 0.5, zorder = 1)
            ax3.set_xlabel('Image instrumental %s mag' % (self.fluxType), fontsize = 10)
            ax3.set_ylabel('N', rotation = 180, fontsize = 10)
    
            # Mag - Zpt
            ax4  = fig.fig.add_axes([0.225, 0.775, 0.675, 0.125], sharex=axis)
            if len(mimgSmag):
                mimgSeb = ax4.errorbar(mimgSmag, (mimgSmag - self.zeroPoint.get(raft, ccd)) - mrefSmag,
                                       yerr = mimgSmerr,
                                       fmt = 'bo', ms = 2, alpha = 0.25, capsize = 0, elinewidth = 0.5)
                mimgSeb[2][0].set_alpha(0.25) # alpha for error bars
    
            if len(mimgGmag):
                mimgGeb = ax4.errorbar(mimgGmag, (mimgGmag - self.zeroPoint.get(raft, ccd)) - mrefGmag,
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
    
            numerator   = (mimgSmag - self.zeroPoint.get(raft, ccd)) - mrefSmag
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
    
