import re
import numpy as num

import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.meas.algorithms as measAlg
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData

from matplotlib.font_manager import FontProperties


class CompletenessQa2(qaAna.QaAnalysis):
    def __init__(self, completenessMagMin, completenessMagMax):
        qaAna.QaAnalysis.__init__(self)
        self.limits = [completenessMagMin, completenessMagMax]
        self.bins   = num.arange(14, 27, 0.5)

    def test(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)

        self.fluxType = fluxType

        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        self.matchListDict = data.getMatchListBySensor(dataId)
        self.ssDict        = data.getSourceSetBySensor(dataId)
        self.sroDict       = data.getRefObjectSetBySensor(dataId)

        self.matchStarObj   = raftCcdData.RaftCcdData(self.detector)
        self.matchGalObj    = raftCcdData.RaftCcdData(self.detector)
        #self.matchStarSrc   = raftCcdData.RaftCcdData(self.detector)
        #self.matchGalSrc    = raftCcdData.RaftCcdData(self.detector)
        self.unmatchCatStar = raftCcdData.RaftCcdData(self.detector)
        self.unmatchCatGal  = raftCcdData.RaftCcdData(self.detector)
        self.unmatchImage   = raftCcdData.RaftCcdData(self.detector)

        self.depth          = raftCcdData.RaftCcdData(self.detector)
        
        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER
        

        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            sroMagsStar = []
            sroMagsGxy = []
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
    
                        star = flags & measAlg.Flags.STAR
    
                        if num.isfinite(mrefmag):
                            if star > 0:
                                sroMagsStar.append(mrefmag)
                            else:
                                sroMagsGxy.append(mrefmag)
            
            self.matchStarObj.set(raftId, ccdId, num.array(sroMagsStar))
            self.matchGalObj.set(raftId, ccdId, num.array(sroMagsGxy))


            if False:
                # All stars/gals through source association; more complete than object assoc
                matchSrcSql  = 'select s.sourceId, sro.refObjectId, sro.%sMag, sro.isStar' % (filterName)
                matchSrcSql += ' from SimRefObject as sro, RefSrcMatch as rsm, Source as s, Science_Ccd_Exposure as sce'
                matchSrcSql += ' where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
                matchSrcSql += ' and (sce.visit = %s)' % (dataId['visit'])
                matchSrcSql += ' and (sce.raftName = "%s")' % (re.sub("R:", "", raftId))
                matchSrcSql += ' and (sce.ccdName = "%s")' % (ccdId[-3:])
                matchSrcSql += ' and (s.sourceId = rsm.sourceId)'
                matchSrcSql += ' and (rsm.refObjectId = sro.refObjectId)'
                matchSrcResults = data.dbInterface.execute(matchSrcSql)
                srcId   = num.array([x[0] for x in matchSrcResults])
                sroId   = num.array([x[1] for x in matchSrcResults])
                sroMags = num.array([x[2] for x in matchSrcResults])
                isStar  = num.array([x[3] for x in matchSrcResults])
                self.matchStarSrc.set(raftId, ccdId, sroMags[num.where(isStar == 1)])
                self.matchGalSrc.set(raftId, ccdId, sroMags[num.where(isStar == 0)])

            
            unmatchImage = []
            if self.ssDict.has_key(key):
                for s in self.ssDict[key]:
                    if not sids.has_key(s.getId()):
                        if self.fluxType == 'psf':
                            f = s.getPsfFlux()
                        else:
                            f = s.getApFlux()
                        if f <= 0.0:
                            continue
                        unmatchImage.append(-2.5*num.log10(f))
            uimgmag = num.array(unmatchImage)
            self.unmatchImage.set(raftId, ccdId, uimgmag)

            # Unmatched reference objects (inserted, but not found ... false negatives)
            unmatchCatStar = []
            unmatchCatGal = []
            if self.sroDict.has_key(key):
                for sro in self.sroDict[key]:
                    if not refIds.has_key(sro.refObjectId):
                        mag = sro.getMag(filterName)
                        if sro.isStar:
                            unmatchCatStar.append(mag)
                        else:
                            unmatchCatGal.append(mag)
            self.unmatchCatStar.set(raftId, ccdId, num.array(unmatchCatStar))
            self.unmatchCatGal.set(raftId, ccdId, num.array(unmatchCatGal))

            # decide which matchStar list to use: matchStarSrc or matchStarObj (Src better?)
            matchStarList = self.matchStarObj.get(raftId, ccdId)
            
            ###########################################################
            histStarSrc     = num.histogram(matchStarList, bins = self.bins)
            maxSrcIdx       = num.argsort(histStarSrc[0])[-1]
            histUnmatchStar = num.histogram(self.unmatchCatStar.get(raftId, ccdId), bins = self.bins)
            
            histRatio = num.zeros(len(histStarSrc[0]))
            w = num.where( histStarSrc[0] + histUnmatchStar[0] != 0)
            histRatio       = histStarSrc[0][w]/(1.0 * (histStarSrc[0][w]+histUnmatchStar[0][w]))


            badDepth = 0.0
            idxLim = None
            # Start at the bin with the most source counts
            for i in range(maxSrcIdx, len(histRatio)):
                if histRatio[i-2] > 0.5 and histRatio[i-1] > 0.5 and histRatio[i] <= 0.5:
                    idxLim = i
                    break
            if idxLim:
                if num.isnan(histStarSrc[1][idxLim-1]) or num.isnan(histStarSrc[1][idxLim]):
                    maxDepth = badDepth
                else:
                    maxDepth = 0.5 * (histStarSrc[1][idxLim-1] + histStarSrc[1][idxLim])
            else:
                maxDepth = badDepth
            self.depth.set(raftId, ccdId, maxDepth)

            areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
            label = "photometric depth "
            comment = "magnitude where completeness drops below 0.5"
            test = testCode.Test(label, maxDepth, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)



            
    def plot(self, data, dataId, showUndefined=False):
        testSet = self.getTestSet(data, dataId)

        # fpa figure
        depths = []
        depthFig = qaFig.FpaQaFigure(data.cameraInfo)
        for raft, ccdDict in depthFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.depth.get(raft, ccd) is None:
                    depth = self.depth.get(raft, ccd)
                    depths.append(depth)
                    depthFig.data[raft][ccd] = depth
                    depthFig.map[raft][ccd] = 'mag=%.2f'%(depth) 

        blue = '#0000ff'
        red  = '#ff0000'
        if len(depths) >= 2:
            vmin = max(num.min(depths), self.limits[0])
            vmax = min(num.max(depths), self.limits[1])
        else:
            vmin = self.limits[0]
            vmax = self.limits[1]
        depthFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=[vmin, vmax],
                            title="Photometric Depth", cmapOver=red, cmapUnder=blue, failLimits=self.limits)
        testSet.addFigure(depthFig, "completenessDepth.png", "Estimate of photometric depth", 
                          navMap=True)

        # Each CCD
        for raft, ccd in self.depth.raftCcdKeys():
            matchStarObjData   = self.matchStarObj.get(raft, ccd)
            matchGalObjData    = self.matchGalObj.get(raft, ccd)
            matchStarSrcData   = num.array([]) #self.matchStarSrc.get(raft, ccd)
            matchGalSrcData    = num.array([]) #self.matchGalSrc.get(raft, ccd)
            unmatchCatStarData = self.unmatchCatStar.get(raft, ccd)
            unmatchCatGalData  = self.unmatchCatGal.get(raft, ccd)
            unmatchImageData   = self.unmatchImage.get(raft, ccd)

            print "Plotting ", ccd

            # Just to get the histogram results
            fig = qaFig.QaFigure()
            sp1 = fig.fig.add_subplot(411)
            sp2 = fig.fig.add_subplot(412, sharex = sp1)
            sp3 = fig.fig.add_subplot(413, sharex = sp1)
            sp4 = fig.fig.add_subplot(414, sharex = sp1)

            if len(matchGalObjData):
                sp1.hist(matchGalObjData, facecolor='r', bins=self.bins, alpha=0.5, label='Galaxies', log=True)
            if len(matchStarObjData):
                sp1.hist(matchStarObjData, facecolor='g', bins=self.bins, alpha=0.5, label='Stars', log=True)

            if len(matchGalSrcData):
                sp2.hist(matchGalSrcData, facecolor='r', bins=self.bins, alpha=0.5, label='Galaxies', log=True)
            if len(matchStarSrcData):
                sp2.hist(matchStarSrcData, facecolor='g', bins=self.bins, alpha=0.5, label='Stars', log=True)

            if len(unmatchCatGalData):
                sp3.hist(unmatchCatGalData, facecolor='r', bins=self.bins, alpha=0.5, label='Galaxies', log=True)
            if len(unmatchCatStarData):
                sp3.hist(unmatchCatStarData, facecolor='g', bins=self.bins, alpha=0.5, label='Stars', log=True)

            if len(unmatchImageData):
                sp4.hist(unmatchImageData, facecolor='b', bins=self.bins, alpha=0.5, label='All', log=True)

            if len(matchGalObjData) or len(matchStarObjData):
                sp1.legend(numpoints=1, prop=FontProperties(size='x-small'), loc = 'upper left')
            if len(unmatchImageData):
                sp4.legend(numpoints=1, prop=FontProperties(size='x-small'), loc = 'upper left')
                
            qaFigUtils.qaSetp(sp1.get_xticklabels()+sp2.get_xticklabels()+sp3.get_xticklabels(), visible=False)

            sp1.text(0.5, 1.0, 'Match to Obj (G:%d S:%s)' % (len(matchGalObjData),
                                                             len(matchStarObjData)),
                     fontsize=8, transform=sp1.transAxes, ha='center')

            sp2.text(0.5, 1.0, 'Match to Src (G:%d S:%s)' % (len(matchGalSrcData),
                                                             len(matchStarSrcData)),
                     fontsize=8, transform=sp2.transAxes, ha='center')

            sp3.text(0.5, 1.0, 'Unmatched Cat (G:%d S:%s)' % (len(unmatchCatGalData),
                                                              len(unmatchCatStarData)),
                     fontsize=8, transform=sp3.transAxes, ha='center')

            sp4.text(0.5, 1.0, 'Unmatched Image (N:%d)' % (len(unmatchImageData)),
                     fontsize=8, transform=sp4.transAxes, ha='center')
            
            sp1.set_ylim(0.75, 999)
            sp2.set_ylim(0.75, 999)
            sp3.set_ylim(0.75, 999)
            sp4.set_ylim(0.75, 999)
            sp4.set_xlabel('Mag')
            qaFigUtils.qaSetp(sp1.get_yticklabels()+sp2.get_yticklabels()+sp3.get_yticklabels()+sp4.get_yticklabels(),
                              fontsize = 8)
            qaFigUtils.qaSetp(sp4.get_xticklabels(), fontsize = 8)
            
            label = data.cameraInfo.getDetectorName(raft, ccd)
            fig.fig.suptitle('%s' % (label), fontsize = 12)
            testSet.addFigure(fig, "completeness2.png", "Photometric detections "+label, areaLabel=label)
