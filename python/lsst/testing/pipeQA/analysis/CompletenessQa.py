import re
import numpy as num

import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.meas.algorithms as measAlg
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData

from matplotlib.font_manager import FontProperties


class CompletenessQa(qaAna.QaAnalysis):
    def __init__(self, completenessMagMin, completenessMagMax):
        qaAna.QaAnalysis.__init__(self)
        self.limits = [completenessMagMin, completenessMagMax]
        self.bins   = num.arange(14, 27, 0.5)

    def free(self):
        del self.detector
        del self.filter
        del self.matchListDictObj
        del self.matchListDictSrc
        del self.ssDict
        del self.sroDict
        del self.matchStarObj
        del self.matchGalObj
        del self.matchStarSrc
        del self.matchGalSrc
        del self.unmatchCatStar
        del self.unmatchCatGal
        del self.unmatchImage
        del self.depth


    def test(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)

        self.fluxType = fluxType

        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        self.matchListDictObj = data.getMatchListBySensor(dataId, useRef='obj')
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.ssDict        = data.getSourceSetBySensor(dataId)
        self.sroDict       = data.getRefObjectSetBySensor(dataId)

        self.matchStarObj   = raftCcdData.RaftCcdData(self.detector)
        self.matchGalObj    = raftCcdData.RaftCcdData(self.detector)
        self.matchStarSrc   = raftCcdData.RaftCcdData(self.detector)
        self.matchGalSrc    = raftCcdData.RaftCcdData(self.detector)
        self.unmatchCatStar = raftCcdData.RaftCcdData(self.detector)
        self.unmatchCatGal  = raftCcdData.RaftCcdData(self.detector)
        self.unmatchImage   = raftCcdData.RaftCcdData(self.detector)

        self.depth          = raftCcdData.RaftCcdData(self.detector)
        
        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER
        

        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            matchStuff = [
                [self.matchListDictObj, self.matchStarObj, self.matchGalObj],
                [self.matchListDictSrc, self.matchStarSrc, self.matchGalSrc],
                ]
            multiples = [0, 0, 0, 0, 0, 0]
            for ms in matchStuff:
                matchListDict, matchStar, matchGal = ms

                sroMagsStar = []
                sroMagsGxy = []
                # Unmatched detections (found by never inserted ... false positives)
                matchSourceIds = {}
                matchSourceRefIds = {}
                if matchListDict.has_key(key):
                    matchList = matchListDict[key]
                    for m in matchList:
                        sref, s, dist = m

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

                        if (fref > 0.0 and f > 0.0 and not flags & badFlags):
                            mrefmag  = -2.5*num.log10(fref)
                            if not matchSourceIds.has_key(s.getId()):
                                matchSourceIds[s.getId()] = 0
                            matchSourceIds[s.getId()] += 1
                            if not matchSourceRefIds.has_key(sref.getId()):
                                matchSourceRefIds[sref.getId()] = 0
                            matchSourceRefIds[sref.getId()] += 1

                            star = flags & measAlg.Flags.STAR

                            if num.isfinite(mrefmag):
                                if star > 0:
                                    sroMagsStar.append(mrefmag)
                                else:
                                    sroMagsGxy.append(mrefmag)

                    matchStar.set(raftId, ccdId, num.array(sroMagsStar))
                    matchGal.set(raftId, ccdId, num.array(sroMagsGxy))


            
            unmatchImage = []
            q = 0
            if self.ssDict.has_key(key):
                for s in self.ssDict[key]:
                    if not matchSourceIds.has_key(s.getId()):
                        if self.fluxType == 'psf':
                            f = s.getPsfFlux()
                        else:
                            f = s.getApFlux()
                        if f <= 0.0:
                            continue
                        unmatchImage.append(-2.5*num.log10(f))
                    else:
                        q += 1

            uimgmag = num.array(unmatchImage)
            
            self.unmatchImage.set(raftId, ccdId, uimgmag)

            # Unmatched reference objects (inserted, but not found ... false negatives)
            unmatchCatStar = []
            unmatchCatGal = []
            if self.sroDict.has_key(key):
                for sro in self.sroDict[key]:
                    if not matchSourceRefIds.has_key(sro.refObjectId):
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
            
            magbins   = 0.5 * (histStarSrc[1][1:] + histStarSrc[1][:-1])
            histRatio = num.zeros(len(histStarSrc[0]))

            w         = num.where((histStarSrc[0] + histUnmatchStar[0]) != 0)
            magbins   = magbins[w]
            histRatio = histStarSrc[0][w]/(1.0 * (histStarSrc[0][w]+histUnmatchStar[0][w]))

            badDepth = 0.0
            idxLim = None
            # Start at the bin with the most source counts
            for i in range(maxSrcIdx, len(histRatio)):
                # Too many failures if there are fluctuations around the 50% region
                # if histRatio[i-2] > 0.5 and histRatio[i-1] > 0.5 and histRatio[i] <= 0.5:
                if histRatio[i-1] > 0.5 and histRatio[i] <= 0.5:
                    idxLim = i
            
            if idxLim:
                if num.isnan(histStarSrc[1][idxLim-1]) or num.isnan(histStarSrc[1][idxLim]):
                    maxDepth = badDepth
                else:
                    maxDepth = magbins[idxLim]
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

        if vmax <= vmin:
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
            matchStarSrcData   = self.matchStarSrc.get(raft, ccd)
            matchGalSrcData    = self.matchGalSrc.get(raft, ccd)
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

            nmss = []
            if len(matchGalSrcData):
                sp2.hist(matchGalSrcData, facecolor='r', bins=self.bins, alpha=0.5, label='Galaxies', log=True)
            if len(matchStarSrcData):
                nmss, bmss, pmss = sp2.hist(matchStarSrcData, facecolor='g', bins=self.bins, alpha=0.5, label='Stars', log=True)

            nuss = []
            w = num.isfinite(unmatchCatGalData)
            if len(unmatchCatGalData):
                sp3.hist(unmatchCatGalData[w], facecolor='r', bins=self.bins, alpha=0.5, label='Galaxies', log=True)
            w = num.isfinite(unmatchCatStarData)
            if len(unmatchCatStarData):
                nuss, buss, puss = sp3.hist(unmatchCatStarData[w], facecolor='g', 
                                            bins=self.bins, alpha=0.5, label='Stars', log=True)

            noss = []
            if len(unmatchImageData):
                noss, boss, poss = sp4.hist(unmatchImageData, facecolor='b', bins=self.bins, alpha=0.5, label='All', log=True)

            # completeness lines
            nmss = num.array(nmss)
            nuss = num.array(nuss)
            noss = num.array(noss)
            magbins = 0.5 * (self.bins[1:] + self.bins[:-1])
            if len(nmss) and len(nuss):
                idx = num.where((nmss + nuss) != 0)
                fracDet = 1.0 * nmss[idx] / (nmss[idx] + nuss[idx])
                sp2x2   = sp2.twinx()
                sp2x2.plot(magbins[idx], fracDet)
                sp2x2.set_ylabel('Match/Tot', fontsize = 8)
                qaFigUtils.qaSetp(sp2x2.get_xticklabels(), visible=False)
                qaFigUtils.qaSetp(sp2x2.get_yticklabels(), fontsize = 6)

                fracuDet = 1.0 * nuss[idx] / (nmss[idx] + nuss[idx])
                sp3x2    = sp3.twinx()
                sp3x2.plot(magbins[idx], fracuDet)
                sp3x2.set_ylabel('UnDet/Tot', fontsize = 8)
                qaFigUtils.qaSetp(sp3x2.get_xticklabels(), visible=False)
                qaFigUtils.qaSetp(sp3x2.get_yticklabels(), fontsize = 6)

            if len(noss) and len(nmss):
                idx = num.where((noss + nmss) != 0)
                fracOrph = 1.0 * noss[idx] / (noss[idx] + nmss[idx])
                sp4x2   = sp4.twinx()
                sp4x2.plot(magbins[idx], fracOrph)
                sp4x2.set_ylabel('Orph/Det', fontsize = 8)
                qaFigUtils.qaSetp(sp4x2.get_xticklabels(), visible=False)
                qaFigUtils.qaSetp(sp4x2.get_yticklabels(), fontsize = 6)
                
                
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
