import re
import numpy as num

import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.meas.algorithms as measAlg
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData

from matplotlib.font_manager import FontProperties

# Until we can make it more robust
hasMinuit = False
try:
    import minuit2
except:
    hasMinuit = False
    

class CompletenessQa(qaAna.QaAnalysis):
    def __init__(self, completenessMagMin, completenessMagMax, **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.limits = [completenessMagMin, completenessMagMax]
        self.bins   = num.arange(14, 27, 0.5)

    def free(self):
        del self.detector
        del self.filter
        del self.matchListDictSrc

        del self.orphan
        del self.matchedStar
        del self.blendedStar
        del self.undetectedStar
        del self.matchedGalaxy
        del self.blendedGalaxy
        del self.undetectedGalaxy
        del self.depth
        if hasMinuit:
            del self.fit

    def limitingMag(self, raftId, ccdId):
        if hasMinuit:
            try:
                return self.limitingMagMinuit(raftId, ccdId)
            except:
                pass

        matchedStar     = num.array(self.matchedStar.get(raftId, ccdId))
        blendedStar     = num.array(self.blendedStar.get(raftId, ccdId))
        undetectedStar  = num.array(self.undetectedStar.get(raftId, ccdId))

        allStars        = num.concatenate((matchedStar, blendedStar, undetectedStar))
        foundStars      = num.concatenate((matchedStar, blendedStar))
        histAll         = num.histogram(allStars, bins=self.bins)
        histFound       = num.histogram(foundStars, bins=self.bins)

        magbins = 0.5 * (histAll[1][1:] + histAll[1][:-1])
        w       = num.where(histAll[0] != 0)
        x       = magbins[w]
        n       = 1.0 * histFound[0][w]
        d       = 1.0 * histAll[0][w]
        y       = n / d  

        for i in num.arange(len(y) - 1, 1, -1):
            if y[i] <= 0.5 and y[i-1] > 0.5:
                return (0.5 - y[i-1]) / (y[i] - y[i-1]) * (x[i] - x[i-1]) + x[i-1]
        return 0.0
        
    def limitingMagMinuit(self, raftId, ccdId):
        # Model is of the form:
        # 0.5 + -1.0 / num.pi * num.arctan(A * x + B)
        
        import minuit2
        matchStarList   = self.matchStarSrc.get(raftId, ccdId)
        unmatchStarList = self.unmatchCatStar.get(raftId, ccdId)
        
        histStarSrc     = num.histogram(matchStarList, bins = self.bins)
        maxSrcIdx       = num.argsort(histStarSrc[0])[-1]
        histUnmatchStar = num.histogram(unmatchStarList, bins = self.bins)

        magbins   = 0.5 * (histStarSrc[1][1:] + histStarSrc[1][:-1])
        histRatio = num.zeros(len(histStarSrc[0]))

        w     = num.where((histStarSrc[0] + histUnmatchStar[0]) != 0)
        x     = magbins[w]

        # approximate
        d     = 1.0 * histStarSrc[0][w]
        u     = 1.0 * histUnmatchStar[0][w]
        dd    = num.sqrt(d)
        n     = d + u
        y     = d / n
        dy    = dd / n

        idx = num.where(dy != 0)
        x   = x[idx]
        y   = y[idx]
        dy  = dy[idx]
              
        def fcn(A, B):
            model  = 0.5 + -1.0 / num.pi * num.arctan(A * x + B)
            chi    = (model - y) / dy
            return num.sum(chi**2)
        
        m = minuit2.Minuit2(fcn)
        m.values['A'] = 1
        m.values['B'] = -10
        m.migrad()

        mx = num.arange(min(x), max(x), 0.1)
        my = 0.5 + -1.0 / num.pi * num.arctan(m.values['A'] * mx + m.values['B'])
        mindx = num.argsort((num.abs(my-0.5)))[0]

        self.fit.set(raftId, ccdId, [m.values['A'], m.values['B']])
        return mx[mindx]

    def test(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)
        self.fluxType = fluxType
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')

        self.orphan           = raftCcdData.RaftCcdVector(self.detector)
        self.matchedStar      = raftCcdData.RaftCcdVector(self.detector)
        self.blendedStar      = raftCcdData.RaftCcdVector(self.detector)
        self.undetectedStar   = raftCcdData.RaftCcdVector(self.detector)
        self.matchedGalaxy    = raftCcdData.RaftCcdVector(self.detector)
        self.blendedGalaxy    = raftCcdData.RaftCcdVector(self.detector)
        self.undetectedGalaxy = raftCcdData.RaftCcdVector(self.detector)
        self.depth            = raftCcdData.RaftCcdVector(self.detector)

        if hasMinuit:
            self.fit = raftCcdData.RaftCcdData(self.detector, initValue=[0.0, 0.0]) 
        
        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER
        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            if self.matchListDictSrc.has_key(key):
                # Detections
                matchSet = [
                    [self.matchListDictSrc[key]['matched'], self.matchedStar, self.matchedGalaxy],
                    [self.matchListDictSrc[key]['blended'], self.blendedStar, self.blendedGalaxy]
                    ]
                for mset in matchSet:
                    mdict, starvec, galvec = mset
                    
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
                        flags = s.getFlagForDetection()
                        if (fref > 0.0 and f > 0.0 and not flags & badFlags):
                            # Use known catalog mag
                            mrefmag  = -2.5*num.log10(fref)
                            star = flags & measAlg.Flags.STAR
                            if num.isfinite(mrefmag):
                                if star > 0:
                                    stars.append(mrefmag)
                                else:
                                    galaxies.append(mrefmag)
                    starvec.set(raftId, ccdId, num.array(stars))
                    galvec.set(raftId, ccdId, num.array(galaxies))
    
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
                        orphans.append(-2.5 * num.log10(f))
                self.orphan.set(raftId, ccdId, num.array(orphans))

                ############ Calculate limiting mag
                
                maxDepth = self.limitingMag(raftId, ccdId)
                self.depth.set(raftId, ccdId, maxDepth)
                
                areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
                label = "photometric depth "
                comment = "magnitude where star completeness drops below 0.5"
                test = testCode.Test(label, maxDepth, self.limits, comment, areaLabel=areaLabel)
                testSet.addTest(test)

    def plot(self, data, dataId, showUndefined = False):
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
        
        if False:
            if len(depths) >= 2:
                vmin = max(num.min(depths), self.limits[0])
                vmax = min(num.max(depths), self.limits[1])
            else:
                vmin = self.limits[0]
                vmax = self.limits[1]

            if vmax <= vmin:
                vmin = self.limits[0]
                vmax = self.limits[1]
        else:
            vmin, vmax = 1.0*self.limits[0], 1.0*self.limits[1]

        depthFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[vmin, vmax],
                            title="Photometric Depth", cmapOver=red, cmapUnder=blue,
                            failLimits=self.limits)
        testSet.addFigure(depthFig, "completenessDepth.png", "Estimate of photometric depth", 
                          navMap=True)

        # Each CCD
        for raft, ccd in self.depth.raftCcdKeys():
            orphan           = self.orphan.get(raft, ccd)
            matchedStar      = self.matchedStar.get(raft, ccd)
            blendedStar      = self.blendedStar.get(raft, ccd)
            undetectedStar   = self.undetectedStar.get(raft, ccd)
            matchedGalaxy    = self.matchedGalaxy.get(raft, ccd)
            blendedGalaxy    = self.blendedGalaxy.get(raft, ccd)
            undetectedGalaxy = self.undetectedGalaxy.get(raft, ccd)

            print "Plotting ", ccd

            fig = qaFig.QaFigure()
            sp1 = fig.fig.add_subplot(211)
            sp2 = fig.fig.add_subplot(212, sharex = sp1)

            # Stacked histogram
            orphanHist         = num.histogram(orphan, bins=self.bins)
            matchedStarHist    = num.histogram(matchedStar, bins=self.bins)
            blendedStarHist    = num.histogram(blendedStar, bins=self.bins)
            undetectedStarHist = num.histogram(undetectedStar, bins=self.bins)
            # For bar, you send the coordinate of the left corner of the bar
            barbins    = orphanHist[1][:-1]
            width      = 1.0 * (orphanHist[1][1] - orphanHist[1][0])
            orphanBar  = sp1.bar(barbins, orphanHist[0], width=width, color='r', alpha = 0.5, label = 'Orphan', capsize = 1)
            bottom     = orphanHist[0]
            matchedBar = sp1.bar(barbins, matchedStarHist[0], width=width, color='g', alpha=0.5, label='Matched',
                                 bottom=bottom, capsize=1)
            bottom    += matchedStarHist[0]
            blendedBar = sp1.bar(barbins, blendedStarHist[0], width=width, color='cyan', alpha=0.5, label='Blended',
                                 bottom=bottom, capsize=1)
            bottom    += blendedStarHist[0]
            unmatBar   = sp1.bar(barbins, undetectedStarHist[0], width=width, color='b', alpha=0.5, label='Unmatched',
                                 bottom=bottom, capsize=1)

            sp1x2           = sp1.twinx()
            allStars        = num.concatenate((matchedStar, blendedStar, undetectedStar))
            foundStars      = num.concatenate((matchedStar, blendedStar))
            histAll         = num.histogram(allStars, bins=self.bins)
            histFound       = num.histogram(foundStars, bins=self.bins)
            
            magbins = 0.5 * (histAll[1][1:] + histAll[1][:-1])
            w       = num.where(histAll[0] != 0)
            x       = magbins[w]
            n       = 1.0 * histFound[0][w]
            d       = 1.0 * histAll[0][w]
            y       = n / d  
            sp1x2.plot(x, y)
            sp1x2.set_ylabel('(Match+Blend)/Tot', fontsize=8)
            sp1x2.axhline(y = 0.5, c='k', linestyle='-', alpha = 0.25)
            sp1x2.axvline(x = self.depth.get(raft, ccd), c='k', linestyle='-', alpha = 0.25)
            qaFigUtils.qaSetp(sp1x2.get_xticklabels(), visible=False)
            qaFigUtils.qaSetp(sp1x2.get_yticklabels(), fontsize = 6)

            sp1.set_ylabel('N Stars', fontsize=10)
            qaFigUtils.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize = 8)

            ##############

            orphanHist        = num.histogram(orphan, bins=self.bins)
            matchedGalHist    = num.histogram(matchedGalaxy, bins=self.bins)
            blendedGalHist    = num.histogram(blendedGalaxy, bins=self.bins)
            undetectedGalHist = num.histogram(undetectedGalaxy, bins=self.bins)
            orphanBar  = sp2.bar(barbins, orphanHist[0], width=width, color='r', alpha = 0.5, label = 'Orphan', capsize = 1, log=False)
            bottom     = orphanHist[0]
            matchedBar = sp2.bar(barbins, matchedGalHist[0], width=width, color='g', alpha=0.5, label='Matched',
                                 bottom=bottom, capsize=1, log=False)
            bottom    += matchedGalHist[0]
            blendedBar = sp2.bar(barbins, blendedGalHist[0], width=width, color='cyan', alpha=0.5, label='Blended',
                                 bottom=bottom, capsize=1, log=False)
            bottom    += blendedGalHist[0]
            unmatBar   = sp2.bar(barbins, undetectedGalHist[0], width=width, color='b', alpha=0.5, label='Unmatched',
                                 bottom=bottom, capsize=1, log=False)
            
            sp2.set_xlabel('Mag', fontsize=10)
            sp2.set_ylabel('N Gals', fontsize=10)
            qaFigUtils.qaSetp(sp2.get_xticklabels()+sp2.get_yticklabels(), fontsize = 8)
            sp2.legend(numpoints = 1, prop=FontProperties(size='x-small'), loc = 'upper left')
            #sp2.set_ylim(0.75, 999)
            #sp2.semilogy()

            sp1.set_xlim(14, 26)
            label = data.cameraInfo.getDetectorName(raft, ccd)
            fig.fig.suptitle('%s Stacked histogram' % (label), fontsize = 12)
            testSet.addFigure(fig, "completeness.png", "Photometric detections "+label, areaLabel=label)
            


###########
###########
###########
###########
###########
###########

            
    def testOld(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)

        self.fluxType = fluxType

        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        self.matchListDictObj = data.getMatchListBySensor(dataId, useRef='obj')
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.ssDict        = data.getSourceSetBySensor(dataId)
        self.sroDict       = data.getRefObjectSetBySensor(dataId)

        self.matchStarObj   = raftCcdData.RaftCcdVector(self.detector)
        self.matchGalObj    = raftCcdData.RaftCcdVector(self.detector)
        self.matchStarSrc   = raftCcdData.RaftCcdVector(self.detector)
        self.matchGalSrc    = raftCcdData.RaftCcdVector(self.detector)
        self.unmatchCatStar = raftCcdData.RaftCcdVector(self.detector)
        self.unmatchCatGal  = raftCcdData.RaftCcdVector(self.detector)
        self.unmatchImage   = raftCcdData.RaftCcdVector(self.detector)

        self.depth          = raftCcdData.RaftCcdData(self.detector)

        if hasMinuit:
            self.fit = raftCcdData.RaftCcdData(self.detector, initValue=[0.0, 0.0]) 
        
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
                    if not matchSourceRefIds.has_key(sro.getId()):
                        mag = sro.getMag(filterName)
                        if sro.getIsStar():
                            unmatchCatStar.append(mag)
                        else:
                            unmatchCatGal.append(mag)
            self.unmatchCatStar.set(raftId, ccdId, num.array(unmatchCatStar))
            self.unmatchCatGal.set(raftId, ccdId, num.array(unmatchCatGal))

            maxDepth = self.limitingMag(raftId, ccdId)
            self.depth.set(raftId, ccdId, maxDepth)

            areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
            label = "photometric depth "
            comment = "magnitude where star completeness drops below 0.5"
            test = testCode.Test(label, maxDepth, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)



            
    def plotOld(self, data, dataId, showUndefined=False):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)

        # fpa figure
        filebase = "completenessDepth"
        depthData, depthMap = testSet.unpickle(filebase, default=[None, None])
        depths = []
        depthFig = qaFig.FpaQaFigure(data.cameraInfo, data=depthData, map=depthMap)

        for raft, ccdDict in depthFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.depth.get(raft, ccd) is None:
                    depth = self.depth.get(raft, ccd)
                    depths.append(depth)
                    depthFig.data[raft][ccd] = depth
                    depthFig.map[raft][ccd] = 'mag=%.2f'%(depth) 

        blue = '#0000ff'
        red  = '#ff0000'

        if False:
            if len(depths) >= 2:
                vmin = max(num.min(depths), self.limits[0])
                vmax = min(num.max(depths), self.limits[1])
            else:
                vmin = self.limits[0]
                vmax = self.limits[1]

            if vmax <= vmin:
                vmin = self.limits[0]
                vmax = self.limits[1]
        else:
            vmin, vmax = 1.0*self.limits[0], 1.0*self.limits[1]
            
        depthFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[vmin, vmax],
                            title="Photometric Depth", cmapOver=red, cmapUnder=blue,
                            failLimits=self.limits)
        testSet.addFigure(depthFig, filebase+".png", "Estimate of photometric depth",  navMap=True)
        testSet.pickle(filebase, [depthFig.data, depthFig.map])

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

                if hasMinuit:
                    A, B = self.fit.get(raft, ccd)
                    curve_x = num.arange(self.bins[0], self.bins[-1], 0.1)
                    curve_y = 0.5 + -1.0 / num.pi * num.arctan(A * curve_x + B)
                    sp2x2.plot(curve_x, curve_y, 'k-', alpha = 0.25)
                else:
                    sp2x2.axhline(y = 0.5, c='k', linestyle='-', alpha = 0.25)
                    sp2x2.axvline(x = self.depth.get(raft, ccd), c='k', linestyle='-', alpha = 0.25)

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

            sp1.set_xlim(14, 26)
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
