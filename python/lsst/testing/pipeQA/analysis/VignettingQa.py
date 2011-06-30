import numpy as num
import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData
import lsst.meas.algorithms as measAlg

class VignettingQa(qaAna.QaAnalysis):
    def __init__(self, maxMedian, maxRms, maxMag):
        qaAna.QaAnalysis.__init__(self)
        self.medLimits = [-1 * maxMedian, maxMedian]
        self.rmsLimits = [0, maxRms]
        self.maxMag    = maxMag

        self.magType1 = "psf"
        self.magType2 = "cat"

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
        del self.detector
        del self.filter
        del self.matchListDictSrc

        del self.dmag
        del self.ids
        del self.radius
        
        del self.medianOffset
        del self.rmsOffset

    def test(self, data, dataId):
        testSet = self.getTestSet(data, dataId)

        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        
        self.dmag    = raftCcdData.RaftCcdVector(self.detector)
        self.ids     = raftCcdData.RaftCcdVector(self.detector)
        self.radius  = raftCcdData.RaftCcdVector(self.detector)

        self.medianOffset = raftCcdData.RaftCcdVector(self.detector)
        self.rmsOffset    = raftCcdData.RaftCcdVector(self.detector)
        
        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER
        
        for key in self.detector.keys():

            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            # We have a trimmed detector
            self.detector[key].setTrimmed(True)
            pixelSize          = self.detector[key].getPixelSize()    # mm
            centerXm, centerYm = self.detector[key].getCenter()       # focal plane mm

            bbox   = self.detector[key].getAllPixels()
            startX = bbox.getBeginX()
            startY = bbox.getBeginY()
            endX   = bbox.getEndX()
            endY   = bbox.getEndY()
            centerXp = 0.5 * (endX - startX)
            centerYp = 0.5 * (endY - startY)
            
            if self.matchListDictSrc.has_key(key):
                mdict    = self.matchListDictSrc[key]['matched']
                for m in mdict:
                    sref, s, dist = m

                    if not sref.getFlagForDetection() & measAlg.Flags.STAR:
                        continue

                    f1 = self._getFlux(self.magType1, s, sref)
                    f2 = self._getFlux(self.magType2, s, sref)

                    flags = s.getFlagForDetection()
                    
                    if (f1 > 0.0 and f2 > 0.0  and not flags & badFlags):
                        m1 = -2.5*num.log10(f1)
                        m2 = -2.5*num.log10(f2)

                        if m2 > self.maxMag:
                            continue

                        if num.isfinite(m1) and num.isfinite(m2):
                            self.dmag.append(raftId, ccdId, m1 - m2)
                            self.ids.append(raftId, ccdId, str(s.getId()))
                            # XY switched
                            xmm     = centerXm + (s.getYAstrom() - centerXp) * pixelSize
                            ymm     = centerYm + (s.getXAstrom() - centerYp) * pixelSize
                            radiusp = num.sqrt(xmm**2 + ymm**2) / pixelSize
                            self.radius.append(raftId, ccdId, radiusp)

                # Calculate stats
                dmags = self.dmag.get(raftId, ccdId)
                med   = num.median(dmags)
                std   = num.std(dmags)
                self.medianOffset.set(raftId, ccdId, med)
                self.rmsOffset.set(raftId, ccdId, std)
                
                areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)

                label = "median offset "
                comment = "median offset from cat mag"
                test = testCode.Test(label, med, self.medLimits, comment, areaLabel=areaLabel)
                testSet.addTest(test)

                label = "stddev offset "
                comment = "stddev of offset from cat mag"
                test = testCode.Test(label, std, self.rmsLimits, comment, areaLabel=areaLabel)
                testSet.addTest(test)

    def plot(self, data, dataId, showUndefined = False):
        testSet = self.getTestSet(data, dataId)

        # fpa figures
        medFig = qaFig.FpaQaFigure(data.cameraInfo)
        for raft, ccdDict in medFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.medianOffset.get(raft, ccd) is None:
                    med = self.medianOffset.get(raft, ccd)
                    medFig.data[raft][ccd] = med
                    medFig.map[raft][ccd] = 'med=%.2f'%(med) 

        stdFig = qaFig.FpaQaFigure(data.cameraInfo)
        for raft, ccdDict in stdFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.rmsOffset.get(raft, ccd) is None:
                    std = self.rmsOffset.get(raft, ccd)
                    stdFig.data[raft][ccd] = std
                    stdFig.map[raft][ccd] = 'stddev=%.2f'%(std) 

        blue = '#0000ff'
        red  = '#ff0000'
        
        medFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=self.medLimits,
                          title="Median offset", cmapOver=red, cmapUnder=blue,
                          failLimits=self.medLimits)
        testSet.addFigure(medFig, "vignettingMedianPhotOffset.png", "Median offset of bright (m < %d) stars as a function of radius" % (self.maxMag), 
                          navMap=True)

        stdFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=self.rmsLimits,
                          title="Stddev offset", cmapOver=red, cmapUnder=blue,
                          failLimits=self.rmsLimits)
        testSet.addFigure(stdFig, "vignettingRmsPhotOffset.png", "Stddev of bright (m < %d) stars as a function of radius" % (self.maxMag), 
                          navMap=True)

        dmagsAll = []
        radiiAll = []
        for raft, ccd in self.dmag.raftCcdKeys():
            dmags = self.dmag.get(raft, ccd)
            radii = self.radius.get(raft, ccd)
            ids   = self.ids.get(raft, ccd)
            
            print "Plotting ", ccd
            fig = qaFig.QaFigure()
            sp1 = fig.fig.add_subplot(111)
            sp1.plot(radii, dmags, 'ro')

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            ddmag = 0.05
            drad  = 0.05 * (max(radii) - min(radii))
            for i in range(len(dmags)):
                info = "nolink:sourceId=%s" % (ids[i])
                area = (radii[i]-drad, dmags[i]-ddmag, radii[i]+drad, dmags[i]+ddmag)
                fig.addMapArea(areaLabel, area, info, axes=sp1)

                dmagsAll.append(dmags[i])
                radiiAll.append(radii[i])
                
            sp1.axhline(y=0, c = 'k', linestyle = ':', alpha = 0.25)
            sp1.axhline(y=num.median(dmags), c = 'b', linestyle = '-')
            sp1.axhspan(ymin = num.median(dmags)-num.std(dmags), ymax = num.median(dmags)+num.std(dmags), fc = 'b', alpha = 0.15)
            sp1x2 = sp1.twinx()
            ylab = sp1x2.set_ylabel('Delta magnitude (%s-%s)' % (self.magType1, self.magType2), fontsize=10)
            ylab.set_rotation(-90)
            sp1.set_xlabel('Dist from focal plane center (pixels)', fontsize=10)
            qaFigUtils.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize = 8)
            qaFigUtils.qaSetp(sp1x2.get_xticklabels()+sp1x2.get_yticklabels(), visible=False)
                        
            label = data.cameraInfo.getDetectorName(raft, ccd)
            testSet.addFigure(fig, "vignetting_dmag.png", "Delta magnitude vs. radius "+label, areaLabel=label)

        fig = qaFig.QaFigure()
        sp1 = fig.fig.add_subplot(111)
        sp1.plot(radiiAll, dmagsAll, 'ro', ms=2, alpha = 0.5)
        sp1.axhline(y=0, c = 'k', linestyle = ':', alpha = 0.25)
        sp1x2 = sp1.twinx()
        ylab = sp1x2.set_ylabel('Delta magnitude (%s-%s)' % (self.magType1, self.magType2), fontsize=10)
        ylab.set_rotation(-90)
        sp1.set_xlabel('Dist from focal plane center (pixels)', fontsize=10)
        qaFigUtils.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize = 8)
        qaFigUtils.qaSetp(sp1x2.get_xticklabels()+sp1x2.get_yticklabels(), visible=False)
        label = "all"
        testSet.addFigure(fig, "vignetting_dmag.png", "Delta magnitude vs. radius "+label, areaLabel=label)

        
                
                
