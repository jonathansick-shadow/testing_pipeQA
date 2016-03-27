import numpy as num
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .QaAnalysisTask import QaAnalysisTask

import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData
import lsst.testing.pipeQA.source as pqaSource
import QaAnalysisUtils as qaAnaUtil
import QaPlotUtils as qaPlotUtil


class VignettingQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str, doc = "Cameras to run VignettingQaTask",
                                  default = ("lsstSim", "cfht", "suprimecam", "hscSim", "sdss", "coadd"))
    maxMedian = pexConfig.Field(dtype = float, doc = "Maximum median magnitude offset", default = 0.02)
    maxRms = pexConfig.Field(dtype = float, doc = "Maximum magnitude offset RMS", default = 0.02)
    maxMag = pexConfig.Field(dtype = float, doc = "Maximum magnitude star to use in VignettingQa test",
                             default = 19.0)


class VignettingQaTask(QaAnalysisTask):
    ConfigClass = VignettingQaConfig
    _DefaultName = "vignettingQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)
        self.medLimits = [-self.config.maxMedian, self.config.maxMedian]
        self.rmsLimits = [0, self.config.maxRms]
        self.maxMag = self.config.maxMag

        self.magType1 = "ap"
        self.magType2 = "cat"

        self.sCatDummy = pqaSource.Catalog()
        self.srefCatDummy = pqaSource.RefCatalog()

        self.description = """
         For each CCD, the difference in aperture and reference catalog magnitudes
         is plotted as a function of radial location from the center of the
         focal plane.  The summary FPA figures show the median offset, as well
         as the standard deviation of this offset, for each chip.
        """

    def _getFlux(self, mType, s, sref):

        # if the source isn't valid, return NaN
        if not hasattr(s, 'getId') or not hasattr(sref, 'getId'):
            return num.NaN

        if mType == "psf":
            return s.getD(self.sCatDummy.PsfFluxKey)
        elif mType == "ap":
            return s.getD(self.sCatDummy.ApFluxKey)
        elif mType == "mod":
            return s.getD(self.sCatDummy.ModelFluxKey)
        elif mType == "cat":
            return sref.getD(self.srefCatDummy.PsfFluxKey)
        elif mType == "inst":
            return s.getD(self.sCatDummy.InstFluxKey)

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
        testSet.addMetadata({"Description": self.description})

        self.detector = data.getDetectorBySensor(dataId)
        self.filter = data.getFilterBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')

        self.dmag = raftCcdData.RaftCcdVector(self.detector)
        self.ids = raftCcdData.RaftCcdVector(self.detector)
        self.radius = raftCcdData.RaftCcdVector(self.detector)

        self.medianOffset = raftCcdData.RaftCcdData(self.detector)
        self.rmsOffset = raftCcdData.RaftCcdData(self.detector)

        #badFlags = pqaSource.INTERP_CENTER | pqaSource.SATUR_CENTER | pqaSource.EDGE

        for key in self.detector.keys():

            if self.detector[key] is None:
                continue

            raftId = self.detector[key].getParent().getId().getName()
            ccdId = self.detector[key].getId().getName()

            # We have a trimmed detector
            self.detector[key].setTrimmed(True)
            pixelSize = self.detector[key].getPixelSize()    # mm
            centerXm, centerYm = self.detector[key].getCenter().getMm()  # focal plane mm

            bbox = self.detector[key].getAllPixels()
            startX = bbox.getBeginX()
            startY = bbox.getBeginY()
            endX = bbox.getEndX()
            endY = bbox.getEndY()
            centerXp = 0.5 * (endX - startX)
            centerYp = 0.5 * (endY - startY)

            if self.matchListDictSrc.has_key(key):
                mdict = self.matchListDictSrc[key]['matched']
                for m in mdict:
                    sref, s, dist = m

                    if s.getD(self.sCatDummy.ExtendednessKey):  # if non-stellar
                        continue

                    f1 = self._getFlux(self.magType1, s, sref)
                    f2 = self._getFlux(self.magType2, s, sref)

                    intcen = s.getD(self.sCatDummy.FlagPixInterpCenKey)
                    satcen = s.getD(self.sCatDummy.FlagPixSaturCenKey)
                    edge = s.getD(self.sCatDummy.FlagPixEdgeKey)

                    if data.cameraInfo.name == 'coadd':
                        flagit = (satcen or edge)  # coadds have excessive area covered by InterpCen flags
                    else:
                        flagit = (intcen or satcen or edge)
                    if (f1 > 0.0 and f2 > 0.0 and not flagit):
                        m1 = -2.5*num.log10(f1)
                        m2 = -2.5*num.log10(f2)

                        if m2 > self.maxMag:
                            continue

                        if num.isfinite(m1) and num.isfinite(m2):
                            self.dmag.append(raftId, ccdId, m1 - m2)
                            self.ids.append(raftId, ccdId, str(s.getId()))

                            if data.cameraInfo.name == 'lsstSim':
                                # XY switched
                                xmm = centerXm + (s.getD(self.sCatDummy.YAstromKey) - centerXp)*pixelSize
                                ymm = centerYm + (s.getD(self.sCatDummy.XAstromKey) - centerYp)*pixelSize
                                radiusp = num.sqrt(xmm**2 + ymm**2) / pixelSize
                            else:
                                # XY not switch, and pixel centers not in mm
                                xmm = centerXm + (s.getD(self.sCatDummy.XAstromKey) - centerXp)
                                ymm = centerYm + (s.getD(self.sCatDummy.YAstromKey) - centerYp)
                                radiusp = num.sqrt(xmm**2 + ymm**2)
                            self.radius.append(raftId, ccdId, radiusp)

                # Calculate stats
                dmags = self.dmag.get(raftId, ccdId)
                if len(dmags) > 0:
                    med = num.median(dmags)
                else:
                    med = 0.0

                std = 0.0
                if len(dmags) > 1:
                    stat = afwMath.makeStatistics(dmags, afwMath.IQRANGE)
                    std = 0.741 * stat.getValue(afwMath.IQRANGE)
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
        testSet.setUseCache(self.useCache)  # cache
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        if self.showFpa:
            # fpa figures
            medFigbase = "vignettingMedianPhotOffset"  # cache
            medFigData, medFigMap = testSet.unpickle(medFigbase, [None, None])  # cache
            medFig = qaFig.FpaQaFigure(data.cameraInfo, data=medFigData, map=medFigMap)  # cache
            for raft, ccdDict in medFig.data.items():
                for ccd, value in ccdDict.items():
                    if self.medianOffset.get(raft, ccd) is not None:
                        med = self.medianOffset.get(raft, ccd)
                        medFig.data[raft][ccd] = med
                        if num.isfinite(med):
                            medFig.map[raft][ccd] = 'med=%.2f'%(med)
                        else:
                            medFig.map[raft][ccd] = 'med=nan'

            stdFigbase = "vignettingRmsPhotOffset"  # cache
            stdFigData, stdFigMap = testSet.unpickle(stdFigbase, [None, None])  # cache
            stdFig = qaFig.FpaQaFigure(data.cameraInfo, data=stdFigData, map=stdFigMap)  # cache
            for raft, ccdDict in stdFig.data.items():
                for ccd, value in ccdDict.items():
                    if self.rmsOffset.get(raft, ccd) is not None:
                        std = self.rmsOffset.get(raft, ccd)
                        stdFig.data[raft][ccd] = std
                        if num.isfinite(std):
                            stdFig.map[raft][ccd] = 'stddev=%.2f'%(std)
                        else:
                            stdFig.map[raft][ccd] = 'stddev=nan'

            testSet.pickle(medFigbase, [medFig.data, medFig.map])  # cache
            testSet.pickle(stdFigbase, [stdFig.data, stdFig.map])  # cache
            blue = '#0000ff'
            red = '#ff0000'

            if not self.delaySummary or isFinalDataId:
                self.log.log(self.log.INFO, "plotting FPAs")
                medFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=self.medLimits,
                                  title="Median offset", cmapOver=red, cmapUnder=blue,
                                  failLimits=self.medLimits)
                testSet.addFigure(medFig, medFigbase+".png",
                                  "Median offset of bright (m<%d) stars versus radius" % (self.maxMag),
                                  navMap=True)
                del medFig
                stdFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=self.rmsLimits,
                                  title="Stddev offset", cmapOver=red, cmapUnder=blue,
                                  failLimits=self.rmsLimits)
                testSet.addFigure(stdFig, stdFigbase+".png",
                                  "Stddev of bright (m < %d) stars as a function of radius" % (self.maxMag),
                                  navMap=True)
                del stdFig
            else:
                del medFig
                del stdFig

        cacheLabel = "vignetting_dmag"  # cache
        shelfData = {}

        # make any individual (ie. per sensor) plots
        for raft, ccd in self.dmag.raftCcdKeys():

            dmags = self.dmag.get(raft, ccd)
            radii = self.radius.get(raft, ccd)
            ids = self.ids.get(raft, ccd)
            med = self.medianOffset.get(raft, ccd)
            std = self.rmsOffset.get(raft, ccd)

            dataDict = {'dmags': dmags, 'radii': radii, 'ids': ids,
                        'offsetStats': [med, std],
                        'magTypes': [self.magType1, self.magType2],
                        'summary': False,
                        }

            self.log.log(self.log.INFO, "plotting %s" % (ccd))
            import VignettingQaPlot as plotModule
            label = data.cameraInfo.getDetectorName(raft, ccd)
            caption = "Delta magnitude vs. radius " + label
            pngFile = cacheLabel+".png"

            if self.lazyPlot.lower() in ['sensor', 'all']:
                testSet.addLazyFigure(dataDict, pngFile, caption,
                                      plotModule, areaLabel=label, plotargs="")
            else:
                testSet.cacheLazyData(dataDict, pngFile, areaLabel=label)
                fig = plotModule.plot(dataDict)
                testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                del fig

        if not self.delaySummary or isFinalDataId:
            self.log.log(self.log.INFO, "plotting Summary figure")

            import VignettingQaPlot as plotModule
            label = 'all'
            caption = "Delta magnitude vs. radius "+label
            pngFile = cacheLabel + ".png"

            if self.lazyPlot in ['all']:
                testSet.addLazyFigure({}, pngFile, caption,
                                      plotModule, areaLabel=label, plotargs="")
            else:
                dataDict, isSummary = qaPlotUtil.unshelveGlob(cacheLabel+"-all.png", testSet=testSet)
                dataDict['summary'] = True
                fig = plotModule.plot(dataDict)
                testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                del fig



