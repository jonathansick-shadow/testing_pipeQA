import re
import numpy as num
import time

import lsst.meas.algorithms as measAlg
import lsst.testing.pipeQA.TestCode as testCode
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .QaAnalysisTask import QaAnalysisTask
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import QaAnalysisUtils as qaAnaUtil
import RaftCcdData as raftCcdData

import lsst.testing.pipeQA.source as pqaSource

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse
import QaPlotUtils as qaPlotUtil


class ZeropointFitQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str,
                                  doc = "Cameras to run ZeropointFitQa",
                                  default = ("lsstSim", "cfht", "sdss", "coadd"))
    offsetMin = pexConfig.Field(dtype = float,
                                doc = "Median offset of stars from zeropoint fit; minimum good value",
                                default = -0.05)
    offsetMax = pexConfig.Field(dtype = float,
                                doc = "Median offset of stars from zeropoint fit; maximum good value",
                                default = +0.05)


class ZeropointFitQaTask(QaAnalysisTask):
    ConfigClass = ZeropointFitQaConfig
    _DefaultName = "zeropointFitQa"

    def __init__(self, figsize=(5.0, 5.0), **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)
        self.figsize = figsize
        self.limits = [self.config.offsetMin, self.config.offsetMax]

        self.sCatDummy = pqaSource.Catalog()
        self.srefCatDummy = pqaSource.RefCatalog()

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
        testSet.addMetadata({"Description": self.description})

        self.fluxType = fluxType
        self.detector = data.getDetectorBySensor(dataId)
        self.filter = data.getFilterBySensor(dataId)
        self.calib = data.getCalibBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')

        # Ignore blends in this analysis
        self.orphan = raftCcdData.RaftCcdVector(self.detector)
        self.matchedStar = raftCcdData.RaftCcdData(self.detector)
        self.undetectedStar = raftCcdData.RaftCcdVector(self.detector)
        self.matchedGalaxy = raftCcdData.RaftCcdData(self.detector)
        self.undetectedGalaxy = raftCcdData.RaftCcdVector(self.detector)

        # Results
        self.zeroPoint = raftCcdData.RaftCcdData(self.detector)
        self.medOffset = raftCcdData.RaftCcdData(self.detector)

        if data.cameraInfo.name == "coadd":
            badFlags = pqaSource.SATUR_CENTER | pqaSource.EDGE  # coadds have excessive area covered by INTERP_CENTER flags
        else:
            badFlags = pqaSource.INTERP_CENTER | pqaSource.SATUR_CENTER | pqaSource.EDGE

        for key in self.detector.keys():
            raftId = self.detector[key].getParent().getId().getName()
            ccdId = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            self.matchedStar.set(raftId, ccdId, {"Refmag": num.array([]),
                                                 "Imgmag": num.array([]),
                                                 "Imgerr": num.array([]), })
            self.matchedGalaxy.set(raftId, ccdId, {"Refmag": num.array([]),
                                                   "Imgmag": num.array([]),
                                                   "Imgerr": num.array([])})
            self.undetectedStar.set(raftId, ccdId, num.array([]))
            self.undetectedGalaxy.set(raftId, ccdId, num.array([]))
            self.orphan.set(raftId, ccdId, num.array([]))

            fmag0 = self.calib[key].getFluxMag0()[0]
            if fmag0 <= 0.0:
                self.zeroPoint.set(raftId, ccdId, NaN)
                continue
            zpt = -2.5*num.log10(fmag0)
            self.zeroPoint.set(raftId, ccdId, zpt)

            if self.matchListDictSrc.has_key(key):
                # Matched
                mdict = self.matchListDictSrc[key]['matched']
                stars = []
                galaxies = []

                for m in mdict:
                    sref, s, dist = m
                    if fluxType == "psf":
                        fref = sref.getD(self.srefCatDummy.PsfFluxKey)
                        f = s.getD(self.sCatDummy.PsfFluxKey)
                        ferr = s.getD(self.sCatDummy.PsfFluxErrKey)
                    else:
                        fref = sref.getD(self.srefCatDummy.PsfFluxKey)
                        f = s.getD(self.sCatDummy.ApFluxKey)
                        ferr = s.getD(self.sCatDummy.ApFluxErrKey)

                    # un-calibrate the magnitudes
                    f *= fmag0

                    intcen = s.getD(self.sCatDummy.FlagPixInterpCenKey)
                    satcen = s.getD(self.sCatDummy.FlagPixSaturCenKey)
                    edge = s.getD(self.sCatDummy.FlagPixEdgeKey)

                    if data.cameraInfo.name == 'coadd':
                        flagit = (satcen or edge)  # coadds have excessive area covered by InterpCen flags
                    else:
                        flagit = (intcen or satcen or edge)

                    if (fref > 0.0 and f > 0.0 and not flagit):
                        mrefmag = -2.5*num.log10(fref)
                        mimgmag = -2.5*num.log10(f)
                        mimgmerr = 2.5 / num.log(10.0) * ferr / f

                        star = 0 if s.getD(self.sCatDummy.ExtendednessKey) else 1

                        if num.isfinite(mrefmag) and num.isfinite(mimgmag):
                            if star:
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
                        f = orphan.getD(self.sCatDummy.PsfFluxKey)
                    else:
                        f = orphan.getD(self.sCatDummy.ApFluxKey)
                    if f > 0.0:
                        # un-calibrate the magnitudes
                        f *= fmag0

                        orphans.append(-2.5 * num.log10(f))

                self.orphan.set(raftId, ccdId, num.array(orphans))

                # Metrics
                offset = num.array(self.matchedStar.get(raftId, ccdId)["Imgmag"])  # make a copy
                offset -= self.zeroPoint.get(raftId, ccdId)
                offset -= self.matchedStar.get(raftId, ccdId)["Refmag"]
                med = num.median(offset)
                self.medOffset.set(raftId, ccdId, med)

                areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
                label = "median offset from zeropoint"
                comment = "Median offset of calibrated stellar magnitude to zeropoint fit"
                test = testCode.Test(label, med, self.limits, comment, areaLabel=areaLabel)
                testSet.addTest(test)

                label = "median zeropoint"
                comment = "Median zeropoint measured for sensor"
                test = testCode.Test(label, zpt, [None, 0], comment, areaLabel=areaLabel)
                testSet.addTest(test)

    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        if self.showFpa:
            # fpa figure
            zpts = []
            zptBase = "zeropoint"
            zptData, zptMap = testSet.unpickle(zptBase, default=[None, None])
            zptFig = qaFig.FpaQaFigure(data.cameraInfo, data=zptData, map=zptMap)

            offsetBase = "medZeropointOffset"
            offsetData, offsetMap = testSet.unpickle(offsetBase, default=[None, None])
            offsetFig = qaFig.FpaQaFigure(data.cameraInfo, data=offsetData, map=offsetMap)

            for raft, ccdDict in zptFig.data.items():
                for ccd, value in ccdDict.items():
                    if self.zeroPoint.get(raft, ccd) is not None:
                        zpt = self.zeroPoint.get(raft, ccd)
                        zpts.append(zpt)
                        zptFig.data[raft][ccd] = zpt
                        zptFig.map[raft][ccd] = 'zpt=%.2f' % (zpt)

                        offset = self.medOffset.get(raft, ccd)
                        offsetFig.data[raft][ccd] = offset
                        offsetFig.map[raft][ccd] = 'offset=%.2f' % (offset)
                    else:
                        if zptFig.data[raft][ccd] is not None:
                            zpts.append(zptFig.data[raft][ccd])

            testSet.pickle(zptBase, [zptFig.data, zptFig.map])
            testSet.pickle(offsetBase, [offsetFig.data, offsetFig.map])

            if not self.delaySummary or isFinalDataId:
                self.log.log(self.log.INFO, "plotting FPAs")

                blue = '#0000ff'
                red = '#ff0000'
                zptFig.makeFigure(showUndefined=showUndefined, cmap="jet",
                                  vlimits=[num.min(zpts)-0.05, num.max(zpts)+0.05],
                                  title="Zeropoint", cmapOver=red, cmapUnder=blue)
                testSet.addFigure(zptFig, zptBase+".png", "Photometric zeropoint", navMap=True)
                del zptFig

                offsetFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=self.limits,
                                     title="Med offset from Zpt Fit", cmapOver=red, failLimits=self.limits,
                                     cmapUnder=blue)
                testSet.addFigure(offsetFig, offsetBase + ".png", "Median offset from photometric zeropoint",
                                  navMap=True)
                del offsetFig
            else:
                del zptFig, offsetFig

        cacheLabel = "zeropointFit"
        shelfData = {}

        # Each CCD
        for raft, ccd in self.zeroPoint.raftCcdKeys():
            zeropt = self.zeroPoint.get(raft, ccd)
            if zeropt == 0.0:
                continue

            self.log.log(self.log.INFO, "Plotting %s" % (ccd))

            # Plot all matched galaxies
            mrefGmag = self.matchedGalaxy.get(raft, ccd)["Refmag"]
            mimgGmag = self.matchedGalaxy.get(raft, ccd)["Imgmag"]
            mimgGmerr = self.matchedGalaxy.get(raft, ccd)["Imgerr"]

            # Plot all matched stars
            mrefSmag = self.matchedStar.get(raft, ccd)["Refmag"]
            mimgSmag = self.matchedStar.get(raft, ccd)["Imgmag"]
            mimgSmerr = self.matchedStar.get(raft, ccd)["Imgerr"]

            urefmag = num.concatenate((self.undetectedStar.get(raft, ccd),
                                       self.undetectedGalaxy.get(raft, ccd)))
            uimgmag = self.orphan.get(raft, ccd)

            label = data.cameraInfo.getDetectorName(raft, ccd)
            dataDict = {
                'mrefGmag': mrefGmag, 'mimgGmag': mimgGmag, 'mimgGmerr': mimgGmerr,
                'mrefSmag': mrefSmag, 'mimgSmag': mimgSmag, 'mimgSmerr': mimgSmerr,
                'urefmag': urefmag, 'uimgmag': uimgmag,
                'zeropt': zeropt, 'title': label, 'figsize': self.figsize,
                'fluxType': self.fluxType,
            }

            import ZeropointFitQaPlot as plotModule
            caption = "Zeropoint fit " + label
            pngFile = cacheLabel + ".png"

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

            label = 'all'
            import ZeropointFitQaPlot as plotModule
            caption = "Zeropoint fit " + label
            pngFile = cacheLabel + ".png"

            if self.lazyPlot.lower() in ['all']:
                testSet.addLazyFigure(dataDict, pngFile, caption,
                                      plotModule, areaLabel=label, plotargs="")
            else:
                dataDict, isSummary = qaPlotUtil.unshelveGlob(cacheLabel+"-all.png", testSet=testSet)
                dataDict['summary'] = True
                fig = plotModule.plot(dataDict)
                testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                del fig

