import sys
import os
import re
import numpy
import lsst.meas.algorithms as measAlg
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

from .QaAnalysisTask import QaAnalysisTask
import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures as qaFig
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil
import QaPlotUtils as qaPlotUtil

import lsst.testing.pipeQA.source as pqaSource

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.font_manager as fm
from matplotlib.collections import LineCollection


class EmptySectorQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str,
                                  doc = "Cameras to run EmptySectorQaTask",
                                  default = ("lsstSim", "hscSim", "suprimecam", "cfht", "sdss", "coadd"))
    maxMissing = pexConfig.Field(dtype = int, doc = "Maximum number of missing CCDs", default = 1)
    nx = pexConfig.Field(dtype = int, doc = "Mesh size in x", default = 4)
    ny = pexConfig.Field(dtype = int, doc = "Mesh size in y", default = 4)


class EmptySectorQaTask(QaAnalysisTask):
    ConfigClass = EmptySectorQaConfig
    _DefaultName = "emptySectorQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)
        self.limits = [0, self.config.maxMissing]
        self.nx = self.config.nx
        self.ny = self.config.ny

        self.sCatDummy = pqaSource.Catalog()
        self.srefCatDummy = pqaSource.RefCatalog()

        self.description = """
         For each CCD, the 1-to-1 matches between the reference catalog and
         sources are plotted as a function of position in the focal plane.
         Portions of each CCD where there are no matches are considered "empty
         sectors" and suggest a problem with the reference catalog matching.
         The summary FPA figure shows the number of empty sectors per CCD.
        """

    def free(self):

        # please free all large data structures here.
        del self.y
        del self.x
        del self.ymat
        del self.xmat
        del self.size
        del self.filter
        del self.detector
        del self.ssDict
        del self.matchListDictSrc

        del self.emptySectors
        del self.emptySectorsMat

    def test(self, data, dataId):

        # get data
        self.ssDict = data.getSourceSetBySensor(dataId)

        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector = data.getDetectorBySensor(dataId)
        self.filter = data.getFilterBySensor(dataId)

        # create containers for data we're interested in
        self.x = raftCcdData.RaftCcdVector(self.detector)
        self.y = raftCcdData.RaftCcdVector(self.detector)
        self.xmat = raftCcdData.RaftCcdVector(self.detector)
        self.ymat = raftCcdData.RaftCcdVector(self.detector)

        # fill containers with values we need for our test
        filter = None
        self.size = raftCcdData.RaftCcdData(self.detector, initValue=[1.0, 1.0])
        for key, ss in self.ssDict.items():
            xKey = self.sCatDummy.XAstromKey
            yKey = self.sCatDummy.YAstromKey

            raft = self.detector[key].getParent().getId().getName()
            ccd = self.detector[key].getId().getName()
            bbox = self.detector[key].getAllPixels(True)
            size = [bbox.getMaxX() - bbox.getMinX(), bbox.getMaxY() - bbox.getMinY()]
            self.size.set(raft, ccd, size)
            filter = self.filter[key].getName()

            for s in ss:
                self.x.append(raft, ccd, s.getD(xKey))
                self.y.append(raft, ccd, s.getD(yKey))

            if self.matchListDictSrc.has_key(key):
                for m in self.matchListDictSrc[key]['matched']:
                    sref, s, dist = m
                    self.xmat.append(raft, ccd, s.getD(xKey))
                    self.ymat.append(raft, ccd, s.getD(yKey))

        # create a testset
        testSet = self.getTestSet(data, dataId)

        # this normally gets set in the plot as that's where the caching happens,
        # here we're stashing the nDetection and nCcd values, so we need to set it early.
        testSet.setUseCache(self.useCache)
        testSet.addMetadata({"Description": self.description})

        # analyse each sensor and put the values in a raftccd container
        self.emptySectors = raftCcdData.RaftCcdData(self.detector, initValue=self.nx*self.ny)
        self.emptySectorsMat = raftCcdData.RaftCcdData(self.detector, initValue=self.nx*self.ny)

        countBase = "countShelf"
        nShelf = testSet.unshelve(countBase)

        for raft, ccd in self.emptySectors.raftCcdKeys():
            x, y = self.x.get(raft, ccd), self.y.get(raft, ccd)
            xmat, ymat = self.xmat.get(raft, ccd), self.ymat.get(raft, ccd)
            xwid, ywid = self.size.get(raft, ccd)

            xlo, ylo = 0.0, 0.0
            if data.cameraInfo.name == 'coadd':
                xlo, ylo, xhi, yhi = x.min(), y.min(), x.max(), y.max()
                xwid, ywid = xhi-xlo, yhi-ylo

            def countEmptySectors(x, y):
                counts = numpy.zeros([self.nx, self.ny])
                for i in range(len(x)):
                    xi, yi = int(self.nx*x[i]/xwid), int(self.ny*y[i]/ywid)
                    if xi >= 0 and xi < self.nx and yi >= 0 and yi < self.ny:
                        counts[xi, yi] += 1
                whereEmpty = numpy.where(counts.flatten() == 0)[0]
                nEmpty = len(whereEmpty)
                return nEmpty

            nEmpty = countEmptySectors(x - xlo, y - ylo)
            nEmptyMat = countEmptySectors(xmat - xlo, ymat - ylo)
            self.emptySectors.set(raft, ccd, nEmpty)
            self.emptySectorsMat.set(raft, ccd, nEmptyMat)

            # add tests for acceptible numpy of empty sectors
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "empty ccd regions"
            comment = "%dx%d (nstar=%d)" % (self.nx, self.ny, len(x))

            test = testCode.Test(label, nEmpty, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)
            test = testCode.Test(label+" (matched)", nEmptyMat, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)

            nShelf[ccd] = len(x)

        testSet.shelve(countBase, nShelf)

        nCcd, nDet = 0, 0
        for k, v in nShelf.items():
            if v > 0:
                nCcd += 1
                nDet += v

        # a bit sketchy adding tests in the plot section, but these are dummies
        # they pass useful numbers through to the display, but don't actually test
        # anything useful
        test = testCode.Test("nDetections", nDet, [1, None], "number of detected sources", areaLabel="all")
        testSet.addTest(test)
        test = testCode.Test("nCcd", nCcd, [1, None], "number of ccds processed", areaLabel="all")
        testSet.addTest(test)

    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        if self.showFpa:
            # make fpa figures - for all detections, and for matched detections
            emptyBase = "emptySectors"
            emptyMatBase = "aa_emptySectorsMat"

            emptyData, emptyMap = testSet.unpickle(emptyBase, [None, None])
            emptyMatData, emptyMatMap = testSet.unpickle(emptyMatBase, [None, None])

            emptyFig = qaFig.FpaQaFigure(data.cameraInfo, data=emptyData, map=emptyMap)
            emptyFigMat = qaFig.FpaQaFigure(data.cameraInfo, data=emptyMatData, map=emptyMatMap)

            for raft, ccdDict in emptyFig.data.items():
                for ccd, value in ccdDict.items():

                    # set values for data[raft][ccd] (color coding)
                    # set values for map[raft][ccd]  (tooltip text)
                    if self.emptySectors.get(raft, ccd) is not None:
                        nEmpty = self.emptySectors.get(raft, ccd)
                        emptyFig.data[raft][ccd] = nEmpty
                        emptyFig.map[raft][ccd] = "%dx%d,empty=%d" % (self.nx, self.ny, nEmpty)

                        nEmptyMat = self.emptySectorsMat.get(raft, ccd)
                        emptyFigMat.data[raft][ccd] = nEmptyMat
                        emptyFigMat.map[raft][ccd] = "%dx%d,empty=%d" % (self.nx, self.ny, nEmptyMat)

            testSet.pickle(emptyBase, [emptyFig.data, emptyFig.map])
            testSet.pickle(emptyMatBase, [emptyFigMat.data, emptyFigMat.map])

            # make the figures and add them to the testSet
            # sample colormaps at: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
            if not self.delaySummary or isFinalDataId:
                self.log.log(self.log.INFO, "plotting FPAs")
                emptyFig.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r",
                                    vlimits=[0, self.nx*self.ny],
                                    title="Empty sectors (%dx%d grid)" % (self.nx, self.ny),
                                    failLimits=self.limits)
                testSet.addFigure(emptyFig, emptyBase+".png",
                                  "Empty Sectors in %dx%d grid." % (self.nx, self.ny), navMap=True)
                del emptyFig

                emptyFigMat.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r",
                                       vlimits=[0, self.nx*self.ny],
                                       title="Empty sectors (matched, %dx%d grid)" % (self.nx, self.ny),
                                       failLimits=self.limits)
                testSet.addFigure(emptyFigMat, emptyMatBase+".png",
                                  "Empty Sectors in %dx%d grid." % (self.nx, self.ny), navMap=True)
                del emptyFigMat
            else:
                del emptyFig, emptyFigMat

        cacheLabel = "pointPositions"
        shelfData = {}

        xlo, xhi, ylo, yhi = 1.e10, -1.e10, 1.e10, -1.e10
        for raft, ccd in data.cameraInfo.raftCcdKeys:
            if data.cameraInfo.name == 'coadd':
                xtmp, ytmp = self.x.get(raft, ccd), self.y.get(raft, ccd)
                if (xtmp is not None) and (ytmp is not None):
                    xxlo, yylo, xxhi, yyhi = xtmp.min(), ytmp.min(), xtmp.max(), ytmp.max()
                else:
                    xxlo, yylo, xxhi, yyhi = xlo, ylo, xhi, yhi
            else:
                xxlo, yylo, xxhi, yyhi = data.cameraInfo.getBbox(raft, ccd)
            if xxlo < xlo:
                xlo = xxlo
            if xxhi > xhi:
                xhi = xxhi
            if yylo < ylo:
                ylo = yylo
            if yyhi > yhi:
                yhi = yyhi

        # make any individual (ie. per sensor) plots
        for raft, ccd in self.emptySectors.raftCcdKeys():

            # get the data we want for this sensor (we stored it here in test() method above)
            x, y = self.x.get(raft, ccd), self.y.get(raft, ccd)
            xmat, ymat = self.xmat.get(raft, ccd), self.ymat.get(raft, ccd)
            xwid, ywid = self.size.get(raft, ccd)

            if data.cameraInfo.name == 'coadd':
                xmin, ymin, xmax, ymax = x.min(), y.min(), x.max(), y.max()
                x -= xmin
                y -= ymin
                xmat -= xmin
                ymat -= ymin
                xxlo, yylo, xxhi, yyhi = xmin, ymin, xmax, ymax
            else:
                xxlo, yylo, xxhi, yyhi = data.cameraInfo.getBbox(raft, ccd)

            dataDict = {'x': x+xxlo, 'y': y+yylo, 'xmat': xmat+xxlo, 'ymat': ymat+yylo,
                        'limits': [0, xwid, 0, ywid],
                        'summary': False, 'alllimits': [xlo, xhi, ylo, yhi],
                        'bbox': [xxlo, xxhi, yylo, yyhi],
                        'nxn': [self.nx, self.ny]}

            self.log.log(self.log.INFO, "plotting %s" % (ccd))
            import EmptySectorQaAnalysisPlot as plotModule
            label = data.cameraInfo.getDetectorName(raft, ccd)
            caption = "Pixel coordinates of all (black) and matched (red) detections." + label
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

            import EmptySectorQaAnalysisPlot as plotModule
            label = 'all'
            caption = "Pixel coordinates of all (black) and matched (red) detections." + label
            pngFile = "pointPositions.png"

            if self.lazyPlot in ['all']:
                testSet.addLazyFigure({}, cacheLabel+".png", caption,
                                      plotModule, areaLabel=label, plotargs="")
            else:
                dataDict, isSummary = qaPlotUtil.unshelveGlob(cacheLabel+"-all.png", testSet=testSet)
                dataDict['summary'] = True
                dataDict['limits'] = dataDict['alllimits']
                fig = plotModule.plot(dataDict)
                testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                del fig


