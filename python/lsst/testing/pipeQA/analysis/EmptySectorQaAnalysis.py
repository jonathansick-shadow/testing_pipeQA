import sys, os, re
import lsst.meas.algorithms        as measAlg
import lsst.testing.pipeQA.figures as qaFig
import numpy

import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.font_manager as fm
from matplotlib.collections import LineCollection


class EmptySectorQaAnalysis(qaAna.QaAnalysis):

    def __init__(self, maxMissing, nx=4, ny=4, **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.limits = [0, maxMissing]
        self.nx = nx
        self.ny = ny

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
        self.ssDict           = data.getSourceSetBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)

        # create containers for data we're interested in
        self.x     = raftCcdData.RaftCcdVector(self.detector)
        self.y     = raftCcdData.RaftCcdVector(self.detector)
        self.xmat  = raftCcdData.RaftCcdVector(self.detector)
        self.ymat  = raftCcdData.RaftCcdVector(self.detector)

        # fill containers with values we need for our test
        filter = None
        self.size = raftCcdData.RaftCcdData(self.detector, initValue=[1.0, 1.0])
        for key, ss in self.ssDict.items():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()
            bbox = self.detector[key].getAllPixels()
            size = [bbox.getMaxX() - bbox.getMinX(), bbox.getMaxY() - bbox.getMinY()]
            self.size.set(raft, ccd, size)
            filter = self.filter[key].getName()
            for s in ss:
                self.x.append(raft, ccd, s.getXAstrom())
                self.y.append(raft, ccd, s.getYAstrom())
            if self.matchListDictSrc.has_key(key):
                for m in self.matchListDictSrc[key]['matched']:
                    sref, s, dist = m
                    self.xmat.append(raft, ccd, s.getXAstrom())
                    self.ymat.append(raft, ccd, s.getYAstrom())

        # create a testset
        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})

        # analyse each sensor and put the values in a raftccd container
        self.emptySectors    = raftCcdData.RaftCcdData(self.detector, initValue=self.nx*self.ny)
        self.emptySectorsMat = raftCcdData.RaftCcdData(self.detector, initValue=self.nx*self.ny)
        for raft, ccd in self.emptySectors.raftCcdKeys():
            x, y       = self.x.get(raft, ccd), self.y.get(raft, ccd)
            xmat, ymat = self.xmat.get(raft, ccd), self.ymat.get(raft, ccd)
            xwid, ywid = self.size.get(raft, ccd)

            def countEmptySectors(x, y):
                counts = numpy.zeros([self.nx, self.ny])
                for i in range(len(x)):
                    xi, yi = int(self.nx*x[i]/xwid), int(self.ny*y[i]/ywid)
                    if xi >= 0 and xi < self.nx and yi >= 0 and yi < self.ny:
                        counts[xi,yi] += 1
                whereEmpty = numpy.where(counts.flatten() == 0)[0]
                nEmpty = len(whereEmpty)
                return nEmpty

            nEmpty = countEmptySectors(x, y)
            nEmptyMat = countEmptySectors(xmat, ymat)
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

            

    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        
        # make fpa figures - for all detections, and for matched detections
        emptyBase = "emptySectors"
        emptyMatBase = "aa_emptySectorsMat"

        emptyData, emptyMap       = testSet.unpickle(emptyBase, [None, None])
        emptyMatData, emptyMatMap = testSet.unpickle(emptyMatBase, [None, None])
        
        emptyFig    = qaFig.FpaQaFigure(data.cameraInfo, data=emptyData, map=emptyMap)
        emptyFigMat = qaFig.FpaQaFigure(data.cameraInfo, data=emptyMatData, map=emptyMatMap)

        for raft, ccdDict in emptyFig.data.items():
            for ccd, value in ccdDict.items():

                # set values for data[raft][ccd] (color coding)
                # set values for map[raft][ccd]  (tooltip text)
                if not self.emptySectors.get(raft, ccd) is None:
                    nEmpty = self.emptySectors.get(raft, ccd)
                    emptyFig.data[raft][ccd] = nEmpty
                    emptyFig.map[raft][ccd] = "%dx%d,empty=%d" % (self.nx, self.ny, nEmpty)
                    
                    nEmptyMat = self.emptySectorsMat.get(raft, ccd)
                    emptyFigMat.data[raft][ccd] = nEmptyMat
                    emptyFigMat.map[raft][ccd] = "%dx%d,empty=%d" % (self.nx, self.ny, nEmptyMat)

        # make the figures and add them to the testSet
        # sample colormaps at: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
        emptyFig.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r", vlimits=[0, self.nx*self.ny],
                            title="Empty sectors (%dx%d grid)" % (self.nx, self.ny),
                            failLimits=self.limits)
        testSet.addFigure(emptyFig, emptyBase+".png", "Empty Sectors in %dx%d grid." % (self.nx, self.ny),
                          navMap=True)
        testSet.pickle(emptyBase, [emptyFig.data, emptyFig.map])

        emptyFigMat.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r", vlimits=[0, self.nx*self.ny],
                               title="Empty sectors (matched, %dx%d grid)" % (self.nx, self.ny),
                               failLimits=self.limits)
        testSet.addFigure(emptyFigMat, emptyMatBase+".png",
                          "Empty Sectors in %dx%d grid." % (self.nx, self.ny), navMap=True)
        testSet.pickle(emptyMatBase, [emptyFigMat.data, emptyFigMat.map])


        # make any individual (ie. per sensor) plots
        figsize = (4.0, 4.0)
        for raft, ccd in self.emptySectors.raftCcdKeys():

            # get the data we want for this sensor (we stored it here in test() method above)
            x, y       = self.x.get(raft, ccd), self.y.get(raft, ccd)
            xmat, ymat = self.xmat.get(raft, ccd), self.ymat.get(raft, ccd)
            xwid, ywid = self.size.get(raft, ccd)

            # handle no-data possibility
            if len(x) == 0:
                x = numpy.array([0.0])
                y = numpy.array([0.0])
            if len(xmat) == 0:
                xmat = numpy.array([0.0])
                ymat = numpy.array([0.0])

            print "plotting ", ccd

            ####################
            # create the plot
            fig = qaFig.QaFigure(size=figsize)
            ax = fig.fig.add_subplot(111)
            ax.plot(x, y, "k.", ms=2.0, label="detected")
            ax.plot(xmat, ymat, "ro", ms=4.0, label="matched",
                    mfc='None', markeredgecolor='r')

            ax.set_xlim([0, xwid])
            ax.set_ylim([0, ywid])
            ax.legend(prop=fm.FontProperties(size ="xx-small"), ncol=2, loc="upper center")
            for tic in ax.get_xticklabels() + ax.get_yticklabels():
                tic.set_size("x-small")

            # show the regions
            for i in range(self.nx):
                xline = (i+1)*xwid/self.nx
                ax.axvline(xline, color="k")
            for i in range(self.ny):
                yline = (i+1)*ywid/self.ny
                ax.axhline(yline, color="k")

            # add map areas to allow mouseover tooltip showing pixel coords
            dx, dy = 20, 20  # on a 4kx4k ccd, < +/-20 pixels is tough to hit with a mouse
            for i in range(len(x)):
                area = x[i]-dx, y[i]-dy, x[i]+dx, y[i]+dy
                fig.addMapArea("no_label_info", area, "nolink:%.1f_%.1f"%(x[i],y[i]))

            # add the plot to the testSet
            areaLabel = emptyFig.getAreaLabel(raft, ccd)
            testSet.addFigure(fig, "pointPositions.png",
                              "Pixel coordinates of all (black) and matched (red) detections.",
                              areaLabel=areaLabel)
            
            
