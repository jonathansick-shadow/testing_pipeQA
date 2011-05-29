import re
import numpy as num

import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData

from matplotlib.font_manager import FontProperties


class CompletenessQa(qaAna.QaAnalysis):
    def __init__(self, completenessMagMin, completenessMagMax):
        qaAna.QaAnalysis.__init__(self)
        self.limits = [completenessMagMin, completenessMagMax]
        self.bins   = num.arange(14, 27, 0.25)

    def test(self, data, dataId, fluxType = "psf"):
        testSet = self.getTestSet(data, dataId)

        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        self.catData       = raftCcdData.RaftCcdData(self.detector)
        self.detData       = raftCcdData.RaftCcdData(self.detector)
        self.depth         = raftCcdData.RaftCcdData(self.detector)

        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            print "Running", ccdId
            
            catsql      = 'select sro.refObjectId, sro.ra, sro.decl, %sMag' % (filterName)
            catsql     += ' from SimRefObject as sro, Science_Ccd_Exposure as sce'
            catsql     += ' where (sce.visit = %s)' % (dataId['visit'])
            catsql     += ' and (sce.raftName = "%s")' % (re.sub("R:", "", raftId))
            catsql     += ' and (sce.ccdName = "%s")' % (ccdId[-3:])
            catsql     += ' and (sro.isStar = 1) '
            catsql     += ' and qserv_ptInSphPoly(sro.ra, sro.decl,'
            catsql     += ' concat_ws(" ", sce.llcRa, sce.llcDecl, sce.lrcRa, sce.lrcDecl, '
            catsql     += ' sce.urcRa, sce.urcDecl, sce.ulcRa, sce.ulcDecl))'
            catresults  = data.dbInterface.execute(catsql)
            self.catData.set(raftId, ccdId, num.array([x[3] for x in catresults]))
            
            detsql      = 'select dnToAbMag(s.%sFlux, sce.fluxMag0)' % (fluxType)
            detsql     += ' from Source as s, Science_Ccd_Exposure as sce' 
            detsql     += ' where ((s.flagForDetection & 0xa01) = 0)'
            detsql     += ' and (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
            detsql     += ' and (sce.visit = %s)' % (dataId['visit'])
            detsql     += ' and (sce.raftName = "%s")' % (re.sub("R:", "", raftId))
            detsql     += ' and (sce.ccdName = "%s")' % (ccdId[-3:])
            detresults  = data.dbInterface.execute(detsql)
            self.detData.set(raftId, ccdId, num.array([x[0] for x in detresults]))

            hist     = num.histogram(self.detData.get(raftId, ccdId), bins = self.bins)
            maxIdx   = num.argsort(hist[0])[-1]
            maxDepth = 0.5 * (hist[1][maxIdx] + hist[1][maxIdx+1]) # len(hist[1]) = len(hist[0]) + 1
            self.depth.set(raftId, ccdId, maxDepth)

            areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
            label = "photometric depth "+ areaLabel
            comment = "magnitude bin with most number of detections"
            test = testCode.Test(label, maxDepth, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)
            
    def plot(self, data, dataId, showUndefined=False):
        testSet = self.getTestSet(data, dataId)

        # fpa figure
        depthFig = qaFig.FpaQaFigure(data.cameraInfo.camera)
        for raft, ccdDict in depthFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.depth.get(raft, ccd) is None:
                    depth = self.depth.get(raft, ccd)
                    depthFig.data[raft][ccd] = depth
                    depthFig.map[raft][ccd] = 'mag=%.2f'%(depth) 

        blue = '#0000ff'
        red  = '#ff0000'
        depthFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=self.limits,
                            title="Photometric Depth", cmapOver=red, cmapUnder=blue, failLimits=self.limits)
        testSet.addFigure(depthFig, "completenessDepth.png", "Estimate of photometric depth", 
                          navMap=True)

        # Each CCD
        for raft, ccd in self.detData.raftCcdKeys():
            catData = self.catData.get(raft, ccd)
            detData = self.detData.get(raft, ccd)

            print "Plotting", ccd

            # Just to get the histogram results
            fig = qaFig.QaFigure()            
            sp1 = fig.fig.add_subplot(111)
            nc, bc, pc = sp1.hist(catData, bins=self.bins, alpha=0.5, label = 'Catalog Stars')
            nd, bd, pd = sp1.hist(detData, bins=self.bins, alpha=0.5, label = 'Detections')
            sp1.legend(numpoints=1, prop=FontProperties(size='small'), loc = 'upper left')

            label = data.cameraInfo.getDetectorName(raft, ccd)
            sp1.set_title('%s' % (label), fontsize = 12)
            testSet.addFigure(fig, "completeness.png", "Photometric detections "+label, areaLabel=label)
