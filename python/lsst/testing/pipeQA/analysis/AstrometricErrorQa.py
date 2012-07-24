import sys, os, re
import lsst.meas.algorithms        as measAlg
import lsst.testing.pipeQA.figures as qaFig
import numpy

import time

import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtil

import lsst.testing.pipeQA.source as pqaSource


class AstrometricErrorQa(qaAna.QaAnalysis):

    def __init__(self, maxErr, **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.limits = [0.0, maxErr]

        self.description = """
         For each CCD, the left figure shows the distance offset between the
         measured centroid of matched objects and the catalog position of these
         objects, represented as an arrow.  The top right panel provides the
         view of these vectors stacked at the position of the reference object,
         with the green circle representing the radius that contains 50% of the
         matches.  The bottom panel provides a histogram of the astrometric
         offsets, with the median indicated.  The summary FPA figure provides
         the median vector (offset and angle) for each chip.
        """
        
    def free(self):
        del self.x
        del self.y
        del self.dDec
        del self.dRa
        del self.filter
        
        del self.detector
        del self.matchListDictSrc

        del self.medErrArcsec
        del self.medThetaRad

        del self.sCatDummy
        del self.srefCatDummy
        
    def test(self, data, dataId):

        # get data
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        
        #self.clusters = data.getSourceClusters(dataId)

        # compute the mean ra, dec for each source cluster

        self.dRa  = raftCcdData.RaftCcdVector(self.detector)
        self.dDec = raftCcdData.RaftCcdVector(self.detector)
        self.x    = raftCcdData.RaftCcdVector(self.detector)
        self.y    = raftCcdData.RaftCcdVector(self.detector)

        self.sCatDummy = pqaSource.Catalog()
        sCatDummy = self.sCatDummy.catalog
        sCatSchema = sCatDummy.getSchema()
        self.srefCatDummy  = pqaSource.RefCatalog()
        srefCatDummy = self.srefCatDummy.catalog
        srefCatSchema = srefCatDummy.getSchema()
        
        xKey      = sCatSchema.find('XAstrom').key
        yKey      = sCatSchema.find('YAstrom').key
        raKey     = sCatSchema.find('Ra').key
        decKey    = sCatSchema.find('Dec').key
        refRaKey  = srefCatSchema.find('Ra').key
        refDecKey = srefCatSchema.find('Dec').key

        filter = None
        for key in self.matchListDictSrc.keys():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()
            filter = self.filter[key].getName()

            matchList = self.matchListDictSrc[key]['matched']
            for m in matchList:
                sref, s, dist = m
                ra, dec, raRef, decRef = \
                    [numpy.radians(x) for x in [s.getD(raKey), s.getD(decKey),
                                                sref.getD(refRaKey), sref.getD(refDecKey)]]
                
                dDec = decRef - dec
                dRa  = (raRef - ra)*abs(numpy.cos(decRef))
                
                if not (s.getD(sCatSchema.find('FlagPixInterpCen').key)):
                    self.dRa.append(raft, ccd, dRa)
                    self.dDec.append(raft, ccd, dDec)
                    self.x.append(raft, ccd, s.getD(xKey))
                    self.y.append(raft, ccd, s.getD(yKey))
                    
                    
        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})

        self.medErrArcsec = raftCcdData.RaftCcdData(self.detector)
        self.medThetaRad  = raftCcdData.RaftCcdData(self.detector)

        for raft,  ccd in self.dRa.raftCcdKeys():
            dRa  = self.dRa.get(raft, ccd)
            dDec = self.dDec.get(raft, ccd)
            
            errArcsec = 206265.0*numpy.sqrt(dRa**2 + dDec**2)
            #errArcsec = 3600.0*numpy.sqrt(dRa**2 + dDec**2)
            thetaRad  = numpy.arctan2(dDec, dRa)

            if len(errArcsec) > 0:
                stat  = afwMath.makeStatistics(errArcsec, afwMath.NPOINT | afwMath.MEDIAN)
                medErrArcsec = stat.getValue(afwMath.MEDIAN)
                stat  = afwMath.makeStatistics(thetaRad, afwMath.NPOINT | afwMath.MEDIAN)
                medThetaRad = stat.getValue(afwMath.MEDIAN)
                n = stat.getValue(afwMath.NPOINT)
            else:
                medErrArcsec = -1.0
                medThetaRad = 0.0
                n = 0

            self.medErrArcsec.set(raft, ccd, medErrArcsec)
            self.medThetaRad.set(raft, ccd, medThetaRad)
            
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "median astrometry error "
            comment = "median sqrt(dRa^2+dDec^2) (arcsec, nstar=%d)" % (n)
            test = testCode.Test(label, medErrArcsec, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)

        
    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)


        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        medAstBase = "medAstError"
        medAstData, medAstMap = testSet.unpickle(medAstBase, default=[None, None])

        # fpa figure
        astFig = qaFig.VectorFpaQaFigure(data.cameraInfo, data=medAstData, map=medAstMap)

        vLen = 5000 # length in pixels for 1 arcsec error vector
        for raft, ccdDict in astFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.medErrArcsec.get(raft, ccd) is None:
                    astErrArcsec = self.medErrArcsec.get(raft, ccd)
                    thetaRad = self.medThetaRad.get(raft, ccd)
                    astFig.data[raft][ccd] = [thetaRad, vLen*astErrArcsec, astErrArcsec]
                    astFig.map[raft][ccd] = "\"/theta=%.2f/%.0f" % (astErrArcsec, numpy.degrees(thetaRad))
                
        testSet.pickle(medAstBase, [astFig.data, astFig.map])
        
        if not self.delaySummary or isFinalDataId:
            print "plotting FPAs"
            astFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 2.0*self.limits[1]],
                              title="Median astrometric error", cmapOver='#ff0000', failLimits=self.limits,
                              cmapUnder="#ff0000")
            testSet.addFigure(astFig, medAstBase+".png", "Median astrometric error",  navMap=True)
            
        del astFig

        prePlot = False

        cacheLabel = "astromError"
        shelfData = {}
        
        for raft, ccd in self.dRa.raftCcdKeys():
            ra = self.dRa.get(raft, ccd)
            dec = self.dDec.get(raft, ccd)
            eLen = 206265.0*numpy.sqrt(ra**2 + dec**2)
            #eLen = 3600.0*numpy.sqrt(ra**2 + dec**2)
            t = numpy.arctan2(dec, ra)
            
            dx = eLen*numpy.cos(t)
            w = numpy.where(numpy.abs(dx) < 10.0)
            
            dx = dx[w]
            dy = (eLen*numpy.sin(t))[w]
            x = (self.x.get(raft, ccd))[w]
            y = (self.y.get(raft, ccd))[w]

            print "plotting ", ccd


            import AstrometricErrorQaPlot
            label = data.cameraInfo.getDetectorName(raft, ccd)
            dataDict = {'x': x, 'y':y, 'dx':dx, 'dy':dy, 'gridVectors':False }
            caption = "Astrometric error" + label
            if prePlot:
                fig = AstrometricErrorQaPlot.plot(dataDict)
                testSet.addFigure(fig, "astromError.png", caption, areaLabel=label)
                del fig
            else:
                testSet.addLazyFigure(dataDict, cacheLabel+".png", caption,
                                      AstrometricErrorQaPlot, areaLabel=label, plotargs="",
                                      pythonpath=os.path.join(os.getenv('TESTING_PIPEQA_DIR'),'python', 'lsst', 'testing', 'pipeQA'))
                                      

            shelfData[label] = [x, y, dx, dy]

        if self.useCache:
            testSet.shelve(cacheLabel, shelfData)

        if not self.delaySummary or isFinalDataId:
            print "plotting Summary figure"

            # unstash the values
            if self.useCache:
                shelfData = testSet.unshelve(cacheLabel)

            xAll  = numpy.array([])
            yAll  = numpy.array([])
            dxAll = numpy.array([])
            dyAll = numpy.array([])
            for k,v in shelfData.items():
                #allCcds.append(k)
                x, y, dx, dy = v
                xAll   = numpy.append(xAll  , x )
                yAll   = numpy.append(yAll  , y )
                dxAll  = numpy.append(dxAll , dx)
                dyAll  = numpy.append(dyAll , dy)

            label = 'all'
            dataDict = {'x': xAll, 'y': yAll, "dx" : dxAll, "dy" : dyAll, 'gridVectors' : True}
            if True: #prePlot:
                allFig = AstrometricErrorQaPlot.plot(dataDict)
                testSet.addFigure(allFig, "astromError.png", "Astrometric error" + label, areaLabel=label)
                del allFig
            else:
                testSet.addLazyFigure(dataDict, cacheLabel+".png", caption,
                                      AstrometricErrorQaPlot, areaLabel=label, plotargs="",
                                      pythonpath=os.path.join(os.getenv('TESTING_PIPEQA_DIR'),'python', 'lsst', 'testing', 'pipeQA'))
                    
            del xAll, yAll, dxAll, dyAll
            
            

