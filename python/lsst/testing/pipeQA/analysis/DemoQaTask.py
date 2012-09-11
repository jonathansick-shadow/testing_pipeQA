import numpy as num

import lsst.pex.config as pexConfig 
import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.source as pqaSource
import RaftCcdData as raftCcdData
from .QaAnalysisTask import QaAnalysisTask
import lsst.testing.pipeQA.figures as qaFig

class DemoQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str, doc = "Cameras to run VignettingQaTask", default = ("lsstSim", "cfht", "suprimecam", "hscSim"))
    minMedian = pexConfig.Field(dtype = float, doc = "Minimum median magnitude", default = 20.)
    maxMedian = pexConfig.Field(dtype = float, doc = "Maximum median magnitude", default = 26.)

class DemoQaTask(QaAnalysisTask):
    ConfigClass = DemoQaConfig
    _DefaultName = "demoQa"
    def __init__(self, **kwargs): 
        QaAnalysisTask.__init__(self, **kwargs)
        self.limits = [self.config.minMedian, self.config.maxMedian]

        self.sCatDummy = pqaSource.Catalog()
        self.srefCatDummy = pqaSource.RefCatalog()

        self.description = """
         Demo pipeQA module.  Calculates median magnitude in chip.
        """

    def free(self):
        del self.detector
        del self.filter
        del self.ssDict

        del self.mags
        del self.medianMag

    def test(self, data, dataId):
        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})

        self.detector     = data.getDetectorBySensor(dataId)
        self.filter       = data.getFilterBySensor(dataId)
        self.ssDict       = data.getSourceSetBySensor(dataId)
        
        self.mags         = raftCcdData.RaftCcdVector(self.detector)
        self.medianMag    = raftCcdData.RaftCcdData(self.detector)
        
        for key, ss in self.ssDict.items():
            if self.detector.has_key(key):
                raft = self.detector[key].getParent().getId().getName()
                ccd  = self.detector[key].getId().getName()
            else:
                continue

            for s in ss:
                flux = s.getD(self.sCatDummy.PsfFluxKey)
                if flux < 0.0:
                    continue
                mag  = -2.5*num.log10(flux)
                if num.isfinite(mag):
                    self.mags.append(raft, ccd, mag)

            mags = self.mags.get(raft, ccd)
            if len(mags) > 0:
                med = num.median(mags)
            else:
                med = 0.0
            self.medianMag.set(raft, ccd, med)

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "AAA median mag"
            comment = "BBB median magnitude in chip"
            test = testCode.Test(label, med, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)

    def plot(self, data, dataId, showUndefined = False):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True
        
        # fpa figures
        medFigbase = "medianMag"
        medFigData, medFigMap = testSet.unpickle(medFigbase, [None, None])
        medFig = qaFig.FpaQaFigure(data.cameraInfo, data=medFigData, map=medFigMap) 
        for raft, ccdDict in medFig.data.items():
            for ccd, value in ccdDict.items():
                med = self.medianMag.get(raft, ccd)
                if not med is None:
                    medFig.data[raft][ccd] = med
                    if num.isfinite(med):
                        medFig.map[raft][ccd] = 'med=%.2f'%(med)
                    else:
                        medFig.map[raft][ccd] = 'med=nan'
        testSet.pickle(medFigbase, [medFig.data, medFig.map]) #cache

        blue = '#0000ff'
        red  = '#ff0000'
        if not self.delaySummary or isFinalDataId:
            self.log.log(self.log.INFO, "plotting FPA figures")
            medFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=self.limits,
                              title="Median mag", cmapOver=red, cmapUnder=blue,
                              failLimits=self.limits)
            testSet.addFigure(medFig, medFigbase+".png", "CCC Median magnitude of stars", navMap=True)
            del medFig

        # Per CCD figures
        cacheLabel = "mediam_mag" #cache
        shelfData = {}
        
        for raft, ccd in self.mags.raftCcdKeys():
            mags = self.mags.get(raft, ccd)
            fig = qaFig.QaFigure(size=(4.0,4.0))
            sp1 = fig.fig.add_subplot(111)
            sp1.hist(mags, bins = num.arange(self.limits[0], self.limits[1], 0.25))
            sp1.axvline(x = self.medianMag.get(raft, ccd), color = "k", ls = "--")
            sp1.set_xlabel("Mag")
            sp1.set_ylabel("N")
            sp1.set_title("%s" % (ccd))
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            testSet.addFigure(fig, "demo_medmag.png", "DDD Median magnitude "+areaLabel, areaLabel=areaLabel)
            shelfData[ccd] = [mags]
            del fig

        if self.useCache:
            testSet.shelve(cacheLabel, shelfData)

        if not self.delaySummary or isFinalDataId:
            self.log.log(self.log.INFO, "plotting Summary figure")
            # unstash the values
            if self.useCache:
                shelfData = testSet.unshelve(cacheLabel)
            
            magsAll = num.array(())
            for k,v in shelfData.items():
                mags = v
                magsAll = num.append(magsAll, mags)
                fig = qaFig.QaFigure(size=(4.0,4.0))
                sp1 = fig.fig.add_subplot(111)
                sp1.hist(mags, bins = num.arange(self.limits[0], self.limits[1], 0.25))
                sp1.axvline(x = num.median(mags), color = "k", ls = "--")
                sp1.set_xlabel("Mag")
                sp1.set_ylabel("N")
                sp1.set_title("All")
                del mags
                label = "all"
                testSet.addFigure(fig, "demo_medmag.png", "EEE Median magnitude "+label, areaLabel=label)


