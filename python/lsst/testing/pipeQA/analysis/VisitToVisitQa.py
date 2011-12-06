import sys, os, re
import lsst.meas.algorithms         as measAlg
import lsst.testing.pipeQA.figures  as qaFig
import numpy                        as num

import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

import matplotlib.cm as cm
import matplotlib.colors as colors

class VisitToVisitQaAnalysis(qaAna.QaAnalysis):

    def __init__(self, matchDatabase, matchVisit, mType, magCut, 
                 deltaMin, deltaMax, rmsMax, 
                 **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.matchDatabase = matchDatabase
        self.matchVisit    = matchVisit
        self.magCut        = magCut
        self.deltaLimits   = [deltaMin, deltaMax]
        self.rmsLimits     = [0.0, rmsMax]
        
        self.xlim          = [14.0, 25.0]

        self.description = """

         For each CCD in the reference visit, its sky boundary in the
         comparison visit is queried for sources.  In the case of the
         the filters being the same, magnitude v. magnitude
         comparisons are made.  If the filters are different, color
         v. magnitude plots are made, and color comparisons made to
         the input catalog.

        """
        def magType(mType):
            if re.search("(psf|PSF)", mType):
                return "psf"
            elif re.search("^ap", mType):
                return "ap"
            elif re.search("^mod", mType):
                return "mod"
            elif re.search("^cat", mType):
                return "cat"
            elif re.search("^inst", mType):
                return "inst"
        self.magType = magType(mType)

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

    def _getFluxErr(self, mType, s, sref):

        # if the source isn't valid, return NaN
        if not hasattr(s, 'getId') or not hasattr(sref, 'getId'):
            return numpy.NaN
        
        if mType=="psf":
            return s.getPsfFluxErr()
        elif mType=="ap":
            return s.getApFluxErr()
        elif mType=="mod":
            return s.getModelFluxErr()
        elif mType=="cat":
            return 0.0
        elif mType=="inst":
            return s.getInstFluxErr()

    def free(self):
        del self.matchListDictSrc
        del self.detector 
        del self.filt 
        del self.visitMatches   

        del self.refMag   
        del self.refColor 
        del self.mag1     
        del self.mag2     
        del self.magErr1  
        del self.magErr2  
        del self.star

    def testSameFilter(self, data, dataId):
        # get data before you switch the plotting options since you need to fill filterCache
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)
        self.filt             = data.getFilterBySensor(dataId)
        self.visitMatches     = data.getVisitMatchesBySensor(self.matchDatabase, self.matchVisit, dataId)
        
        km  = self.matchListDictSrc.keys()
        kv  = self.visitMatches.keys()
        if len(km) == 0:
            self.filtm0 = None
            self.filtv0 = None
            return False
        elif len(kv) == 0:
            km0 = km[0]
            self.filtm0 = self.filt[km0].getName()
            self.filtv0 = self.filtm0
            return False
        else:
            km0 = km[0]
            self.filtm0 = self.filt[km0].getName()

            kv0 = kv[0]
            if len(self.visitMatches[kv0]) == 0:
                self.filtv0 = self.filtm0
                return False
            else:
                self.filtv0 = self.visitMatches[kv0][0][2].getName()
                return self.filtm0 == self.filtv0

    def test(self, data, dataId):
        if self.testSameFilter(data, dataId):
            # magnitude comparisons
            self.testM(data, dataId)
        else:
            # color comparisons
            #self.testC(data, dataId)
            self.testM(data, dataId)


    def testM(self, data, dataId):
        # create containers for data we're interested in
        self.refMag           = raftCcdData.RaftCcdVector(self.detector)
        self.refColor         = raftCcdData.RaftCcdVector(self.detector)
        self.mag1             = raftCcdData.RaftCcdVector(self.detector)
        self.mag2             = raftCcdData.RaftCcdVector(self.detector)
        self.magErr1          = raftCcdData.RaftCcdVector(self.detector)
        self.magErr2          = raftCcdData.RaftCcdVector(self.detector)
        self.star             = raftCcdData.RaftCcdVector(self.detector)

        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER | measAlg.Flags.EDGE

        for key in self.matchListDictSrc.keys():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()
            filt = self.filt[key].getName()

            srcMatchList = self.matchListDictSrc[key]['matched']
            visMatchList = self.visitMatches[key]
            
            # List of reference object ids
            srcObjIds = num.array([x[0].getId() for x in srcMatchList])
            matObjIds = num.array([x[0].getId() for x in visMatchList])


            for i in range(len(srcObjIds)):
                srcObjId = srcObjIds[i]                   # ref id for this visit
                idx = num.where(matObjIds == srcObjId)[0] # ref id for other visit

                # only take 1-to-1 matches
                if len(idx) != 1:
                    continue

                sref1 = srcMatchList[i][0]
                srcv1 = srcMatchList[i][1]

                sref2 = visMatchList[idx[0]][0]
                srcv2 = visMatchList[idx[0]][1]

                f1  = self._getFlux(self.magType, srcv1, sref1)
                f2  = self._getFlux(self.magType, srcv2, sref2)
                df1 = self._getFluxErr(self.magType, srcv1, sref1)
                df2 = self._getFluxErr(self.magType, srcv2, sref2)

                fref1 = self._getFlux("cat", srcv1, sref1)
                fref2 = self._getFlux("cat", srcv2, sref2)

                flags1 = srcv1.getFlagForDetection()
                flags2 = srcv2.getFlagForDetection()
                
                if (f1 > 0.0 and f2 > 0.0) and (not flags1 & badFlags) and (not flags2 & badFlags):
                    m1  = -2.5 * num.log10(f1) 
                    m2  = -2.5 * num.log10(f2) 
                    dm1 = 2.5 / num.log(10.0) * df1 / f1
                    dm2 = 2.5 / num.log(10.0) * df2 / f2

                    M1  = -2.5 * num.log10(fref1) 
                    M2  = -2.5 * num.log10(fref2) 

                    star1 = flags1 & measAlg.Flags.STAR
                    
                    if num.isfinite(m1) and num.isfinite(m2) and num.isfinite(M1) and num.isfinite(M2):
                        self.refMag.append(raft, ccd, M1)
                        self.refColor.append(raft, ccd, M1-M2)
                        self.mag1.append(raft, ccd, m1)
                        self.mag2.append(raft, ccd, m2)
                        self.magErr1.append(raft, ccd, dm1)
                        self.magErr2.append(raft, ccd, dm2)
                        self.star.append(raft, ccd, star1)

        testSet = self.getTestSet(data, dataId, label=self.magType)
        testSet.addMetadata('magType', self.magType)
        testSet.addMetadata({"Description": self.description})
        self.means   = raftCcdData.RaftCcdData(self.detector)
        self.medians = raftCcdData.RaftCcdData(self.detector)
        self.stds    = raftCcdData.RaftCcdData(self.detector)
        self.derrs   = raftCcdData.RaftCcdData(self.detector)        

        for raft, ccd in self.mag1.raftCcdKeys():
            m1  = self.mag1.get(raft, ccd)
            m2  = self.mag2.get(raft, ccd)
            dm1 = self.magErr1.get(raft, ccd)
            dm2 = self.magErr2.get(raft, ccd)
            M   = self.refMag.get(raft, ccd)
            C   = self.refColor.get(raft, ccd)
            star= self.star.get(raft, ccd)

            idxS = num.where( (M > 10) & (M < self.magCut) & (star > 0) )[0]
            idxG = num.where( (M > 10) & (M < self.magCut) & (star == 0) )[0]
            
            dmS  = m1[idxS] - m2[idxS] - C[idxS]
            ddmS = num.sqrt(dm1[idxS]**2 + dm2[idxS]**2)

            dmG  = m1[idxG] - m2[idxG] - C[idxG]
            ddmG = num.sqrt(dm1[idxG]**2 + dm2[idxG]**2)

            if len(dmS) > 0:
                stat   = afwMath.makeStatistics(dmS, afwMath.NPOINT | afwMath.MEANCLIP |
                                                afwMath.STDEVCLIP | afwMath.MEDIAN)
                mean   = stat.getValue(afwMath.MEANCLIP)
                median = stat.getValue(afwMath.MEDIAN)
                std    = stat.getValue(afwMath.STDEVCLIP)
                npts   = stat.getValue(afwMath.NPOINT)

                # Common
                tag = self.magType
                areaLabel = data.cameraInfo.getDetectorName(raft, ccd)

                # MEAN
                self.means.set(raft, ccd, mean)
                label = "mean "+tag 
                comment = "mean "+tag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmS), npts)
                testSet.addTest( testCode.Test(label, mean, self.deltaLimits, comment, areaLabel=areaLabel))

                # MEDIAN
                self.medians.set(raft, ccd, median)
                label = "median "+tag 
                comment = "median "+tag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmS), npts)
                testSet.addTest( testCode.Test(label, median, self.deltaLimits, comment, areaLabel=areaLabel))

                # STD
                self.stds.set(raft, ccd, std)
                label = "stdev "+tag 
                comment = "stdev "+tag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmS), npts)
                testSet.addTest( testCode.Test(label, std, self.rmsLimits, comment, areaLabel=areaLabel))
                


    def plot(self, data, dataId, showUndefined=False):
        if self.testSameFilter(data, dataId):
            # magnitude comparisons
            self.plotM(data, dataId, showUndefined=False) 
        else:
            # color comparisons
            #self.plotC(data, dataId, showUndefined=False)
            self.plotM(data, dataId, showUndefined=False)


    def plotM(self, data, dataId, showUndefined=False):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        
        isFinalDataId = False
        if len(data.brokenDataIdList) > 0 and data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        # fpa figure
        tag  = "$\Delta m_{\mathrm{"+self.magType.upper()+"}}$"
        dtag = "d"+self.magType
        wtag = self.magType
        meanFilebase = "mean" + wtag
        stdFilebase  = "std" + wtag
        meanData, meanMap   = testSet.unpickle(meanFilebase, default=[None, None])
        stdData, stdMap     = testSet.unpickle(stdFilebase, default=[None, None])

        meanFig  = qaFig.FpaQaFigure(data.cameraInfo, data=meanData, map=meanMap)
        stdFig   = qaFig.FpaQaFigure(data.cameraInfo, data=stdData, map=stdMap)

        for raft, ccd in self.mag1.raftCcdKeys():
            meanFig.data[raft][ccd] = self.means.get(raft, ccd)
            stdFig.data[raft][ccd] = self.stds.get(raft, ccd)
            
        testSet.pickle(meanFilebase, [meanFig.data, meanFig.map])
        testSet.pickle(stdFilebase, [stdFig.data, stdFig.map])

        blue, red = '#0000ff', '#ff0000'
        if not self.delaySummary or isFinalDataId:
            print "plotting FPAs"
            meanFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[-0.03, 0.03],
                               title="Mean "+tag, cmapOver=red, cmapUnder=blue, failLimits=self.deltaLimits)
            testSet.addFigure(meanFig, "f01"+meanFilebase+".png",
                              "mean "+dtag+" mag   (brighter than %.1f)" % (self.magCut), navMap=True)
            del meanFig
            
            stdFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 0.03],
                              title="Stdev "+tag, cmapOver=red, failLimits=self.rmsLimits)
            testSet.addFigure(stdFig, "f02"+stdFilebase+".png",
                              "stdev "+dtag+" mag  (brighter than %.1f)" % (self.magCut), navMap=True)
            del stdFig

        else:
            del meanFig
            del stdFig


        # Per CCD figures
        figbase = "vvdiff_"+dtag
        shelfData = {}
        for raft, ccd in self.mag1.raftCcdKeys():
            m1  = self.mag1.get(raft, ccd)
            m2  = self.mag2.get(raft, ccd)
            dm1 = self.magErr1.get(raft, ccd)
            dm2 = self.magErr2.get(raft, ccd)
            M   = self.refMag.get(raft, ccd)
            C   = self.refColor.get(raft, ccd)
            star= self.star.get(raft, ccd)

            idxSB = num.where( (M > 10) & (M < self.magCut) & (star > 0) )[0]
            idxSF = num.where( (M > 10) & (M >= self.magCut) & (star > 0) )[0]

            dmSB  = m1[idxSB] - m2[idxSB] - C[idxSB]
            ddmSB = num.sqrt(dm1[idxSB]**2 + dm2[idxSB]**2)

            dmSF  = m1[idxSF] - m2[idxSF] - C[idxSF]
            ddmSF = num.sqrt(dm1[idxSF]**2 + dm2[idxSF]**2)

            fig = qaFig.QaFigure()
            sp1 = fig.fig.add_subplot(111)
            if len(dmSB):
                sp1.errorbar(M[idxSB], dmSB, yerr=ddmSB, fmt='ro', ms=3)
            if len(dmSF):
                sp1.errorbar(M[idxSF], dmSF, yerr=ddmSF, fmt='ko', ms=3)

            sp1.set_xlabel("${\mathrm{M_{%s; cat}}}$" % self.filtm0, fontsize=12)
            if self.filtm0 == self.filtv0:
                lab = "${\mathrm{(%s_{%s} - %s_{%s})_{%s} - (%s - %s)_{cat}}}$" % (self.filtm0,
                                                                                   dataId['visit'],
                                                                                   self.filtv0,
                                                                                   self.matchVisit,
                                                                                   self.magType,
                                                                                   self.filtm0,
                                                                                   self.filtv0)
            else:
                lab = "${\mathrm{(%s_{%s} - %s_{%s})_{%s} - (%s - %s)_{cat}}}$" % (self.filtm0,
                                                                                   dataId['visit'],
                                                                                   self.filtv0,
                                                                                   self.matchVisit,
                                                                                   self.magType,
                                                                                   self.filtm0,
                                                                                   self.filtv0)
            sp1.set_title(lab, fontsize=12)
            qaFigUtils.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize=8)
            sp1.set_ylim(-0.5, 0.5)
            sp1.set_xlim(self.xlim[0], self.xlim[1])

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            statBlurb = "Points used for statistics shown in red."
            testSet.addFigure(fig, figbase+".png",
                              dtag+" vs. "+self.magType + "(stars). "+statBlurb,
                              areaLabel=areaLabel)
