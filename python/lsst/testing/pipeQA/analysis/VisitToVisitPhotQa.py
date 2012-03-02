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
from matplotlib.font_manager import FontProperties

class VisitToVisitPhotQaAnalysis(qaAna.QaAnalysis):

    def __init__(self, database, visits, mType, magCut, 
                 deltaMin, deltaMax, rmsMax, 
                 **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.database      = database
        self.visits        = visits
        self.magCut        = magCut
        self.deltaLimits   = [deltaMin, deltaMax]
        self.rmsLimits     = [0.0, rmsMax]
        
        self.maglim        = [14.0, 25.0]
        self.colorlim      = [-0.98, 2.48]

        self.alloc()
        
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

    def alloc(self):
        # Filter name
        self.ownFilt       = None

        # Data caches
        self.visitFilters  = {}
        self.visitMatches  = {} 

        # Inputs to tests
        self.mag         = {}
        self.magErr      = {}
        self.refMag      = {}
        self.visitMag    = {}
        self.visitMagErr = {}
        self.visitRefMag = {}
        self.refId       = {}
        self.star        = {}

        # Results of tests
        self.meanDmags   = {}
        self.medianDmags = {}
        self.stdDmags    = {}
        
    def free(self):
        del self.matchListDictSrc
        del self.detector 
        del self.visitMatches   
        del self.visitFilters

        del self.mag
        del self.magErr
        del self.refMag   
        del self.visitMag
        del self.visitMagErr 
        del self.visitRefMag
        del self.refId
        del self.star

        # realloc what needs to be reused
        self.alloc()

    def test(self, data, dataId):
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)

        if type(self.ownFilt) == type(None):
            filt = data.getFilterBySensor(dataId)
            self.ownFilt = filt.values()[0]

        if self.visits == None:
            visitList = [dataId['visit'],]
        else:
            visitList = self.visits

        for visit in visitList:
            self.visitMatches[visit] = data.getVisitMatchesBySensor(self.database, visit, dataId)
            kvs = self.visitMatches[visit].keys()
            self.visitFilters[visit] = None
            if len(kvs) > 0:
                for kv in kvs:
                    if len(self.visitMatches[visit][kv]) > 0:
                        self.visitFilters[visit] = self.visitMatches[visit][kv][0][2] 
                        break

            # create containers for data we're interested in
            self.mag[visit]         = raftCcdData.RaftCcdVector(self.detector)
            self.magErr[visit]      = raftCcdData.RaftCcdVector(self.detector)
            self.refMag[visit]      = raftCcdData.RaftCcdVector(self.detector)
            self.visitMag[visit]    = raftCcdData.RaftCcdVector(self.detector)
            self.visitMagErr[visit] = raftCcdData.RaftCcdVector(self.detector)
            self.visitRefMag[visit] = raftCcdData.RaftCcdVector(self.detector)
            self.refId[visit]       = raftCcdData.RaftCcdVector(self.detector)
            self.star[visit]        = raftCcdData.RaftCcdVector(self.detector)
    
            badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER | measAlg.Flags.EDGE

            for key in self.matchListDictSrc.keys():
                raft = self.detector[key].getParent().getId().getName()
                ccd  = self.detector[key].getId().getName()
    
                # All detections matched to ref objects
                srcMatchList = self.matchListDictSrc[key]['matched']
                # All visit dets in footprint of original image
                visMatchList = self.visitMatches[visit][key]
                
                # List of reference object ids
                srcObjIds = num.array([x[0].getId() for x in srcMatchList])
                visObjIds = num.array([x[0].getId() for x in visMatchList])

                # Iterate over all object ids
                for i in range(len(srcObjIds)):
                    srcObjId = srcObjIds[i]                   # ref id for this visit
                    idx = num.where(visObjIds == srcObjId)[0] # ref id for other visit
    
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
                            self.mag[visit].append(raft, ccd, m1)
                            self.magErr[visit].append(raft, ccd, dm1)

                            self.visitMag[visit].append(raft, ccd, m2)
                            self.visitMagErr[visit].append(raft, ccd, dm2)

                            self.refMag[visit].append(raft, ccd, M1)
                            self.visitRefMag[visit].append(raft, ccd, M2)

                            self.refId[visit].append(raft, ccd, srcObjId)
                            self.star[visit].append(raft, ccd, star1)

            # TMI
            #testLabel = "%s_%s_%s" % (self.database, visit, self.magType)

            testLabel = "%s_%s" % (visit, self.magType)
            testSet = self.getTestSet(data, dataId, label=testLabel)
            testSet.addMetadata("magType", "%s mags from %s" % (self.magType, self.database))
            testSet.addMetadata({"Description": self.description})

            self.meanDmags[visit]   = raftCcdData.RaftCcdData(self.detector)
            self.medianDmags[visit] = raftCcdData.RaftCcdData(self.detector)
            self.stdDmags[visit]    = raftCcdData.RaftCcdData(self.detector)
    
            for raft, ccd in self.mag[visit].raftCcdKeys():
                m1  = self.mag[visit].get(raft, ccd)
                m2  = self.visitMag[visit].get(raft, ccd)
                dm1 = self.magErr[visit].get(raft, ccd)
                dm2 = self.visitMagErr[visit].get(raft, ccd)
                M1  = self.refMag[visit].get(raft, ccd)
                M2  = self.visitRefMag[visit].get(raft, ccd)
                star= self.star[visit].get(raft, ccd)
    
                idxS = num.where( (M1 > 10) & (M1 < self.magCut) & (star > 0) )[0]
                #idxG = num.where( (M1 > 10) & (M1 < self.magCut) & (star == 0) )[0]
                
                dmS  = m1[idxS] - m2[idxS] - (M1[idxS] - M2[idxS])
                ddmS = num.sqrt(dm1[idxS]**2 + dm2[idxS]**2)
    
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
                    self.meanDmags[visit].set(raft, ccd, mean)
                    label = "mean "+tag 
                    comment = "mean "+tag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmS), npts)
                    testSet.addTest( testCode.Test(label, mean, self.deltaLimits, comment, areaLabel=areaLabel))
    
                    # MEDIAN
                    self.medianDmags[visit].set(raft, ccd, median)
                    label = "median "+tag 
                    comment = "median "+tag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmS), npts)
                    testSet.addTest( testCode.Test(label, median, self.deltaLimits, comment, areaLabel=areaLabel))
    
                    # STD
                    self.stdDmags[visit].set(raft, ccd, std)
                    label = "stdev "+tag 
                    comment = "stdev "+tag+" (mag lt %.1f, nstar/clip=%d/%d)" % (self.magCut, len(dmS), npts)
                    testSet.addTest( testCode.Test(label, std, self.rmsLimits, comment, areaLabel=areaLabel))
                


    def plot(self, data, dataId, showUndefined=False):
        if self.visits == None:
            visitList = [dataId['visit'],]
        else:
            visitList = self.visits

        # In all cases plot a delta-mag plot
        for i in range(len(visitList)):
            visiti = visitList[i]
            self.visitMatches[visiti] = data.getVisitMatchesBySensor(self.database, visiti, dataId)
            self.plotdM(data, dataId, visiti, showUndefined)
            
            # In some cases plot color-magnitude diagram
            if self.ownFilt.getName() != self.visitFilters[visiti].getName():
                self.plotCmd(data, dataId, visiti, showUndefined)

            # In some cases plot color-color diagrams
            for j in range(i+1, len(visitList)):
                visitj = visitList[j]
                self.visitMatches[visitj] = data.getVisitMatchesBySensor(self.database, visitj, dataId)
                if self.ownFilt.getName() != self.visitFilters[visiti].getName() and \
                        self.visitFilters[visiti].getName() != self.visitFilters[visitj].getName():
                    self.plotCcd(data, dataId, visiti, visitj, showUndefined)

            
    def plotdM(self, data, dataId, visit, showUndefined=False):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        
        isFinalDataId = False
        if len(data.brokenDataIdList) > 0 and data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        # fpa figure
        title        = "$\Delta m_{\mathrm{"+self.magType.upper()+"}}$"
        cachetag     = "%s_%s_%s" % (self.database, visit, self.magType)
        nametag      = "d"+self.magType
        meanFilebase = "mean" + cachetag
        stdFilebase  = "std" + cachetag
        meanData, meanMap   = testSet.unpickle(meanFilebase, default=[None, None])
        stdData, stdMap     = testSet.unpickle(stdFilebase, default=[None, None])

        meanFig  = qaFig.FpaQaFigure(data.cameraInfo, data=meanData, map=meanMap)
        stdFig   = qaFig.FpaQaFigure(data.cameraInfo, data=stdData, map=stdMap)

        for raft, ccd in self.mag[visit].raftCcdKeys():
            meanFig.data[raft][ccd] = self.meanDmags[visit].get(raft, ccd)
            stdFig.data[raft][ccd] = self.stdDmags[visit].get(raft, ccd)
            meanFig.map[raft][ccd] = "mean=%.4f" % (self.meanDmags[visit].get(raft, ccd))
            stdFig.map[raft][ccd] = "std=%.4f" % (self.stdDmags[visit].get(raft, ccd))
            
        testSet.pickle(meanFilebase, [meanFig.data, meanFig.map])
        testSet.pickle(stdFilebase, [stdFig.data, stdFig.map])

        blue, red = '#0000ff', '#ff0000'
        if not self.delaySummary or isFinalDataId:
            print "plotting FPAs"
            meanFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[-0.03, 0.03],
                               title="Mean "+title, cmapOver=red, cmapUnder=blue, failLimits=self.deltaLimits)
            testSet.addFigure(meanFig, "f01"+meanFilebase+".png",
                              "mean "+nametag+" mag   (brighter than %.1f)" % (self.magCut), navMap=True)
            del meanFig
            
            stdFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 0.03],
                              title="Stdev "+title, cmapOver=red, failLimits=self.rmsLimits)
            testSet.addFigure(stdFig, "f02"+stdFilebase+".png",
                              "stdev "+nametag+" mag  (brighter than %.1f)" % (self.magCut), navMap=True)
            del stdFig

        else:
            del meanFig
            del stdFig


        # Per CCD figures
        figbase = "vvdiff_"+nametag
        toggle  = "Diff_"+visit
        shelfData = {}
        for raft, ccd in self.mag[visit].raftCcdKeys():
            m1  = self.mag[visit].get(raft, ccd)
            m2  = self.visitMag[visit].get(raft, ccd)
            dm1 = self.magErr[visit].get(raft, ccd)
            dm2 = self.visitMagErr[visit].get(raft, ccd)
            M1  = self.refMag[visit].get(raft, ccd)
            M2  = self.visitRefMag[visit].get(raft, ccd)
            star= self.star[visit].get(raft, ccd)

            idxSB = num.where( (M1 > 10) & (M1 < self.magCut) & (star > 0) )[0]
            idxSF = num.where( (M1 > 10) & (M1 >= self.magCut) & (star > 0) )[0]

            dmSB  = m1[idxSB] - m2[idxSB] - (M1[idxSB] - M2[idxSB])
            ddmSB = num.sqrt(dm1[idxSB]**2 + dm2[idxSB]**2)

            dmSF  = m1[idxSF] - m2[idxSF] - (M1[idxSF] - M2[idxSF])
            ddmSF = num.sqrt(dm1[idxSF]**2 + dm2[idxSF]**2)

            fig = qaFig.QaFigure()
            sp1 = fig.fig.add_subplot(111)
            if len(dmSB):
                sp1.errorbar(M1[idxSB], dmSB, yerr=ddmSB, fmt='ro', ms=3)
            if len(dmSF):
                sp1.errorbar(M1[idxSF], dmSF, yerr=ddmSF, fmt='ko', ms=3)

            sp1.set_xlabel("${\mathrm{M_{%s; cat}}}$" % self.ownFilt.getName(), fontsize=12)
            if self.ownFilt.getName() == self.visitFilters[visit].getName():
                lab = "${\mathrm{(%s_{%s} - %s_{%s})_{%s}}}$" % (self.ownFilt.getName(),
                                                                 dataId['visit'],
                                                                 self.visitFilters[visit].getName(),
                                                                 visit,
                                                                 self.magType)
            else:
                lab = "${\mathrm{(%s_{%s} - %s_{%s})_{%s} - (%s - %s)_{cat}}}$" % (self.ownFilt.getName(),
                                                                                   dataId['visit'],
                                                                                   self.visitFilters[visit].getName(),
                                                                                   visit,
                                                                                   self.magType,
                                                                                   self.ownFilt.getName(),
                                                                                   self.visitFilters[visit].getName())
            sp1.set_title(lab, fontsize=12)
            qaFigUtils.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize=8)
            sp1.set_ylim(-0.5, 0.5)
            sp1.set_xlim(self.maglim[0], self.maglim[1])

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            statBlurb = "Points used for statistics shown in red."
            testSet.addFigure(fig, figbase+".png",
                              nametag+" vs. "+self.magType + "(stars). "+statBlurb,
                              areaLabel=areaLabel, toggle=toggle)

    def plotCmd(self, data, dataId, visit, showUndefined=False):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        
        isFinalDataId = False
        if len(data.brokenDataIdList) > 0 and data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        # Per CCD figures
        nametag = "cmd_"+self.magType
        figbase = "vvcmd_"+self.magType
        toggle  = "Cmd_"+visit
        shelfData = {}
        for raft, ccd in self.mag[visit].raftCcdKeys():
            m1  = self.mag[visit].get(raft, ccd)
            m2  = self.visitMag[visit].get(raft, ccd)
            dm1 = self.magErr[visit].get(raft, ccd)
            dm2 = self.visitMagErr[visit].get(raft, ccd)
            M1  = self.refMag[visit].get(raft, ccd)
            M2  = self.visitRefMag[visit].get(raft, ccd)
            star= self.star[visit].get(raft, ccd)

            idxS = num.where( (star >  0) )[0]
            idxG = num.where( (star == 0) )[0]

            order = ['u', 'g', 'r', 'i', 'z', 'y']
            idx1  = order.index(self.ownFilt.getName())
            idx2  = order.index(self.visitFilters[visit].getName())
            if (idx1 < idx2):
                firstD  = m1
                secondD = m2
                firstR  = M1
                secondR = M2
                lab1    = self.ownFilt.getName()
                lab2    = self.visitFilters[visit].getName()
                v1      = dataId['visit']
                v2      = visit
            else:
                firstD  = m2
                secondD = m1
                firstR  = M2
                secondR = M1
                lab1    = self.visitFilters[visit].getName()
                lab2    = self.ownFilt.getName()
                v1      = visit
                v2      = dataId['visit']
                
            fig = self.panelPlot(firstD[idxS]-secondD[idxS], secondD[idxS],
                                 firstD[idxG]-secondD[idxG], secondD[idxG],
                                 firstR[idxS]-secondR[idxS], secondR[idxS],
                                 firstR[idxG]-secondR[idxG], secondR[idxG],
                                 '${\mathrm{(%s - %s)_{%s}}}$' % (lab1, lab2, self.magType),
                                 '${\mathrm{%s_{%s}}}$' % (lab2, self.magType),
                                 self.colorlim[0], self.colorlim[1],
                                 self.maglim[1], self.maglim[0])
                                 
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            statBlurb = "Color magnitude diagram; filters %s vs. (%s - %s); visits %s vs. (%s - %s)" % (
                lab2, lab1, lab2, v2, v1, v2)
            testSet.addFigure(fig, figbase+".png",
                              nametag+". "+statBlurb,
                              areaLabel=areaLabel, toggle=toggle)

    def plotCcd(self, data, dataId, visitA, visitB, showUndefined=False):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        
        isFinalDataId = False
        if len(data.brokenDataIdList) > 0 and data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        # Per CCD figures
        nametag = "ccd_"+self.magType
        figbase = "vvccd_"+self.magType
        toggle  = "Ccd_"+visitA+"_"+visitB
        shelfData = {}
        for raft, ccd in self.mag[visitA].raftCcdKeys():
            idA = self.refId[visitA].get(raft, ccd)
            idB = self.refId[visitB].get(raft, ccd)
            idX = num.intersect1d(idA, idB)
            sliceA, sliceB = num.where(num.equal.outer(idA, idB) == True)
            
            m1A  = self.mag[visitA].get(raft, ccd)[sliceA]
            m2A  = self.visitMag[visitA].get(raft, ccd)[sliceA]
            dm1A = self.magErr[visitA].get(raft, ccd)[sliceA]
            dm2A = self.visitMagErr[visitA].get(raft, ccd)[sliceA]
            M1A  = self.refMag[visitA].get(raft, ccd)[sliceA]
            M2A  = self.visitRefMag[visitA].get(raft, ccd)[sliceA]
            starA= self.star[visitA].get(raft, ccd)[sliceA]

            m1B  = self.mag[visitB].get(raft, ccd)[sliceB]
            m2B  = self.visitMag[visitB].get(raft, ccd)[sliceB]
            dm1B = self.magErr[visitB].get(raft, ccd)[sliceB]
            dm2B = self.visitMagErr[visitB].get(raft, ccd)[sliceB]
            M1B  = self.refMag[visitB].get(raft, ccd)[sliceB]
            M2B  = self.visitRefMag[visitB].get(raft, ccd)[sliceB]
            starB= self.star[visitB].get(raft, ccd)[sliceB]

            idxS = num.where( (starA >  0) )[0]
            idxG = num.where( (starA == 0) )[0]

            order = ['u', 'g', 'r', 'i', 'z', 'y']
            f1    = self.ownFilt.getName()
            f2    = self.visitFilters[visitA].getName()
            f3    = self.visitFilters[visitB].getName()
            idx1  = order.index(f1)
            idx2  = order.index(f2)
            idx3  = order.index(f3)

            slist = [ [idx1, m1A, M1A, f1, dataId['visit']],
                      [idx2, m2A, M2A, f2, visitA],
                      [idx3, m2B, M2B, f3, visitB]
                      ]
            slist.sort()

            fig = self.panelPlot(slist[0][1][idxS] - slist[1][1][idxS],  # xdataS
                                 slist[1][1][idxS] - slist[2][1][idxS],  # ydataS
                                 slist[0][1][idxG] - slist[1][1][idxG],  # xdataG
                                 slist[1][1][idxG] - slist[2][1][idxG],  # ydataG
                                 slist[0][2][idxS] - slist[1][2][idxS],  # xmodelS
                                 slist[1][2][idxS] - slist[2][2][idxS],  # ymodelS
                                 slist[0][2][idxG] - slist[1][2][idxG],  # xmodelG
                                 slist[1][2][idxG] - slist[2][2][idxG],  # ymodelG
                                 '${\mathrm{(%s - %s)_{%s}}}$' % (slist[0][3], slist[1][3], self.magType), # xlabel
                                 '${\mathrm{(%s - %s)_{%s}}}$' % (slist[1][3], slist[2][3], self.magType), # ylabel
                                 self.colorlim[0], self.colorlim[1], 
                                 self.colorlim[0], self.colorlim[1])

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            statBlurb = "Color color diagram; filters (%s - %s) vs (%s - %s); visits (%s - %s) vs (%s - %s)" % (
                slist[1][3], slist[2][3], slist[0][3], slist[1][3],
                slist[1][4], slist[2][4], slist[0][4], slist[1][4])                
            testSet.addFigure(fig, figbase+".png",
                              nametag+". "+statBlurb,
                              areaLabel=areaLabel, toggle=toggle)


    def panelPlot(self, xdataS, ydataS, xdataG, ydataG, xmodelS, ymodelS, xmodelG, ymodelG, xlabel, ylabel, xmin, xmax, ymin, ymax):
        fig = qaFig.QaFigure(size = (6.5, 5.0))
        fig.fig.subplots_adjust(wspace=0.0, hspace=0.0)

        sp1 = fig.fig.add_subplot(221)
        sp2 = fig.fig.add_subplot(222, sharex=sp1, sharey=sp1)
        sp3 = fig.fig.add_subplot(223, sharex=sp1, sharey=sp1)
        sp4 = fig.fig.add_subplot(224, sharex=sp1, sharey=sp1)

        sp1.plot(xdataS,  ydataS,  'ro', ms=3, alpha=0.25)
        sp2.plot(xmodelS, ymodelS, 'ro', ms=3, alpha=0.25)
        sp3.plot(xdataG,  ydataG,  'bs', ms=3, alpha=0.25)
        sp4.plot(xmodelG, ymodelG, 'bs', ms=3, alpha=0.25)

        qaFigUtils.qaSetp(sp1.get_xticklabels()+sp2.get_xticklabels(), visible=False)
        qaFigUtils.qaSetp(sp2.get_yticklabels()+sp4.get_yticklabels(), visible=False)

        qaFigUtils.qaSetp(sp3.get_xticklabels()+sp4.get_xticklabels(), fontsize=8)
        qaFigUtils.qaSetp(sp1.get_yticklabels()+sp3.get_yticklabels(), fontsize=8)

        sp1.set_title('Data', fontsize=12)
        sp2.set_title('Ref Cat', fontsize=12)
        sp3.set_xlabel(xlabel, fontsize=12)
        sp4.set_xlabel(xlabel, fontsize=12)
        sp1.set_ylabel(ylabel, fontsize=12)
        sp3.set_ylabel(ylabel, fontsize=12)

        # right side
        twinsp2 = sp2.twinx()
        twinsp2.set_ylabel('Star', fontsize=12, rotation=-90)
        twinsp4 = sp4.twinx()
        twinsp4.set_ylabel('Gal', fontsize=12, rotation=-90)
        qaFigUtils.qaSetp(twinsp2.get_xticklabels()+twinsp2.get_yticklabels(), visible=False)
        qaFigUtils.qaSetp(twinsp4.get_xticklabels()+twinsp4.get_yticklabels(), visible=False)

        #fig.fig.text(0.025, 0.5, ylabel, fontsize=10, rotation=90)
        sp1.set_xlim(xmin, xmax)
        sp1.set_ylim(ymin, ymax)
        return fig
