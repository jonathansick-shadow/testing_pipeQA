import sys, os, re
import lsst.meas.algorithms         as measAlg
import lsst.testing.pipeQA.figures  as qaFig
import numpy                        as num

import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

import matplotlib.cm as cm
import matplotlib.colors as colors

class VisitToVisitQaAnalysis(qaAna.QaAnalysis):

    def __init__(self, matchDatabase, matchVisit, mType, **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.matchDatabase = matchDatabase
        self.matchVisit    = matchVisit

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
        del self.ssDict
        del self.matchListDictSrc
        del self.detector 
        del self.filter 
        del self.visitMatches   

    def test(self, data, dataId):
        # get data
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        self.visitMatches     = data.getVisitMatchesBySensor(self.matchDatabase, self.matchVisit, dataId)
        
        # create containers for data we're interested in
        self.refMag           = raftCcdData.RaftCcdVector(self.detector)
        self.refColor         = raftCcdData.RaftCcdVector(self.detector)
        self.mag1             = raftCcdData.RaftCcdVector(self.detector)
        self.mag2             = raftCcdData.RaftCcdVector(self.detector)
        self.dmag1            = raftCcdData.RaftCcdVector(self.detector)
        self.dmag2            = raftCcdData.RaftCcdVector(self.detector)

        for key in self.matchListDictSrc.keys():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()
            filter = self.filter[key].getName()

            srcMatchList = self.matchListDictSrc[key]['matched']
            visMatchList = self.visitMatches[key]
            
            # List of objectIds
            srcObjIds = num.array([x[0].getId() for x in srcMatchList])
            matObjIds = num.array([x[0].getId() for x in visMatchList])

            print srcObjIds
            print matObjIds

            for i in range(len(srcObjIds)):
                srcObjId = srcObjIds[i]
                idx = num.where(matObjIds = srcObjId)
                print srcObjId, len(idx)
                for j in idx:
                    ss1 = srcObjId[2]
                    ss2 = matObjIds[1][j]
                    
                    
                    f1 = self._getFlux(self.magType1, s, s)
                    f2 = self._getFlux(self.magType2, s, s)
                    df1 = self._getFluxErr(self.magType1, s, s)
                    df2 = self._getFluxErr(self.magType2, s, s)
                    flags = s.getFlagForDetection()
                    
                    if ((f1 > 0.0 and f2 > 0.0) and not flags & badFlags):
                        m1 = -2.5*numpy.log10(f1) #self.calib[key].getMagnitude(f1)
                        m2 = -2.5*numpy.log10(f2) #self.calib[key].getMagnitude(f2)
                        dm1 = 2.5 / numpy.log(10.0) * df1 / f1
                        dm2 = 2.5 / numpy.log(10.0) * df2 / f2

                        star = flags & measAlg.Flags.STAR

                    self.refMag.append(raft, ccd, )
                    self.refColor.append(raft, ccd, )
                    self.mag1.append(raft, ccd, )
                    self.mag2.append(raft, ccd, )
                    self.dmag1.append(raft, ccd, )
                    self.dmag2.append(raft, ccd, )
 
    def plot(self, data, dataId):
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        
