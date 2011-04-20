import sys, os
import lsst.meas.algorithms        as measAlg
import lsst.testing.pipeQA.figures as qaFig
import numpy

import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

class PsfVsApertureQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
	qaAna.QaAnalysis.__init__(self)

    def test(self, data, dataId):
	
	# get data
	self.ssDict        = data.getSourceSetBySensor(dataId)
	self.detector      = data.getDetectorBySensor(dataId)
	
	self.diff = raftCcdData.RaftCcdVector(self.detector)

	
	for key, ss in self.ssDict.items():
	    raft = self.detector[key].getParent().getId().getName()
	    ccd  = self.detector[key].getId().getName()
	    
	    qaAnaUtil.isStar(ss)  # sets the 'STAR' flag
	    for s in ss:
		psfFlux = s.getPsfFlux()
		apFlux  = s.getApFlux()
		psfMag  = 2.5*numpy.log10(psfFlux)
		apMag   = 2.5*numpy.log10(apFlux)
		if s.getFlagForDetection() & measAlg.Flags.STAR:
		    self.diff.append(raft, ccd, psfMag - apMag)
		    
	group = dataId['visit']
	testSet = self.getTestSet(group)

	# get the mean difference and add tests
	self.means = raftCcdData.RaftCcdData(self.detector)
	for raft, ccd, mean in self.diff.listKeysAndValues('mean'):
	    self.means.set(raft, ccd, mean)
	    label = "mean psf_vs_aper " + ccd
	    comment = "mean psf-ap magnitude difference (nstar = %d)" % (len(self.diff.get(raft, ccd)))
	    testSet.addTest( testCode.Test(label, mean, [-0.02, 0.02], comment) )

	# get the mean difference and add tests
	self.medians = raftCcdData.RaftCcdData(self.detector)
	for raft, ccd, median in self.diff.listKeysAndValues('median'):
	    self.medians.set(raft, ccd, median)
	    label = "median psf_vs_aper " + ccd
	    comment = "median psf-ap magnitude difference (nstar = %d)" % (len(self.diff.get(raft, ccd)))
	    testSet.addTest( testCode.Test(label, median, [-0.02, 0.02], comment) )
	    
	# get the stdev and add tests
	self.stds  = raftCcdData.RaftCcdData(self.detector)
	for raft, ccd, std in self.diff.listKeysAndValues('std'):
	    self.stds.set(raft, ccd, std)
	    label = "stdev psf_vs_aper " + ccd
	    comment = "stdev of psf-ap magnitude difference (nstar = %d)" % (len(self.diff.get(raft, ccd)))
	    testSet.addTest( testCode.Test(label, std, [0.0, 0.1], comment) )

    def plot(self, data, dataId, showUndefined=False):

	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	for raft, ccdDict in fig.data.items():
	    for ccd, value in ccdDict.items():
		fig.data[raft][ccd] = self.means.get(raft, ccd)

	fig.makeFigure(showUndefined=showUndefined)

	group = dataId['visit']
	testSet = self.getTestSet(group)
	testSet.addFigure(fig, "psfMinusAperture.png", "Psf - aper magnitude differences.")
	

