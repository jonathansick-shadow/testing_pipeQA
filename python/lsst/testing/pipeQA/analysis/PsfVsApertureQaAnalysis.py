import sys, os
import lsst.meas.algorithms        as measAlg
import lsst.testing.pipeQA.figures as qaFig
import numpy

import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.font_manager as fm

class PsfVsApertureQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
	qaAna.QaAnalysis.__init__(self)

    def test(self, data, dataId):
	
	# get data
	self.ssDict        = data.getSourceSetBySensor(dataId)
	self.detector      = data.getDetectorBySensor(dataId)
	
	self.diff = raftCcdData.RaftCcdVector(self.detector)
	self.mag  = raftCcdData.RaftCcdVector(self.detector)
	
	for key, ss in self.ssDict.items():
	    raft = self.detector[key].getParent().getId().getName()
	    ccd  = self.detector[key].getId().getName()
	    
	    qaAnaUtil.isStar(ss)  # sets the 'STAR' flag
	    for s in ss:
		psfFlux = s.getPsfFlux()
		apFlux  = s.getApFlux()
		psfMag  = -2.5*numpy.log10(psfFlux)
		apMag   = -2.5*numpy.log10(apFlux)
		if s.getFlagForDetection() & measAlg.Flags.STAR:
		    self.diff.append(raft, ccd, psfMag - apMag)
		    self.mag.append(raft, ccd, apMag)
		    
	group = dataId['visit']
	testSet = self.getTestSet(group)

	nLowest = 100
	
	# get the mean difference and add tests
	self.means = raftCcdData.RaftCcdData(self.detector)
	for raft, ccd, mean in self.diff.listKeysAndValues('mean', nLowest=nLowest):
	    self.means.set(raft, ccd, mean)
	    label = "mean psf_vs_aper " + ccd
	    comment = "mean psf-ap magnitude difference (nstar = %d)" % (len(self.diff.get(raft, ccd)))
	    testSet.addTest( testCode.Test(label, mean, [-0.02, 0.02], comment) )

	# get the mean difference and add tests
	self.medians = raftCcdData.RaftCcdData(self.detector)
	for raft, ccd, median in self.diff.listKeysAndValues('median', nLowest=nLowest):
	    self.medians.set(raft, ccd, median)
	    label = "median psf_vs_aper " + ccd
	    comment = "median psf-ap magnitude difference (nstar = %d)" % (len(self.diff.get(raft, ccd)))
	    testSet.addTest( testCode.Test(label, median, [-0.02, 0.02], comment) )
	    
	# get the stdev and add tests
	self.stds  = raftCcdData.RaftCcdData(self.detector)
	for raft, ccd, std in self.diff.listKeysAndValues('std', nLowest=nLowest):
	    self.stds.set(raft, ccd, std)
	    label = "stdev psf_vs_aper " + ccd
	    comment = "stdev of psf-ap magnitude difference (nstar = %d)" % (len(self.diff.get(raft, ccd)))
	    testSet.addTest( testCode.Test(label, std, [0.0, 0.1], comment) )

    def plot(self, data, dataId, showUndefined=False):

	group = dataId['visit']
	testSet = self.getTestSet(group)

	# fpa figure
	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	for raft, ccdDict in fig.data.items():
	    for ccd, value in ccdDict.items():
		fig.data[raft][ccd] = self.means.get(raft, ccd)

	fig.makeFigure(showUndefined=showUndefined)
	testSet.addFigure(fig, "psfMinusAperture.png", "Psf - aper magnitude differences.")
	

	# dmag vs mag
	fig = qaFig.QaFig()
	
	ax = fig.fig.add_subplot(111)

	nKeys = len(self.mag.raftCcdKeys())
	norm = colors.Normalize(vmin=0, vmax=nKeys)
	sm = cm.ScalarMappable(norm, cmap=cm.jet)

	i = 0
	xmin, xmax = 0.0, 0.0
	for raft, ccd in self.mag.raftCcdKeys():
	    mag  = self.mag.get(raft, ccd)
	    diff = self.diff.get(raft, ccd)
	    min, max = mag.min(), mag.max()
	    if min < xmin: xmin = min
	    if max < xmax: xmax = max
	    
	    size = 1.0
	    color = sm.to_rgba(i)
	    ax.scatter(mag, diff, size, color=color, label=ccd)
	    i += 1
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, 0.75*box.width, box.height])
	fp = fm.FontProperties(size="small")
	ax.legend(prop=fp, loc=(1.05, 0.0))
	ax.set_xlim([xmin, xmax])
	ax.set_ylim([-2.0, 2.0])

	testSet.addFigure(fig.fig, "diff_psf-ap.png", "psf-ap vs. ap")
