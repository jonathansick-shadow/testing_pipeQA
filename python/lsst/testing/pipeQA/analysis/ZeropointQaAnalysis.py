import lsst.testing.pipeQA.figures as qaFig
import numpy

import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData

class ZeropointQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
	qaAna.QaAnalysis.__init__(self)


    def test(self, data, dataId):

	# get data
	self.detector = data.getDetectorBySensor(dataId)
	self.calib    = data.getCalibBySensor(dataId)

	# create a container for Raft/Ccd data
	self.data = raftCcdData.RaftCcdData(self.detector, None)

	# get the zeropoints and put them in the raftccd container
	for key, detector in self.detector.items():
	    raft = detector.getParent().getId().getName()
	    ccd = detector.getId().getName()
	    fluxMag0, fluxMag0Err = self.calib[key].getFluxMag0()
	    zp = 2.5*numpy.log10(fluxMag0)
	    self.data.set(raft, ccd, zp)

	# get stats on the zeropoints
	self.mean = self.data.summarize('mean')
	self.std  = self.data.summarize('std')

	# go through all the values and test each one
	group = dataId['visit']
	testSet = self.getTestSet(group)
	for raftName, ccdName, zp in self.data.listKeysAndValues():
	    label = "zeropoint: "+ccdName
	    value = (zp - self.mean)/self.std
	    limits = [-3.0, 3.0]
	    comment = "stdev with respect other ccds."
	    t = testCode.Test(label, value, limits, comment)
	    testSet.addTest(t)
	    self.data.set(raftName, ccdName, value)

		
    def plot(self, data, dataId, showUndefined=False):

	# create an fpaFigure and put in any data we have
	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	for raft, ccdDict in fig.data.items():
	    for ccd, value in ccdDict.items():
		fig.data[raft][ccd] = self.data.get(raft, ccd)

	# create the figure
	fig.makeFigure(showUndefined=showUndefined, vlimits=[-3.0, 3.0], cmap="gray")

	# put the figure in the appropriate testSet
	group = dataId['visit']
	testSet = self.getTestSet(group)
	
	testSet.addFigure(fig, "zeropoint.png", "Zeropoint caption")

