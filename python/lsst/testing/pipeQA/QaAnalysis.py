import QaFigure as qaFig
import numpy

import TestCode as testCode

class QaAnalysis(object):

    def __init__(self):
	self.testSets = {}

    def test(self):
	return []

    def plot(self):
	return []

    def getTestSet(self, group):
	label = self.__class__.__name__
	if not self.testSets.has_key(group):
	    self.testSets[group] = testCode.TestSet(label, group=group)
	return self.testSets[group]


class SourceBoundsQaAnalysis(QaAnalysis):

    def __init__(self):
	QaAnalysis.__init__(self)

    def test(self, data, dataId):

	# get data
	self.ss = data.getSourceSet(dataId)

	# check the numbers for each visit
	label = "SourceBoundsQaAnalysis"
	value = 0
	limits = [-1, 1]
	comment = "dummy"

	group = dataId['visit']
	testSet = self.getTestSet(group)
	t = testCode.Test(label, value, limits, comment)
	testSet.addTest(t)

    def plot(self, data, dataId, showUndefined=False):

	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	data = fig.data

	fig.makeFigure(showUndefined=showUndefined)

	group = dataId['visit']
	testSet = self.getTestSet(group)
	testSet.addFigure(fig, "sourceBounds.png", "test caption")



class ZeropointQaAnalysis(QaAnalysis):

    def __init__(self):
	QaAnalysis.__init__(self)


    def test(self, data, dataId):

	# get data
	self.detector = data.getDetectorBySensor(dataId)
	self.calib = data.getCalibBySensor(dataId)

	group = dataId['visit']
	testSet = self.getTestSet(group)


	self.values = []
	self.data = {}
	for key, detector in self.detector.items():
	    raft = detector.getParent().getId().getName()
	    ccd = detector.getId().getName()
	    fluxMag0, fluxMag0Err = self.calib[key].getFluxMag0()
	    zp = 2.5*numpy.log10(fluxMag0)

	    if not self.data.has_key(raft):
		self.data[raft] = {}
	    self.data[raft][ccd] = zp
	    self.values.append(zp)

	self.zp = numpy.array(self.values)

	self.mean = self.zp.mean()
	self.std = self.zp.std() #numpy.sqrt(numpy.variance(self.zp))

	for raftName, ccdDict in self.data.items():
	    for ccdName, zp in ccdDict.items():
		label = "zeropoint: "+ccdName
		value = (zp - self.mean)/self.std
		limits = [-3.0, 3.0]
		comment = "stdev with respect other ccds."
		t = testCode.Test(label, value, limits, comment)
		testSet.addTest(t)

		self.data[raftName][ccdName] = value
		
    def plot(self, data, dataId, showUndefined=False):

	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	for raft, ccdDict in fig.data.items():
	    for ccd, value in ccdDict.items():
		if self.data.has_key(raft) and self.data[raft].has_key(ccd):
		    fig.data[raft][ccd] = self.data[raft][ccd]
		else:
		    fig.data[raft][ccd] = None

	fig.makeFigure(showUndefined=showUndefined, vlimits=[-3.0, 3.0], cmap="gray")

	group = dataId['visit']
	testSet = self.getTestSet(group)
	
	testSet.addFigure(fig, "zeropoint.png", "Zeropoint caption")
