import QaFigure as qaFig
import numpy

import TestCode as testCode

class QaAnalysis(object):

    def __init__(self, testSet=None, label=""):
	if not testSet is None:
	    self.testSet = testSet
	else:
	    self.testSet = testCode.TestSet(label=label)

    def test(self):
	return []

    def plot(self):
	return []



class SourceBoundsQaAnalysis(QaAnalysis):

    def __init__(self, testSet=None):
	QaAnalysis.__init__(self, testSet, label="SourceBounds")

    def test(self, data, dataId):

	# get data
	self.ss = data.getSourceSet(dataId)

	# check the numbers for each visit
	label = "SourceBoundsQaAnalysis"
	value = 0
	limits = [-1, 1]
	comment = "dummy"
	t = testCode.Test(label, value, limits, comment)

	self.testSet.addTest(t)
	
    def plot(self, data, dataId):

	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	data = fig.data

	fig.makeFigure()

	self.testSet.addFigure(fig, "sourceBounds.png", "test caption")


class ZeropointQaAnalysis(QaAnalysis):

    def __init__(self, testSet=None):
	QaAnalysis.__init__(self, testSet, label="Zeropoint")


    def test(self, data, dataId):

	# get data
	self.detector = data.getDetectorBySensor(dataId)
	self.calib = data.getCalibBySensor(dataId)

	
	# check the numbers for each visit
	label = "Zeropoint"
	value = 0
	limits = [-1, 1]
	comment = "dummy"
	t = testCode.Test(label, value, limits, comment)

	self.testSet.addTest(t)
	
    def plot(self, data, dataId):

	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)

	for key, detector in self.detector.items():
	    raft = detector.getParent().getId().getName()
	    ccd = detector.getId().getName()
	    fluxMag0, fluxMag0Err = self.calib[key].getFluxMag0()
	    fig.data[raft][ccd] = 2.5*numpy.log10(fluxMag0)
	    print raft, ccd

	fig.makeFigure()

	self.testSet.addFigure(fig, "zeropoint.png", "Zeropoint caption")
