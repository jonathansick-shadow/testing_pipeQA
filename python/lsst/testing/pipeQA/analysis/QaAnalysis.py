import lsst.testing.pipeQA.figures as qaFig
import numpy

import lsst.testing.pipeQA.TestCode as testCode

class QaAnalysis(object):

    def __init__(self, testLabel=None):
	self.testSets = {}
	self.testLabel = testLabel

    def test(self):
	return []

    def plot(self):
	return []

    def getTestSet(self, group, label=None):
	if not label is None:
	    label = self.__class__.__name__ + "."+ label
	else:
	    label = self.__class__.__name__
	    
	if not self.testSets.has_key(group):
	    self.testSets[group] = testCode.TestSet(label, group=group)
	return self.testSets[group]

    def __str__(self):
	testLabel = ""
	if not self.testLabel is None:
	    testLabel = "."+self.testLabel
	return self.__class__.__name__ + testLabel
    

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



