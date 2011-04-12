
class QaAnalysis(object):

    def __init__(self):
	pass

    def test(self):
	return []

    def plot(self):
	return []


import Test as qaTest
import TestSet as qaTestSet
import QaFigure as qaFig
import numpy

class SourceBoundsQaAnalysis(QaAnalysis):

    def __init__(self):
	QaAnalysis.__init__(self, testSet=None)

	if not testSet is None:
	    self.testSet = testSet
	else:
	    self.testSet = qaTestSet.TestSet(label='sourceBounds')


    def test(self, data, dataId):

	# get data
	self.ss = data.getSourceSet(dataId)

	# check the numbers for each visit
	label = "SourceBoundsQaAnalysis"
	value = 0
	limits = [-1, 1]
	comment = "dummy"
	t = qaTest.Test(label, value, limits, comment)

	self.testSet.addTest(t)
	
	return [t]


    def plot(self, data, dataId):

	fig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	data = fig.data

	fig.makeFigure()

	for k, v in data.items():
	    print k, v
	    
	return [fig]
