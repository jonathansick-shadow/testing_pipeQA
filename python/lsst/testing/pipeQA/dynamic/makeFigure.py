import sys, os, re

def makeFigure(data, dataId, dynaFig):
    moduleName = "lsst.testing.pipeQA.dynamic."+dynaFig
    __import__(moduleName)
    module = sys.modules[moduleName]
    module.plot(data, dataId)
