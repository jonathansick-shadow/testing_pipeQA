#!/usr/bin/env python
# Original filename: bin/pipeQa.py
#
# Author: 
# Email: 
# Date: Mon 2011-04-11 13:11:01
# 
# Summary: 
# 

import sys
import re
import optparse
import os
import datetime

import lsst.testing.pipeQA as pipeQA

#############################################################
#
# Main body of code
#
#############################################################

def main():

    ts = pipeQA.TestSet()

    for camInfo in pipeQA.getCameraInfoAvailable():
	print "trying camera: ", camInfo.name
	camera = camInfo.camera
	fig = pipeQA.FpaQaFigure(camera)
	fig.makeFigure()
	ts.addFigure(fig, camInfo.name, "FPA for camera: "+camInfo.name)


if __name__ == '__main__':
    main()
    
