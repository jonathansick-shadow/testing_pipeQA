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

    ts = pipeQA.TestSet(group="debug")

    for camInfo in pipeQA.getCameraInfoAvailable():
	print "trying camera: ", camInfo.name
	camera = camInfo.camera
	fig = pipeQA.FpaQaFigure(camera)

	# set values to range 0 to 1
	i = 0
	for raft,ccdDict in fig.data.items():
	    for ccd, value in ccdDict.items():
		fig.data[raft][ccd] = 1.0*i
		i += 1
	
	fig.makeFigure(size=(1024, 1024), showUndefined=True)
	ts.addFigure(fig, camInfo.name+".png", "FPA for camera: "+camInfo.name)



if __name__ == '__main__':
    main()
    
