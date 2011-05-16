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

import numpy

import lsst.testing.pipeQA         as pipeQA
import lsst.testing.pipeQA.figures as qaFig

#############################################################
#
# Main body of code
#
#############################################################

# This tutorial assumes you've read through dispQaTutorial1.py !

def main():


    # create a TestSet
    ts = pipeQA.TestSet(group="tutorial", label="fpa-howto")

    # Adding metadata
    ts.addMetadata("time-run", datetime.datetime.now().strftime("%a %Y-%m-%d %H:%M:%S"))


    ###############################
    # Adding an FPA figure

    # To demo multiple TestSets, let's add this figure to a different TestSet
    tsFpa = pipeQA.TestSet(group="tutorial", label="FPA-figure-demo")

    # get a cameraInfo object for the camera you're interested in
    # - available: LsstSimCameraInfo, CfhtCameraInfo, HscCameraInfo, SuprimecamCamerainfo
    camInfo = pipeQA.LsstSimCameraInfo()

    # pull out the camera, and pass it to the FPA figure constructor
    camera = camInfo.camera
    sfig = qaFig.FpaQaFigure(camera)

    # Fill the sfig.data attribute with the values to display.
    # This will color-code ccds in the focal plane according to these values
    # ... just sequentially for this example ... in range 0 --> 1
    i = 0
    data = sfig.data
    for raft in sorted(data.keys()):
	ccdDict = data[raft]
	for ccd, value in ccdDict.items():
	    data[raft][ccd] = 1.0*i
	    i += 1

    # now tell sfig to make the matplotlib figure, and add the fpaFigure to the TestSet
    sfig.makeFigure(showUndefined=True)
    tsFpa.addFigure(sfig, camInfo.name+".png", "FPA for camera: "+camInfo.name)



if __name__ == '__main__':
    main()
    
