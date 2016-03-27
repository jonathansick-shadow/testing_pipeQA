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

import lsst.testing.pipeQA as pipeQA
import lsst.testing.pipeQA.figures as qaFig

#############################################################
#
# Main body of code
#
#############################################################

# This tutorial assumes you've read through dispQaTutorial1.py and dispQaTutorial2.py !


def main():

    ##############################
    # Adding a plain FPA figure

    # To demo multiple TestSets, let's add this figure to a different TestSet
    tsFpa = pipeQA.TestSet(group="tutorial", label="plain-FPA-figure-demo")

    # get a cameraInfo object for the camera you're interested in
    # - available: LsstSimCameraInfo, CfhtCameraInfo, HscCameraInfo, SuprimecamCamerainfo
    camInfo = pipeQA.LsstSimCameraInfo()

    # pass it to the FPA figure constructor
    fpaFig = qaFig.FpaQaFigure(camInfo)

    # Fill the fpaFig.data attribute with the values to display.
    # This will color-code ccds in the focal plane according to these values
    data = fpaFig.data
    for raft in sorted(data.keys()):
        ccdDict = data[raft]
        for ccd, value in ccdDict.items():
            data[raft][ccd] = numpy.random.normal(0.0, 1.0)

    # now tell fpaFig to make the matplotlib figure, and add the fpaFigure to the TestSet
    vlimits = [-1.0, 1.0]  # the colormap limits
    cmap = "copper"        # the name of a matplotlib.cm colormap ('jet' by default)

    # color to use if below/above the cmap range, defaults are min/max of cmap
    cmapUnder, cmapOver = "#0000ff", "#ff0000"

    fpaFig.makeFigure(vlimits=vlimits, cmap="copper", cmapUnder=cmapUnder, cmapOver=cmapOver)
    tsFpa.addFigure(fpaFig, camInfo.name+".png", "FPA for camera: "+camInfo.name)


if __name__ == '__main__':
    main()

