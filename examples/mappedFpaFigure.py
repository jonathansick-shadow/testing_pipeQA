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
    # Adding a mapped FPA figure
    # - The mapped fpa will appear on the right 'nav' panel
    # - when you click on a ccd, you'll select data to display for that ccd

    # let's add this to a different TestSet
    tsMapFpa = pipeQA.TestSet(group="tutorial", label="mapped-FPA-figure-demo")

    camInfo = pipeQA.LsstSimCameraInfo()
    fpaFig = qaFig.FpaQaFigure(camInfo)

    # Fill the fpaFig.data attribute with the values to display.
    # This will color-code ccds in the focal plane according to these values
    # - this time, assume we only have 1 raft (the remaining sensors will be greyed-out)
    raftName = "2,2"
    data = fpaFig.data
    map = fpaFig.map
    for raft in sorted(data.keys()):
        if not re.search(raftName, raft):
            continue
        ccdDict = data[raft]
        for ccd, value in ccdDict.items():
            data[raft][ccd] = numpy.random.normal(0.0, 1.0)

            # this time, add a mouse-over string (no spaces) for the map
            map[raft][ccd] = "data-from-" + ccd

            # make a plot for this ccd
            qafig = qaFig.QaFigure()
            fig = qafig.fig
            ax = fig.add_subplot(111)
            x = 2.0*numpy.pi*numpy.arange(100)/100
            ax.plot(x, numpy.sin(x)+numpy.random.normal(0.0, 0.1, len(x)), 'r-')
            ax.set_ylim([-1.0, 1.0])

            areaLabel = fpaFig.getAreaLabel(raft, ccd)  # areaLabels are handled by FpaQaFigure
            tsMapFpa.addFigure(qafig, "sine-wave.png", "Sine wave from " + ccd, areaLabel=areaLabel)

    # now tell fpaFig to make the matplotlib figure, and add the fpaFigure to the TestSet
    fpaFig.makeFigure(vlimits=[-1.0, 1.0], cmap="gray", cmapUnder="#0000ff", cmapOver="#ff0000")
    tsMapFpa.addFigure(fpaFig, camInfo.name+".png", "FPA for camera: "+camInfo.name, navMap=True)


if __name__ == '__main__':
    main()

