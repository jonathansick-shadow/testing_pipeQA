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
    ts = pipeQA.TestSet(group="tutorial", label="clickplot-howto")

    # Adding metadata
    ts.addMetadata("time-run", datetime.datetime.now().strftime("%a %Y-%m-%d %H:%M:%S"))


    # make some fake data, two clusters of stars in a 1024x1024 image
    nx, ny = 1024, 1024
    n = [100, 50]
    size = [200, 125]
    xcen = [300, 700]
    ycen = [300, 600]
    magPeak  = [24.0, 20.0]
    dmag = [2.0, 1.0]

    x = numpy.array([])
    y = numpy.array([])
    mag = numpy.array([])


    # We'll make an plot showing x,y positions of the stars
    # We'll make histograms of each cluster, and add map 'areas' to the x,y plot
    #     ... clicking in the 'area' boxes will show the histogram for the cluster we clicked on

    qafigNav = qaFig.QaFig() # this will be the x,y plot
    
    for i in range(len(n)):

	# add a map area around the cluster of points

	# create a unique id 'areaLabel' for this area
	# areaLabel must appear in:
	#  -- any Test label (it can be a prefix/suffix etc)
	#  -- and any QaFigure name which is to be displayed when this map area is clicked (prefix/suffix ok)
	areaLabel = "cluster%04d" % (i)

	# define the lower left and upper right corners of the clickable region
	x0, y0, x1, y1 = xcen[i]-size[i], ycen[i]-size[i], xcen[i]+size[i], ycen[i]+size[i]
	area      = [x0, y0, x1, y1]

	# add text to be displayed on mouse-over
	areaInfo  = areaLabel + " n=%d"%(n[i])

	qafigNav.addMapArea(areaLabel, area, areaInfo)


	# make fake x,y,mag coords and make the histogram plot
	r      = size[i]*numpy.random.uniform(0.0, 1.0, n[i])**2
	theta  = numpy.random.uniform(0.0, 2.0*numpy.pi, n[i])
	
	x      = numpy.append(x, (xcen[i] - r*numpy.cos(theta)))
	y      = numpy.append(y, (ycen[i] - r*numpy.sin(theta)))
	magTmp = numpy.random.normal(magPeak[i], dmag[i], n[i])
	mag    = numpy.append(mag, magTmp)

	qafig = qaFig.QaFig()
	fig = qafig.fig
	ax = fig.add_subplot(111)
	ax.hist(magTmp)
	filebase = "hist_" + areaLabel  # notice areaLabel is added as a part of the filename
	ts.addFigure(qafig, filebase+".png", "Histogram of stars in cluster %d" % (i))


	# add a test - notice the 'areaLabel' is part of the Test label
	ts.addTest("count_" + areaLabel, n[i], [0, None], "Verify > 0 stars in cluster.")


    # make a histogram of *all* stars
    # - this will be displayed if no 'area' region is active, or if you click 'all'
    #  --> you can skip this if you like ... then no image is displayed unless a map region is selected
    qafigAll = qaFig.QaFig()
    fig = qafigAll.fig
    ax = fig.add_subplot(111)
    ax.hist(mag)
    ts.addFigure(qafigAll, "hist_all.png", "Histogram of all stars")

    # make the x,y plot which will contain the clickable map areas
    fig = qafigNav.fig
    ax = fig.add_subplot(111)
    ax.plot(x, y, "r.")
    ax.set_xlim([0, nx])
    ax.set_ylim([0, ny])

    ts.addFigure(qafigNav, "stars%d.png" % (i), "Stars in cluster %d" % (i), saveMap=True, navMap=True)



if __name__ == '__main__':
    main()
    
