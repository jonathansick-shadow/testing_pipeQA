#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Test to verify quality of PSF photometry on test frames
"""
import os
import unittest
import lsst.testing.pipeQA as pipeQA
import numpy

import lsst.afw.detection as afwDet

import matplotlib
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas

#import lsst.meas.extensions.shapeHSM.hsmLib as shapeHSM


def testPipetteAllMappers():
    """Verify pipette runs to completion with all supported mappers."""


    doLsstSim    = 1#True
    doCfht       = 1#True
    doSuprimecam = 1#True
    doHscSim     = 1#True

    
    pr = pipeQA.PipeRunner()

    
    #########################
    # create the testdata objects
    
    ## LsstSim
    if doLsstSim:
	tdLsstSim = pipeQA.makeTestData("imsimTestData001",
					visit='85501867', snap='0', raft='1,1', sensor='1,1',
					verifyChecksum=False, outDir='local',
					astrometryNetData="imsim_20100625")
	pr.addTestData(tdLsstSim)
    

    ## megacam
    if doCfht:
	tdCfht = pipeQA.makeTestData("cfhtTestData001", visit="788033", ccd='17',
				      verifyChecksum=False, outDir='local',
				      astrometryNetData="hsc-dc2-2011-02-27plus")
	pr.addTestData(tdCfht)


    ## Suprimecam
    visit = '101412' # ZR
    #visit = '108504' # I+
    if doSuprimecam:
	tdScSim = pipeQA.makeTestData("suprimeTestData001", visit=visit, ccd='2', 
				      verifyChecksum=False, outDir='local',
				      astrometryNetData="hsc-dc1-2010-08-04.1-starsonly")
	pr.addTestData(tdScSim)

    ## HscSim
    if doHscSim:
	tdHscSim = pipeQA.makeTestData("hscsimTestData001", visit='200', ccd='50', 
				       verifyChecksum=False, outDir='local',
				       astrometryNetData="hsc-dc2-2011-02-27plus")
	pr.addTestData(tdHscSim)
	


				     

    ##########################
    # run the pipe
    hsmConfig = os.path.join(os.getenv('MEAS_EXTENSIONS_SHAPEHSM_DIR'), "policy", "hsmShape.paf")
    pr.run(force=False, overrideConfig=[])



    ##########################
    # Test the outputs exist ... mostly sane
    # - note that this is not intended to verify quality per se, just to
    #   check that pipette still runs.
    
    ts = pipeQA.TestSet(mainDisplayFull=True)
    
    ts.importLogs(pr.getLogFiles())
    ts.importEupsSetups(pr.getEupsSetupFiles())
    ts.importExceptionDict(pr.getUncaughtExceptionDict())


    ##########################
    # get the data we want from the piperunner and perform a test
    sourceSets = {
	"lsstSim"    : pr.getSourceSet(visit='85501867', raft='1,1'),
	"hscSim"     : pr.getSourceSet(visit='200', ccd='50'),
	"suprimecam" : pr.getSourceSet(visit=visit, ccd='2'),
	"cfht"       : pr.getSourceSet(visit='788033', ccd='17'),
	}

    minDetections = 100
    for label, ss in sourceSets.items():
	n       = len(ss)
	limits  = [minDetections, None]
	comment = "verify no. detections > %d"%(minDetections)
	ts.addTest(label, n, limits, comment)



    
if __name__ == "__main__":
    testPipetteAllMappers()

