#!/usr/bin/env python
#
# Original filename: bin/pipeQa.py
#
# Author: 
# Email: 
# Date: Mon 2011-04-11 13:11:01
# 
# Summary: 
# 
"""
%prog [options] dataset
  dataset = Directory in TESTBED_PATH, full path to data, or database name
"""

import sys
import re
import optparse
import os
import datetime
import copy

import lsst.testing.pipeQA as pipeQA
import lsst.testing.pipeQA.analysis     as qaAnalysis


#############################################################
#
# Main body of code
#
#############################################################

def main(dataset, dataIdInput, rerun=None):

    data = pipeQA.makeQaData(dataset, rerun=rerun)

    if data.cameraInfo.name == 'lsstSim' and dataIdInput.has_key('ccd'):
	dataIdInput['sensor'] = dataIdInput['ccd']
	del dataIdInput['ccd']

    # take what we need for this camera, ignore the rest
    dataId = {}
    for name in data.dataIdNames:
	if dataIdInput.has_key(name):
	    dataId[name] = dataIdInput[name]

    # if they requested a key that doesn't exist for this camera ... throw
    for k, v in dataIdInput.items():
	if (not k in data.dataIdNames) and (v != '.*'):
	    raise Exception("Key "+k+" not available for this dataset (camera="+data.cameraInfo.name+")")


    # split by visit for now
    visits = data.getVisits(dataId)

    magCut = 20.0
    analysisList = [
	#qaAnalysis.SourceBoundsQaAnalysis(),
	#qaAnalysis.ZeropointQaAnalysis(),
	qaAnalysis.AstrometricErrorQaAnalysis(),
	qaAnalysis.PhotCompareQaAnalysis("psf", "cat", cut=magCut),
	qaAnalysis.PhotCompareQaAnalysis("psf", "ap",  cut=magCut),
	qaAnalysis.PhotCompareQaAnalysis("psf", "mod", cut=magCut),
	qaAnalysis.PsfEllipticityQaAnalysis(),
	]

    for visit in visits:
	for a in analysisList:
	    print "Running " + str(a), "  visit:", visit
	    dataIdVisit = copy.copy(dataId)
	    dataIdVisit['visit'] = visit
	    a.test(data, dataIdVisit)
	    a.plot(data, dataIdVisit, showUndefined=False)
	    a.free()
	data.clearCache()


#############################################################
# end
#############################################################




if __name__ == '__main__':
    
    ########################################################################
    # command line arguments and options
    ########################################################################
    
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-v", "--visit", default="-1",
		      help="Specify visit as regex. Use neg. number for last 'n' visits. (default=%default)")
    parser.add_option("-c", "--ccd", default=".*",
		      help="Specify ccd as regex (default=%default)")
    parser.add_option("-r", "--raft", default=".*",
		      help="Specify raft as regex (default=%default)")
    parser.add_option("-s", "--snap", default=".*",
		      help="Specify snap as regex (default=%default)")
    parser.add_option("-C", "--camera", default="L",
		      help="(L)sst, (C)fht, (H)sc, (S)uprime")
    parser.add_option("-R", "--rerun", default=None,
		      help="Rerun to analyse - only valid for hsc/suprimecam (default=%default)")
    
    opts, args = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)

    dataId = {
	'visit': opts.visit,
	'ccd': opts.ccd,
	'raft': opts.raft,
	'snap': opts.snap,
	}
    rerun = opts.rerun
    dataset, = args


    if False:
	if opts.camera=='L':
	    dataset = 'buildbot_DC3b_u_weekly_production_trunk_2011_0402_143716_science'
	    dataset = 'buildbot_weekly_latest'
	    #dataset = "update"
	    dataId = {'visit':'85501858', 'snap':'0', 'raft':'2,2', 'sensor':'.*'}
	    dataId = {'visit':'8.*', 'snap':'0', 'raft':'.*', 'sensor':'.*'}
	    #dataId = {'visit':'855.*', 'snap':'0', 'raft':'.*', 'sensor':'.*'}
	    dataId = {'visit':'857064441', 'snap':'0', 'raft':'2,2', 'sensor':'.*'}
	    #dataId = {'visit':'855018581', 'snap':'0', 'raft':'.*', 'sensor':'1,1'}

	elif opts.camera=='H':
	    dataset = 'hscsimTestData002'
	    dataId['visit'] = '.*'
	    #dataset = "Subaru"
	    #dataId = {'visit':'.*', 'ccd':'.*'}
	    #dataId = {'visit':'216', 'ccd':'.*'}
	    #rerun = "price"
	elif opts.camera=='C':
	    dataset = 'cfhtTestData002'
	    dataId['visit'] = '.*'
	elif opts.camera=="S":
	    dataset = 'suprimeTestData002'
	    dataId['visit'] = '.*'

    
    main(dataset, dataId, rerun)
