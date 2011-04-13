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

import lsst.testing.pipeQA as pipeQA

#############################################################
#
# Main body of code
#
#############################################################

def main(dataset, dataIdInput):

    data = pipeQA.makeQaData(dataset)

    # take what we need for this camera, ignore the rest
    dataId = {}
    for name in data.dataIdNames:
	if dataIdInput.has_key(name):
	    dataId[name] = dataIdInput[name]

    # if they requested a key that doesn't exist for this camera ... throw
    for k, v in dataIdInput.items():
	if (not k in data.dataIdNames) and (v != '.*'):
	    raise Exception("Key "+k+" not available for this dataset (camera="+data.cameraInfo.name+")")


    analysisList = [
	pipeQA.SourceBoundsQaAnalysis(),
	pipeQA.ZeropointQaAnalysis(),
	]

    for a in analysisList:
	a.test(data, dataId)
	a.plot(data, dataId)
	
	


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
    dataset, = args
    
    #dataset = 'buildbot_DC3b_u_weekly_production_trunk_2011_0402_143716_science'
    #dataId = {'visit':'85501858', 'snap':'0', 'raft':'2,2', 'sensor':'0,0'}

    dataset = 'hscsimTestData002'
    #dataset = 'cfhtTestData002'
    #dataset = 'suprimeTestData002'
    dataId['visit'] = '.*'

    main(dataset, dataId)
