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

def main(dataset, dataIdInput, rerun=None, testRegex=".*"):

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
        #qaAnalysis.ZeropointQaAnalysis(),
        qaAnalysis.EmptySectorQaAnalysis(4, 4),
        qaAnalysis.AstrometricErrorQaAnalysis(),
        qaAnalysis.PhotCompareQaAnalysis("psf", "cat", cut=magCut),
        qaAnalysis.PhotCompareQaAnalysis("psf", "ap",  cut=magCut),
        qaAnalysis.PhotCompareQaAnalysis("psf", "mod", cut=magCut),
        qaAnalysis.PsfEllipticityQaAnalysis(),
        ]

    for visit in visits:
        for a in analysisList:
            
            test = str(a)
            if not re.search(testRegex, test):
                continue
            
            print "Running " + test + "  visit:" + str(visit)
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
    parser.add_option("-R", "--rerun", default=None,
                      help="Rerun to analyse - only valid for hsc/suprimecam (default=%default)")
    parser.add_option("-s", "--snap", default=".*",
                      help="Specify snap as regex (default=%default)")
    parser.add_option("-t", "--test", default=".*",
                      help="Regex specifying which QaAnalysis to run (default=%default)")
    
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
    
    main(dataset, dataId, rerun, opts.test)
