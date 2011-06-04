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
import lsst.pex.policy as pexPolicy


#############################################################
#
# Main body of code
#
#############################################################

def main(dataset, dataIdInput, rerun=None, testRegex=".*", camera=None):

    data = pipeQA.makeQaData(dataset, rerun=rerun, retrievalType=camera)

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

    # Policy for which plots to make and value of quality metrics
    policyDictName = "PipeQa.paf"
    policyFile = pexPolicy.DefaultPolicyFile("testing_pipeQA", policyDictName, "policy")
    policy  = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    
    analysisList = []
#    if policy.get("doZptQa"):
#        zptMin = policy.get("zptQaMetricMin")
#        zptMax = policy.get("zptQaMetricMax")
#        analysisList.append(qaAnalysis.ZeropointQaAnalysis(zptMin, zptMax))
    if policy.get("doZptFitQa"):
        offsetMin = policy.get("zptFitQaOffsetMin")
        offsetMax = policy.get("zptFitQaOffsetMax")
        analysisList.append(qaAnalysis.ZeropointFitQa(offsetMin, offsetMax))
#    if policy.get("doEmptySectorQa"):
#        maxMissing = policy.get("emptySectorMaxMissing")
#        analysisList.append(qaAnalysis.EmptySectorQaAnalysis(maxMissing, nx = 4, ny = 4))
#    if policy.get("doAstromQa"):
#        analysisList.append(qaAnalysis.AstrometricErrorQaAnalysis(policy.get("astromQaMaxErr")))
#    if policy.get("doPhotCompareQa"):
#        magCut   = policy.get("photCompareMagCut")
#        deltaMin = policy.get("photCompareDeltaMin")
#        deltaMax = policy.get("photCompareDeltaMax")
#        rmsMax   = policy.get("photCompareRmsMax")
#        slopeMin = policy.get("photCompareSlopeMin")
#        slopeMax = policy.get("photCompareSlopeMax")
#        for types in policy.getStringArray("photCompareTypes"):
#            cmp1, cmp2 = types.split()
#            analysisList.append(qaAnalysis.PhotCompareQaAnalysis(cmp1, cmp2, magCut, deltaMin, deltaMax,
#                                                                 rmsMax, slopeMin, slopeMax))
#    if policy.get("doPsfEllipQa"):
#        analysisList.append(qaAnalysis.PsfEllipticityQaAnalysis(policy.get("psfEllipMax")))
#    if policy.get("doCompleteQa"):
#        analysisList.append(qaAnalysis.CompletenessQa(policy.get("completeMinMag"),
#                                                      policy.get("completeMaxMag")))
#    if policy.get("doCompleteQa2"):
#        analysisList.append(qaAnalysis.CompletenessQa2(policy.get("completeMinMag"),
#                                                       policy.get("completeMaxMag")))
       
    for visit in visits:
        for a in analysisList:
            
            test = str(a)
            if not re.search(testRegex, test):
                continue
            
            print "Running " + test + "  visit:" + str(visit)
            dataIdVisit = copy.copy(dataId)
            dataIdVisit['visit'] = visit
            try:
                a.test(data, dataIdVisit)
            except Exception, e:
                print 'WARNING:', e
                pass
            else:
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
    parser.add_option("-C", "--camera", default=None,
		      help="Specify a camera and override auto-detection (default=%default)")
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
    
    main(dataset, dataId, rerun, opts.test, opts.camera)
