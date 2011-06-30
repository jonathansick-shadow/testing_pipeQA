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
import traceback
import copy
import time

import lsst.testing.pipeQA as pipeQA
import lsst.testing.pipeQA.analysis     as qaAnalysis
import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Trace

import numpy

def getMemUsageThisPid(size="rss"):
    """Generalization; memory sizes: rss, rsz, vsz."""
    #return 1.0
    return int(os.popen('ps -p %d -o %s | tail -1' % (os.getpid(), size)).read())

#############################################################
#
# Main body of code
#
#############################################################

def main(dataset, dataIdInput, rerun=None, testRegex=".*", camera=None,
         exceptExit=False, keep=False, wwwCache=True, breakBy='visit'):

    if exceptExit:
        numpy.seterr(all="raise")
    
    data = pipeQA.makeQaData(dataset, rerun=rerun, retrievalType=camera)

    ccdConvention = 'ccd'
    if data.cameraInfo.name == 'lsstSim':
        ccdConvention = 'sensor'
        if dataIdInput.has_key('ccd'):
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


    # Policy for which plots to make and value of quality metrics
    policyDictName = "PipeQa.paf"
    policyFile = pexPolicy.DefaultPolicyFile("testing_pipeQA", policyDictName, "policy")
    policy  = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)
    
    analysisList = []
    if data.cameraInfo.name in policy.getStringArray("doZptQa"):
        zptMin = policy.get("zptQaMetricMin")
        zptMax = policy.get("zptQaMetricMax")
        analysisList.append(qaAnalysis.ZeropointQaAnalysis(zptMin, zptMax, useCache=keep, wwwCache=wwwCache))
    if data.cameraInfo.name in policy.getStringArray("doZptFitQa"):
        offsetMin = policy.get("zptFitQaOffsetMin")
        offsetMax = policy.get("zptFitQaOffsetMax")
        analysisList.append(qaAnalysis.ZeropointFitQa(offsetMin, offsetMax, useCache=keep, wwwCache=wwwCache))
    if data.cameraInfo.name in policy.getStringArray("doEmptySectorQa"):
        maxMissing = policy.get("emptySectorMaxMissing")
        analysisList.append(qaAnalysis.EmptySectorQaAnalysis(maxMissing, nx = 4, ny = 4, useCache=keep,
                                                             wwwCache=wwwCache))
    if data.cameraInfo.name in policy.getStringArray("doAstromQa"):
        analysisList.append(qaAnalysis.AstrometricErrorQaAnalysis(policy.get("astromQaMaxErr"),
                                                                  useCache=keep, wwwCache=wwwCache))
    if data.cameraInfo.name in policy.getStringArray("doPhotCompareQa"):
        magCut   = policy.get("photCompareMagCut")
        deltaMin = policy.get("photCompareDeltaMin")
        deltaMax = policy.get("photCompareDeltaMax")
        rmsMax   = policy.get("photCompareRmsMax")
        slopeMin = policy.get("photCompareSlopeMinSigma")
        slopeMax = policy.get("photCompareSlopeMaxSigma")
        for types in policy.getStringArray("photCompareTypes"):
            cmp1, cmp2 = types.split()
            analysisList.append(qaAnalysis.PhotCompareQaAnalysis(cmp1, cmp2, magCut, deltaMin, deltaMax,
                                                                 rmsMax, slopeMin, slopeMax, useCache=keep,
                                                                 wwwCache=wwwCache))
    if data.cameraInfo.name in policy.getStringArray("doPsfShapeQa"):
        analysisList.append(qaAnalysis.PsfShapeQaAnalysis(policy.get("psfEllipMax"),
                                                          policy.get("psfFwhmMax"), useCache=keep,
                                                          wwwCache=wwwCache))
    if data.cameraInfo.name in policy.getStringArray("doCompleteQa"):
        analysisList.append(qaAnalysis.CompletenessQa(policy.get("completeMinMag"),
                                                      policy.get("completeMaxMag"), useCache=keep,
                                                      wwwCache=wwwCache))
    if data.cameraInfo.name in policy.getStringArray("doVignettingQa"):
        analysisList.append(qaAnalysis.VignettingQa(policy.get("vigMaxMedian"),
                                                    policy.get("vigMagRms"),
                                                    policy.get("vigMaxMag")))


    useFp = open("runtimePerformance.dat", 'w')
    useFp.write("#%-11s %-24s %-32s %10s  %16s\n" %
                ("timestamp", "dataId", "testname", "t-elapsed", "resident-memory-kb-Mb"))
    
    # split by visit
    visits = data.getVisits(dataId)
    for visit in visits:

        visit_t0 = time.time()
        
        testset = pipeQA.TestSet(group="", label="QA-failures", wwwCache=wwwCache)

        dataIdVisit = copy.copy(dataId)
        dataIdVisit['visit'] = visit

        # now break up the run into eg. rafts or ccds
        #  ... if we only run one raft or ccd at a time, we use less memory
        brokenDownDataIdList = data.breakDataId(dataIdVisit, breakBy)
        
        for thisDataId in brokenDownDataIdList:

            for a in analysisList:

                test_t0 = time.time()
                test = str(a)
                if not re.search(testRegex, test):
                    continue

                print "Running " + test + "  visit:" + str(visit)


                # For debugging, it's useful to exit on failure, and get
                # the full traceback
                memory = 0
                if exceptExit:
                    a.test(data, thisDataId)
                    a.plot(data, thisDataId)
                    memory = getMemUsageThisPid()
                    a.free()

                # otherwise, we want to continue gracefully
                else:
                    try:
                        a.test(data, thisDataId)
                    except Exception, e:
                        exc_type, exc_value, exc_traceback = sys.exc_info()
                        s = traceback.format_exception(exc_type, exc_value,
                                                       exc_traceback)
                        label = "visit_%s_analysis_%s" % (visit, test)
                        print "Warning: Exception in QA processing of visit:%s, analysis:%s" % (visit, test)
                        #print "       :", "".join(s)
                        testset.addTest(label, 1, [0, 0], "QA exception thrown", backtrace="".join(s))
                    else:
                        a.plot(data, thisDataId)
                        memory = getMemUsageThisPid()
                        a.free()


                test_tf = time.time()
                tstamp = time.mktime(datetime.datetime.now().timetuple())
                idstamp = "v"+str(visit)+"r"+thisDataId['raft']+"s"+thisDataId[ccdConvention]
                useFp.write("%-12.1f %-24s %-32s %9.2fs %7d %7.2f\n" %
                            (tstamp, str(idstamp), test, test_tf-test_t0, memory, memory/1024.0))
                useFp.flush()

            # we're now done this dataId ... can clear the cache
            data.clearCache()
            
    useFp.close()


#############################################################
# end
#############################################################




if __name__ == '__main__':
    
    ########################################################################
    # command line arguments and options
    ########################################################################
    
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-b", "--breakBy", default="visit",
                      help="Break the run by 'visit','raft', or 'ccd' (default=%default)")
    parser.add_option("-v", "--visit", default="-1",
                      help="Specify visit as regex. Use neg. number for last 'n' visits. (default=%default)")
    parser.add_option("-c", "--ccd", default=".*",
                      help="Specify ccd as regex (default=%default)")
    parser.add_option("-C", "--camera", default=None,
		      help="Specify a camera and override auto-detection (default=%default)")
    parser.add_option("-e", "--exceptExit", default=False, action='store_true',
                      help="Don't capture exceptions, fail and exit (default=%default)")
    parser.add_option("-k", "--keep", default=False, action="store_true",
                      help="Keep existing outputs (default=%default)")
    parser.add_option("-r", "--raft", default=".*",
                      help="Specify raft as regex (default=%default)")
    parser.add_option("-R", "--rerun", default=None,
                      help="Rerun to analyse - only valid for hsc/suprimecam (default=%default)")
    parser.add_option("-s", "--snap", default=".*",
                      help="Specify snap as regex (default=%default)")
    parser.add_option("-t", "--test", default=".*",
                      help="Regex specifying which QaAnalysis to run (default=%default)")
    parser.add_option("-V", "--verbosity", default=1,
                      help="Trace level for lsst.testing.pipeQA")

    parser.add_option("--noWwwCache", default=False, action="store_true",
                      help="Disable caching of pass/fail (needed to run in parallel) (default=%default)")
    
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

    Trace.setVerbosity('lsst.testing.pipeQA', int(opts.verbosity))

    if not re.search("^(visit|raft|ccd)$", opts.breakBy):
        print "breakBy (-b) must be 'visit', 'raft', or 'ccd'"
        sys.exit()


    wwwCache = not opts.noWwwCache

    main(dataset, dataId, rerun=rerun,
         testRegex=opts.test,          camera=opts.camera,
         exceptExit=opts.exceptExit,   keep=opts.keep,      wwwCache=wwwCache,
         breakBy=opts.breakBy)
        
