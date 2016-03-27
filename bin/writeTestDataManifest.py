#!/usr/bin/env python
#
# Original filename: bin/writeTestDataManifest.py
#
# Author: Steve Bickerton
# Email:
# Date: Fri 2011-01-07 13:19:55
#
# Summary:
#
"""
%prog [options] arg
"""

import sys
import re
import optparse
import os
import datetime
import hashlib
import lsst.testing.pipeQA as pipeQA

#############################################################
#
# Main body of code
#
#############################################################


def main(testdataDir, hashtype):
    manifest = pipeQA.Manifest(testdataDir)
    manifest.create(hashtype)
    manifest.write()


#############################################################
# end
#############################################################

if __name__ == '__main__':

    ########################################################################
    # command line arguments and options
    ########################################################################

    parser = optparse.OptionParser(usage=__doc__)
    hashtypes = ",".join(pipeQA.hashtypesDefined())
    parser.add_option("-t", "--hashtype", default="crc32",
                      help="Hashtype to use. (default=%default) [options="+hashtypes + "]")
    opts, args = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)

    testdata, = args
    testbedPath = os.getenv('TESTBED_PATH')

    # maybe we're in the testbed directory
    if os.path.exists(testdata):
        main(os.path.join(os.getcwd(), testdata), opts.hashtype)

    # maybe it's specified in the environment
    elif testbedPath is not None:
        testbedDirs = testbedPath.split(":")
        found = False
        for testbedDir in testbedDirs:
            testdataDir = os.path.join(testbedDir, testdata)
            if os.path.exists(testdataDir):
                main(testdataDir, opts.hashtype)
                found = True
        if not found:
            print "test data set "+testdata+" not found in TESTBED_PATH."
            sys.exit()

    else:
        print "Test data set " + testdata + " not found in current working directory,"
        print "and TESTBED_PATH is not defined.  Can't create manifest.  Exiting."
        sys.exit()
