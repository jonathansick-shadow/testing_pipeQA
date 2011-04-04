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
import os
import sys
from optparse import OptionParser

import lsst.testing.pipeQA as pipeQA

import lsst.pex.logging as pexLog
from lsst.pex.logging import Trace
pexLog.Trace_setVerbosity("lsst.testing.pipeQA", 2)

import eups
simdir        = eups.productDir("obs_lsstSim")
cameraGeomPaf = os.path.join(simdir, "description", "Full_STA_geom.paf")

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-R', '--registry', dest='registry',
                      default=None,
                      help='Name of registry for inputs (e.g. /lsst2/imsim-TuesdayRuns/imSim/registry.sqlite3)')
    parser.add_option('-i', '--input', dest='inRoot',
                      default = None,
                      help='Data input root for butler (e.g. /lsst3/weekly/datarel-runs/test_trunk_prod3_2011_0324)')
    parser.add_option('-v', '--visit', type='int', dest='visit', action='append', type='int', default=[])
    parser.add_option('-r', '--raft', dest='raft', action='append', default=[])
    parser.add_option('-s', '--sensor', dest='sensor', action='append', default=[])
    parser.add_option('-o', '--output', dest='outRoot', default='pipeQA', help='Output directory prefix')

    parser.add_option('--centshift', dest='docentshift', action='store_true', default=False,
                      help='Make plot of centroid shifts between snaps?')

    (opt, args) = parser.parse_args()
    registry = opt.registry
    inRoot   = opt.inRoot
    if (registry == None) or (inRoot == None):
        Trace("lsst.testing.pipeQA.testDbQueries", 1, "Error: registry (-R) and input root dir (-i) required")
        parser.print_help()
        sys.exit(1)
    butler = pipeQA.PipeQaUtils.getInputButler(inRoot, registry)
        
    outRoot     = os.path.join(opt.outRoot)
    if not os.path.isdir(outRoot):
        Trace("lsst.testing.pipeQA.testDbQueries", 1, "Making output dir: %s" % (outRoot))
        os.makedirs(outRoot)

    if opt.docentshift:
        if len(opt.visit):
            visits = [opt.visit,]
        else:
            visits = butler.queryMetadata('raw', 'visit')

        cfig = pipeQA.CentroidFpaFigure(cameraGeomPaf)
        for visitId in visits:
            cfig.retrieveDataViaButler(butler, visitId)
            cfig.makeFigure(doLabel = True)
            cfig.saveFigure(os.path.join(outRoot, "cshift_%d.png" % (visitId)))
            
            

        
                                              
