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
from optparse import OptionParser

import lsst.testing.pipeQA as pipeQA

import lsst.pex.logging as pexLog
from lsst.pex.logging import Trace
pexLog.Trace_setVerbosity("lsst.testing.pipeQA", 1)

import eups
simdir        = eups.productDir("obs_lsstSim")
cameraGeomPaf = os.path.join(simdir, "description", "Full_STA_geom.paf")

if __name__ == '__main__':


    parser = OptionParser()
    parser.add_option('-D', '--database', dest='database',
                      default='rplante_DC3b_u_weeklytest_2011_0218_science',
                      help='Name of database for queries')
    parser.add_option('-v', '--visit', type='int', dest='visit', action='append', type='int', default=[])
    parser.add_option('-r', '--raft', dest='raft', action='append', default=[])
    parser.add_option('-s', '--sensor', dest='sensor', action='append', default=[])
    parser.add_option('-o', '--output', dest='outRoot', default='pipeQA', help='Output directory')
    parser.add_option('--photrms', dest='dophotrms', action='store_true', default=False,
                      help='Make photometric RMS plot?')
    parser.add_option('--zptfpa', dest='dozptfpa', action='store_true', default=False,
                      help='Make FPA plot of zeropoint')
    parser.add_option('--zptfit', dest='dozptfit', action='store_true', default=False,
                      help='Make photometric zeropoint fit plot?')
    
    (opt, args) = parser.parse_args()
    database    = opt.database
    outRoot     = opt.outRoot
    if not os.path.isdir(outRoot):
        Trace("lsst.testing.pipeQA.testDbQueries", 1, "Making output dir: %s" % (outRoot))
        os.makedirs(outRoot)

    # Moderates database interface
    dbId        = pipeQA.DatabaseIdentity(database)
    dbInterface = pipeQA.LsstSimDbInterface(dbId)

    if opt.dophotrms:
        sql = 'select distinct(filterName) from Science_Ccd_Exposure'
        results     = dbInterface.execute(sql)
        if len(results) == 0:
            print 'No filter data, skipping...'
        else:
            prmsfig = pipeQA.PhotometricRmsFigure()
            for filter in results:
                prmsfig.retrieveData(database, filter[0])
                prmsfig.makeFigure()
                prmsfig.saveFigure(os.path.join(outRoot, "photRms_%s.png" % (filter)))
            
    if opt.dozptfpa:
        if len(opt.visit) == 0:
            sql      = 'select distinct(visit) from Science_Ccd_Exposure'
            results  = dbInterface.execute(sql)
            visitIds = [x[0] for x in results]
        else:
            visitIds = opt.visit

        zptfpafig = pipeQA.ZeropointFpaFigure(cameraGeomPaf)
        for visitId in visitIds:
            zptfpafig.retrieveData(database, visitId)
            zptfpafig.makeFigure(doLabel = True)
            zptfpafig.saveFigure(os.path.join(outRoot, "zptFPA_%d.png" % (visitId)))

    if opt.dozptfit:
        htmlf     = pipeQA.HtmlFormatter()
        fptfitfig = pipeQA.ZeropointFitFigure()

        # Do 1 sensor only
        if len(opt.visit) == 1 and len(opt.raft) == 1 and len(opt.sensor) == 1:
            visitId = opt.visit[0]
            raft    = opt.raft[0]
            ccd     = opt.sensor[0]
            
            sql        = 'select distinct(filterName) from Science_Ccd_Exposure where visit = %d' % (visitId)
            results    = dbInterface.execute(sql) # need mag for reference catalog query
            filterName = results[0][0]

            fptfitfig.retrieveData(database, visitId, filterName, raft, ccd)
            fptfitfig.makeFigure()
            
            prefix = 'zptFit_%d' % (visitId)
            outdir = os.path.join(outRoot, prefix)
            fptfitfig.saveFigure(htmlf.generateFileName(os.path.join(outdir, prefix), raft, ccd))
            
        else:
            # Do an entire focal plane, HTML and all
            if len(opt.visit) == 0:
                sql      = 'select distinct(visit) from Science_Ccd_Exposure'
                results  = dbInterface.execute(sql)
                visitIds = [x[0] for x in results]
            else:
                visitIds = opt.visit

            for visitId in visitIds:
                sql2 = 'select distinct(filterName) from Science_Ccd_Exposure where visit = %d' % (visitId)
                results2 = dbInterface.execute(sql2) # need mag for reference catalog query
                filterName = results2[0][0]

                prefix = 'zptFit_%d' % (visitId)
                outdir = os.path.join(outRoot, prefix)
                if not os.path.isdir(outdir):
                    Trace("lsst.testing.pipeQA.testDbQueries", 1, "Making output dir: %s" % (outdir))
                    os.makedirs(outdir)
                outhtml = open(outdir+'.html', 'w')
                htmlf.generateHtml(outhtml, os.path.join(prefix, prefix))
                outhtml.close()

                sql3     = 'select raftName, ccdName from Science_Ccd_Exposure where visit = %s' % (visitId)
                results3 = dbInterface.execute(sql3)
                for raftccd in results3:
                    raft, ccd = raftccd

                    fptfitfig.retrieveData(database, visitId, filterName, raft, ccd)
                    fptfitfig.makeFigure()
                    fptfitfig.saveFigure(htmlf.generateFileName(os.path.join(outdir, prefix), raft, ccd))
                
  
        
