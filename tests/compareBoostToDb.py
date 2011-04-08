import sys

from optparse import OptionParser
import lsst.testing.pipeQA as pipeQA

from lsst.pex.logging import Trace
import lsst.pex.logging as pexLog
pexLog.Trace_setVerbosity("lsst.testing.pipeQA", 2)

filterMap = {0:"u", 1:"g", 2:"r", 3:"i", 4:"z"}

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-D', '--database', dest='database',
                      default=None,
                      help='Name of database for queries (e.g. rplante_DC3b_u_weeklytest_2011_0218_science')
    parser.add_option('-R', '--registry', dest='registry',
                      default=None,
                      help='Name of registry for inputs (e.g. /lsst2/imsim-TuesdayRuns/imSim/registry.sqlite3)')
    parser.add_option('-i', '--input', dest='inRoot',
                      default = None,
                      help='Data input root for butler (e.g. /lsst3/weekly/datarel-runs/test_trunk_prod3_2011_0324)')

    parser.add_option('-v', '--visit', type='int', dest='visit', action='append', type='int', default=[])
    parser.add_option('-r', '--raft', dest='raft', action='append', default=[])
    parser.add_option('-s', '--sensor', dest='sensor', action='append', default=[])

    (opt, args) = parser.parse_args()
    registry = opt.registry
    inRoot   = opt.inRoot
    if (registry == None) or (inRoot == None):
        Trace("lsst.testing.pipeQA.compareBoostToDb", 1, "Error: registry (-R) and input root dir (-i) required")
        parser.print_help()
        sys.exit(1)

    database = opt.database
    if (database == None):
        Trace("lsst.testing.pipeQA.compareBoostToDb", 1, "Error: database (-D) required")
        parser.print_help()
        sys.exit(1)
    
    # database
    dbId        = pipeQA.DatabaseIdentity(database)
    dbInterface = pipeQA.LsstSimDbInterface(dbId)

    # boost 
    butler = pipeQA.PipeQaUtils.getInputButler(inRoot, registry)

    # all allowed/requested combos
    allkeys = pipeQA.PipeQaUtils.getAllKeysOpt(opt, butler)
    for key in allkeys:
        print 'key', key
        if not butler.datasetExists('src', key):
            Trace("lsst.testing.pipeQA.compareBoostToDb", 1, "Error: butler does not recognize %s" % (key))

        psv = butler.get('src', key)
        for source in psv.getSources():
            sourceId = source.getId()

            sql = 'select psfFlux,psfFluxSigma,apFlux,apFluxSigma,filterId from Source where sourceId=%d' % (sourceId)
            results = dbInterface.execute(sql)

            # Sources that cannot be associated into Objects by SourceAssoc get put into the BadSource table.
            # We can probably safely ignore them here
            if len(results) == 0:
                continue

            # Values from the DB
            psfFluxDb, psfFluxSigmaDb, apFluxDb, apFluxSigmaDb, filterId = results[0]
            filterName = filterMap[filterId]

            # Values from Boost
            psfFluxSrc      = source.getPsfFlux()
            psfFluxSigmaSrc = source.getPsfFluxErr()
            apFluxSrc       = source.getApFlux()
            apFluxSigmaSrc  = source.getApFluxErr()

            # OK, these are all 1.0, basically
            print '%d %.3f %.3f %.3f %.3f' % (sourceId,
                                              psfFluxDb/psfFluxSrc, psfFluxSigmaDb/psfFluxSigmaSrc,
                                              apFluxDb/apFluxSrc, apFluxSigmaDb/apFluxSigmaSrc),

            # Get ref mag
            sql  = 'select s.objectId, '                                                         # 0: objectId
            sql += ' sro.%sMag,' % (filterName)                                                  # 1: catMag
            sql += ' dnToAbMag(s.apFlux, sce.fluxMag0),'                                         # 2: apMag
            sql += ' dnToAbMagSigma(s.apFlux, s.apFluxSigma, sce.fluxMag0, sce.fluxMag0Sigma),'  # 3: apMagSigma
            sql += ' dnToAbMag(s.psfFlux, sce.fluxMag0),'                                        # 4: psfMag
            sql += ' dnToAbMagSigma(s.psfFlux, s.psfFluxSigma, sce.fluxMag0, sce.fluxMag0Sigma)' # 5: psfMagSigma
            sql += ' from Source as s, Science_Ccd_Exposure as sce,'
            sql += ' RefObjMatch as rom, SimRefObject as sro'
            sql += ' where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
            sql += ' and (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
            sql += ' and (sro.isStar = 1) and (sro.isVar = 0)'
            sql += ' and (s.sourceId = %d)' % (sourceId)
            results = dbInterface.execute(sql)
            #print sql

            # Things not matched to ref catalog as non-variable stars we can also ignore here
            if len(results) == 0:
                print
                continue

            # Problems?
            objectId, catMag, apMag, apMagSigma, psfMag, psfMagSigma = results[0]
            if apMag == None or psfMag == None:
                print
                continue
            
            print '%.3f %.3f %.3f' % (catMag, (apMag - catMag), (psfMag - catMag))
