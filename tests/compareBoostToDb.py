from optparse import OptionParser
import lsst.testing.pipeQA as pipeQA

import lsst.daf.persistence as dafPersist
import lsst.daf.base as dafBase
import lsst.afw.detection as afwDet

def sourceset_read_boost(boost):
    print '# Reading', boost, 
    loc = dafPersist.LogicalLocation(boost)
    storageList = dafPersist.StorageList()
    additionalData = dafBase.PropertySet()
    persistence = dafPersist.Persistence.getPersistence(pexPolicy.Policy())
    storageList.append(persistence.getRetrieveStorage("BoostStorage", loc))
    psvptr = persistence.unsafeRetrieve("PersistableSourceVector", storageList, additionalData)
    psv = afwDet.PersistableSourceVector.swigConvert(psvptr)
    ss = psv.getSources()
    print 'with', ss.size(), 'sources'
    return ss

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
        Trace("lsst.testing.pipeQA.testDbQueries", 1, "Error: registry (-R) and input root dir (-i) required")
        parser.print_help()
        sys.exit(1)

    database = opt.database
    if (database == None):
        Trace("lsst.testing.pipeQA.testDbQueries", 1, "Error: database (-D) required")
        parser.print_help()
        sys.exit(1)
    
    # database
    dbId        = pipeQA.DatabaseIdentity(database)
    dbInterface = pipeQA.LsstSimDbInterface(dbId)

    # boost 
    butler = pipeQA.PipeQaUtils.getInputButler(inRoot, registry)

    # all combos
    allkeys = pipeQA.PipeQaUtils.getAllKeysOpt(opt, butler)
    print allkeys
