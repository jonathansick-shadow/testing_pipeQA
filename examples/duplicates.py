#!/usr/bin/env python
#
import sys
import os
import re
from lsst.testing.pipeQA.DatabaseQuery import LsstSimDbInterface, DatabaseIdentity

#############################################################
# Main body of code
#############################################################


def main():

    label = "buildbot_PT1_2_u_wp_trunk_2011_0601_171743"
    dbId = DatabaseIdentity(label)
    dbInterface = LsstSimDbInterface(dbId)

    sql_source = "select sce.visit, sce.raftName, sce.ccdName,s.objectId, "\
                 "  s.ra, s.decl, s.xAstrom, s.yAstrom, s.flagForDetection  "\
                 "    from Source as s, Science_Ccd_Exposure as sce  "\
                 "    where (s.scienceCcdExposureId = sce.scienceCcdExposureId) "\
                 "      and (sce.visit = 887252941) "\
                 "      and (sce.raftName = '2,2') "\
                 "      and (sce.ccdName = '1,1')"
    source_results = dbInterface.execute(sql_source)

    sql_match = "select sce.visit, sce.raftName, sce.ccdName, sro.refObjectId, s.objectId,"\
                "  s.ra,s.decl,s.xAstrom,s.yAstrom,s.flagForDetection  "\
                "    from Source as s, Science_Ccd_Exposure as sce, "\
                "      RefSrcMatch as rom, SimRefObject as sro "\
                "    where (s.scienceCcdExposureId = sce.scienceCcdExposureId) "\
                "      and (s.sourceId = rom.sourceId) "\
                "      and (rom.refObjectId = sro.refObjectId) "\
                "      and (rom.closestToSrc = 1)    "\
                "      and (sce.visit = 887252941) "\
                "      and (sce.raftName = '2,2') "\
                "      and (sce.ccdName = '1,1')"
    match_results = dbInterface.execute(sql_match)

    nObjIds = {}
    nRefIds = {}
    for row in match_results:
        refId, objId = row[3], row[4]
        if not nObjIds.has_key(objId):
            nObjIds[objId] = 0
        nObjIds[objId] += 1
        if not nRefIds.has_key(refId):
            nRefIds[refId] = 0
        nRefIds[refId] += 1
        if nRefIds[refId] > 1 and False:
            print "duplicate RefId: ", refId

    orphans = 0
    nSource = {}
    for row in source_results:
        objId = row[3]
        if not nObjIds.has_key(objId):
            orphans += 1
        if not nSource.has_key(objId):
            nSource[objId] = 0
        nSource[objId] += 1

    print "--> sources: ", len(source_results)
    print "unique source objIds: ", len(nSource.keys())
    print ""
    print "--> matches: ", len(match_results)
    print "unique match objIds: ", len(nObjIds.keys())
    print "unique match refIds: ", len(nRefIds.keys())
    print ""
    print "orphans: ", orphans


if __name__ == '__main__':
    main()
