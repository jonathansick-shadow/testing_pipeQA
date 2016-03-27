import sys
import os
import re
import copy
import time
import numpy

import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.cameraGeom as cameraGeom

import CameraInfo as qaCamInfo

from DatabaseQuery import LsstSimDbInterface, DatabaseIdentity
from QaData import QaData

import QaDataUtils as qaDataUtils
import simRefObject as simRefObj
import source as pqaSource


class Source(object):

    def __init__(self):
        self.sourceId = None


class MatchedSource(object):

    def __init__(self):
        self.simRefObjectId = None
        self.u = Source()
        self.g = Source()
        self.r = Source()
        self.i = Source()
        self.z = Source()
        self.y = Source()


class Timer(object):

    def __init__(self):
        self.t0 = {}
        self.t = {}
        self.labels = []

    def start(self, label):
        if not self.t.has_key(label):
            self.labels.append(label)
            self.t[label] = 0.0
        self.t0[label] = time.time()

    def stop(self, label):
        self.t[label] += time.time() - self.t0[label]

    def write(self):
        for label in self.labels:
            print label, "%.2f" % (self.t[label])


#########################################################################
#
#
#
#########################################################################
class DbQaData(QaData):
    # Qa__init__(self, label, rerun, dataInfo):

    def __init__(self, database, rerun, cameraInfo, **kwargs):
        """
        @param database The name of the database to connect to
        @param rerun The data rerun to use
        @param cameraInfo A cameraInfo object describing the camera for these data
        """
        QaData.__init__(self, database, rerun, cameraInfo)
        self.dbId = DatabaseIdentity(self.label)
        self.dbInterface = LsstSimDbInterface(self.dbId)

        self.refStr = {'obj': ('Obj', 'object'), 'src': ('Src', 'source')}

        self.coaddTable = kwargs.get('coaddTable', 'goodSeeing')

        coaddTables = ['goodSeeing', 'chiSquared', 'deep', 'keith']
        if not self.coaddTable in coaddTables:
            raise ValueError, "coaddTable must be on of: %s" % (", ".join(coaddTables))

        if self.coaddTable == 'chiSquared' and cameraInfo.name == 'coadd':
            cameraInfo.setFilterless()

        self.useForced = kwargs.get('useForced', False)
        forced = ''
        if self.useForced:
            forced = 'Forced'

        cTabUpper = self.coaddTable[0].title() + self.coaddTable[1:]

        if type(cameraInfo) == qaCamInfo.CoaddCameraInfo:
            # We need to get the valid dataIds
            sql = "select distinct tract,patch,filterName from %sCoadd where tract=%s" % (cTabUpper, kwargs[
                                                                                          'tract'])
            results = self.dbInterface.execute(sql)
            dataIds = [dict(tract=x[0], patch=map(int, x[1].split(",")), filter=x[2]) for x in results]
            cameraInfo.skyMapToCamera(dataIds)

        # define the tables to use for this camera
        self.sceTables = {
            'lsstSim': 'Science_Ccd_Exposure',
            'sdss': 'Science_Ccd_Exposure',
            'coadd': '%sCoadd' % (cTabUpper),
        }

        self.sceIds = {
            'lsstSim': 'scienceCcdExposureId',
            'sdss': 'scienceCcdExposureId',
            'coadd': '%sCoaddId' % (self.coaddTable),
        }

        self.sTables = {
            'lsstSim': 'Source',
            'sdss': '%s%sSource' % ((cTabUpper, forced) if self.useForced else ('', '')),
            'coadd': '%s%sSource' % (cTabUpper, forced),
        }

        self.sIds = {
            'lsstSim': 'sourceId',
            'sdss': self.coaddTable+'SourceId' if self.useForced else 'sourceId',
            'coadd': '%sSourceId' % (self.coaddTable),
        }

        self.romTables = {
            'lsstSim': 'Ref%sMatch',
            'sdss': 'Ref'+cTabUpper+'%sMatch' if self.useForced else 'Ref%sMatch',
            'coadd': 'Ref%sSrcMatch' % (cTabUpper),
        }
        self.sceReplacements = {
            'lsstSim': {},  # 'fwhm' : 'fwhm',         'scienceCcdExposureId' : 'scienceCcdExposureId' },
            'sdss': {},  # 'fwhm' : 'fwhm',         'scienceCcdExposureId' : 'scienceCcdExposureId' },
            'coadd': {'fwhm': 'measuredFwhm', 'scienceCcdExposureId': '%sCoaddId' % (cTabUpper)},
        }
        self.yMagColumn = {
            'lsstSim': True,
            'sdss': False,
            'coadd': False,  # assuming coadds are SDSS coadds ... will require update eventually.
        }

        defaultCamera = 'lsstSim'
        self.sceTable = self.sceTables.get(cameraInfo.name, defaultCamera)
        self.sceId = self.sceIds.get(cameraInfo.name, defaultCamera)
        self.sTable = self.sTables.get(cameraInfo.name, defaultCamera)
        self.sId = self.sIds.get(cameraInfo.name, defaultCamera)
        self.romTable = self.romTables.get(cameraInfo.name, defaultCamera)
        self.sceReplace = self.sceReplacements.get(cameraInfo.name, defaultCamera)
        self.haveYmag = self.yMagColumn.get(cameraInfo.name, defaultCamera)

    def initCache(self):

        QaData.initCache(self)
        # need to intialize these differently than base class
        # ... Db has 'object' and 'source' matching to be cached
        self.matchListCache = {'obj': {}, 'src': {}}
        self.matchQueryCache = {'obj': {}, 'src': {}}

    def calibFluxError(self, f, df, f0, df0):
        val = numpy.NaN
        try:
            return (df/f + df0/f0)*f/f0
        except FloatingPointError:
            val = numpy.NaN
        return val

    def getMatchListBySensor(self, dataIdRegex, useRef='src'):
        """Get a dict of all SourceMatches matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        calib = self.getCalibBySensor(dataIdRegex)
        sourcesDict = self.getSourceSetBySensor(dataIdRegex)
        refObjectsDict = self.getRefObjectSetBySensor(dataIdRegex)

        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        setMethods = [x for x in qaDataUtils.getSourceSetAccessors()]
        selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames()]
        selectStr = ", ".join(selectList)

        sql = 'select sce.filterId, sce.filterName from '+self.sceTable+' as sce'
        sql += ' where '
        haveAllKeys = True

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
        ]
        nDataId = len(sceNames)

        idWhereList = []
        for keyNames in sceNames:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                likeEqual = self._sqlLikeEqual(sqlName, dataIdRegex[key])
                idWhereList.append(likeEqual)
            else:
                haveAllKeys = False
        idWhere = " and ".join(idWhereList)

        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query

        if not re.search("\%", idWhere) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy, defineFully=True)
            if self.matchListCache[useRef].has_key(key):
                return {key: self.matchListCache[useRef][key]}

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.matchQueryCache[useRef].has_key(dataIdStr):
            matchListDict = {}
            # get only the ones that match the request
            for key, matchList in self.matchListCache[useRef].items():
                if re.search(dataIdStr, key):
                    matchListDict[key] = matchList
            return matchListDict

        sql += idWhere
        result = self.dbInterface.execute(sql)
        filterId, filterName = result[0]

        romTable = self.romTable
        # if the romTable has a %s in it, sprintf the refStr (Obj or Src)
        # otherwise assume it's hardcoded (ie. for RefGoodSeeingSourceMatch table)
        if re.search('%s', romTable):
            romTable = self.romTable % (self.refStr[useRef][0])

        useIndex = ''
        if self.cameraInfo.name == 'sdss':  # and self.useForced:
            useIndex = 'use index()'

        # this will have to be updated for the different dataIdNames when non-lsst cameras get used.

        sql  = 'select ' + ", ".join(zip(*sceNames)[1])+', sro.%sMag, sro.ra, sro.decl, sro.isStar, sro.refObjectId, s.%s, '%(filterName, self.sId) \
            + ' rom.n%sMatches,' % (self.refStr[useRef][0]) \
            + selectStr \
            + '  from '+self.sTable+' as s, '+self.sceTable+' as sce %s,' % (useIndex) \
            + '    '+romTable+' as rom, RefObject as sro' \
            + '  where (s.'+self.sceId+' = sce.'+self.sceId+')' \
            + '    and (s.'+self.sId+' = rom.'+self.sId+') and (rom.refObjectId = sro.refObjectId)'

        if useRef == 'obj':
            sql += '    and (s.objectID is not NULL) '
        sql += '    and '+idWhere

        self.printStartLoad("Loading MatchList (" + self.refStr[useRef][1] + ") for: " + dataIdStr + "...")

        # run the query
        results = self.dbInterface.execute(sql)

        self.sqlCache['match'][dataIdStr] = sql

        # parse results and put them in a sourceSet
        multiplicity = {}
        matchListDict = {}
        for row in results:

            nFields = 7 + nDataId

            mag, ra, dec, isStar, refObjId, srcId, nMatches = row[nDataId:nFields]
            dataIdTmp = {}
            for j in range(nDataId):
                idName = sceNames[j][0]
                dataIdTmp[idName] = row[j]

            key = self._dataIdToString(dataIdTmp, defineFully=True)
            self.dataIdLookup[key] = dataIdTmp

            if not matchListDict.has_key(key):
                refCatObj = pqaSource.RefCatalog()
                refCat = refCatObj.catalog
                catObj = pqaSource.Catalog()
                cat = catObj.catalog

                matchListDict[key] = []

                refRaKey = refCatObj.keyDict['Ra']
                refDecKey = refCatObj.keyDict['Dec']
                refPsfKey = refCatObj.keyDict['PsfFlux']
                refApKey = refCatObj.keyDict['ApFlux']
                refModKey = refCatObj.keyDict['ModelFlux']
                refInstKey = refCatObj.keyDict['InstFlux']

                psfKey = catObj.keyDict['PsfFlux']
                apKey = catObj.keyDict['ApFlux']
                modKey = catObj.keyDict['ModelFlux']
                instKey = catObj.keyDict['InstFlux']

                psfErrKey = catObj.keyDict['PsfFluxErr']
                apErrKey = catObj.keyDict['ApFluxErr']
                modErrKey = catObj.keyDict['ModelFluxErr']
                instErrKey = catObj.keyDict['InstFluxErr']

            matchList = matchListDict[key]

            # reference objects
            sref = refCat.addNew()

            sref.setId(refObjId)
            sref.setD(refRaKey, ra)
            sref.setD(refDecKey, dec)

            # clip at -30
            if mag < -30:
                mag = -30
            flux = 10**(-mag/2.5)

            sref.setD(refPsfKey, flux)
            sref.setD(refApKey, flux)
            sref.setD(refModKey, flux)
            sref.setD(refInstKey, flux)

            # sources
            s = cat.addNew()
            s.setId(srcId)
            s.setD(catObj.keyDict['Extendedness'], 0.0 if isStar else 1.0)

            i = 0
            for value in row[nFields:]:
                if value is not None:
                    setKey = catObj.setKeys[i]
                    if isinstance(value, str):
                        value = 1.0 if ord(value) else 0.0
                    s.setD(setKey, value)
                i += 1

            #sref.setFlagForDetection(sss.getFlagForDetection() | pqaSource.STAR)

            fmag0, fmag0Err = calib[key].getFluxMag0()

            # fluxes
            s.setD(psfKey, s.getD(psfKey)/fmag0)
            s.setD(apKey, s.getD(apKey)/fmag0)
            s.setD(modKey, s.getD(modKey)/fmag0)
            s.setD(instKey, s.getD(instKey)/fmag0)

            # flux errors
            psfFluxErr = qaDataUtils.calibFluxError(s.getD(psfKey), s.getD(psfErrKey),
                                                    fmag0, fmag0Err)
            s.setD(psfErrKey, psfFluxErr)

            apFluxErr = qaDataUtils.calibFluxError(s.getD(psfKey), s.getD(apErrKey),
                                                   fmag0, fmag0Err)
            s.setD(apErrKey, apFluxErr)

            modFluxErr = qaDataUtils.calibFluxError(s.getD(modKey), s.getD(modErrKey),
                                                    fmag0, fmag0Err)
            s.setD(modErrKey, modFluxErr)

            instFluxErr = qaDataUtils.calibFluxError(s.getD(instKey), s.getD(instErrKey),
                                                     fmag0, fmag0Err)
            s.setD(instErrKey, instFluxErr)

            dist = 0.0

            matchList.append([sref, s, dist])
            multiplicity[s.getId()] = nMatches

        ######
        ######
        ######
        # Determine which are orphans, blends, straight matches, and non-detections
        typeDict = {}
        for key in matchListDict.keys():
            matchList = matchListDict[key]

            sources = sourcesDict[key]
            if refObjectsDict.has_key(key):
                refObjects = refObjectsDict[key]
            else:
                refObjects = simRefObj.SimRefObjectSet()  # an empty set

            typeDict[key] = {}

            refIds = []
            for ro in refObjects:
                refIds.append(ro.getId())

            srcIds = []
            for so in sources:
                srcIds.append(so.getId())

            matRef = []
            matSrc = []
            for ma in matchList:
                matRef.append(ma[0].getId())
                matSrc.append(ma[1].getId())

            refIds = set(refIds)
            srcIds = set(srcIds)
            matRef = set(matRef)
            matSrc = set(matSrc)

            undetectedIds = refIds - matRef
            orphanIds = srcIds - matSrc
            matchedIds = srcIds & matSrc   # Does not know about duplicates
            # print 'Undet, orphan, matched:', len(undetectedIds), len(orphanIds), len(matchedIds)

            undetected = []
            orphans = []
            matched = []
            blended = []
            for ro in refObjects:
                if ro.getId() in undetectedIds:
                    undetected.append(ro)
            matchListById = dict([(m[1].getId(), m) for m in matchList])
            matchIdsSet = set(matchListById.keys())
            for so in sources:
                soid = so.getId()
                if soid in orphanIds:
                    orphans.append(so)
                if soid in matchedIds and soid in matchIdsSet:
                    if multiplicity[soid] == 1:
                        matched.append(matchListById[soid])
                    else:
                        #print -2.5*numpy.log10(so.getD(psfKey)), multiplicity[soid]
                        blended.append(matchListById[soid])

            self.printMidLoad('\n        %s: Undet, orphan, matched, blended = %d %d %d %d' % (
                key, len(undetected), len(orphans), len(matched), len(blended))
            )

            typeDict[key]['orphan'] = orphans
            typeDict[key]['matched'] = matched
            typeDict[key]['blended'] = blended
            typeDict[key]['undetected'] = undetected

            typeDict[key]['sql'] = sql

            # cache it
            self.matchListCache[useRef][key] = typeDict[key]

            # Determine which are orphans, blends, straight matches, and non-detections
            ######
            ######
            ######

        # cache it
        self.matchQueryCache[useRef][dataIdStr] = True

        self.printStopLoad()

        return typeDict

    def getSourceSetBySensor(self, dataIdRegex):
        """Get a dict of all Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        setMethods = [x for x in qaDataUtils.getSourceSetAccessors()]
        selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames()]
        selectStr = ", ".join(selectList)

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
        ]

        # this will have to be updated for the different dataIdNames when non-lsst cameras get used.
        #sql  = 'select '+", ".join(zip(*sceNames)[1])+', s.'+self.sId+', '+selectStr
        #sql += '  from '+self.sTable+' as s, '+self.sceTable+' as sce'
        #sql += '  where (s.'+self.sceId+' = sce.'+self.sceId+')'
        # above query too slow on DC_W13_Stripe82, rewrote using subquery below
        # Note that qserv will not support subqueries
        sql  = 'select '+", ".join(zip(*sceNames)[1])+', s.'+self.sId+', '+selectStr \
            + '  from '+self.sTable+' as s, '+'(select ' + ", ".join([w.replace('sce.', '') for w in zip(*sceNames)[1]]) \
            + ', '+self.sceId+' from '+self.sceTable

        haveAllKeys = True

        whereList = []
        for keyNames in sceNames:
            key, sqlName = [w.replace('sce.', '') for w in keyNames]
            if dataIdRegex.has_key(key):
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
        if whereList:
            sql += " where %s" % (" and ".join(whereList))
        else:
            haveAllKeys = False

        sql += ') as sce  where (s.'+self.sceId+' = sce.'+self.sceId+')'

        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query

        if not re.search("\%", sql) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy, defineFully=True)
            if self.sourceSetCache.has_key(key):
                return {key: self.sourceSetCache[key]}

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.queryCache.has_key(dataIdStr):
            ssDict = {}
            # get only the ones that match the request
            for key, ss in self.sourceSetCache.items():
                if re.search(dataIdStr, key):
                    ssDict[key] = ss
            return ssDict

        self.queryCache[dataIdStr] = True

        self.printStartLoad("Loading SourceSets for: " + dataIdStr + "...")

        # run the query
        results = self.dbInterface.execute(sql)
        self.sqlCache['src'][dataIdStr] = sql
        calib = self.getCalibBySensor(dataIdRegex)

        # parse results and put them in a sourceSet
        ssDict = {}
        for k in calib.keys():
            catObj = pqaSource.Catalog()
            ssDict[k] = catObj.catalog

            psfKey = catObj.keyDict['PsfFlux']
            apKey = catObj.keyDict['ApFlux']
            modKey = catObj.keyDict['ModelFlux']
            instKey = catObj.keyDict['InstFlux']

            psfErrKey = catObj.keyDict['PsfFluxErr']
            apErrKey = catObj.keyDict['ApFluxErr']
            modErrKey = catObj.keyDict['ModelFluxErr']
            instErrKey = catObj.keyDict['InstFluxErr']

        for row in results:

            # get the values for the dataId
            i = 0
            dataIdTmp = {}
            for idName, dbName in sceNames:
                dataIdTmp[idName] = row[i]
                i += 1
            sid = row[i]
            nIdKeys = i+1

            key = self._dataIdToString(dataIdTmp, defineFully=True)
            self.dataIdLookup[key] = dataIdTmp

            s = ssDict[key].addNew()

            s.setId(sid)

            i = 0
            for value in row[nIdKeys:]:
                if value is not None:
                    setKey = catObj.setKeys[i]
                    # print value, type(value)
                    if isinstance(value, str) and len(value) == 1:
                        value = 1.0 if ord(value) else 0.0
                    s.setD(setKey, value)
                    # print catObj.setNames[i], value
                i += 1

            # calibrate it
            fmag0, fmag0Err = calib[key].getFluxMag0()

            if (fmag0 == 0.0):
                continue

            pf, af, sf = s.getD(psfKey), s.getD(apKey), s.getD(modKey)

            # if pf > 1000.0:
            #    print dataIdTmp
            #    print "%8.2f %8.2f %8.2f    %8.4f %8.4f" % (pf, af, sf, pf/af - 1.0, pf/sf - 1.0)

            # fluxes
            s.setD(psfKey, s.getD(psfKey)/fmag0)
            s.setD(apKey, s.getD(apKey)/fmag0)
            s.setD(modKey, s.getD(modKey)/fmag0)
            s.setD(instKey, s.getD(instKey)/fmag0)

            # flux errors
            psfFluxErr = qaDataUtils.calibFluxError(s.getD(psfKey), s.getD(psfErrKey),
                                                    fmag0, fmag0Err)
            s.setD(psfErrKey, psfFluxErr)

            apFluxErr = qaDataUtils.calibFluxError(s.getD(psfKey), s.getD(apErrKey),
                                                   fmag0, fmag0Err)
            s.setD(apErrKey, apFluxErr)

            modFluxErr = qaDataUtils.calibFluxError(s.getD(modKey), s.getD(modErrKey),
                                                    fmag0, fmag0Err)
            s.setD(modErrKey, modFluxErr)

            instFluxErr = qaDataUtils.calibFluxError(s.getD(instKey), s.getD(instErrKey),
                                                     fmag0, fmag0Err)
            s.setD(instErrKey, instFluxErr)

        # cache it
        for k, ss in ssDict.items():
            self.sourceSetCache[k] = ssDict[k]

        self.printStopLoad()

        return ssDict

    def getDataIdsFromRegex(self, dataIdRegex):

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
        ]

        haveAllKeys = True
        sqlDataId = []
        for keyNames in sceNames:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        sqlDataId = " and ".join(sqlDataId)

        sql  = "select "+", ".join(zip(*sceNames)[1]) \
            + "  from "+self.sceTable+" as sce " \
            + "  where " + sqlDataId

        dataIdList = []
        results = self.dbInterface.execute(sql)
        nIds = len(sceNames)
        for r in results:
            dataId = {}
            i = 0
            for idName, dbName in sceNames:
                dataId[idName] = r[i]
                i += 1

            dataIdList.append(dataId)

        return dataIdList

    def getVisitMatchesBySensor(self, matchDatabase, matchVisit, dataIdRegex):
        """ Get a dict of all Catalog Sources matching dataId, but
        within another Science_Ccd_Exposure's polygon"""

        # If the dataIdEntry is identical to an earlier query, we must already have all the data
        # E.g. visit862826551-snap.*-raft.*-sensor.*
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)

        if self.visitMatchQueryCache.has_key(matchDatabase):
            if self.visitMatchQueryCache[matchDatabase].has_key(matchVisit):
                vmqCache = self.visitMatchQueryCache[matchDatabase][matchVisit]

                if vmqCache.has_key(dataIdStr):
                    vmCache = self.visitMatchCache[matchDatabase][matchVisit]
                    vmDict = {}
                    for key, ss in vmCache.items():
                        if re.search(dataIdStr, key):
                            vmDict[key] = ss
                    return vmDict

        # Load each of the dataIds
        dataIdList = self.getDataIdsFromRegex(dataIdRegex)

        # Set up the outputs
        calib = self.getCalibBySensor(dataIdRegex)
        vmDict = {}
        for k in calib.keys():
            vmDict[k] = []

        for dataIdEntry in dataIdList:
            visit, raft, sensor = dataIdEntry['visit'], dataIdEntry['raft'], dataIdEntry['sensor']
            # E.g. visit862826551-snap0-raft30-sensor20
            dataIdEntryStr = self._dataIdToString(dataIdEntry, defineFully=True)

            haveAllKeys = True
            sqlDataId = []
            for keyNames in [['visit', 'sce.visit'], ['raft', 'sce.raftName'], ['sensor', 'sce.ccdName']]:
                key, sqlName = keyNames
                if dataIdEntry.has_key(key):
                    sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdEntry[key]))
                else:
                    haveAllKeys = False
            sqlDataId = " and ".join(sqlDataId)

            # Poly comes from our own database
            sql1  = 'SELECT poly FROM '+self.sceTable+' as sce ' \
                + 'WHERE %s ' % (sqlDataId) \
                + 'INTO @poly;'

            sql2 = 'CALL scisql.scisql_s2CPolyRegion(@poly, 20);'

            # Selection of source matches from the comparison database
            self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)
            setMethods = ["set"+x for x in qaDataUtils.getSourceSetAccessors()]
            selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames()]
            selectStr = ", ".join(selectList)
            sql3  = 'SELECT sce.visit, sce.raftName, sce.ccdName, sce.filterName, ' \
                + ' sce.fluxMag0, sce.fluxMag0Sigma,'                               \
                + '   CASE WHEN sce.filterId = 0 THEN sro.uMag' \
                + '        WHEN sce.filterId = 1 THEN sro.gMag' \
                + '        WHEN sce.filterId = 2 THEN sro.rMag' \
                + '        WHEN sce.filterId = 3 THEN sro.iMag' \
                + '        WHEN sce.filterId = 4 THEN sro.zMag' \
                + '        WHEN sce.filterId = 5 THEN sro.yMag' \
                + '   END as mag,'                              \
                + ' sro.ra, sro.decl, sro.isStar, sro.refObjectId,' \
                + selectStr \
                + ' FROM %s.'+self.sTable+' AS s USE INDEX FOR JOIN(IDX_htmId20)' % (matchDatabase) \
                + ' INNER JOIN %s.'+self.sceTable+' AS sce ' % (matchDatabase) \
                + ' ON (s.'+self.sceId+' = sce.'+self.sceId+') AND (sce.visit = %s)' % (matchVisit) \
                + '   INNER JOIN %s.'+self.romTable+' AS rsm ON (s.'+self.sId+' = rsm.'+self.sId+')' % (matchDatabase) \
                + '   INNER JOIN %s.RefObject AS sro ON (sro.refObjectId = rsm.refObjectId)'  % (matchDatabase) \
                + '   INNER JOIN scisql.Region AS reg ON (s.htmId20 BETWEEN reg.htmMin AND reg.htmMax) ' \
                + 'WHERE scisql_s2PtInCPoly(s.ra, s.decl, @poly) = 1;'

            # if not re.search("\%", sql1) and haveAllKeys:
            #    dataIdCopy = copy.copy(dataIdEntry)
            #    dataIdCopy['snap'] = "0"
            #    key = self._dataIdToString(dataIdCopy)
            #    if self.visitMatchCache.has_key(key):
            #        vmDict[key] = self.visitMatchCache[key]
            #        continue

            self.printStartLoad("Loading DatasetMatches for: " + dataIdEntryStr + "...")
            self.dbInterface.execute(sql1)
            self.dbInterface.execute(sql2)
            results = self.dbInterface.execute(sql3)

            self.printMidLoad("Found %d matches..." % (len(results)))

            for row in results:
                s = pqaSource.Source()
                qaDataUtils.setSourceBlobsNone(s)
                sref = pqaSource.RefSource()
                qaDataUtils.setSourceBlobsNone(sref)

                nValues = 11
                mvisit, mraft, mccd, mfilt, fmag0, fmag0Err, mag, ra, dec, isStar, refObjId = row[:nValues]
                filt = afwImage.Filter(mfilt, True)

                sref.setId(refObjId)
                sref.setRa(ra)
                sref.setDec(dec)
                flux = 10**(-mag/2.5)
                sref.setPsfFlux(flux)
                sref.setApFlux(flux)
                sref.setModelFlux(flux)
                sref.setInstFlux(flux)

                i = 0
                for value in row[nValues:]:
                    method = getattr(s, setMethods[i])
                    if value is not None:
                        method(value)
                    i += 1

                for sss in [s, sref]:
                    if isStar == 1:
                        sss.setFlagForDetection(sss.getFlagForDetection() | pqaSource.STAR)
                    else:
                        sss.setFlagForDetection(sss.getFlagForDetection() & ~pqaSource.STAR)

                # fluxes
                s.setPsfFlux(s.getPsfFlux()/fmag0)
                s.setApFlux(s.getApFlux()/fmag0)
                s.setModelFlux(s.getModelFlux()/fmag0)
                s.setInstFlux(s.getInstFlux()/fmag0)

                # flux errors
                psfFluxErr = qaDataUtils.calibFluxError(s.getPsfFlux(), s.getPsfFluxErr(), fmag0, fmag0Err)
                s.setPsfFluxErr(psfFluxErr)

                apFluxErr = qaDataUtils.calibFluxError(s.getApFlux(), s.getApFluxErr(), fmag0, fmag0Err)
                s.setApFluxErr(apFluxErr)

                modFluxErr = qaDataUtils.calibFluxError(
                    s.getModelFlux(), s.getModelFluxErr(), fmag0, fmag0Err)
                s.setModelFluxErr(modFluxErr)

                instFluxErr = qaDataUtils.calibFluxError(s.getInstFlux(), s.getInstFluxErr(), fmag0, fmag0Err)
                s.setInstFluxErr(instFluxErr)

                vm = vmDict[dataIdEntryStr]
                vm.append([sref, s, filt])

            self.printStopLoad()

        # cache it
        if not self.visitMatchQueryCache.has_key(matchDatabase):
            self.visitMatchQueryCache[matchDatabase] = {}
            self.visitMatchCache[matchDatabase] = {}

        if not self.visitMatchQueryCache[matchDatabase].has_key(matchVisit):
            self.visitMatchQueryCache[matchDatabase][matchVisit] = {}
            self.visitMatchCache[matchDatabase][matchVisit] = {}

        self.visitMatchQueryCache[matchDatabase][matchVisit][dataIdStr] = True
        for k, ss in vmDict.items():
            self.visitMatchCache[matchDatabase][matchVisit][k] = ss

        return vmDict

    def getRefObjectSetBySensor(self, dataIdRegex):
        """Get a dict of all Catalog Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
        ]

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        sroFields = simRefObj.fields
        if not self.haveYmag:
            sroFields = [x for x in sroFields if x != 'yMag']
        sroFieldStr = ", ".join(["sro."+field for field in sroFields])

        # if the dataIdEntry is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.refObjectQueryCache.has_key(dataIdStr):
            # get only the ones that match the request
            sroDict = {}
            for key, sro in self.refObjectCache.items():
                if re.search(dataIdStr, key):
                    sroDict[key] = sro
            return sroDict

        # get a list of matching dataIds
        dataIdList = self.getDataIdsFromRegex(dataIdRegex)

        # Load each of the dataIds
        sroDict = {}
        for dataIdEntry in dataIdList:

            dataIdEntryStr = self._dataIdToString(dataIdEntry, defineFully=True)

            haveAllKeys = True
            sqlDataId = []
            for keyNames in sceNames:
                key, sqlName = keyNames
                if dataIdEntry.has_key(key):
                    sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdEntry[key]))
                else:
                    haveAllKeys = False
            sqlDataId = " and ".join(sqlDataId)

            # if there are no regexes (ie. actual wildcard expressions),
            #  we can check the cache, otherwise must run the query

            if not re.search("\%", sqlDataId) and haveAllKeys:
                dataIdCopy = copy.copy(dataIdEntry)

                key = self._dataIdToString(dataIdCopy, defineFully=True)
                if self.refObjectCache.has_key(key):
                    sroDict[key] = self.refObjectCache[key]
                    continue

            self.printStartLoad("Loading RefObjects for: " + dataIdEntryStr + "...")

            # usePoly is currently hardwired to use scisql polygon search function, but migration to
            # qserv database will require use of a different function.  The code inside the else
            # statement below is an example of the qserv implementation.  Once we migrate to
            # qserv, will remove or disable the usePoly block
            usePoly = True
            if usePoly:
                sql  = 'SELECT scisql_s2CPolyToBin(' \
                    + '   sce.corner1Ra, sce.corner1Decl, ' \
                    + '   sce.corner2Ra, sce.corner2Decl, ' \
                    + '   sce.corner3Ra, sce.corner3Decl, ' \
                    + '   sce.corner4Ra, sce.corner4Decl) ' \
                    + 'FROM '+self.sceTable+' as sce ' \
                    + 'WHERE %s ' % (sqlDataId) \
                    + 'INTO @poly; '
                self.dbInterface.execute(sql)

                sql2 = 'SELECT %s ' % (sroFieldStr) \
                    + 'FROM ' \
                    + '    RefObject AS sro ' \
                    + 'WHERE ' \
                    + '    (scisql_s2PtInCPoly(sro.ra, sro.decl, @poly) = 1) '
                results = self.dbInterface.execute(sql2)

            else:
                sql  = 'SELECT ' \
                    + '   sce.corner1Ra, sce.corner1Decl, ' \
                    + '   sce.corner2Ra, sce.corner2Decl, ' \
                    + '   sce.corner3Ra, sce.corner3Decl, ' \
                    + '   sce.corner4Ra, sce.corner4Decl  ' \
                    + 'FROM '+self.sceTable+' as sce ' \
                    + 'WHERE %s ' % (sqlDataId)

                corners = self.dbInterface.execute(sql)
                raLL, decLL, raLR, decLR, raUR, decUR, raUL, decUL = corners[0]

                sql2 = 'SELECT %s ' % (sroFieldStr) \
                    + 'FROM ' \
                    + '    RefObject AS sro ' \
                    + 'WHERE ' \
                    + "    qserv_areaspec_poly(sro.ra, sro.decl, %f, %f,   %f, %f,   %f, %f,   %f, %f);" % (
                    raLL, decLL, raLR, decLR, raUR, decUR, raUL, decUL)
                results = self.dbInterface.execute(sql2)

            # parse results and put them in a sourceSet
            raftName, ccdName = self.cameraInfo.getRaftAndSensorNames(dataIdEntry)
            #visit, raft, sensor = dataIdEntry['visit'], dataIdEntry['raft'], dataIdEntry['sensor']
            wcs = self.getWcsBySensor(dataIdEntry)[dataIdEntryStr]
            #raftName = self.cameraInfo.raftKeyToName(raft)
            #ccdName = self.cameraInfo.ccdKeyToName(sensor)
            bbox = self.cameraInfo.getBbox(raftName, ccdName)

            for row in results:
                sroStuff = list(row[:])
                if not self.haveYmag:
                    sroStuff.append(0.0)  # dummy yMag

                # ignore things near the edge
                # ... they wouldn't be detected, and we should know about them
                ra, dec = sroStuff[2], sroStuff[3]
                x, y = wcs.skyToPixel(afwCoord.Coord(afwGeom.PointD(ra, dec)))
                if qaDataUtils.atEdge(bbox, x, y):
                    continue

                dataIdTmp = dataIdEntry  # {'visit':str(visit), 'raft':raft, 'sensor':sensor, 'snap':'0'}
                key = self._dataIdToString(dataIdTmp, defineFully=True)
                self.dataIdLookup[key] = dataIdTmp

                if not sroDict.has_key(key):
                    sroDict[key] = simRefObj.SimRefObjectSet()
                sros = sroDict[key]
                sros.push_back(simRefObj.SimRefObject(*sroStuff))

            self.refObjectQueryCache[dataIdStr] = True

            self.printStopLoad()

        # cache it
        for k, sro in sroDict.items():
            self.refObjectCache[k] = sroDict[k]

        return sroDict

    def getVisits(self, dataIdRegex):
        """ Return explicit visits matching for a dataIdRegex.

        @param dataIdRegex dataId dict containing regular expressions of data to retrieve.
        """

        columnsTmp = zip(*self.cameraInfo.dataInfo)[0]
        visitLike = zip(*self.cameraInfo.dataInfo)[1]
        dbNames = []
        columns = []
        for i in range(len(columnsTmp)):
            c = columnsTmp[i]
            n = visitLike[i]
            if not re.search("snap", c) and n > 0:
                columns.append(self.cameraInfo.dataIdDbNames[c])
                dbNames.append([c, self.cameraInfo.dataIdDbNames[c]])

        sql = "select distinct "+", ".join(columns)+" from "+self.sceTable
        sql += "   where "
        haveAllKeys = True

        whereList = []
        for keyNames in dbNames:
            key, sqlName = keyNames
            if re.search("snap", key):
                continue
            if dataIdRegex.has_key(key):
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        sql += " and ".join(whereList)

        results = self.dbInterface.execute(sql)

        visits = []
        for r in results:
            dataId = {}
            for i in range(len(columns)):
                dataId[columns[i]] = str(r[i])
            visits.append(self.cameraInfo.dataIdCameraToStandard(dataId)['visit'])

        return sorted(set(visits))

    def breakDataId(self, dataIdRegex, breakBy):
        """Take a dataId with regexes and return a list of dataId regexes
        which break the dataId by raft, or ccd.

        @param dataId    ... to be broken
        @param breakBy   'visit', 'raft', or 'ccd'
        """

        if not re.search("(visit|raft|ccd)", breakBy):
            raise Exception("breakBy must be 'visit','raft', or 'ccd'")

        if breakBy == 'visit':
            return [dataIdRegex]

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
        ]

        sql = "select "+", ".join(zip(*sceNames)[1])+" from "+self.sceTable+" as sce"
        sql += "   where "
        whereList = []
        for keyNames in sceNames:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
        sql += " and ".join(whereList)

        results = self.dbInterface.execute(sql)

        dataIdDict = {}
        #ccdConvention = 'ccd'
        # if not dataIdRegex.has_key('ccd'):
        #    ccdConvention = 'sensor'

        for r in results:

            i = 0
            thisDataId = {}
            for idName, dbName in sceNames:
                thisDataId[idName] = r[i]
                i += 1

            if breakBy == 'raft':
                # handle lsst/hsc different naming conventions
                ccd = dataIdRegex[self.cameraInfo.dataIdTranslationMap['sensor']]

            key = self._dataIdToString(thisDataId, defineFully=True)
            dataIdDict[key] = thisDataId

        # store the list of broken dataIds
        self.brokenDataIdList = []
        for key in sorted(dataIdDict.keys()):
            self.brokenDataIdList.append(dataIdDict[key])

        return copy.copy(self.brokenDataIdList)

    def loadCalexp(self, dataIdRegex):
        """Load the calexp data for data matching dataIdRegex.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceDataIdNames = [
            #[x[0], "sce."+x[1]]
            x for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
        ]

        selectList = ["sce."+x for x in qaDataUtils.getSceDbNames(sceDataIdNames, self.sceReplace)]
        selectStr = ", ".join(selectList)

        sql  = 'select '+selectStr \
            + '  from '+self.sceTable+' as sce' \
            + '  where '

        haveAllKeys = True

        whereList = []
        for keyNames in [[x[0], "sce."+x[1]] for x in sceDataIdNames]:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        sql += " and ".join(whereList)

        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query
        if not re.search("\%", sql) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy, defineFully=True)
            if self.calexpQueryCache.has_key(key):
                return

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.calexpQueryCache.has_key(dataIdStr) and self.calexpQueryCache[dataIdStr]:
            return

        self.printStartLoad("Loading Calexp for: " + dataIdStr + "...")

        # run the query
        results = self.dbInterface.execute(sql)

        for row in results:

            rowDict = dict(zip(qaDataUtils.getSceDbNames(sceDataIdNames, self.sceReplace), row))

            for k, v in self.sceReplace.items():
                if rowDict.has_key(v):
                    rowDict[k] = rowDict[v]
                    del rowDict[v]

            dataIdTmp = {}
            for idName, dbName in sceDataIdNames:
                dataIdTmp[idName] = rowDict[dbName]

            #visit, raft, sensor = rowDict['visit'], rowDict['raftName'], rowDict['ccdName']
            #dataIdTmp = {'visit':visit, 'raft':raft, 'sensor':sensor, 'snap':'0'}
            key = self._dataIdToString(dataIdTmp, defineFully=True)
            self.dataIdLookup[key] = dataIdTmp

            # print rowDict
            if not self.wcsCache.has_key(key):
                crval = afwCoord.Coord(afwGeom.PointD(rowDict['crval1'], rowDict['crval2']))
                crpix = afwGeom.PointD(rowDict['crpix1'], rowDict['crpix2'])
                cd11, cd12, cd21, cd22 = rowDict['cd1_1'], rowDict[
                    'cd1_2'], rowDict['cd2_1'], rowDict['cd2_2']
                wcs = afwImage.makeWcs(crval, crpix, cd11, cd12, cd21, cd22)
                self.wcsCache[key] = wcs

            if not self.detectorCache.has_key(key):
                raftName, ccdName = self.cameraInfo.getRaftAndSensorNames(dataIdTmp)
                if self.cameraInfo.detectors.has_key(ccdName):
                    self.detectorCache[key] = self.cameraInfo.detectors[ccdName]  # ccdDetector
                if self.cameraInfo.detectors.has_key(raftName):
                    self.raftDetectorCache[key] = self.cameraInfo.detectors[raftName]

            if not self.filterCache.has_key(key):
                filt = afwImage.Filter(rowDict['filterName'], True)
                self.filterCache[key] = filt

            if not self.calibCache.has_key(key):
                calib = afwImage.Calib()
                calib.setFluxMag0(rowDict['fluxMag0'], rowDict['fluxMag0Sigma'])
                self.calibCache[key] = calib

            self.calexpCache[key] = rowDict
            self.calexpQueryCache[key] = True

        self.calexpQueryCache[dataIdStr] = True

        self.printStopLoad()

    def getCalexpEntryBySensor(self, cache, dataIdRegex):
        """Fill and return the dict for a specified calexp cache.

        @param cache The cache dictionary to return
        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # get the datasets corresponding to the request
        self.loadCalexp(dataIdRegex)
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        entryDict = {}
        for dataKey in cache.keys():
            if re.search(dataIdStr, dataKey):
                entryDict[dataKey] = cache[dataKey]
        return entryDict

    def getSourceSet(self, dataIdRegex):
        """Get a SourceSet of all Sources matching dataId.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        ssDict = self.getSourceSetBySensor(dataIdRegex)
        ssReturn = []
        for key, ss in ssDict.items():
            ssReturn += ss

        return ssReturn

    ########################################
    #
    # utility method to generate a string for an sql 'where' clause
    # given 'field' and 'regex':
    #    - if regex is just a value - return "field = regex"
    #    - if regex is a regex - sub '.*' and '?' to '%' and return "field like 'regex'"
    ########################################
    def _sqlLikeEqual(self, field, regex):
        """Utility to convert a dataId regex to an sql 'where' clause.
        """

        regex = str(regex)

        clause = ""
        # just a number
        if re.search('^\d+$', regex):
            clause += field + " = %s" % (regex)
        # comma-sep numbers
        elif re.search('^[\d,]+$', regex):
            clause += field + " = '%s'" % (regex)
        # if it's a filter
        elif re.search('^[ugrizy]$', regex):
            clause += field + " = '%s'" % (regex)
        # .*  ?  % followed/preceeding by comma-sep numbers
        elif re.search('^(\.\*|\?|\%)[\d,]*$', regex) or re.search('^[\d,]*(\.\*|\?|\%)$', regex):
            regexSql = re.sub("(\.\*|\?)", "%", regex)
            clause += field + " like '%s'" % (regexSql)
        else:
            raise Exception("Regex for SQL can only use '.*', '?', or '%' at beginning or end of string (" +
                            field+"="+regex+")")

        return "("+clause+")"


###################################################
# Factory for dbQaData
# - curently only lsstSim is available by database, so this is a trivial factory
###################################################
def makeDbQaData(label, rerun=None, camera=None, **kwargs):
    """Factory for a DbQaData object.

    @param database The name of the database to connect to
    @param rerun The data rerun to use
    """

    cameraInfos = {
        #       "cfht": qaCamInfo.CfhtCameraInfo(), # XXX CFHT camera geometry is currently broken following #1767
        "hsc": qaCamInfo.HscCameraInfo(),
        "suprimecam": qaCamInfo.SuprimecamCameraInfo(),
        "suprimecam-old": qaCamInfo.SuprimecamCameraInfo(True),
        "sdss": qaCamInfo.SdssCameraInfo(),
        "coadd": qaCamInfo.CoaddCameraInfo(**kwargs),
        "lsstSim": qaCamInfo.LsstSimCameraInfo(),
    }

    cameraToUse = None
    if camera is not None:
        cameraToUse = cameraInfos[camera]
    else:
        cameraToUse = cameraInfos['lsstSim']

    return DbQaData(label, rerun, cameraToUse, **kwargs)
