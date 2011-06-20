import sys, os, re, copy

import lsst.afw.detection               as afwDet
import lsst.afw.image                   as afwImage
import lsst.meas.algorithms             as measAlg
import lsst.afw.geom                    as afwGeom
import lsst.afw.cameraGeom              as cameraGeom

import CameraInfo                       as qaCamInfo

from DatabaseQuery import LsstSimDbInterface, DatabaseIdentity
from QaData        import QaData

import QaDataUtils as qaDataUtils
import simRefObject as simRefObj

def mem(size="rss"):
    """Generalization; memory sizes: rss, rsz, vsz."""
    return int(os.popen('ps -p %d -o %s | tail -1' % (os.getpid(), size)).read())

#########################################################################
#
#
#
#########################################################################
class DbQaData(QaData):
    #Qa__init__(self, label, rerun, dataInfo):

    def __init__(self, database, rerun, cameraInfo):
        """
        @param database The name of the database to connect to
        @param rerun The data rerun to use
        @param cameraInfo A cameraInfo object describing the camera for these data
        """
        QaData.__init__(self, database, rerun, cameraInfo)
        self.dbId        = DatabaseIdentity(self.label)
        self.dbInterface = LsstSimDbInterface(self.dbId)

        self.refStr = {'obj' : ('Obj', 'object'), 'src' : ('Src', 'source') }

        

    def initCache(self):

        QaData.initCache(self)
        # need to intialize these differently than base class
        # ... Db has 'object' and 'source' matching to be cached
        self.matchListCache = { 'obj': {}, 'src': {} }
        self.matchQueryCache = { 'obj' : {}, 'src': {} }
        

    def getMatchListBySensor(self, dataIdRegex, useRef='src'):
        """Get a dict of all SourceMatches matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        setMethods = ["set"+x for x in qaDataUtils.getSourceSetAccessors()]
        selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames()]
        selectStr = ",".join(selectList)
        
        sql  = 'select sce.filterId, sce.filterName from Science_Ccd_Exposure as sce'
        sql += ' where '
        haveAllKeys = True

        idWhereList = []
        for keyNames in [['visit', 'sce.visit'], ['raft', 'sce.raftName'], ['sensor', 'sce.ccdName']]:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                idWhereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        idWhere = " and ".join(idWhereList)

        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query
        
        if not re.search("\%", idWhere) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy)
            if self.matchListCache[useRef].has_key(key):
                return {key : self.matchListCache[useRef][key]}

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex)
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

        
        # this will have to be updated for the different dataIdNames when non-lsst cameras get used.
        sql  = 'select sce.visit, sce.raftName, sce.ccdName, sro.%sMag, sro.ra, sro.decl, sro.isStar, sro.refObjectId, '%(filterName)
        sql += selectStr
        sql += '  from Source as s, Science_Ccd_Exposure as sce,'
        sql += '    Ref%sMatch as rom, SimRefObject as sro' % (self.refStr[useRef][0])
        sql += '  where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
        sql += '    and (s.%sId = rom.%sId) and (rom.refObjectId = sro.refObjectId)' % \
               (self.refStr[useRef][1], self.refStr[useRef][1])
        if useRef == 'obj':
            sql += '    and (s.objectID is not NULL) '
            #sql += '    and (rom.nRefMatches = 1) '
            sql += '    and (rom.nObjMatches = 1) '
        else:
            sql += '    and (rom.nSrcMatches = 1) '
            #sql += '    and (rom.nRefMatches = 1) '
        sql += '    and '+idWhere
        

        self.printStartLoad("Loading MatchList ("+ self.refStr[useRef][1]  +") for: " + dataIdStr + "...")


        # run the query
        results  = self.dbInterface.execute(sql)

        # parse results and put them in a sourceSet
        matchListDict = {}
        for row in results:
            s = afwDet.Source()
            qaDataUtils.setSourceBlobsNone(s)
            sref = afwDet.Source()
            qaDataUtils.setSourceBlobsNone(sref)
            
            visit, raft, sensor, mag, ra, dec, isStar, refObjId = row[0:8]
            dataIdTmp = {'visit':str(visit), 'raft':raft, 'sensor':sensor, 'snap':'0'}
            key = self._dataIdToString(dataIdTmp)
            self.dataIdLookup[key] = dataIdTmp

            if not matchListDict.has_key(key):
                matchListDict[key] = []
            matchList = matchListDict[key]

            sref.setId(refObjId)
            sref.setRa(ra)
            sref.setDec(dec)
            flux = 10**(-mag/2.5)
            sref.setPsfFlux(flux)
            sref.setApFlux(flux)
            sref.setModelFlux(flux)
            
            i = 0
            for value in row[8:]:
                method = getattr(s, setMethods[i])
                if not value is None:
                    method(value)
                i += 1

            for sss in [s, sref]:
                if isStar == 1:
                    sss.setFlagForDetection(sss.getFlagForDetection() | measAlg.Flags.STAR)
                else:
                    sss.setFlagForDetection(sss.getFlagForDetection() & ~measAlg.Flags.STAR)

            # calibrate it
            calib = self.getCalibBySensor(dataIdTmp)
            fmag0, fmag0Err = calib[key].getFluxMag0()
            s.setPsfFlux(s.getPsfFlux()/fmag0)
            s.setApFlux(s.getApFlux()/fmag0)
            s.setModelFlux(s.getModelFlux()/fmag0)

            dist = 0.0
            matchList.append([sref, s, dist])
            
        # cache it
        self.matchQueryCache[useRef][dataIdStr] = True
        for k, matchList in matchListDict.items():
            self.matchListCache[useRef][k] = matchListDict[k]

        self.printStopLoad()
        
        return matchListDict



    def getSourceSetBySensor(self, dataIdRegex):
        """Get a dict of all Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        setMethods = ["set"+x for x in qaDataUtils.getSourceSetAccessors()]
        selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames()]
        selectStr = ",".join(selectList)

        
        # this will have to be updated for the different dataIdNames when non-lsst cameras get used.
        sql  = 'select sce.visit, sce.raftName, sce.ccdName,'+selectStr
        sql += '  from Source as s, Science_Ccd_Exposure as sce'
        sql += '  where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
        haveAllKeys = True

        for keyNames in [['visit', 'sce.visit'], ['raft', 'sce.raftName'], ['sensor', 'sce.ccdName']]:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                sql += '    and '+self._sqlLikeEqual(sqlName, dataIdRegex[key])
            else:
                haveAllKeys = False
        
        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query
        
        if not re.search("\%", sql) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy)
            if self.sourceSetCache.has_key(key):
                return {key : self.sourceSetCache[key]}

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex)
        if self.queryCache.has_key(dataIdStr):
            ssDict = {}
            # get only the ones that match the request
            for key, ss in self.sourceSetCache.items():
                if re.search(dataIdStr, key):
                    ssDict[key] = ss
            return ssDict

        self.printStartLoad("Loading SourceSets for: " + dataIdStr + "...")

        # run the query
        results  = self.dbInterface.execute(sql)

        # parse results and put them in a sourceSet
        ssDict = {}
        for row in results:
            s = afwDet.Source()
            qaDataUtils.setSourceBlobsNone(s)

            visit, raft, sensor = row[0:3]
            dataIdTmp = {'visit':str(visit), 'raft':raft, 'sensor':sensor, 'snap':'0'}
            key = self._dataIdToString(dataIdTmp)
            self.dataIdLookup[key] = dataIdTmp

            if not ssDict.has_key(key):
                ssDict[key] = []
            ss = ssDict[key]
                
            i = 0
            for value in row[3:]:
                method = getattr(s, setMethods[i])
                if not value is None:
                    method(value)
                i += 1


            # calibrate it
            calib = self.getCalibBySensor(dataIdTmp)
            fmag0, fmag0Err = calib[key].getFluxMag0()
            s.setPsfFlux(s.getPsfFlux()/fmag0)
            s.setApFlux(s.getApFlux()/fmag0)
            s.setModelFlux(s.getModelFlux()/fmag0)
            
            ss.append(s)


        # set the STAR flag if we have matched objects
        #if self.matchQueryCache.has_key(dataIdStr) and self.matchQueryCache[dataIdStr]:

        # Orphans defined by RefSrcMatch
        matchListDict = self.getMatchListBySensor(dataIdRegex, useRef='src')
        for k, matchList in matchListDict.items():
            index = {}
            for m in matchList:
                sref, s, dist = m
                sid = s.getId()
                isStar = s.getFlagForDetection() & measAlg.Flags.STAR
                index[sid] = isStar
                
            for s in ssDict[k]:
                if index.has_key(s.getId()):
                    s.setFlagForDetection(s.getFlagForDetection() | index[s.getId()])


        # cache it
        self.queryCache[dataIdStr] = True
        for k, ss in ssDict.items():
            self.sourceSetCache[k] = ssDict[k]
        
        self.printStopLoad()

        return ssDict



    def getRefObjectSetBySensor(self, dataIdRegex):
        """Get a dict of all Catalog Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        sroFields = simRefObj.SimRefObject.fields
        sroFieldStr = ",".join(["sro."+field for field in sroFields])

        oldWay = False

        # if the dataIdEntry is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex)
        if self.refObjectQueryCache.has_key(dataIdStr):
            # get only the ones that match the request
            sroDict = {}
            for key, sro in self.refObjectCache.items():
                if re.search(dataIdStr, key):
                    sroDict[key] = sro
            return sroDict


        # get a list of matching dataIds 
        if oldWay:
            dataIdList = [dataIdRegex]
        else:
            haveAllKeys = True
            sqlDataId = []
            for keyNames in [['visit', 'sce.visit'], ['raft', 'sce.raftName'], ['sensor', 'sce.ccdName']]:
                key, sqlName = keyNames
                if dataIdRegex.has_key(key):
                    sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
                else:
                    haveAllKeys = False
            sqlDataId = " and ".join(sqlDataId)

            sql  = "select sce.visit, sce.raftName, sce.ccdName"
            sql += "  from Science_Ccd_Exposure as sce "
            sql += "  where " + sqlDataId

            dataIdList = []
            results  = self.dbInterface.execute(sql)
            for r in results:
                visit, raft, ccd = map(str, r[0:3])
                dataIdList.append({'visit':visit, 'raft':raft, 'sensor':ccd, 'snap':'0'})


        # Load each of the dataIds
        sroDict = {}
        for dataIdEntry in dataIdList:

            dataIdEntryStr = self._dataIdToString(dataIdEntry)

            haveAllKeys = True
            sqlDataId = []
            for keyNames in [['visit', 'sce.visit'], ['raft', 'sce.raftName'], ['sensor', 'sce.ccdName']]:
                key, sqlName = keyNames
                if dataIdEntry.has_key(key):
                    sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdEntry[key]))
                else:
                    haveAllKeys = False
            sqlDataId = " and ".join(sqlDataId)
            
            # this will have to be updated for the different dataIdNames when non-lsst cameras get used.
            if oldWay:
                sql  = 'select sce.visit, sce.raftName, sce.ccdName, %s' % (sroFieldStr)
                sql += '  from SimRefObject as sro, RefSrcMatch as rsm, Science_Ccd_Exposure as sce'
                sql += '  where (qserv_ptInSphPoly(sro.ra, sro.decl,'
                sql += '          concat_ws(" ", sce.llcRa, sce.llcDecl, sce.lrcRa, sce.lrcDecl, '
                sql += '          sce.urcRa, sce.urcDecl, sce.ulcRa, sce.ulcDecl)) = 1) '
                sql += '        and (rsm.refObjectId = sro.refObjectId) '
                #sql += '        ans (rsm.nSrcMatches = 0) '
                sql += '        and ' + sqlDataId

            else:

                sql  = 'SELECT concat_ws(" ", '
                sql += '   sce.llcRa, sce.llcDecl, '
                sql += '   sce.lrcRa, sce.lrcDecl, '
                sql += '   sce.urcRa, sce.urcDecl, '
                sql += '   sce.ulcRa, sce.ulcDecl) '
                sql += 'FROM Science_Ccd_Exposure as sce '
                sql += 'WHERE %s ' % (sqlDataId)
                #sq += '   (sce.visit = 887252941) AND'
                #sq += '   (sce.raftName = \'2,2\') AND'
                #sq += '   (sce.ccdName = \'1,1\');'
                sql += 'INTO @poly; '

                sql2 = 'SELECT %s ' % (sroFieldStr)
                sql2 += 'FROM '
                sql2 += '    SimRefObject AS sro, '
                sql2 += '    RefSrcMatch AS rsm '
                sql2 += 'WHERE '
                sql2 += '    (qserv_ptInSphPoly(sro.ra, sro.decl, @poly) = 1) AND '
                #sql2 += '    (rsm.nSrcMatches = 0) and '
                sql2 += '    (rsm.refObjectId = sro.refObjectId); '


            # if there are no regexes (ie. actual wildcard expressions),
            #  we can check the cache, otherwise must run the query

            if not re.search("\%", sql) and haveAllKeys:
                dataIdCopy = copy.copy(dataIdEntry)
                dataIdCopy['snap'] = "0"
                key = self._dataIdToString(dataIdCopy)
                if self.refObjectCache.has_key(key):
                    sroDict[key] = self.refObjectCache[key]
                    continue

                        
            self.printStartLoad("Loading RefObjects for: " + dataIdEntryStr + "...")

            # run the query
            if oldWay:
                results  = self.dbInterface.execute(sql)
            else:
                self.dbInterface.execute(sql)
                results = self.dbInterface.execute(sql2)

            # parse results and put them in a sourceSet
            for row in results:
                if oldWay:
                    visit, raft, sensor = row[0:3]
                    sroStuff = row[3:]
                else:
                    visit, raft, sensor = dataIdEntry['visit'], dataIdEntry['raft'], dataIdEntry['sensor']
                    sroStuff = row[:] #row[3:]
                    
                dataIdTmp = {'visit':str(visit), 'raft':raft, 'sensor':sensor, 'snap':'0'}
                key = self._dataIdToString(dataIdTmp)
                self.dataIdLookup[key] = dataIdTmp

                if not sroDict.has_key(key):
                    sroDict[key] = []
                sros = sroDict[key]
                sros.append(simRefObj.SimRefObject(sroStuff))

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

        sql = "select distinct visit from Science_Ccd_Exposure"
        sql += "   where "
        haveAllKeys = True

        whereList = []
        for keyNames in [['visit', 'visit'], ['raft', 'raftName'], ['sensor', 'ccdName']]:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        sql += " and ".join(whereList)

        results = self.dbInterface.execute(sql)
        visits = map(str, zip(*results)[0])
        return visits


    def loadCalexp(self, dataIdRegex):
        """Load the calexp data for data matching dataIdRegex.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        selectList = ["sce."+x for x in qaDataUtils.getSceDbNames()]
        selectStr = ",".join(selectList)

        sql  = 'select '+selectStr
        sql += '  from Science_Ccd_Exposure as sce'
        sql += '  where '

        haveAllKeys = True

        whereList = []
        for keyNames in [['visit', 'sce.visit'], ['raft', 'sce.raftName'], ['sensor', 'sce.ccdName']]:
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
            key = self._dataIdToString(dataIdCopy)
            if self.calexpQueryCache.has_key(key):
                return

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex)
        if self.calexpQueryCache.has_key(dataIdStr) and self.calexpQueryCache[dataIdStr]:
            return

        self.printStartLoad("Loading Calexp for: " + dataIdStr + "...")

        # run the query
        results  = self.dbInterface.execute(sql)

        for row in results:

            rowDict = dict(zip(qaDataUtils.getSceDbNames(), row))

            visit, raft, sensor = rowDict['visit'], rowDict['raftName'], rowDict['ccdName']
            dataIdTmp = {'visit':visit, 'raft':raft, 'sensor':sensor, 'snap':'0'}
            key = self._dataIdToString(dataIdTmp)
            self.dataIdLookup[key] = dataIdTmp
            
            #print rowDict
            if not self.wcsCache.has_key(key):
                crval = afwGeom.PointD(rowDict['crval1'], rowDict['crval2'])
                crpix = afwGeom.PointD(rowDict['crpix1'], rowDict['crpix2'])
                cd11, cd12, cd21, cd22 = rowDict['cd1_1'], rowDict['cd1_2'], rowDict['cd2_1'], rowDict['cd2_2']
                wcs = afwImage.createWcs(crval, crpix, cd11, cd12, cd21, cd22)
                self.wcsCache[key] = wcs

            if not self.detectorCache.has_key(key):
                raftName = "R:"+rowDict['raftName']
                ccdName = raftName + " S:"+rowDict['ccdName']
                #raftId = cameraGeom.Id(rowDict['raft']) #cameraGeom.Id(raftName)
                #ccdId = cameraGeom.Id(rowDict['ccd']) #cameraGeom.Id(ccdName)
                #ccdDetector = cameraGeom.Detector(ccdId)
                #raftDetector = cameraGeom.Detector(raftId)
                #ccdDetector.setParent(raftDetector)
                #self.raftDetectorCache[key] = self.cameraInfo.camera.findDetector(raftId) #raftDetector
                self.detectorCache[key] = self.cameraInfo.detectors[ccdName] #ccdDetector
                self.raftDetectorCache[key] = self.cameraInfo.detectors[raftName]

            if not self.filterCache.has_key(key):
                filter = afwImage.Filter(rowDict['filterName'], True)
                self.filterCache[key] = filter
            
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
        dataIdStr = self._dataIdToString(dataIdRegex
                                         )
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

        clause = ""
        if re.search('^\d+$', regex):
            clause += field + " = %s" % (regex)
        elif re.search('^[\d,]+$', regex):
            clause += field + " = '%s'" % (regex)
        elif re.search('^(\.\*|\?|\%)[\d,]*$', regex) or re.search('^[\d,]*(\.\*|\?|\%)$', regex):
            regexSql = re.sub("(\.\*|\?)", "%", regex)
            clause += field + " like '%s'" % (regexSql)
        else:
            raise Exception("Regex for SQL can only use '.*', '?', or '%' at beginning or end of string ("+
                            field+"="+regex+")")

        return "("+clause+")"



###################################################
# Factory for dbQaData
# - curently only lsstSim is available by database, so this is a trivial factory
###################################################
def makeDbQaData(label, rerun=None, **kwargs):
    """Factory for a DbQaData object.
    
    @param database The name of the database to connect to
    @param rerun The data rerun to use
    """
    return DbQaData(label, rerun, qaCamInfo.LsstSimCameraInfo())


