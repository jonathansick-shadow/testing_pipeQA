import sys, os, re, copy

import lsst.afw.detection               as afwDet
import lsst.afw.image                   as afwImage
import lsst.afw.geom                    as afwGeom
import lsst.afw.cameraGeom              as cameraGeom

import CameraInfo                       as qaCamInfo

from DatabaseQuery import LsstSimDbInterface, DatabaseIdentity
from QaData        import QaData

import QaDataUtils as qaDataUtils

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
        QaData.__init__(self, database, rerun, cameraInfo)
        self.dbId        = DatabaseIdentity(self.label)
        self.dbInterface = LsstSimDbInterface(self.dbId)


    def getMatchListBySensor(self, dataIdRegex):
	
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
	    if self.matchListCache.has_key(key):
		return {key : self.matchListCache[key]}

	# if the dataIdRegex is identical to an earlier query, we must already have all the data
	dataIdStr = self._dataIdToString(dataIdRegex)
	if self.matchQueryCache.has_key(dataIdStr):
	    matchListDict = {}
	    # get only the ones that match the request
	    for key, matchList in self.matchListCache.items():
		if re.search(dataIdStr, key):
		    matchListDict[key] = matchList
	    return matchListDict

	sql += idWhere
	result = self.dbInterface.execute(sql)
	filterId, filterName = result[0]

	
        # this will have to be updated for the different dataIdNames when non-lsst cameras get used.
        sql  = 'select sce.visit, sce.raftName, sce.ccdName, sro.%sMag, sro.ra, sro.decl, '%(filterName)
	sql += selectStr
        sql += '  from Source as s, Science_Ccd_Exposure as sce,'
        sql += '    RefObjMatch as rom, SimRefObject as sro'
        sql += '  where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
        sql += '    and (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
        sql += '    and (sro.isStar = 1) and (sro.isVar = 0)'
        sql += '    and (s.filterId = %d) and ((s.flagForDetection & 0xa01) = 0)' % (filterId)
        sql += '    and s.objectID is not NULL'        
	sql += '    and '+idWhere
	

	# run the query
        results  = self.dbInterface.execute(sql)

        # parse results and put them in a sourceSet
        matchListDict = {}
        for row in results:
            s = afwDet.Source()
	    qaDataUtils.setSourceBlobsNone(s)
	    sref = afwDet.Source()
	    qaDataUtils.setSourceBlobsNone(sref)
	    
	    visit, raft, sensor, mag, ra, dec = row[0:6]
	    dataIdTmp = {'visit':str(visit), 'raft':raft, 'sensor':sensor, 'snap':'0'}
	    key = self._dataIdToString(dataIdTmp)
	    self.dataIdLookup[key] = dataIdTmp

	    if not matchListDict.has_key(key):
		matchListDict[key] = []
	    matchList = matchListDict[key]

	    sref.setRa(ra)
	    sref.setDec(dec)
	    flux = 10**(-mag/2.5)
	    sref.setPsfFlux(flux)
	    sref.setApFlux(flux)
	    sref.setModelFlux(flux)
	    
            i = 0
	    for value in row[6:]:
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

	    dist = 0.0
            matchList.append([sref, s, dist])

	# cache it
	self.matchQueryCache[dataIdStr] = True
	for k, matchList in matchListDict.items():
	    self.matchListCache[k] = matchListDict[k]
	
	return matchListDict



    def getSourceSetBySensor(self, dataIdRegex):

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


	# cache it
	self.queryCache[dataIdStr] = True
	for k, ss in ssDict.items():
	    self.sourceSetCache[k] = ssDict[k]
	
	return ssDict



    def getVisits(self, dataIdRegex):
	""" Return explicit visits matching for a dataIdRegex."""

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
	    if self.calexpCache.has_key(key):
		return

	# if the dataIdRegex is identical to an earlier query, we must already have all the data
	#dataIdStr = self._dataIdToString(dataIdRegex)
	#if self.queryCache.has_key(dataIdStr):
	#    return

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
		raftId = cameraGeom.Id(raftName)
		ccdId = cameraGeom.Id(ccdName)
		ccdDetector = cameraGeom.Detector(ccdId)
		raftDetector = cameraGeom.Detector(raftId)
		ccdDetector.setParent(raftDetector)
		self.raftDetectorCache[key] = raftDetector
		self.detectorCache[key] = ccdDetector

	    if not self.filterCache.has_key(key):
		filter = afwImage.Filter(rowDict['filterName'], True)
		self.filterCache[key] = filter
	    
	    if not self.calibCache.has_key(key):
		calib = afwImage.Calib()
		calib.setFluxMag0(rowDict['fluxMag0'], rowDict['fluxMag0Sigma'])
		self.calibCache[key] = calib

	    self.calexpCache[key] = True



    def getCalexpEntryBySensor(self, cache, dataIdRegex):

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
	""" """

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
    return DbQaData(label, rerun, qaCamInfo.LsstSimCameraInfo())


