import sys, os, glob, re, stat, copy
import traceback

import numpy

import sqlite

import eups
import lsst.pex.policy                  as pexPolicy
import lsst.pex.logging                 as pexLog
import lsst.daf.persistence             as dafPersist
import lsst.afw.detection               as afwDet

from .DatabaseQuery import LsstSimDbInterface, DatabaseIdentity

import Manifest   as manifest
import CameraInfo as qaCamInfo


import QaDataUtils as qaDataUtils


#######################################################################
#
#
#
#######################################################################
class QaData(object):
    """ """

    #######################################################################
    #
    #######################################################################
    def __init__(self, label, rerun, cameraInfo):
        self.label = label
        self.rerun = rerun
	self.cameraInfo = cameraInfo
        self.dataInfo = self.cameraInfo.dataInfo
        
        self.dataIdNames   = []
        self.dataIdDiscrim = []
        for array in self.dataInfo:
            dataIdName, dataIdDiscrim = array
            self.dataIdNames.append(dataIdName)
            self.dataIdDiscrim.append(dataIdDiscrim)

	# cache the dataId requests
	# they may contain regexes ... we won't know if our cached sourceSets have
	#   all the entries that match unless we redo the query
	#   But, if we've already done the identical query, we know we have everything
	self.queryCache = {}
	self.columnQueryCache = {}

        # cache source sets to avoid reloading the same thing
        self.sourceSetCache = {}
	self.sourceSetColumnCache = {}
	
        # cache calexp to avoid reloading
        self.calexpCache = {}


    def getSourceSetColumnsBySensor(self, dataIdRegex, accessors):

	dataIdStr = self._dataIdToString(dataIdRegex)

	# if they want a specific dataId and we have it ....
	ssTDict = {}
	if self.sourceSetColumnCache.has_key(dataIdStr):
	    ssTDict[dataIdStr] = {}
	    haveIt = True
	    for accessor in accessors:
		if self.sourceSetColumnCache[dataIdStr].has_key(accessor):
		    ssTDict[dataIdStr][accessor] = self.sourceSetColumnCache[dataIdStr][accessor]
		else:
		    haveIt = False
		    break
	    if haveIt:
		return ssTDict

	# if they've made this exact query before ... we must have it!
	if False:
	    ssTDict = {}
	    if self.columnQueryCache.has_key(dataIdStr+"-"+accessor):
		for k, ssDict in self.sourceSetColumnCache():
		    if re.search(dataIdStr, k):
			if not ssTDict.has_key(k):
			    ssTDict[k] = {}
			ssTDict[k][accessor] = ssDict[accessor]
		return ssTDict

	# get it and put the transpose in the cache
	# return a copy of what we cached.
	ssDict = self.getSourceSetBySensor(dataIdRegex)
	ssTDict = {}
	for k, ss in ssDict.items():
	    if not self.sourceSetColumnCache.has_key(k):
		self.sourceSetColumnCache[k] = {}
	    ssTDict[k] = {}
	    for accessor in accessors:
		tmp = numpy.array([])
		for s in ss:
		    method = getattr(s, "get"+accessor)
		    value = method()
		    tmp = numpy.append(tmp, value)
		self.sourceSetColumnCache[k][accessor] = tmp
		ssTDict[k][accessor] = tmp

	#self.transposeQueryCache[dataIdStr+'-'+accessor] = True
	
	return ssTDict

    

    def verifyDataIdKeys(self, dataIdKeys, raiseOnFailure=True):
        missingKeys = []
        for key in dataIdKeys:
            if not key in self.dataIdNames:
                missingKeys.append(key)
        if len(missingKeys) > 0:
            keyStr = ",".join(missingKeys)
            if raiseOnFailure:
                raise Exception("Key(s): "+keyStr+" not valid for camera "+str(self.cameraInfo))
            return False
        return True


    #######################################################################
    #
    #######################################################################
    def _dataTupleToString(self, dataTuple):
        """Represent a dataTuple as a string."""
        
        s = []
        for i in range(len(self.dataIdNames)):
            name = self.dataIdNames[i]
            value = re.sub("[,]", "", str(dataTuple[i]))
            s.append(name + value)
        return "-".join(s)

    #######################################################################
    # utility to convert a data tuple to a dictionary using dataId keys
    #######################################################################
    def _dataTupleToDataId(self, dataTuple):
        """ """
        dataId = {}
        for i in range(len(self.dataIdNames)):
            dataIdName = self.dataIdNames[i]
            dataId[dataIdName] = dataTuple[i]
        return dataId

    #######################################################################
    # utility to convert a dataId dictionary to a tuple
    #######################################################################
    def _dataIdToDataTuple(self, dataId):
	""" """
	# if snap isn't specified, we'll add it.
	dataIdCopy = copy.copy(dataId)
	if not dataIdCopy.has_key('snap'):
	    dataIdCopy['snap'] = '0'
	    
	dataList = []
	for i in range(len(self.dataIdNames)):
	    dataIdName = self.dataIdNames[i]
	    if dataIdCopy.has_key(dataIdName):
		dataList.append(dataIdCopy[dataIdName])
	    else:
		raise Exception("key: "+dataIdName+" not in dataId: "+ str(dataIdCopy))
	return tuple(dataList)

    def _dataIdToString(self, dataId):
	return self._dataTupleToString(self._dataIdToDataTuple(dataId))
    
    
    def getSourceSet(self, dataId):
        pass

    
#######################################################################
#
#
#
#######################################################################
class ButlerQaData(QaData):
    """ """

    #######################################################################
    #
    #######################################################################
    def __init__(self, label, rerun, cameraInfo, dataDir, **kwargs):
        """
        keyword args:
        haveManifest = boolean, verify files in dataDir are present according to manifest
        verifyChecksum = boolean, verify files in dataDir have correct checksum as listed in manifest
        """
        QaData.__init__(self, label, rerun, cameraInfo)
        self.rerun = rerun
        self.dataDir = dataDir


        ###############################################
        # handle keyword args
        ###############################################
        self.kwargs      = kwargs
	self.dataId         = self.kwargs.get('dataId', {})
        self.haveManifest   = self.kwargs.get('haveManifest', False)
        self.verifyChecksum = self.kwargs.get('verifyChecksum', False)

        ###############################################
        # check the manifest, if requested
        # haveManifest = True is a bit slowish
        # verifyChecksum = True is quite slow
        manifest.verifyManifest(self.dataDir, verifyExists=self.haveManifest,
				verifyChecksum=self.verifyChecksum)


        # This (dataId fetching) needs a better design, but will require butler/mapper change, I think.
        #
        # these obscure things refer to the names assigned to levels in the data hierarchy
        # eg. for lsstSim:   dataInfo  = [['visit',1], ['snap', 0], ['raft',0], ['sensor',0]]
        # a level is considered a discriminator if it represents different pictures of the same thing
        # ... so the same object may appear in multiple 'visits', but not on multiple 'sensors'
        # dataInfo is passed in from the derived class as it's specific to each mapper
        
        dataIdRegexDict = {}
        for array in self.dataInfo:
            dataIdName, dataIdDiscrim = array

            # if the user requested eg. visit=1234.*
            # pull that out of kwargs and put it in dataIdRegexDict
            if self.dataId.has_key(dataIdName):
                dataIdRegexDict[dataIdName] = self.dataId[dataIdName]
                

        #######################################
        # get butler
        self.outMapper = self.cameraInfo.getMapper(self.dataDir, rerun=self.rerun)
        self.outButler = dafPersist.ButlerFactory(mapper=self.outMapper).create()

        
        ####################################################
        # make a list of the frames we're asked to care about

        # get all the available raw inputs
        self.availableDataTuples = self.outButler.queryMetadata('raw', self.dataIdNames,
                                                                format=self.dataIdNames)

        # of the data available, get a list of the ones the user actually wants us
        #  to run.  A bit sketchy here ... kwargs contains non-idname info as well.
        self.dataTuples = self._regexMatchDataIds(dataIdRegexDict, self.availableDataTuples)


                

    #######################################################################
    #
    #######################################################################
    def getSourceSetBySensor(self, dataIdRegex):
        """Get sources for requested data as one sourceSet."""
        
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)

        # get the datasets corresponding to the request
        ssDict = {}
        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            dataKey = self._dataTupleToString(dataTuple)
            
            if self.sourceSetCache.has_key(dataKey):
                ssDict[dataKey] = copy.copy(self.sourceSetCache[dataKey])
                continue

            # make sure we actually have the output file
            isWritten = self.outButler.datasetExists('src', dataId)
            if isWritten:
                persistableSourceVector = self.outButler.get('src', dataId)
                sourceSetTmp = persistableSourceVector.getSources()

                if self.outButler.datasetExists('calexp', dataId):
                    postIsrCcd = self.outButler.get('calexp', dataId)
                    calib = postIsrCcd.getCalib()
                    
                    fmag0, fmag0err = calib.getFluxMag0()
                    for s in sourceSetTmp:
                        apFlux  = s.getApFlux()
                        psfFlux = s.getPsfFlux()
                        s.setApFlux(apFlux/fmag0)
                        s.setPsfFlux(psfFlux/fmag0)

                self.sourceSetCache[dataKey] = sourceSetTmp
                ssDict[dataKey] = copy.copy(sourceSetTmp)
            else:
                print str(dataTuple) + " output file missing.  Skipping."
                
        return ssDict

    def getSourceSet(self, dataIdRegex):

	ssDict = self.getSourceSetBySensor(dataIdRegex)
	ssReturn = []
	for key, ss in ssDict.items():
	    ssReturn += ss
	return ssReturn


    #######################################################################
    # utility to go through a list of data Tuples and return
    #  the ones which match regexes for the corresponding data type
    # so user can say eg. raft='0,\d', visit='855.*', etc
    #######################################################################
    def _regexMatchDataIds(self, dataIdRegexDict, availableDataTuples):
        """ """

        # go through the list of what's available, and compare to what we're asked for
        # Put matches in a list of tuples, eg. [(vis1,sna1,raf1,sen1),(vis2,sna2,raf2,sen2)] 
        dataTuples = []
        for dataTuple in availableDataTuples:

            # start true and fail if any dataId keys fail ... eg. 'visit' doesn't match
            match = True
            for i in range(len(self.dataIdNames)):
                dataIdName = self.dataIdNames[i]   # eg. 'visit', 'sensor', etc
                regexForThisId = dataIdRegexDict.get(dataIdName, '.*') # default to '.*' or 'anything'
                dataId = dataTuple[i]

                # if it doesn't match, this frame isn't to be run.
                if not re.search(regexForThisId,  str(dataId)):
                    match = False

            if match:
                dataTuples.append(dataTuple)
                
        return dataTuples
                

    



#######################################################################
#
#
#
#######################################################################
def makeButlerQaData(label, rerun=None, **kwargs):
        
    testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)

    # make sure LsstSim is last in the list (its 'verifyRegistries()' will pass for all cameras)
    cameraInfos = [
        qaCamInfo.CfhtCameraInfo(),
        qaCamInfo.HscCameraInfo(),
        qaCamInfo.SuprimecamCameraInfo(),
        qaCamInfo.LsstSimCameraInfo(),
        ]

    # try each camera ***in-order***
    # note that LsstSim will look like all other as it uses the same registry for data and calib
    # ... must test it last
    cameraToUse = None
    for cameraInfo in cameraInfos:
        # if the mapper couldn't be found, we can't use this camera
        hasMapper = not cameraInfo.mapperClass is None
        validReg = cameraInfo.verifyRegistries(testdataDir)
        if hasMapper and validReg:
            cameraToUse = cameraInfo
            break

    if cameraToUse is None:
        raise Exception("Can't find registries usable with any mappers.")
    else:
        if rerun is None:
            rerun = cameraToUse.getDefaultRerun()
        return ButlerQaData(label, rerun, cameraToUse, testdataDir, **kwargs)




    
#########################################################################
#
#
#
#########################################################################
class DbQaData(QaData):
    #Qa__init__(self, label, rerun, dataInfo):
    #ButlerQaData.__init__(self, label, rerun, cameraInfo, dataDir, **kwargs):

    def __init__(self, database, rerun, cameraInfo):
        QaData.__init__(self, database, rerun, cameraInfo)
        self.dbId        = DatabaseIdentity(self.label)
        self.dbInterface = LsstSimDbInterface(self.dbId)


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

	    visit, raft, sensor = row[0:3]
	    dataIdTmp = {'visit':visit, 'raft':raft, 'sensor':sensor, 'snap':'0'}
	    key = self._dataIdToString(dataIdTmp)
	    if not ssDict.has_key(key):
		ssDict[key] = []
	    ss = ssDict[key]
		
            i = 0
	    for value in row[3:]:
                method = getattr(s, setMethods[i])
                if not value is None:
                    method(value)
                i += 1
                
            ss.append(s)

	# cache it
	self.queryCache[dataIdStr] = True
	for k, ss in ssDict.items():
	    self.sourceSetCache[k] = ssDict[k]
	
	return ssDict



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


###################################################
# Factory for QaData
###################################################
def makeQaData(label, rerun=None, retrievalType=None, **kwargs):
    """ """
    
    if retrievalType is None:
	
	# if there's only one possibility, use that
	
	# see if there's a testbed directory called 'label'
	validButler = False
	testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label, raiseOnFailure=False)
	if (not testbedDir is None) and (not testdataDir is None):
	    validButler = True

	# see if we can connect to a database with name 'label'
	# NOTE: must update if/when non-lsst databases get used
	validDb = True
	try:
	    dbInterface = LsstSimDbInterface(DatabaseIdentity(label))
	except Exception, e:
	    validDb = False

	if validButler and not validDb:
	    retrievalType = 'butler'
	if validDb and not validButler:
	    retrievalType = 'db'
	if validDb and validButler:
	    raise Exception("The label "+label+" is present as both a testbed directory and a database."\
			    "Please specify retrievalType='butler', or retrievalType='db'.")
	if not validDb and not validButler:
	    raise Exception("Unable to find "+label+" as a testbed directory or a database.")

	    
    if re.search("[Bb]utler", retrievalType):
	return makeButlerQaData(label, rerun, **kwargs)

    if re.search("^([Dd][Bb]|[Dd]atabase)$", retrievalType):
	return makeDbQaData(label, rerun, **kwargs)

