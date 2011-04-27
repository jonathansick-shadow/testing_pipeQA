import sys, os, re, copy
import numpy

import lsst.ap.cluster as apCluster

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

	self.matchQueryCache = {}
	self.matchListCache = {}
	
        # cache calexp info, but not the MaskedImage ... it's too big.
        self.calexpCache = {}
	self.wcsCache = {}
	self.detectorCache = {}
	self.raftDetectorCache = {}
	self.filterCache = {}
	self.calibCache = {}

	# store the explicit dataId (ie. no regexes) for each key used in a cache
	self.dataIdLookup = {}


    def getDataName(self):
	return self.label+" rerun="+str(self.rerun)
	

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




    def getWcsBySensor(self, dataIdRegex):
	return self.getCalexpEntryBySensor(self.wcsCache, dataIdRegex)
    def getDetectorBySensor(self, dataIdRegex):
	return self.getCalexpEntryBySensor(self.detectorCache, dataIdRegex)
    def getFilterBySensor(self, dataIdRegex):
	return self.getCalexpEntryBySensor(self.filterCache, dataIdRegex)
    def getCalibBySensor(self, dataIdRegex):
	return self.getCalexpEntryBySensor(self.calibCache, dataIdRegex)


    

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


    def getSourceClusters(self, dataIdRegex,
			  epsilonArcsec = 0.5,
			  minNeighbors = 1,
			  pointsPerLeaf = 100,
			  leafExtentThresholdArcsec = 0.5):
	
	policy = pexPolicy.Policy()
	policy.set("epsilonArcsec", epsilonArcsec)
	policy.set("minNeighbors", minNeighbors)
	policy.set("pointsPerLeaf", pointsPerLeaf)
	policy.set("leafExtentThresholdArcsec", leafExtentThresholdArcsec)

	sources = self.getSourceSet(dataIdRegex)
	sourceClusters = apCluster.cluster(sources, policy)
	
	return sourceClusters

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
    

    def keyStringsToVisits(self, keyList):
	visits = {}
	for key in keyList:
	    visit = str(self.dataIdLookup[key]['visit'])
	    visits[visit] = True
	return visits.keys()
    
	
    def getSourceSet(self, dataId):
        pass

    
