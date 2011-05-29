import os, re
import eups

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils

####################################################################
#
# CameraInfo base class
#
####################################################################


class CameraInfo(object):
    
    def __init__(self, name, dataInfo, mapperClass=None, camera=None):
        """
        @param name         A name for this camera
        @param dataInfo     List of dataId [name,distinguisher] pairs (distinguisher now depricated)
        @param mapperClass  Mapper class for this camera
        @param camera       LSST camera object for this camera
        """
        
        self.name        = name
        self.dataInfo    = dataInfo
        self.mapperClass = mapperClass
        self.camera      = camera
        
        self.detectors   = {}
        self.nSensor     = 0
        for r in self.camera:
            raft = cameraGeom.cast_Raft(r)
	    raftName = raft.getId().getName().strip()
            self.detectors[raftName] = raft
            for c in raft:
                ccd = cameraGeom.cast_Ccd(c)
		ccdName = ccd.getId().getName().strip()
                self.detectors[ccdName] = ccd
                self.nSensor += 1

    def getDetectorName(self, raft, ccd):
        ccdId = self.detectors[ccd].getId()
        name = re.sub("\s+", "_", ccdId.getName())
        serial = "%04d" % (ccdId.getSerial())
        return name + "--" + serial
                
    def getRoots(self, data, calib, output):
        """Store data directories in a dictionary
        
        @param data    Input directory containing registry.sqlite file
        @param calib   Input calibration directory containing calibRegistry.sqlite file
        @param output  Output directory
        """
        return {'data': data, 'calib': calib, 'output': output}

    def getRegistries(self, baseDir):
        """Given a directory, get the registry files for this data set.
        
        @param baseDir  Directory wherein lives the registry
        """
        roots = self.getRoots(baseDir)
        registry      = os.path.join(roots['data'], "registry.sqlite3")
        calibRegistry = os.path.join(roots['calib'], "calibRegistry.sqlite3")
        return registry, calibRegistry

        
    def verifyRegistries(self, baseDir):
        """Verify that registry.sqlite files exist in the specified directory

        @param baseDir  Directory to check for registries.
        """
        
        registry, calibRegistry = self.getRegistries(baseDir)
        haveReg = os.path.exists(registry)
        haveCalib = os.path.exists(calibRegistry)
        #print self.name, registry, calibRegistry, haveReg, haveCalib
        return haveReg and haveCalib


    def __str__(self):
        return name
    
    # some butler's use a 'rerun', others don't ... base class should default to None
    def getDefaultRerun(self):
        """Get our rerun."""
        return None

    def getSensorCount(self):
        """Get the number of sensors (ie. ccds) for this camera."""
        return self.nSensor
                


####################################################################
####################################################################    




####################################################################
#
# LsstSimCameraInfo class
#
####################################################################
class LsstSimCameraInfo(CameraInfo):
    
    def __init__(self):
        """ """
        try:
            import lsst.obs.lsstSim           as obsLsst
            mapper = obsLsst.LsstSimMapper
        except Exception, e:
            print "Failed to import lsst.obs.lsstSim", e
            mapper = None
        dataInfo       = [['visit',1], ['snap', 0], ['raft',0], ['sensor',0]]

        #simdir        = eups.productDir("obs_lsstSim")
        simdir        = os.environ['OBS_LSSTSIM_DIR']
        cameraGeomPaf = os.path.join(simdir, "description", "Full_STA_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        
        CameraInfo.__init__(self, "lsstSim", dataInfo, mapper, camera)

        self.doLabel = False


    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = baseDir
        if not output is None:
            baseOut = output
        return CameraInfo.getRoots(self, baseDir, baseDir, baseOut)

    def verifyRegistries(self, baseDir):
        """Verify that registry.sqlite files exist in the specified directory

        @param baseDir  Directory to check for registries.
        """
        roots = self.getRoots(baseDir)
        registry = os.path.join(roots['data'], "registry.sqlite3")
        #calibRegistry = os.path.join(roots['data'], "registry.sqlite3")
        return os.path.exists(registry)

    def getMapper(self, baseDir, rerun=None):
        """Get a mapper for data in specified directory

        @param baseDir  Directory where the registry files are to be found.
        @param rerun    The rerun of the data we want.
        """
        roots = self.getRoots(baseDir)
        registry, calibRegistry = self.getRegistries(baseDir)

        return self.mapperClass(root=roots['output'], calibRoot=roots['calib'], registry=registry)

####################################################################
#
# CfhtCameraInfo class
#
####################################################################
class CfhtCameraInfo(CameraInfo):
    def __init__(self):
        try:
            import lsst.obs.cfht              as obsCfht
            mapper = obsCfht.CfhtMapper
        except Exception, e:
            print "Failed to import lsst.obs.cfht", e
            mapper = None
        dataInfo       = [['visit',1], ['ccd', 0]]
        
        #simdir        = eups.productDir("obs_cfht")
        simdir         = os.environ['OBS_CFHT_DIR']
        cameraGeomPaf = os.path.join(simdir, "megacam", "description", "Full_Megacam_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)

        CameraInfo.__init__(self, "cfht", dataInfo, mapper, camera)

        self.doLabel = True
        
    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = baseDir
        if not output is None:
            baseOut = output
        return CameraInfo.getRoots(self, baseDir, os.path.join(baseDir, "calib"), baseDir)


    def getMapper(self, baseDir, rerun=None):
        """Get a mapper for data in specified directory
        
        @param baseDir  Directory where the registry files are to be found.
        @param rerun    The rerun of the data we want
        """
        roots = self.getRoots(baseDir)
        registry, calibRegistry = self.getRegistries(baseDir)
        return self.mapperClass(root=roots['output'], calibRoot=roots['calib'], registry=registry)


####################################################################
#
# HscCameraInfo class
#
####################################################################
class HscCameraInfo(CameraInfo):
    
    def __init__(self):
        """ """
        try:
            import lsst.obs.hscSim            as obsHsc
            mapper = obsHsc.HscSimMapper
        except Exception, e:
            print "Failed to import lsst.obs.hscSim", e
            mapper = None
        dataInfo       = [['visit',1], ['ccd', 0]]

        #simdir        = eups.productDir("obs_subaru")
        simdir         = os.environ['OBS_SUBARU_DIR']
        cameraGeomPaf = os.path.join(simdir, "hscSim", "description", "hscSim_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)

        CameraInfo.__init__(self, "hscSim", dataInfo, mapper, camera)

        self.doLabel = False
        
    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = os.path.join(baseDir, "HSC")
        if not output is None:
            baseOut = output
        return CameraInfo.getRoots(self, os.path.join(baseDir, "HSC"), os.path.join(baseDir, "CALIB"), baseOut)

    def getDefaultRerun(self):
        return "pipeQA"
    

    def getMapper(self, baseDir, rerun=None):
        """Get a mapper for data in specified directory

        @param baseDir  Directory where the registry files are to be found.
        @param rerun    The rerun of the data we want
        """
       
        roots = self.getRoots(baseDir)
        registry, calibRegistry = self.getRegistries(baseDir)
        return self.mapperClass(rerun, root=roots['output'], calibRoot=roots['calib'], registry=registry)


####################################################################
#
# SuprimecamCameraInfo class
#
####################################################################
class SuprimecamCameraInfo(CameraInfo):
    def __init__(self):
        try:
            import lsst.obs.suprimecam        as obsSuprimecam
            mapper = obsSuprimecam.SuprimecamMapper
        except Exception, e:
            print "Failed to import lsst.obs.suprimecam", e
            mapper = None
        dataInfo       = [['visit',1], ['ccd', 0]]

        #simdir        = eups.productDir("obs_subaru")
        simdir         = os.environ['OBS_SUBARU_DIR']
        cameraGeomPaf = os.path.join(simdir, "suprimecam", "description", "Full_Suprimecam_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)

        CameraInfo.__init__(self, "suprimecam", dataInfo, mapper, camera)

        self.doLabel = True
        
    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = os.path.join(baseDir, "SUPA")
        if not output is None:
            baseOut = output
        data = os.path.join(baseDir, "SUPA")
        calib = os.path.join(baseDir, "SUPA", "CALIB")
        return CameraInfo.getRoots(self, data, calib, baseOut)

    def getDefaultRerun(self):
        return "pipeQA"
    

    def getMapper(self, baseDir, rerun=None):
        """Get a mapper for data in specified directory

        @param baseDir  Directory where the registry files are to be found.
        @param rerun    The rerun of the data we want
        """
       
        roots = self.getRoots(baseDir)
        registry, calibRegistry = self.getRegistries(baseDir)
        return self.mapperClass(rerun, root=roots['output'], calibRoot=roots['calib'], registry=registry)
    



def getCameraInfoAvailable():
    """Get a list of available CameraInfo objects."""
    
    available = []

    def tryLoad(cameraInfo):
        haveCam = True
        try:
            ci = cameraInfo()
        except Exception, e:
            haveCam = False
        return haveCam

    all = [
        LsstSimCameraInfo,
        CfhtCameraInfo,
        HscCameraInfo,
        SuprimecamCameraInfo,
        ]

    for camInfo in all:
        if tryLoad(LsstSimCameraInfo):
            available.append(camInfo())

    return available
    
