import os
import eups

import lsst.afw.cameraGeom.utils as cameraGeomUtils

####################################################################
#
# CameraInfo base class
#
####################################################################


class CameraInfo(object):
    
    def __init__(self, name, dataInfo, mapperClass=None, camera=None):
        self.name        = name
        self.dataInfo    = dataInfo
        self.mapperClass = mapperClass
	self.camera      = camera
        

    def getRoots(self, data, calib, output):
        return {'data': data, 'calib': calib, 'output': output}

    def getRegistries(self, baseDir):
        roots = self.getRoots(baseDir)
        registry      = os.path.join(roots['data'], "registry.sqlite3")
        calibRegistry = os.path.join(roots['calib'], "calibRegistry.sqlite3")
        return registry, calibRegistry
        
    def verifyRegistries(self, baseDir):
        registry, calibRegistry = self.getRegistries(baseDir)
        return os.path.exists(registry) and os.path.exists(calibRegistry)

    def getMapper(self, baseDir, rerun=None):
        roots = self.getRoots(baseDir)
        registry, calibRegistry = self.getRegistries(baseDir)
        if rerun is None:
            return self.mapperClass(root=roots['output'], calibRoot=roots['calib'], registry=registry)
        else:
            return self.mapperClass(rerun, root=roots['output'], calibRoot=roots['calib'], registry=registry)

    def __str__(self):
        return name
    
    # some butler's use a 'rerun', others don't ... base class should default to None
    def getDefaultRerun(self):
        return None




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

	simdir        = eups.productDir("obs_lsstSim")
	cameraGeomPaf = os.path.join(simdir, "description", "Full_STA_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)
	
        CameraInfo.__init__(self, "lsstSim", dataInfo, mapper, camera)


    def getRoots(self, baseDir, output=None):
        baseOut = baseDir
        if not output is None:
            baseOut = output
        return CameraInfo.getRoots(self, baseDir, baseDir, baseOut)

    def verifyRegistries(self, baseDir):
        roots = self.getRoots(baseDir)
        registry = os.path.join(roots['data'], "registry.sqlite3")
        #calibRegistry = os.path.join(roots['data'], "registry.sqlite3")
        return os.path.exists(registry)



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
	
	simdir        = eups.productDir("obs_cfht")
	cameraGeomPaf = os.path.join(simdir, "megacam", "description", "Full_Megacam_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)

        CameraInfo.__init__(self, "cfht", dataInfo, mapper, camera)

        
    def getRoots(self, baseDir, output=None):
        baseOut = baseDir
        if not output is None:
            baseOut = output
        return CameraInfo.getRoots(self, baseDir, os.path.join(baseDir, "calib"), baseDir)


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

	simdir        = eups.productDir("obs_subaru")
	cameraGeomPaf = os.path.join(simdir, "hscSim", "description", "hscSim_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)

        CameraInfo.__init__(self, "hscSim", dataInfo, mapper, camera)
        
    def getRoots(self, baseDir, output=None):
        baseOut = os.path.join(baseDir, "HSC")
        if not output is None:
            baseOut = output
        return CameraInfo.getRoots(self, os.path.join(baseDir, "HSC"), os.path.join(baseDir, "CALIB"), baseOut)

    def getDefaultRerun(self):
        return "pipeQA"
    

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

	simdir        = eups.productDir("obs_subaru")
	cameraGeomPaf = os.path.join(simdir, "suprimecam", "description", "Full_Suprimecam_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)

        CameraInfo.__init__(self, "suprimecam", dataInfo, mapper, camera)

    def getRoots(self, baseDir, output=None):
        baseOut = os.path.join(baseDir, "SUPA")
        if not output is None:
            baseOut = output
        data = os.path.join(baseDir, "SUPA")
        calib = os.path.join(baseDir, "SUPA", "CALIB")
        return CameraInfo.getRoots(self, data, calib, baseOut)

    def getDefaultRerun(self):
        return "pipeQA"
    
    
