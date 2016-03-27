import os
import re
import copy
import eups
import numpy

import lsst.daf.persistence as dafPersist
import lsst.afw.geom as afwGeom
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

        self.name = name
        self.dataInfo = dataInfo
        self.mapperClass = mapperClass
        self.camera = camera

        self.rafts = {}
        self.sensors = {}
        self.detectors = {}
        self.nSensor = 0

        self.raftCcdKeys = []

        if self.camera is None:
            return

        self.rawName = "raw"

        for r in self.camera:
            raft = cameraGeom.cast_Raft(r)
            raftName = raft.getId().getName().strip()
            self.detectors[raftName] = raft
            self.rafts[raftName] = raft
            for c in raft:
                ccd = cameraGeom.cast_Ccd(c)
                ccdName = ccd.getId().getName().strip()
                self.detectors[ccdName] = ccd
                self.sensors[ccdName] = ccd
                self.nSensor += 1
                self.raftCcdKeys.append([raftName, ccdName])

        self.dataIdTranslationMap = {
            # standard  : camera
            'visit': 'visit',
            'snap': 'snap',
            'raft': 'raft',
            'sensor': 'sensor',
        }
        self.dataIdDbNames = {
            'visit': 'visit',
            'raft': 'raftName',
            'sensor': 'ccdName',
            'snap': 'snap',
        }

    def standardToDbName(self, name):
        return self.dataIdDbNames[self.dataIdTranslationMap[name]]

    def dataIdCameraToStandard(self, dataIdIn):
        """Put this camera dataId in standard visit,raft,sensor format"""

        dataId = copy.copy(dataIdIn)

        transMap = self.dataIdTranslationMap

        for sKey, cKey in transMap.items():

            # it may be a simple string
            if isinstance(cKey, str):
                if dataId.has_key(cKey) and cKey != sKey:
                    dataId[sKey] = dataId[cKey]
                    del dataId[cKey]

            # or it may be a list e.g. visit might be coded ('-' sep.) to map to run-frame
            elif isinstance(cKey, list):
                haveKeys = True
                values = []
                for c in cKey:
                    if not dataId.has_key(c):
                        haveKeys = False
                    else:
                        values.append(dataId[c])

                if haveKeys:
                    dataId[sKey] = "-".join(map(str, values))
                    for c in cKey:
                        del dataId[c]
                else:
                    raise KeyError, "Can't map "+",".join(cKey)+" to "+sKey+". "+str(dataId)

        return dataId

    def dataIdStandardToCamera(self, dataIdIn):
        """ Convert an input dataId (visit,raft,sensor) to use key,value pairs for this specific camera. """

        dataId = copy.copy(dataIdIn)
        # check all keys in the translation map
        transMap = self.dataIdTranslationMap
        for inKey, outKey in transMap.items():

            # if we have this key, and it's has to be remapped ...
            if dataId.has_key(inKey) and inKey != outKey:

                # it may be a simple string
                if isinstance(outKey, str):
                    dataId[outKey] = dataId[inKey]

                # or it may be a list e.g. visit might be coded ('-' sep.) to map to run-frame
                elif isinstance(outKey, list):

                    # value for inKey must be '-' joined or '.*'
                    if dataId[inKey] == '.*':
                        dataId[inKey] = '-'.join(['.*']*len(outKey))
                    if not re.search('-', dataId[inKey]):
                        raise ValueError, "Combined keys must be dash '-' separated."
                    inValues = dataId[inKey].split('-')

                    # inValues must have the same number of values as outKey
                    if len(inValues) != len(outKey):
                        raise IndexError, "Can't map %s.  %d keys != %d values" % (
                            inKey, len(outKey), len(inValues))
                    for i in range(len(outKey)):
                        dataId[outKey[i]] = inValues[i]

                del dataId[inKey]
        return dataId

    def setDataId(self, dataId, key, value):
        if dataId.has_key(key):
            dataId[key] = value
        else:
            cKey = self.dataIdTranslationMap[key]
            if isinstance(cKey, str):
                dataId[cKey] = value
            elif isinstance(cKey, list):
                # make sure value can be split
                if not re.search("-", value):
                    raise ValueError, "Cannot split split "+value+" into "+str(cKey)

                values = value.split("-")
                if len(values) != len(cKey):
                    raise IndexError, "Number of values and keys don't match: " + str(cKey)+" "+str(values)

                for i in range(len(cKey)):
                    dataId[cKey[i]] = values[i]

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
        registry = os.path.join(roots['data'], "registry.sqlite3")
        calibRegistry = os.path.join(roots['calib'], "calibRegistry.sqlite3")
        return registry, calibRegistry

    def verifyRegistries(self, baseDir):
        """Verify that registry.sqlite files exist in the specified directory

        @param baseDir  Directory to check for registries.
        """

        registry, calibRegistry = self.getRegistries(baseDir)
        haveReg = os.path.exists(registry)
        haveCalib = os.path.exists(calibRegistry)
        # print self.name, registry, calibRegistry, haveReg, haveCalib
        return haveReg and haveCalib

    def __str__(self):
        return self.name

    # some butler's use a 'rerun', others don't ... base class should default to None
    def getDefaultRerun(self):
        """Get our rerun."""
        return None

    def getSensorCount(self):
        """Get the number of sensors (ie. ccds) for this camera."""
        return self.nSensor

    def raftKeyToName(self, raft):
        for r in self.rafts.keys():
            if re.search(raft, r):
                return r
        return None

    def ccdKeyToName(self, ccd):
        for c in self.sensors.keys():
            if re.search("^R:\d,\d S:\d,\d$", c):
                if re.search("^R:\d,\d S:"+ccd, c):
                    return c
            else:
                if re.search(ccd, c):
                    return c
        return None

    def getRaftAndSensorNames(self, dataId):
        raftName = "R:"+str(dataId[self.dataIdTranslationMap['raft']])
        ccdName = raftName + " S:"+str(dataId[self.dataIdTranslationMap['sensor']])
        return raftName, ccdName

    def getSensorName(self, raft, ccd):
        ccdName = raftName + " S:"+str(rowDict[self.cameraInfo.standardToDbName('sensor')])
        pass

    def getBbox(self, raftName, ccdName):

        rxc, ryc = 0.0, 0.0
        if self.rafts.has_key(raftName):
            raft = self.rafts[raftName]
            # NOTE: all ccd coords are w.r.t. the *center* of the raft, not its LLC
            rxc = raft.getCenterPixel().getX()
            ryc = raft.getCenterPixel().getY()

        ccd = self.sensors[ccdName]
        cxc = ccd.getCenterPixel().getX()
        cyc = ccd.getCenterPixel().getY()
        orient = ccd.getOrientation()
        nQuart = ccd.getOrientation().getNQuarter()
        yaw = orient.getYaw()

        cbbox = ccd.getAllPixels(True)
        cwidth = cbbox.getMaxX() - cbbox.getMinX()
        cheight = cbbox.getMaxY() - cbbox.getMinY()
        if abs(yaw.asRadians() - numpy.pi/2.0) < 1.0e-3:  # nQuart == 1 or nQuart == 3:
            ctmp = cwidth
            cwidth = cheight
            cheight = ctmp

        cx0 = rxc + cxc - cwidth/2
        cy0 = ryc + cyc - cheight/2
        cx1 = cx0 + cwidth
        cy1 = cy0 + cheight

        return [cx0, cy0, cx1, cy1]


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
            import lsst.obs.lsstSim as obsLsst
            mapper = obsLsst.LsstSimMapper
        except Exception, e:
            print "Failed to import lsst.obs.lsstSim", e
            mapper = None
        dataInfo = [['visit', 1], ['snap', 0], ['raft', 0], ['sensor', 0]]

        #simdir        = eups.productDir("obs_lsstSim")
        if os.environ.has_key('OBS_LSSTSIM_DIR'):
            simdir = os.environ['OBS_LSSTSIM_DIR']
            cameraGeomPaf = os.path.join(simdir, "description", "Full_STA_geom.paf")
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera = None

        CameraInfo.__init__(self, "lsstSim", dataInfo, mapper, camera)

        self.doLabel = False

        self.dataIdTranslationMap = {
            # input  : return
            'visit': 'visit',
            'snap': 'snap',
            'raft': 'raft',
            'sensor': 'sensor',
        }

    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = baseDir
        if output is not None:
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
            import lsst.obs.cfht as obsCfht
            mapper = obsCfht.CfhtMapper
        except Exception, e:
            print "Failed to import lsst.obs.cfht", e
            mapper = None
        dataInfo = [['visit', 1], ['ccd', 0]]

        #simdir        = eups.productDir("obs_cfht")
        if os.environ.has_key('OBS_CFHT_DIR'):
            simdir = os.environ['OBS_CFHT_DIR']
            cameraGeomPaf = os.path.join(simdir, "megacam", "description", "Full_Megacam_geom.paf")
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera = None

        CameraInfo.__init__(self, "cfht", dataInfo, mapper, camera)
        self.doLabel = True

    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = baseDir
        if output is not None:
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
            import lsst.obs.hscSim as obsHsc
            mapper = obsHsc.HscSimMapper
        except Exception, e:
            print "Failed to import lsst.obs.hscSim", e
            print "  Did you setup obs_subaru?"
            mapper = None
        dataInfo = [['visit', 1], ['ccd', 0]]

        #simdir        = eups.productDir("obs_subaru")
        if os.environ.has_key('OBS_SUBARU_DIR'):
            simdir = os.environ['OBS_SUBARU_DIR']
            cameraGeomPaf = os.path.join(simdir, "hscSim", "description", "hscSim_geom.paf")
            if not os.path.exists(cameraGeomPaf):
                cameraGeomPaf = os.path.join(simdir, "hscSim", "hscSim_geom.paf")
                if not os.path.exists(cameraGeomPaf):
                    raise Exception("Unable to find cameraGeom Policy file: %s" % (cameraGeomPaf))
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera = None

        CameraInfo.__init__(self, "hscSim", dataInfo, mapper, camera)

        self.doLabel = False

    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = os.path.join(baseDir, "HSC")
        if output is not None:
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
        return self.mapperClass(rerun=rerun, root=roots['output'], calibRoot=roots['calib'],
                                registry=registry)


####################################################################
#
# SuprimecamCameraInfo class
#
####################################################################
class SuprimecamCameraInfo(CameraInfo):

    def __init__(self, mit=False):
        try:
            import lsst.obs.suprimecam as obsSuprimecam
            mapper = obsSuprimecam.SuprimecamMapper
        except Exception, e:
            print "Failed to import lsst.obs.suprimecam", e
            mapper = None
        dataInfo = [['visit', 1], ['ccd', 0]]

        #simdir        = eups.productDir("obs_subaru")
        if os.environ.has_key('OBS_SUBARU_DIR'):
            simdir = os.environ['OBS_SUBARU_DIR']
            cameraGeomPaf = os.path.join(simdir, "suprimecam", "description", "Full_Suprimecam_geom.paf")
            if not os.path.exists(cameraGeomPaf):
                cameraGeomPaf = os.path.join(simdir, "suprimecam", "Full_Suprimecam_geom.paf")
                if not os.path.exists(cameraGeomPaf):
                    raise Exception("Unable to find cameraGeom Policy file: %s" % (cameraGeomPaf))
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera = None

        CameraInfo.__init__(self, "suprimecam", dataInfo, mapper, camera)

        self.doLabel = True
        self.mit = mit

    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = os.path.join(baseDir, "SUPA")
        if output is not None:
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
        return self.mapperClass(rerun=rerun, mit=self.mit, root=roots['output'],
                                calibRoot=roots['calib'], registry=registry)


####################################################################
#
# SdssCameraInfo class
#
####################################################################
class SdssCameraInfo(CameraInfo):

    def __init__(self):
        try:
            import lsst.obs.sdss as obsSdss
            mapper = obsSdss.SdssMapper
        except Exception, e:
            print "Failed to import lsst.obs.sdss", e
            mapper = None

        messingWithNames = True
        if messingWithNames:
            dataInfo = [['run', 1], ['filter', 0], ['field', 1], ['camcol', 0]]
        else:
            dataInfo = [['run', 1], ['band', 0], ['frame', 1], ['camcol', 0]]

        #simdir        = eups.productDir("obs_subaru")
        if os.environ.has_key('OBS_SDSS_DIR'):
            simdir = os.environ['OBS_SDSS_DIR']
            #cameraGeomPaf = os.path.join(simdir, "sdss", "description", "Full_Suprimecam_geom.paf")
            # if not os.path.exists(cameraGeomPaf):
            #    cameraGeomPaf = os.path.join(simdir, "sdss", "Full_Sdss_geom.paf")
            #    if not os.path.exists(cameraGeomPaf):
            #        raise Exception("Unable to find cameraGeom Policy file: %s" % (cameraGeomPaf))
            #cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            #camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)
            camera = obsSdss.makeCamera.makeCamera(name='SDSS')
        else:
            camera = None

        CameraInfo.__init__(self, "sdss", dataInfo, mapper, camera)
        self.rawName = "fpC"

        self.doLabel = True

        if messingWithNames:
            self.dataIdTranslationMap = {
                'visit': ['run', 'field'],
                'raft': 'camcol',
                'sensor': 'filter',
            }

            self.dataIdDbNames = {
                'run': 'run',
                'field': 'field',
                'camcol': 'camcol',
                'filter': 'filterName',
            }
        else:
            self.dataIdTranslationMap = {
                'visit': ['run', 'frame'],
                'raft': 'camcol',
                'sensor': 'band',
            }

            self.dataIdDbNames = {
                'run': 'run',
                'frame': 'frame',
                'camcol': 'camcol',
                'band': 'filterName',
            }

    def getRaftAndSensorNames(self, dataId):
        raftName = str(dataId[self.dataIdTranslationMap['raft']])
        ccdName = str(dataId[self.dataIdTranslationMap['sensor']]) + raftName
        return raftName, ccdName

    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = baseDir
        if output is not None:
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

    def getDefaultRerun(self):
        return "pipeQA"

    def getMapper(self, baseDir, rerun=None):
        """Get a mapper for data in specified directory

        @param baseDir  Directory where the registry files are to be found.
        @param rerun    The rerun of the data we want
        """

        roots = self.getRoots(baseDir)
        registry, calibRegistry = self.getRegistries(baseDir)
        return self.mapperClass(root=roots['output'],
                                calibRoot=roots['calib'], registry=registry)


####################################################################
#
# SdssCameraInfo class
#
####################################################################
class CoaddCameraInfo(CameraInfo):

    def __init__(self, **kwargs):
        try:
            # Use this until we can determine the mapper from the archive
            import lsst.obs.lsstSim as obsLsst
            mapper = obsLsst.LsstSimMapper(root = kwargs["skymapRep"])
            butlerFactory = dafPersist.ButlerFactory(mapper = mapper)
            butler = butlerFactory.create()
            self.skyMap = butler.get(datasetType=kwargs["coaddTable"] + "Coadd_skyMap")
        except Exception, e:
            print "Failed to import lsst.obs.lsstSim", e
            mapper = None

        dataInfo = [['tract', 1], ['patch', 1], ['filterName', 1]]
        camera = []  # Empty until we get an idea of the skymap footprint

        CameraInfo.__init__(self, "coadd", dataInfo, mapper, camera)

        self.doLabel = False

        self.dataIdTranslationMap = {
            'visit': ['tract', 'filterName'],
            'raft': None,
            'sensor': 'patch',
        }

        self.dataIdDbNames = {
            'patch': 'patch',
            'tract': 'tract',
            'filterName': 'filterName',
        }

    def skyMapToCamera(self, dataIds):
        """Make a minimal camera based on the skymap; ONLY DESIGNED TO WORK WITH 1 TRACT (sorry, future developer)"""
        tracts = set(x["tract"] for x in dataIds)
        assert(len(tracts) == 1)
        self.tract = tracts.pop()

        cols = set([x["patch"][0] for x in dataIds])  # x
        rows = set([x["patch"][1] for x in dataIds])  # y
        col0 = min(cols)
        row0 = min(rows)
        ncols = max(cols) - col0 + 1
        nrows = max(rows) - row0 + 1

        origBBox = afwGeom.Box2I()

        raftId = cameraGeom.Id(0, str(self.tract))
        raft = cameraGeom.Raft(raftId, ncols, nrows)
        for nId, dataId in enumerate(dataIds):
            patch = self.skyMap[self.tract].getPatchInfo(dataId["patch"])
            detId = cameraGeom.Id(nId, "%d-%d,%d" % (self.tract, patch.getIndex()[0], patch.getIndex()[1]))
            bbox = patch.getInnerBBox()  # outer bbox overlaps, inner doesn't
            origBBox.include(bbox)       # capture the full extent of the skymap
            # need the bbox to be w.r.t. the raft coord system; i.e. 0
            bbox.shift(-afwGeom.Extent2I(bbox.getBegin()))
            ccd = cameraGeomUtils.makeDefaultCcd(bbox, detId=detId)

            col = patch.getIndex()[0] - col0
            row = patch.getIndex()[1] - row0

            xc = bbox.getBeginX() + 0.5 * bbox.getWidth()
            yc = bbox.getBeginY() + 0.5 * bbox.getHeight()

            raft.addDetector(afwGeom.Point2I(col, row),
                             cameraGeom.FpPoint(afwGeom.Point2D(xc, yc)),
                             cameraGeom.Orientation(),
                             ccd)

        raftBbox = raft.getAllPixels()
        xc = origBBox.getBeginX() + 0.5 * origBBox.getWidth()
        yc = origBBox.getBeginY() + 0.5 * origBBox.getHeight()
        cameraId = cameraGeom.Id(0, "Skymap")
        camera = cameraGeom.Camera(cameraId, 1, 1)
        camera.addDetector(afwGeom.Point2I(0, 0),
                           cameraGeom.FpPoint(afwGeom.Point2D(xc, yc)),
                           cameraGeom.Orientation(), raft)

        self.camera = camera
        # Now, redo the assignments in __init__
        for r in self.camera:
            raft = cameraGeom.cast_Raft(r)
            raftName = raft.getId().getName().strip()
            self.detectors[raftName] = raft
            self.rafts[raftName] = raft
            for c in raft:
                ccd = cameraGeom.cast_Ccd(c)
                ccdName = ccd.getId().getName().strip()
                self.detectors[ccdName] = ccd
                self.sensors[ccdName] = ccd
                self.nSensor += 1
                self.raftCcdKeys.append([raftName, ccdName])

    def setFilterless(self):
        self.dataIdTranslationMap['visit'] = ['tract', 'patch']
        del self.dataIdDbNames['filterName']
        self.dataInfo = self.dataInfo[0:2]

    def getRaftAndSensorNames(self, dataId):
        ccdName = str(dataId['tract']) + '-' + str(dataId['patch'])
        return None, ccdName

    def getRoots(self, baseDir, output=None):
        """Get data directories in a dictionary

        @param baseDir The base directory where the registries can be found.
        """
        baseOut = baseDir
        if output is not None:
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

    def getDefaultRerun(self):
        return "pipeQA"

    def getMapper(self, baseDir, rerun=None):
        """Get a mapper for data in specified directory

        @param baseDir  Directory where the registry files are to be found.
        @param rerun    The rerun of the data we want
        """

        roots = self.getRoots(baseDir)
        registry, calibRegistry = self.getRegistries(baseDir)
        return self.mapperClass(root=roots['output'],
                                calibRoot=roots['calib'], registry=registry)


def getCameraInfoAvailable():
    """Get a list of available CameraInfo objects."""

    available = []

    def tryLoad(cameraInfo):
        haveCam = True
        ci = cameraInfo()
        if ci.camera is None:
            haveCam = False
        return haveCam

    all = [
        SdssCameraInfo,
        CoaddCameraInfo,
        LsstSimCameraInfo,
        # CfhtCameraInfo,
        HscCameraInfo,
        SuprimecamCameraInfo,
    ]

    for camInfo in all:
        if tryLoad(camInfo):
            available.append(camInfo())

    return available

