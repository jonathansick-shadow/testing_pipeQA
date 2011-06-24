import numpy

useSource = False
useRefSource = False
useSource = True
useRefSource = True

if useRefSource:
    #from lsst.afw.detection import Source as RefSource
    import pipeQaLib as pipeQaLib
    RefSource = pipeQaLib.Source
    RefSourceSet = pipeQaLib.SourceSet
else:
    class RefSource(object):

        def __init__(self):
            self.ints = numpy.zeros(2, dtype=numpy.int)
            self.radec = numpy.zeros(2)
            self.flux  = numpy.zeros(4, dtype=numpy.float32)
            #self._id               = 0            
            #self._flagForDetection = 0x0
            #self._ra               = 0.0
            #self._dec              = 0.0
            #self._psfFlux          = 0.0
            #self._apFlux           = 0.0
            #self._modelFlux        = 0.0

        def setPhotometry(self, val):
            pass
        def setAstrometry(self, val):
            pass
        def setShape(self, val):
            pass

        def getId(self):      return self.ints[0]
        def setId(self, val): self.ints[0] = val
        def getFlagForDetection(self):  return self.ints[1]
        def setFlagForDetection(self, val): self.ints[1] = val

        def getRa(self):        return self.radec[0]
        def setRa(self, val):   self.radec[0] = val
        def getDec(self):       return self.radec[1]
        def setDec(self, val):  self.radec[1] = val

        def getPsfFlux(self):      return self.flux[0]
        def setPsfFlux(self, val): self.flux[0] = val
        def getApFlux(self):       return self.flux[1]
        def setApFlux(self, val):  self.flux[1] = val
        def getInstFlux(self):     return self.flux[2]
        def setInstFlux(self, val): self.flux[2] = val
        def getModelFlux(self):    return self.flux[3]
        def setModelFlux(self, val): self.flux[3] = val



if useSource:
    import pipeQaLib as pipeQaLib
    Source = pipeQaLib.Source
    SourceSet = pipeQaLib.SourceSet
    #from lsst.afw.detection import Source
    
else:
    class Source(object):

        def __init__(self):
            self.ints              = numpy.zeros(2, dtype=int)
            self.radec             = numpy.zeros(2)
            self.astrom            = numpy.zeros(2, dtype=numpy.float32)
            self.flux              = numpy.zeros(4, dtype=numpy.float32)
            self.fluxErr           = numpy.zeros(4, dtype=numpy.float32)
            self.izz               = numpy.zeros(3, dtype=numpy.float32)
            
            #self._id               = 0  
            #self._ra               = 0.0
            #self._dec              = 0.0
            #self._xAstrom          = 0.0
            #self._yAstrom          = 0.0
            #self._psfFlux          = 0.0
            #self._psfFluxErr       = 0.0
            #self._apFlux           = 0.0
            #self._apFluxErr        = 0.0
            #self._instFlux         = 0.0

            #self._instFluxErr      = 0.0
            #self._modelFlux        = 0.0
            #self._modelFluxErr     = 0.0
            #self._ixx              = 0.0
            #self._ixxErr           = 0.0
            #self._iyy              = 0.0
            #self._iyyErr           = 0.0
            #self._ixy              = 0.0
            #self._ixyErr           = 0.0
            #self._psfIxx           = 0.0

            #self._psfIxxErr        = 0.0
            #self._psfIyy           = 0.0
            #self._psfIyyErr        = 0.0
            #self._psfIxy           = 0.0
            #self._psfIxyErr        = 0.0
            #self._resolution       = 0.0
            #self._e1               = 0.0
            #self._e1Err            = 0.0
            #self._e2               = 0.0
            #self._e2Err            = 0.0

            #self._shear1           = 0.0
            #self._shear1Err        = 0.0
            #self._shear2           = 0.0
            #self._shear2Err        = 0.0
            #self._sigma            = 0.0
            #self._sigmaErr         = 0.0
            ##self._shapeStatus    = 0.0  
            #self._flagForDetection = 0x0


        def setPhotometry(self, val): pass
        def setAstrometry(self, val): pass
        def setShape(self, val):      pass
        
        def getId(self):            return self.ints[0]
        def setId(self, val):       self.ints[0] = val
        def getFlagForDetection(self):      return self.ints[1] 
        def setFlagForDetection(self, val): self.ints[1] = val


        def getRa(self):            return self.radec[0]
        def setRa(self, val):       self.radec[0] = val
        def getDec(self):           return self.radec[1]
        def setDec(self, val):      self.radec[1] = val

        def getXAstrom(self):       return self.astrom[0]
        def setXAstrom(self, val):  self.astrom[0] = val
        def getYAstrom(self):       return self.astrom[1]
        def setYAstrom(self, val):  self.astrom[1] = val

        def getPsfFlux(self):       return self.flux[0]
        def setPsfFlux(self, val):  self.flux[0] = val
        def getPsfFluxErr(self):      return self.fluxErr[0]
        def setPsfFluxErr(self, val): self.fluxErr[0] = val
        def getApFlux(self):        return self.flux[1]
        def setApFlux(self, val):   self.flux[1] = val
        def getApFluxErr(self):       return self.fluxErr[1]
        def setApFluxErr(self, val):  self.fluxErr[1] = val
        def getInstFlux(self):      return self.flux[2]
        def setInstFlux(self, val): self.flux[2] = val
        #def getInstFluxErr(self):     return self.fluxErr[2]
        #def setInstFluxErr(self, val): self.fluxErr[2] = val
        def getModelFlux(self):      return self.flux[3]
        def setModelFlux(self, val): self.flux[3] = val
        #def getModelFluxErr(self):   return self.flux[4]
        #def setModelFluxErr(self, val):  self.flux[4] = val

        def getIxx(self):             return self.izz[0]
        def setIxx(self, val):        self.izz[0] = val
        #def getIxxErr(self):          return self.izz[0]Err
        #def setIxxErr(self, val):     self.izz[0]Err = val
        def getIyy(self):             return self.izz[1]
        def setIyy(self, val):        self.izz[1] = val
        #def getIyyErr(self):          return self.izz[1]Err
        #def setIyyErr(self, val):     self.izz[1]Err = val
        def getIxy(self):              return self.izz[2]
        def setIxy(self, val):         self.izz[2] = val
        #def getIxyErr(self):           return self.izz[2]Err
        #def setIxyErr(self, val):      self.izz[2]Err = val


        def getPsfIxx(self): return 0.0
        def setPsfIxx(self, val): pass
        def getPsfIxxErr(self): return 0.0
        def setPsfIxxErr(self, val): pass
        def getPsfIyy(self): return 0.0
        def setPsfIyy(self, val): pass
        def getPsfIyyErr(self): return 0.0
        def setPsfIyyErr(self, val): pass
        def getPsfIxy(self): return 0.0
        def setPsfIxy(self, val): pass
        def getPsfIxyErr(self): return 0.0
        def setPsfIxyErr(self, val): pass

        def getResolution(self): return 0.0
        def setResolution(self, val): pass

        def getE1(self): return 0.0   
        def setE1(self, val): pass
        def getE1Err(self): return 0.0
        def setE1Err(self, val): pass
        def getE2(self): return 0.0   
        def setE2(self, val): pass
        def getE2Err(self): return 0.0
        def setE2Err(self, val): pass

        def getShear1(self): return 0.0
        def setShear1(self, val): pass
        def getShear1Err(self): return 0.0
        def setShear1Err(self, val): pass
        def getShear2(self): return 0.0
        def setShear2(self, val): pass
        def getShear2Err(self): return 0.0
        def setShear2Err(self, val): pass

        def getSigma(self): return 0.0
        def setSigma(self, val): pass
        def getSigmaErr(self): return 0.0
        def setSigmaErr(self, val): pass


#         def getPsfIxx(self):             return self._psfIxx
#         def setPsfIxx(self, val):             self._psfIxx = val
#         def getPsfIxxErr(self):             return self._psfIxxErr
#         def setPsfIxxErr(self, val):             self._psfIxxErr = val
#         def getPsfIyy(self):             return self._psfIyy
#         def setPsfIyy(self, val):             self._psfIyy = val
#         def getPsfIyyErr(self):             return self._psfIyyErr
#         def setPsfIyyErr(self, val):             self._psfIyyErr = val
#         def getPsfIxy(self):             return self._psfIxy
#         def setPsfIxy(self, val):             self._psfIxy = val
#         def getPsfIxyErr(self):             return self._psfIxyErr
#         def setPsfIxyErr(self, val):             self._psfIxyErr = val

#         def getResolution(self):             return self._resolution
#         def setResolution(self, val):             self._resolution = val

#         def getE1(self):                return self._e1
#         def setE1(self, val):             self._e1 = val
#         def getE1Err(self):             return self._e1Err
#         def setE1Err(self, val):             self._e1Err = val
#         def getE2(self):                return self._e2
#         def setE2(self, val):             self._e2 = val
#         def getE2Err(self):             return self._e2Err
#         def setE2Err(self, val):             self._e2Err = val

#         def getShear1(self):             return self._shear1
#         def setShear1(self, val):             self._shear1 = val
#         def getShear1Err(self):             return self._shear1Err
#         def setShear1Err(self, val):             self._shear1Err = val
#         def getShear2(self):             return self._shear2
#         def setShear2(self, val):             self._shear2 = val
#         def getShear2Err(self):             return self._shear2Err
#         def setShear2Err(self, val):             self._shear2Err = val

#         def getSigma(self):             return self._sigma
#         def setSigma(self, val):             self._sigma = val
#         def getSigmaErr(self):             return self._sigmaErr
#         def setSigmaErr(self, val):             self._sigmaErr = val



