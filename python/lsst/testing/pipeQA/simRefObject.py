import numpy

class SimRefObject(object):

    fields = [
        "refObjectId",    
        "isStar",         
        #"varClass",       
        #"ra",             
        #"decl",           
        #"gLat",           
        #"gLon",           
        #"sedName",        
        "uMag",           
        "gMag",           
        "rMag",           
        "iMag",           
        "zMag",           
        "yMag",           
        #"muRa",           
        #"muDecl",         
        #"parallax",       
        #"vRad",           
        #"redshift",       
        #"semiMajorBulge", 
        #"semiMinorBulge", 
        #"semiMajorDisk",  
        #"semiMinorDisk",  
        #"uCov",           
        #"gCov",           
        #"rCov",           
        #"iCov",           
        #"zCov",           
        #"yCov",           
        ]

    def __init__(self, *args):

        i = 0
        for field in SimRefObject.fields:
            if len(args) == 1:
                values = args[0]
                value = values[i]
            else:
                value = numpy.NaN
            setattr(self, field, value)
            i += 1

    def setMag(self, magNew, filter):
        mag = getattr(self, filter+"Mag")
        mag = magNew

    def setFlux(self, fluxNew, filter):
        mag = getattr(self, filter+"Mag")
        if fluxNew > 0 and not numpy.isnan(fluxNew):
            mag = -2.5*numpy.log10(fluxNew)
        else:
            mag = numpy.NaN
        
    def getMag(self, filter):
        mag = getattr(self, filter+"Mag")
        return mag

    def getFlux(self, filter):
        mag = getattr(self, filter+"Mag")
        return 10.0**(-0.4*mag)

    #def getCov(self, filter):
    #    cov = getattr(self, filter+"Cov")
    #    return cov
    
