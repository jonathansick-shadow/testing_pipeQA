import sys, os, re, copy
import numpy

import lsst.afw.math        as afwMath
import lsst.meas.algorithms as measAlg

# NOTE: please replace this with (s.getFlagForDetection() & measAlg.Flags.STAR)
#       when we eventually start setting it.
def isStarMoment(ss):
    """Quick and dirty isStar() based on comparison to psf ixx/yy/xy"""

    vixx, vixy, viyy = [], [], []
    for s0 in ss:

        s = s0
        if isinstance(s, list):
            sref, s, d = s0
            
        ixx, ixy, iyy    = s.getIxx(), s.getIyy(), s.getIxy()

        vixx.append(ixx)
        vixy.append(ixy)
        viyy.append(iyy)

    sxx = afwMath.makeStatistics(vixx, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    sxy = afwMath.makeStatistics(vixy, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    syy = afwMath.makeStatistics(viyy, afwMath.MEANCLIP | afwMath.STDEVCLIP)

    for s0 in ss:

        s = s0
        if isinstance(s, list):
            sref, s, d = s0
            
        xxOk = (sxx.getValue(afwMath.MEANCLIP) - s.getIxx())/sxx.getValue(afwMath.STDEVCLIP) < 3.0
        xyOk = (sxy.getValue(afwMath.MEANCLIP) - s.getIxy())/sxx.getValue(afwMath.STDEVCLIP) < 3.0
        yyOk = (syy.getValue(afwMath.MEANCLIP) - s.getIyy())/sxx.getValue(afwMath.STDEVCLIP) < 3.0

        if xxOk and xyOk and yyOk:
            s.setFlagForDetection(s.getFlagForDetection() | measAlg.Flags.STAR)



def isStarDeltaMag(ss):

    for s0 in ss:

        s = s0
        if isinstance(s, list):
            sref, s, d = s0
        f_psf, f_mod = s.getPsfFlux(), s.getModelFlux()
        if f_psf > 0 and f_mod > 0:
            m_psf, m_mod = -2.5*numpy.log10(f_psf), -2.5*numpy.log10(f_mod)

            if abs(m_psf - m_mod) < 0.1:
                s.setFlagForDetection(s.getFlagForDetection() | measAlg.Flags.STAR)


def isStar(ss):
    return isStarDeltaMag(ss)



def robustPolyFit(x, y, order, sigma=3.0, niter=3):

    xNew, yNew = copy.copy(x), copy.copy(y)

    for i in range(niter):

        ret = numpy.polyfit(xNew, yNew, order)
        p = ret[0:2]
        residuals = yNew - numpy.polyval(p, xNew)

        if i == 0:
            mean = numpy.median(residuals)
        else:
            mean = numpy.mean(residuals)
            
        std = numpy.std(residuals)
        
        w = numpy.where( (numpy.abs(residuals - mean)/std) < sigma )
        xNew = xNew[w]
        yNew = yNew[w]
        
    return p
