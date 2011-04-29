import sys, os, re
import numpy

import lsst.afw.math        as afwMath
import lsst.meas.algorithms as measAlg

# NOTE: please replace this with (s.getFlagForDetection() & measAlg.Flags.STAR)
#       when we eventually start setting it.
def isStar(ss):
    """Quick and dirty isStar() based on comparison to psf ixx/yy/xy"""

    vixx, vixy, viyy = [], [], []
    for s in ss:
	ixx, ixy, iyy    = s.getIxx(), s.getIyy(), s.getIxy()

	vixx.append(ixx)
	vixy.append(ixy)
	viyy.append(iyy)

    sxx = afwMath.makeStatistics(vixx, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    sxy = afwMath.makeStatistics(vixy, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    syy = afwMath.makeStatistics(viyy, afwMath.MEANCLIP | afwMath.STDEVCLIP)

    for s in ss:
	xxOk = (sxx.getValue(afwMath.MEANCLIP) - s.getIxx())/sxx.getValue(afwMath.STDEVCLIP) < 3.0
	xyOk = (sxy.getValue(afwMath.MEANCLIP) - s.getIxy())/sxx.getValue(afwMath.STDEVCLIP) < 3.0
	yyOk = (syy.getValue(afwMath.MEANCLIP) - s.getIyy())/sxx.getValue(afwMath.STDEVCLIP) < 3.0

	if xxOk and xyOk and yyOk:
	    s.setFlagForDetection(s.getFlagForDetection() | measAlg.Flags.STAR)


