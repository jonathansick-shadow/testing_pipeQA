import sys, os, re
import lsst.meas.algorithms        as measAlg
import lsst.testing.pipeQA.figures as qaFig
import numpy

import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.font_manager as fm
from  matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection

class AstrometricErrorQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
	qaAna.QaAnalysis.__init__(self)

    def free(self):
	del self.x
	del self.y
	del self.dDec
	del self.dRa
	del self.filter
	
	del self.detector
	del self.matchListDict

	del self.medErrArcsec
	del self.medThetaRad
	
    def test(self, data, dataId):
	
	# get data
	self.matchListDict = data.getMatchListBySensor(dataId)
	self.detector      = data.getDetectorBySensor(dataId)
	self.filter        = data.getFilterBySensor(dataId)
	
	#self.clusters = data.getSourceClusters(dataId)

	# compute the mean ra, dec for each source cluster

	self.dRa  = raftCcdData.RaftCcdVector(self.detector)
	self.dDec = raftCcdData.RaftCcdVector(self.detector)
	self.x    = raftCcdData.RaftCcdVector(self.detector)
	self.y    = raftCcdData.RaftCcdVector(self.detector)

	filter = None
	for key, matchList in self.matchListDict.items():
	    raft = self.detector[key].getParent().getId().getName()
	    ccd  = self.detector[key].getId().getName()
	    filter = self.filter[key].getName()

	    for m in matchList:
		sref, s, dist = m
		ra, dec, raRef, decRef = \
		    [x*numpy.pi/180.0 for x in [s.getRa(), s.getDec(), sref.getRa(), sref.getDec()]]
		
		dDec = decRef - dec
		dRa  = (raRef - ra)*abs(numpy.cos(decRef))

		if not (s.getFlagForDetection() & measAlg.Flags.INTERP_CENTER ):
		    self.dRa.append(raft, ccd, dRa)
		    self.dDec.append(raft, ccd, dDec)
		    self.x.append(raft, ccd, s.getXAstrom())
		    self.y.append(raft, ccd, s.getYAstrom())
		    
		    
	testSet = self.getTestSet(data, dataId)
	testSet.addMetadata('dataset', data.getDataName())
	testSet.addMetadata('visit', dataId['visit'])
	testSet.addMetadata('filter', filter)

	self.medErrArcsec = raftCcdData.RaftCcdData(self.detector)
	self.medThetaRad  = raftCcdData.RaftCcdData(self.detector)
	
	self.maxErr = 0.5 # arcsec
	for raft,  ccd in self.dRa.raftCcdKeys():
	    dRa  = self.dRa.get(raft, ccd)
	    dDec = self.dDec.get(raft, ccd)
	    
	    errArcsec = 206265.0*numpy.sqrt(dRa**2 + dDec**2)
	    #errArcsec = 3600.0*numpy.sqrt(dRa**2 + dDec**2)
	    thetaRad  = numpy.arctan2(dDec, dRa)

	    if len(errArcsec) > 0:
		stat  = afwMath.makeStatistics(errArcsec, afwMath.NPOINT | afwMath.MEDIAN)
		medErrArcsec = stat.getValue(afwMath.MEDIAN)
		stat  = afwMath.makeStatistics(thetaRad, afwMath.NPOINT | afwMath.MEDIAN)
		medThetaRad = stat.getValue(afwMath.MEDIAN)
		n = stat.getValue(afwMath.NPOINT)
	    else:
		medErrArcsec = 10.0*self.maxErr
		medThetaRad = 0.0
		n = 0

	    self.medErrArcsec.set(raft, ccd, medErrArcsec)
	    self.medThetaRad.set(raft, ccd, medThetaRad)
	    
	    label = "median astrometry error "+re.sub("\s+", "_", ccd)
	    comment = "median sqrt(dRa^2+dDec^2) (arcsec, nstar=%d)" % (n)
	    testSet.addTest( testCode.Test(label, medErrArcsec, [0.0, self.maxErr], comment) )


    def plot(self, data, dataId, showUndefined=False):

	testSet = self.getTestSet(data, dataId)

	# fpa figure
	astFig = qaFig.VectorFpaQaFigure(data.cameraInfo.camera)
	vLen = 5000 # length in pixels for 1 arcsec error vector
	for raft, ccdDict in astFig.data.items():
	    for ccd, value in ccdDict.items():
		if not self.medErrArcsec.get(raft, ccd) is None:
		    astErrArcsec = self.medErrArcsec.get(raft, ccd)
		    thetaRad = self.medThetaRad.get(raft, ccd)
		    astFig.data[raft][ccd] = [thetaRad, vLen*astErrArcsec, astErrArcsec]
		    astFig.map[raft][ccd] = "\"/theta=%.2f/%.0f" % (astErrArcsec, (180/numpy.pi)*thetaRad)
		
	astFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, self.maxErr],
			  title="Median astrometric error", cmapOver='#ff0000')
	testSet.addFigure(astFig, "medAstError.png", "Median astrometric error", 
			  saveMap=True, navMap=True)

	#
	figsize = (6.5, 3.25)
	conv = colors.ColorConverter()
	black = conv.to_rgb('k')

	
	i = 0
	for raft, ccd in self.dRa.raftCcdKeys():
	    ra = self.dRa.get(raft, ccd)
	    dec = self.dDec.get(raft, ccd)
	    eLen = 206265.0*numpy.sqrt(ra**2 + dec**2)
	    #eLen = 3600.0*numpy.sqrt(ra**2 + dec**2)
	    t = numpy.arctan2(dec, ra)
	    
	    dx = eLen*numpy.cos(t)
	    dy = eLen*numpy.sin(t)
	    x = self.x.get(raft, ccd)
	    y = self.y.get(raft, ccd)

	    # if there's no data, dump a single point at 0,0
	    if not len(x) > 0:
		x = numpy.array([0.0])
		y = numpy.array([0.0])
		dx = numpy.array([0.0])
		dy = numpy.array([0.0])
		
	    # round up to nearest 1024 for limits
	    xmax, ymax = x.max(), y.max()
	    xlim = [0, 1024*int(xmax/1024.0 + 0.5)]
	    ylim = [0, 1024*int(ymax/1024.0 + 0.5)]

	    r = numpy.sqrt(dx**2 + dy**2)
	    rmax = r.max()
	    
	    print "plotting ", ccd
	    
	    fig = qaFig.QaFig(size=figsize)
	    fig.fig.subplots_adjust(left=0.1)

	    ################
	    # ccd view
	    ax = fig.fig.add_subplot(121)
	    box = ax.get_position()
	    ax.set_position([box.x0, box.y0 + 0.1*box.height, box.width, 0.9*box.height])

	    ax.scatter(x, y, 0.5, color='r')
	    q = ax.quiver(x, y, dx, dy, color='k', scale=5.0, angles='xy')
	    ax.quiverkey(q, 0.9, -0.2, 1.0, "1 arcsec", coordinates='axes',
			 fontproperties={'size':"small"})

	    ax.xaxis.set_major_locator(MaxNLocator(8))
	    ax.yaxis.set_major_locator(MaxNLocator(8))
	    #ax.set_title(label)
	    ax.set_xlabel("x [pixels]")
	    ax.set_ylabel("y [pixels]")
	    ax.set_xlim(xlim)
	    ax.set_ylim(ylim)
	    for tic in ax.get_xticklabels() + ax.get_yticklabels():
		tic.set_size("x-small")

	    ################
	    # rose view
	    ax0 = fig.fig.add_subplot(122)
	    box = ax0.get_position()
	    ax0.set_position([box.x0, box.y0 + 0.1*box.height, box.width, 0.9*box.height])
	    ax = ax0.twinx()

	    # this is much faster than calling plot() in a loop, and quiver() scale length buggy
	    z = numpy.zeros(len(dx))
	    xy2 = zip(dx, dy)
	    xy1 = zip(z, z)
	    lines = zip(xy1, xy2)
	    p = LineCollection(lines, colors=black*len(lines))
	    ax.add_collection(p)
	    ax.scatter(dx, dy, s=1.0, color='r')

	    ax.set_xlabel("dRa [arcsec]")
	    ax.set_ylabel("dDec [arcsec]")
	    ax.set_xlim([-rmax, rmax])
	    ax.set_ylim([-rmax, rmax])
	    for tic in ax.get_xticklabels() + ax.get_yticklabels() + ax0.get_xticklabels():
		tic.set_size("x-small")

	    ax0.set_yticklabels([])


	    label = re.sub("\s+", "_", ccd)
	    testSet.addFigure(fig, "astromError_"+label+".png", "Astrometric error "+label)
	    
	    i += 1

