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

class AstrometricErrorQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
	qaAna.QaAnalysis.__init__(self)

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
		s, sref = m
		#print "%.10f %.10f %.10f %.10f" % (s.getRa(), sref.getRa(), s.getDec(), sref.getDec())
		dDec = sref.getDec() - s.getDec()
		dRa  = (sref.getRa() - s.getRa())*abs(numpy.cos(sref.getDec()))
		
		if not (s.getFlagForDetection() & measAlg.Flags.INTERP_CENTER ):
		    self.dRa.append(raft, ccd, dRa)
		    self.dDec.append(raft, ccd, dDec)
		    self.x.append(raft, ccd, s.getXAstrom())
		    self.y.append(raft, ccd, s.getYAstrom())
		    
		    
	group = dataId['visit']
	testSet = self.getTestSet(group)
	testSet.addMetadata('dataset', data.getDataName())
	testSet.addMetadata('visit', dataId['visit'])
	testSet.addMetadata('filter', filter)

	self.medErrArcsec = raftCcdData.RaftCcdData(self.detector)
	self.medThetaRad  = raftCcdData.RaftCcdData(self.detector)
	
	lim = 1.0 # arcsec
	for raft,  ccd in self.dRa.raftCcdKeys():
	    dRa  = self.dRa.get(raft, ccd)
	    dDec = self.dDec.get(raft, ccd)
	    
	    errArcsec = 206265.0*numpy.sqrt(dRa**2 + dDec**2)
	    thetaRad  = numpy.arctan2(dDec, dRa)
	    
	    stat  = afwMath.makeStatistics(errArcsec, afwMath.NPOINT | afwMath.MEDIAN)
	    medErrArcsec = stat.getValue(afwMath.MEDIAN)
	    stat  = afwMath.makeStatistics(thetaRad, afwMath.NPOINT | afwMath.MEDIAN)
	    medThetaRad = stat.getValue(afwMath.MEDIAN)
	    n = stat.getValue(afwMath.NPOINT)

	    self.medErrArcsec.set(raft, ccd, medErrArcsec)
	    self.medThetaRad.set(raft, ccd, medThetaRad)
	    
	    label = "median astrometry error "+re.sub("\s+", "_", ccd)
	    comment = "median sqrt(dRa^2+dDec^2) (arcsec, nstar=%d)" % (n)
	    testSet.addTest( testCode.Test(label, medErrArcsec, [0.0, lim], comment) )


    def plot(self, data, dataId, showUndefined=False):

	group = dataId['visit']
	testSet = self.getTestSet(group)

	# fpa figure
	astFig = qaFig.VectorFpaQaFigure(data.cameraInfo.camera)
	vLen = 100
	for raft, ccdDict in astFig.data.items():
	    for ccd, value in ccdDict.items():
		if not self.medErrArcsec.get(raft, ccd) is None:
		    astErrArcsec = self.medErrArcsec.get(raft, ccd)
		    thetaRad = self.medThetaRad.get(raft, ccd)
		    astFig.data[raft][ccd] = [thetaRad, vLen*astErrArcsec, astErrArcsec]
		    astFig.map[raft][ccd] = "\"/theta=%.2f/%.0f" % (astErrArcsec, (180/numpy.pi)*thetaRad)
		
	astFig.makeFigure(showUndefined=showUndefined, cmap="YlOrRd", vlimits=[0.0, 1.0],
			  title="Median astrometric error")
	testSet.addFigure(astFig, "medAstError.png", "Median astrometric error", 
			  saveMap=True, navMap=True)

	#
	figsize = (6.5, 3.25)
	
	i = 0
	for raft, ccd in self.dRa.raftCcdKeys():
	    ra = self.dRa.get(raft, ccd)
	    dec = self.dDec.get(raft, ccd)
	    eLen = 206265.0*numpy.sqrt(ra**2 + dec**2)
	    t = numpy.arctan2(dec, ra)
	    
	    dx = eLen*numpy.cos(t)
	    dy = eLen*numpy.sin(t)
	    x = self.x.get(raft, ccd)
	    y = self.y.get(raft, ccd)

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
	    q = ax.quiver(x, y, dx, dy, color='k', scale=100.0, angles='xy')
	    ax.quiverkey(q, 0.9, -0.2, 10.0, "10 arcsec", coordinates='axes',
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
	    
	    z = numpy.zeros(len(dx))
	    for i in range(len(dx)):
		ax.plot([0.0, dx[i]], [0.0, dy[i]], '-r')
	    ax.scatter(dx, dy, s=1.0, color='k')

	    #ax.set_title("Astrometric displacements")
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

