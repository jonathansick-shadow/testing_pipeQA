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

class PsfVsApertureQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
	qaAna.QaAnalysis.__init__(self)

    def test(self, data, dataId):
	
	# get data
	self.ssDict        = data.getSourceSetBySensor(dataId)
	self.detector      = data.getDetectorBySensor(dataId)
	self.calib         = data.getCalibBySensor(dataId)
	self.filter        = data.getFilterBySensor(dataId)
	
	self.diff = raftCcdData.RaftCcdVector(self.detector)
	self.mag  = raftCcdData.RaftCcdVector(self.detector)

	filter = None
	for key, ss in self.ssDict.items():
	    raft = self.detector[key].getParent().getId().getName()
	    ccd  = self.detector[key].getId().getName()

	    filter = self.filter[key].getName()
	    
	    qaAnaUtil.isStar(ss)  # sets the 'STAR' flag
	    for s in ss:
		f1 = s.getPsfFlux()
		f2 = s.getApFlux()
		if ((f1 > 0.0 and f2 > 0.0) and
		    not (s.getFlagForDetection() & measAlg.Flags.INTERP_CENTER )):
		    #(s.getFlagForDetection() & measAlg.Flags.STAR)):
		    
		    m1 = -2.5*numpy.log10(f1) #self.calib[key].getMagnitude(f1)
		    m2 = -2.5*numpy.log10(f2) #self.calib[key].getMagnitude(f2)

		    self.diff.append(raft, ccd, m1 - m2)
		    self.mag.append(raft, ccd, m2)
		    
	group = dataId['visit']
	testSet = self.getTestSet(group)
	testSet.addMetadata('dataset', data.getDataName())
	testSet.addMetadata('visit', dataId['visit'])
	testSet.addMetadata('filter', filter)

	nLowest = 100
	
	self.means = raftCcdData.RaftCcdData(self.detector)
	self.medians = raftCcdData.RaftCcdData(self.detector)
	self.stds  = raftCcdData.RaftCcdData(self.detector)
	for raft,  ccd in self.mag.raftCcdKeys():
	    dmag = self.diff.get(raft, ccd)
	    mag = self.mag.get(raft, ccd)
	    w = numpy.where((mag > 10) & (mag < 20))
	    dmag = dmag[w]
	    
	    stat = afwMath.makeStatistics(dmag, afwMath.NPOINT | afwMath.MEANCLIP |
					  afwMath.STDEVCLIP | afwMath.MEDIAN)
	    mean = stat.getValue(afwMath.MEANCLIP)
	    median = stat.getValue(afwMath.MEDIAN)
	    std = stat.getValue(afwMath.STDEVCLIP)
	    n = stat.getValue(afwMath.NPOINT)

	    self.means.set(raft, ccd, mean)
	    label = "mean psf_vs_aper " + re.sub("\s+", "_", ccd)
	    comment = "mean psf-ap (mag lt 20, nstar/clip=%d/%d)" % (len(dmag),n)
	    testSet.addTest( testCode.Test(label, mean, [-0.02, 0.02], comment) )

	    self.medians.set(raft, ccd, median)
	    label = "median psf_vs_aper "+re.sub("\s+", "_", ccd)
	    comment = "median psf-ap (mag lt 20, nstar/clip=%d/%d)" % (len(dmag), n)
	    testSet.addTest( testCode.Test(label, median, [-0.02, 0.02], comment) )

	    self.stds.set(raft, ccd, std)
	    label = "stdev psf_vs_aper " + re.sub("\s+", "_", ccd)
	    comment = "stdev of psf-ap (mag lt 20, nstar/clip=%d/%d)" % (len(dmag), n)
	    testSet.addTest( testCode.Test(label, std, [0.0, 0.02], comment) )
		


    def plot(self, data, dataId, showUndefined=False):

	group = dataId['visit']
	testSet = self.getTestSet(group)

	# fpa figure
	meanFig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	stdFig = qaFig.FpaQaFigure(data.cameraInfo.camera)
	for raft, ccdDict in meanFig.data.items():
	    for ccd, value in ccdDict.items():
		meanFig.data[raft][ccd] = self.means.get(raft, ccd)
		stdFig.data[raft][ccd] = self.stds.get(raft, ccd)
		if not self.means.get(raft, ccd) is None:
		    meanFig.map[raft][ccd] = "mean=%.4f" % (self.means.get(raft, ccd))
		    stdFig.map[raft][ccd] = "std=%.4f" % (self.stds.get(raft, ccd))
		
	meanFig.makeFigure(showUndefined=showUndefined, cmap="Spectral_r", vlimits=[-0.02, 0.02],
			   title="Mean m$_{PSF}$-m$_{Ap}$")
	testSet.addFigure(meanFig, "meanPsfMinusAperture.png", "mean Psf - aper mag",
			  saveMap=True, navMap=True)
	stdFig.makeFigure(showUndefined=showUndefined, cmap="YlOrRd", vlimits=[0.0, 0.03],
			  title="Stdev m$_{PSF}$-m$_{Ap}$")
	testSet.addFigure(stdFig, "stdPsfMinusAperture.png", "stdev Psf - aper mag",
			  saveMap=True, navMap=True)
	

	# dmag vs mag
	figsize = (4.0, 4.0)
	fig0 = qaFig.QaFig(size=figsize)
	fig0.fig.subplots_adjust(left=0.195)
	ax0 = fig0.fig.add_subplot(111)
	
	nKeys = len(self.mag.raftCcdKeys())
	norm = colors.Normalize(vmin=0, vmax=nKeys)
	sm = cm.ScalarMappable(norm, cmap=cm.jet)

	xlim = [14.0, 25.0]
	ylim = [-0.4, 0.4]
	
	i = 0
	xmin, xmax = 1.0e99, -1.0e99
	for raft, ccd in self.mag.raftCcdKeys():
	    mag  = self.mag.get(raft, ccd)
	    diff = self.diff.get(raft, ccd)

	    print "plotting ", ccd
	    
	    # if there's no data, dump a single point at xmax, 0.0
	    if not len(mag) > 0:
		mag = numpy.array([xmax])
		diff = numpy.array([0.0])

	    min, max = mag.min(), mag.max()
	    if min < xmin: xmin = min
	    if max > xmax: xmax = max
	    
	    size = 1.0


	    fig = qaFig.QaFig(size=figsize)
	    fig.fig.subplots_adjust(left=0.195)
	    ax = fig.fig.add_subplot(111)
	    ax.scatter(mag, diff, size, color='k', label=ccd)

	    #box = ax.get_position()
	    #ax.set_position([box.x0, box.y0, 0.75*box.width, box.height])
	    #fp = fm.FontProperties(size="small")
	    #ax.legend(prop=fp, loc=(1.05, 0.0))
	    ax.set_title("m$_{PSF}$ - m$_{Ap}$  versus  m$_{Ap}$")
	    ax.set_xlabel("m$_{Ap}$")
	    ax.set_ylabel("m$_{PSF}$ - m$_{Ap}$")
	    ax.set_xlim(xlim)
	    ax.set_ylim(ylim)

	    label = re.sub("\s+", "_", ccd)
	    testSet.addFigure(fig, "diff_psf-ap_"+label+".png", "psf-ap vs. ap")
	    
	    color = sm.to_rgba(i)
	    ax0.scatter(mag, diff, size, color=color, label=ccd)
	    ax0.set_xlim(xlim)
	    ax0.set_ylim(ylim)
	    
	    dmag = 0.1
	    ddiff = 0.02
	    for j in range(len(mag)):
		area = (mag[j]-dmag, diff[j]-ddiff, mag[j]+dmag, diff[j]+ddiff)
		fig0.addMapArea(label, area, "%.3f_%.3f"% (mag[j], diff[j]))
	    i += 1



	#box = ax0.get_position()
	#ax0.set_position([box.x0, box.y0, 0.75*box.width, box.height])
	#fp = fm.FontProperties(size="small")
	#ax0.legend(prop=fp, loc=(1.05, 0.0))
	ax0.set_title("m$_{PSF}$ - m$_{Ap}$  versus  m$_{Ap}$")
	ax0.set_xlabel("m$_{Ap}$")
	ax0.set_ylabel("m$_{PSF}$ - m$_{Ap}$")
	#ax0.set_xlim(xlim)
	#ax0.set_ylim(ylim)

	    
	testSet.addFigure(fig0, "diff_psf-ap_all.png", "psf-ap vs. ap", saveMap=True)

