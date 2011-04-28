import os, sys, re
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import numpy
import numpy.ma as numpyMa

#import pylab
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse
from matplotlib import cm
from matplotlib import colors

import QaFigureUtils as qaFigUtils


class QaFig(object):

    def __init__(self, size=(4.0, 4.0), dpi=100): # (512, 512), DPI=100):
        self.fig         = figure.Figure(figsize=size)
	self.fig.set_dpi(dpi)
	self.canvas      = FigCanvas(self.fig)
	self.map         = {}
        #self.fig.set_size_inches(size[0] / DPI, size[1] / DPI)
	self.mapAreas    = []
	
    def reset(self):
        self.fig.clf()
	

    def validate(self):
	pass
    
    def makeFigure(self):
        # override
        pass

    def getFigure(self):
	return self.fig

    def savefig(self, path, **kwargs):
	self.fig.savefig(path, dpi=self.fig.get_dpi(), **kwargs)

    def savemap(self, path):
	mapList = self.getMapInfo()
	if len(mapList) > 0:
	    fp = open(path, 'w')

	    # don't include overplotted map areas
	    n = 100
	    xpmax, ypmax = self.fig.transFigure.transform((1.0, 1.0))
	    haveLookup = numpy.zeros([n, n])
	    for array in mapList:

		label, x0, y0, x1, y1, info = array
		ix = int(n*0.5*(x0 + x1)/xpmax)
		iy = int(n*0.5*(y0 + y1)/ypmax)
		ixOk = ix >= 0 and ix < n
		iyOk = iy >= 0 and iy < n
		if ixOk and iyOk and haveLookup[ix,iy] == 0:
		    fp.write("%s %d %d %d %d %s\n" % (label, x0, y0, x1, y1, info))
		    haveLookup[ix,iy] = 1
		    
	    fp.close()


    def addMapArea(self, label, area, info, axes=None):
	if axes is None:
	    axes = self.fig.gca()

	x0, y0, x1, y1 = area
	tr = self.fig.transFigure.transform((1.0, 1.0))
	xpmax, ypmax = tr
	xy1 = axes.transData.transform((x0, y0))
	xy2 = axes.transData.transform((x1, y1))
	left, bottom = xy1
	right, top = xy2
	self.mapAreas.append([label, left, ypmax-top, right, ypmax-bottom, info])

    
    def getMapInfo(self):
	return self.mapAreas



class FpaQaFigure(QaFig):

    def __init__(self, camera, data=None):
        QaFig.__init__(self)
	self.camera = camera
        self.centers, self.rectangles, self.raftBoundaries, self.ccdBoundaries = \
		      qaFigUtils.cameraToRectangles(self.camera)

        # To be filled in by child class
	if not data is None:
	    if not self.validate():
		raise Exception("Data did not pass validation.")

	self.data = {}
	self.reset()

	self.reset(data=self.map)
	
    def reset(self, data=None):
	if data is None:
	    data = self.data
	    
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            data[rlabel] = {}
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                data[rlabel][clabel] = None
		
    def validate(self):
        # Since we establish the structure of data in __init__, it
        # should always be valid.  Unless someone mucks with self.data
        # directly...
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            if not self.data.has_key(rlabel):
                return False
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                if not self.data[rlabel].has_key(clabel):
                    return False
        return True



    def getMapInfo(self):

	mapList = []
	axes = self.fig.gca()

	tr = self.fig.transFigure.transform((1.0, 1.0))
	xpmax, ypmax = tr
	
	for r in self.camera:
	    raft = cameraGeom.cast_Raft(r)
	    rlabel = raft.getId().getName()
	    for c in raft:
		ccd = cameraGeom.cast_Ccd(c)
		clabel = ccd.getId().getName()
		info = self.map[rlabel][clabel]

		if ((not info is None) and
		    self.ccdBoundaries.has_key(clabel) and
		    (not self.ccdBoundaries[clabel] is None)):
		    bound = self.ccdBoundaries[clabel]
		    x0, x1 = bound[0]
		    y0, y1 = bound[1]
		    xy1 = axes.transData.transform((x0, y0))
		    xy2 = axes.transData.transform((x1, y1))
		    left, bottom = xy1
		    right, top = xy2
		    label = re.sub("\s+", "_", clabel)
		    mapList.append([label, left, ypmax-top, right, ypmax-bottom, info])

	return mapList


    def makeFigure(self, 
                   borderPix = 100,
                   boundaryColors = 'r', doLabel = False, showUndefined=False,
		   vlimits=None, cmap="jet", title=None):

	self.fig.subplots_adjust(left=0.175, right=0.95)
	
        sp     = self.fig.gca()

        values = []  # needs to be synchronized with self.rectangles
	patches = []
	allValues = []
	missingCcds = {}
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
		value = self.data[rlabel][clabel]
		allValues.append(value)
		#if (not value is None) or (showUndefined):
		if value is None:
		    value = numpy.NaN
		    missingCcds[clabel] = self.ccdBoundaries[clabel]
		values.append(value)
		patches.append(self.rectangles[clabel])


	if not vlimits is None:
	    norm = colors.Normalize(vmin=vlimits[0], vmax=vlimits[1], clip=False)
	else:
	    norm = colors.Normalize()
	    
	if len(patches) == 0:
	    patches = self.rectangles.values()
	    values = allValues

	cmap = getattr(cm, cmap)
	cmap.set_bad('k', 0.2)
        p = PatchCollection(patches, norm=norm, cmap=cmap)
	value_array = numpy.array(values)
	masked_value_array = numpyMa.masked_where(numpy.isnan(value_array), value_array)
        p.set_array(masked_value_array)
        cb = self.fig.colorbar(p)
        sp.add_collection(p)

        for b in self.raftBoundaries:
            sp.plot(b[0], b[1], '%s-' % (boundaryColors), lw=3)
	for b in self.ccdBoundaries.values():
	    x0, x1 = b[0]
	    y0, y1 = b[1]
	    x = [x0, x0, x1, x1, x0]
	    y = [y0, y1, y1, y0, y0]
	    sp.plot(numpy.array(x), numpy.array(y), 'k-', lw=1.0)
	for b in missingCcds.values():
	    x0, x1 = b[0]
	    y0, y1 = b[1]
	    x = [x0, x1, x0, x1]
	    y = [y0, y1, y1, y0]
	    sp.plot(numpy.array(x), numpy.array(y), 'k-', lw=0.5)


        if doLabel:
            for r in self.rectangles.values():
                label = r.get_label()
                bbox  = r.get_bbox()
                xplot = 0.5 * (bbox.x0 + bbox.x1)
                yplot = bbox.y1 - size[1]//2
                sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 6, weight = 'bold')

	if not title is None:
	    sp.set_title(title)
        sp.set_xlabel("Focal Plane X", fontsize = 10, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 10, weight = 'bold')

	for tic in cb.ax.get_yticklabels():
	    tic.set_size("x-small")
	for tic in sp.get_xticklabels():
	    tic.set_size("x-small")
	for tic in sp.get_yticklabels():
	    tic.set_size("x-small")
	    tic.set_rotation(45)
	    
        xmin = +1e10
        ymin = +1e10
        xmax = -1e10
        ymax = -1e10
        for r in self.rectangles.values():
            bbox  = r.get_bbox()
            
            if (bbox.x0 < xmin):
                xmin = bbox.x0
            if (bbox.x1 > xmax):
                xmax = bbox.x1
            if (bbox.y0 < ymin):
                ymin = bbox.y0
            if (bbox.y1 > ymax):
                ymax = bbox.y1
        sp.set_xlim((xmin - borderPix, xmax + borderPix))
        sp.set_ylim((ymin - borderPix, ymax + borderPix))
	#self.fig.subplotparams.update(left=0.15)

class VectorFpaQaFigure(FpaQaFigure):

    def __init__(self, camera, data=None):
        FpaQaFigure.__init__(self, camera)


    def makeFigure(self, 
                   borderPix = 100,
                   raftBoundColors = 'r', doLabel = False, showUndefined=False,
		   vlimits=None, cmap="jet", title=None):

	self.fig.subplots_adjust(left=0.175, right=0.95)
        sp     = self.fig.gca()
        colorValues = []  # needs to be synchronized with self.rectangles
	patches = []
	allValues = []

	radiansWrtX = {}
	lenInPix = {}
	colorScalar = {}
	haveColors = False
	
	for r in self.camera:
	    raft   = cameraGeom.cast_Raft(r)
	    rlabel = raft.getId().getName()
	    for c in raft:
		ccd    = cameraGeom.cast_Ccd(c)
		clabel = ccd.getId().getName()
		values = self.data[rlabel][clabel]

		if isinstance(values, list):
		    if len(values) == 3:
			radiansWrtXtmp, lenInPixtmp, colorScalartmp = values
		    elif len(values) == 2:
			radiansWrtXtmp, lenInPixtmp = values
			colorScalartmp = None
		    else:
			raise Exception("values for Vector must be float or [radians, lenInPix, [colorFloat]].")
		else:
		    if not values is None:
			values = float(values)
		    radiansWrtXtmp, lenInPixtmp, colorScalartmp = values, 1500.0, None
		
		allValues.append(radiansWrtXtmp)
		if (not radiansWrtXtmp is None): # or (showUndefined):
		    if not colorScalartmp is None:
			colorValues.append(colorScalartmp)
			patches.append(self.rectangles[clabel])
			colorScalar[clabel] = colorScalartmp
			haveColors = True
		    radiansWrtX[clabel] = radiansWrtXtmp
		    lenInPix[clabel]    = lenInPixtmp
		else:
		    colorValues.append(numpy.NaN)
		    patches.append(self.rectangles[clabel])


	if not vlimits is None:
	    norm = colors.Normalize(vmin=vlimits[0], vmax=vlimits[1], clip=False)
	else:
	    norm = colors.Normalize()

	if len(patches) > 0:
	    #patches = self.rectangles.values()
	    #colorValues = allValues

	    cmap = getattr(cm, cmap)
	    cmap.set_bad('k', 0.2)
	    #cmap.set_over('r', 0.8)
	    #cmap.set_under('b', 0.8)
	    p = PatchCollection(patches, norm=norm, cmap=cmap)
	    value_array = numpy.array(colorValues)
	    masked_value_array = numpyMa.masked_where(numpy.isnan(value_array), value_array)
	    p.set_array(masked_value_array)
	    if haveColors:
		cb = self.fig.colorbar(p)
	    sp.add_collection(p)


	for label, angle in radiansWrtX.items():
	    xc, yc = self.centers[label]
	    arrowLen = lenInPix[label]
	    x = xc - 0.5*arrowLen*numpy.cos(angle)
	    y = yc - 0.5*arrowLen*numpy.sin(angle)
	    dx = arrowLen*numpy.cos(angle)
	    dy = arrowLen*numpy.sin(angle)
	    sp.arrow(x, y, dx, dy) #, ec="k", lw=3)

        for b in self.raftBoundaries:
            sp.plot(b[0], b[1], '%s-' % (raftBoundColors), lw=3)
        for b in self.ccdBoundaries.values():
	    x0, x1 = b[0]
	    y0, y1 = b[1]
	    x = numpy.array([x0, x0, x1, x1, x0])
	    y = numpy.array([y0, y1, y1, y0, y0])
            sp.plot(x, y, 'k-', lw=1)

        if doLabel:
            for r in self.rectangles.values():
                label = r.get_label()
                bbox  = r.get_bbox()
                xplot = 0.5 * (bbox.x0 + bbox.x1)
                yplot = bbox.y1 - size[1]//2
                sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 8, weight = 'bold')

	if not title is None:
	    sp.set_title(title)
        sp.set_xlabel("Focal Plane X", fontsize = 10, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 10, weight = 'bold')

	if len(patches) > 0 and haveColors:
	    for tic in cb.ax.get_yticklabels():
		tic.set_size("x-small")
	for tic in sp.get_xticklabels():
	    tic.set_size("x-small")
	for tic in sp.get_yticklabels():
	    tic.set_size("x-small")
	    tic.set_rotation(45)

        xmin = +1e10
        ymin = +1e10
        xmax = -1e10
        ymax = -1e10
        for r in self.rectangles.values():
            bbox  = r.get_bbox()
            
            if (bbox.x0 < xmin):
                xmin = bbox.x0
            if (bbox.x1 > xmax):
                xmax = bbox.x1
            if (bbox.y0 < ymin):
                ymin = bbox.y0
            if (bbox.y1 > ymax):
                ymax = bbox.y1
        sp.set_xlim((xmin - borderPix, xmax + borderPix))
        sp.set_ylim((ymin - borderPix, ymax + borderPix))
    
