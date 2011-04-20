import os, sys
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import numpy

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

    def __init__(self):
        self.fig         = figure.Figure()
	self.canvas      = FigCanvas(self.fig)

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
	self.fig.savefig(path, **kwargs)



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
	else:
	    self.data = {}
	    self.reset()
	
    def reset(self):
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            self.data[rlabel] = {}
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                self.data[rlabel][clabel] = None
		
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


    def makeFigure(self, 
                   DPI = 100., size = (1024, 1024), borderPix = 100,
                   boundaryColors = 'r', doLabel = False, showUndefined=False,
		   vlimits=None, cmap="jet"):

        self.fig.set_size_inches(size[0] / DPI, size[1] / DPI)
        
        sp     = self.fig.gca()
        values = []  # needs to be synchronized with self.rectangles
	patches = []
	allValues = []
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
		value = self.data[rlabel][clabel]
		allValues.append(value)
		if (not value is None) or (showUndefined):
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
	cmap.set_over('r', 0.8)
	cmap.set_under('b', 0.8)
        p = PatchCollection(patches, norm=norm, cmap=cmap)
        p.set_array(numpy.array(values))
        cb = self.fig.colorbar(p)
        sp.add_collection(p)

        for b in self.raftBoundaries:
            sp.plot(b[0], b[1], '%s-' % (boundaryColors), lw=3)
        for b in self.ccdBoundaries:
            sp.plot(b[0], b[1], 'k-', lw=1)

        if doLabel:
            for r in self.rectangles.values():
                label = r.get_label()
                bbox  = r.get_bbox()
                xplot = 0.5 * (bbox.x0 + bbox.x1)
                yplot = bbox.y1 - size[1]//2
                sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 8, weight = 'bold')

        sp.set_xlabel("Focal Plane X", fontsize = 20, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 20, weight = 'bold')

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
    


class VectorFpaQaFigure(FpaQaFigure):

    def __init__(self, camera, data=None):
        FpaQaFigure.__init__(self, camera)


    def makeFigure(self, 
                   DPI = 100., size = (1024, 1024), borderPix = 100,
                   raftBoundColors = 'r', doLabel = False, showUndefined=False,
		   vlimits=None, cmap="jet"):

        self.fig.set_size_inches(size[0] / DPI, size[1] / DPI)
        
        sp     = self.fig.gca()
        colorValues = []  # needs to be synchronized with self.rectangles
	patches = []
	allValues = []

	radiansWrtX = {}
	lenInPix = {}
	colorScalar = {}
	    
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
		    radiansWrtXtmp, lenInPixtmp, colorScalartmp = float(values), 1500.0, None
		
		allValues.append(radiansWrtXtmp)
		if (not radiansWrtXtmp is None) or (showUndefined):
		    if not colorScalartmp is None:
			colorValues.append(colorScalartmp)
			patches.append(self.rectangles[clabel])
			colorScalar[clabel] = colorScalartmp
		    radiansWrtX[clabel] = radiansWrtXtmp
		    lenInPix[clabel]    = lenInPixtmp
		    

	if not vlimits is None:
	    norm = colors.Normalize(vmin=vlimits[0], vmax=vlimits[1], clip=False)
	else:
	    norm = colors.Normalize()

	if len(patches) > 0:
	    #patches = self.rectangles.values()
	    #colorValues = allValues

	    cmap = getattr(cm, cmap)
	    cmap.set_over('r', 0.8)
	    cmap.set_under('b', 0.8)
	    p = PatchCollection(patches, norm=norm, cmap=cmap)
	    p.set_array(numpy.array(colorValues))
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
        for b in self.ccdBoundaries:
            sp.plot(b[0], b[1], 'k-', lw=1)

        if doLabel:
            for r in self.rectangles.values():
                label = r.get_label()
                bbox  = r.get_bbox()
                xplot = 0.5 * (bbox.x0 + bbox.x1)
                yplot = bbox.y1 - size[1]//2
                sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 8, weight = 'bold')

        sp.set_xlabel("Focal Plane X", fontsize = 20, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 20, weight = 'bold')

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
    
