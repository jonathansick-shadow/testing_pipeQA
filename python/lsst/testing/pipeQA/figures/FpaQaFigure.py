import os, sys, re
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import numpy
import numpy.ma as numpyMa

#import pylab
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse
from matplotlib import cm
from matplotlib import colors

import QaFigureUtils as qaFigUtils
from QaFigure import QaFigure

class FpaQaFigure(QaFigure):

    def __init__(self, camera, data=None, map=None):
        """
        @param camera     The device whose focal plane area we're representing.
        @param data       Data to display
        @param map        Map areas to use.
        """
        
        QaFigure.__init__(self)
        self.camera = camera
        self.centers, self.rectangles, self.raftBoundaries, self.ccdBoundaries = \
                      qaFigUtils.cameraToRectangles(self.camera)

        # To be filled in by child class
        if not data is None:
            if not self.validate():
                raise Exception("Data did not pass validation.")
            self.data = data
        else:
            self.data = {}
            self.reset()

        if not map is None:
            self.map = map
        else:
            self.reset(data=self.map)

        
    def reset(self, data=None):
        """Set all values in data dictionary to None."""
        
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


    def getAreaLabel(self, raft, ccd):
        """Get the area label to use for this raft,ccd."""
        return re.sub("\s+", "_", ccd)


    def setMapInfo(self):
        """Establish map areas for any defined CCDs"""
        
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
                    label = self.getAreaLabel(rlabel, clabel)
                    self.addMapArea(label, [x0, y0, x1, y1], info)
                    


    def plotRaftBoundaries(self, sp, boundaryColors):
        for b in self.raftBoundaries:
            sp.plot(b[0], b[1], '%s-' % (boundaryColors), lw=3)
    def plotCcdBoundaries(self, sp):
        for b in self.ccdBoundaries.values():
            x0, x1 = b[0]
            y0, y1 = b[1]
            x = [x0, x0, x1, x1, x0]
            y = [y0, y1, y1, y0, y0]
            sp.plot(numpy.array(x), numpy.array(y), 'k-', lw=1.0)
    def markMissingCcds(self, sp, missingCcds):
        for b in missingCcds.values():
            x0, x1 = b[0]
            y0, y1 = b[1]
            x = [x0, x1, x0, x1]
            y = [y0, y1, y1, y0]
            sp.plot(numpy.array(x), numpy.array(y), 'k-', lw=0.5)
    def markFailedCcds(self, sp, failedCcds):
        for label, failSign in failedCcds.items():
            b = self.ccdBoundaries[label]
            x0, x1 = b[0]
            y0, y1 = b[1]
            x = 0.5*(x0+x1)
            y = 0.5*(y0+y1)
            text = "F"
            sp.text(x, y, text, horizontalalignment="center", verticalalignment="center",
                    fontsize=8, weight='bold')


    def labelSensors(self, sp):
        for r in self.rectangles.values():
            label = r.get_label()
            bbox  = r.get_bbox()
            xplot = 0.5 * (bbox.x0 + bbox.x1)
            yplot = bbox.y1 - size[1]//2
            sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 8, weight = 'bold')
        
    def adjustTickLabels(self, sp, cb):
        for tic in cb.ax.get_yticklabels():
            tic.set_size("x-small")
        for tic in sp.get_xticklabels():
            tic.set_size("x-small")
            tic.set_rotation(22)
        for tic in sp.get_yticklabels():
            tic.set_size("x-small")
            tic.set_rotation(45)
        
    def getFpaLimits(self):
        x, y = [], []
        for r in self.rectangles.values():
            bbox  = r.get_bbox()
            x += [bbox.x0, bbox.x1]
            y += [bbox.y0, bbox.y1]
        x, y = numpy.array(x), numpy.array(y)
        return x.min(), y.min(), x.max(), y.max()
            



    def makeFigure(self, 
                   borderPix = 100,
                   boundaryColors = 'r', doLabel = False, showUndefined=False,
                   vlimits=None, cmap="jet", title=None,
                   cmapOver=None, cmapUnder=None,
                   failLimits=None, failColor='k',
                   ):
        """Make the figure.

        @param borderPix        width of border in pixels
        @param boundaryColors   matplotlib color specifier for border
        @param doLabel          add ccd labels to the figure
        @param showUndefined    show sensors even if their value is None
        @param vlimits          [low, high] limits for colormap
        @param cmap             string name of an matplotlib.cm colormap
        @param title            title of the figure
        @param cmapOver         Color to use if above cmap high limit.
        @param cmapUnder        Color to use if below cmap low limit.
        @param failLimits       Limits to mark failure.
        @param failColor        Color to use to mark failed sensors.
        """


        if failLimits is None:
            failLimits = vlimits
            
        self.fig.subplots_adjust(left=0.175, right=0.95, bottom=0.15)
        
        sp     = self.fig.gca()

        values = []  # needs to be synchronized with self.rectangles
        patches = []
        allValues = []
        missingCcds = {}
        failedCcds = {}
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
                if value < failLimits[0]:
                    failedCcds[clabel] = -1
                if value > failLimits[1]:
                    failedCcds[clabel] = 1
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
        if not cmapOver is None:
            cmap.set_over(cmapOver, 1.0)
        if not cmapUnder is None:
            cmap.set_under(cmapUnder, 1.0)
        p = PatchCollection(patches, norm=norm, cmap=cmap)
        value_array = numpy.array(values)
        masked_value_array = numpyMa.masked_where(numpy.isnan(value_array), value_array)
        p.set_array(masked_value_array)
        cb = self.fig.colorbar(p)
        sp.add_collection(p)

        self.plotRaftBoundaries(sp, boundaryColors)
        self.plotCcdBoundaries(sp)
        self.markMissingCcds(sp, missingCcds)
        self.markFailedCcds(sp, failedCcds)
        if doLabel:
            self.labelSensors(sp)

        if not title is None:
            sp.set_title(title)
        sp.set_xlabel("Focal Plane X", fontsize = 10, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 10, weight = 'bold')

        self.adjustTickLabels(sp, cb)

        x0, y0, x1, y1 = self.getFpaLimits()
        sp.set_xlim((x0 - borderPix, x1 + borderPix))
        sp.set_ylim((y0 - borderPix, y1 + borderPix))

        self.setMapInfo()
        


class VectorFpaQaFigure(FpaQaFigure):

    def __init__(self, camera, data=None):
        FpaQaFigure.__init__(self, camera)


    def makeFigure(self, 
                   borderPix = 100,
                   boundaryColors = 'r', doLabel = False, showUndefined=False,
                   vlimits=None, cmap="jet", title=None,
                   cmapOver=None, cmapUnder=None,
                   failLimits=None, failColor='k',
                   ):
        """Make the figure.

        @param borderPix        width of border in pixels
        @param boundaryColors   matplotlib color specifier for border
        @param doLabel          add ccd labels to the figure
        @param showUndefined    show sensors even if their value is None
        @param vlimits          [low, high] limits for colormap
        @param cmap             string name of an matplotlib.cm colormap
        @param title            title of the figure
        @param cmapOver         Color to use if above cmap high limit.
        @param cmapUnder        Color to use if below cmap low limit.
        @param failLimits       Limits to mark failure.
        @param failColor        Color to use to mark failed sensors.
        """


        self.fig.subplots_adjust(left=0.175, right=0.95, bottom=0.15)
        sp     = self.fig.gca()
        colorValues = []  # needs to be synchronized with self.rectangles
        patches = []
        allValues = []

        radiansWrtX = {}
        lenInPix = {}
        colorScalar = {}
        haveColors = False
        
        missingCcds = {}
        failedCcds = {}
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
                        colorScalartmp = 0.0
                    else:
                        raise Exception("values for Vector must be float or [radians, lenInPix, [colorFloat]].")
                else:
                    if not values is None:
                        values = float(values)
                    radiansWrtXtmp, lenInPixtmp, colorScalartmp = values, 1500.0, None

                #print clabel, radiansWrtXtmp, lenInPixtmp, colorScalartmp
                lenInPix[clabel]    = lenInPixtmp
                allValues.append(radiansWrtXtmp)
                if (not radiansWrtXtmp is None): # or (showUndefined):
                    if not colorScalartmp is None:
                        colorValues.append(colorScalartmp)
                        patches.append(self.rectangles[clabel])
                        colorScalar[clabel] = colorScalartmp
                        haveColors = True
                        if colorScalartmp < failLimits[0]:
                            failedCcds[clabel] = -1
                        if colorScalartmp > failLimits[1]:
                            failedCcds[clabel] = 1
                    else:
                        colorValues.append(numpy.NaN)
                        patches.append(self.rectangles[clabel])
                        missingCcds[clabel] = self.ccdBoundaries[clabel]
                    radiansWrtX[clabel] = radiansWrtXtmp

                else:
                    colorValues.append(numpy.NaN)
                    patches.append(self.rectangles[clabel])
                    missingCcds[clabel] = self.ccdBoundaries[clabel]
                    

        if not vlimits is None:
            norm = colors.Normalize(vmin=vlimits[0], vmax=vlimits[1], clip=False)
        else:
            norm = colors.Normalize()

        if len(patches) > 0:

            cmap = getattr(cm, cmap)
            cmap.set_bad('k', 0.2)
            if not cmapOver is None:
                cmap.set_over(cmapOver, 1.0)
            if not cmapUnder is None:
                cmap.set_under(cmapUnder, 1.0)

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


        self.plotRaftBoundaries(sp, boundaryColors)
        self.plotCcdBoundaries(sp)
        self.markMissingCcds(sp, missingCcds)
        self.markFailedCcds(sp, failedCcds)
        if doLabel:
            self.labelSensors(sp)

        if not title is None:
            sp.set_title(title)
        sp.set_xlabel("Focal Plane X", fontsize = 10, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 10, weight = 'bold')

        self.adjustTickLabels(sp, cb)

        x0, y0, x1, y1 = self.getFpaLimits()
        sp.set_xlim((x0 - borderPix, x1 + borderPix))
        sp.set_ylim((y0 - borderPix, y1 + borderPix))

        self.setMapInfo()

