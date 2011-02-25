from .DatabaseQuery import LsstSimDbInterface

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils

import numpy as num
import pylab
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

class HtmlFormatter:
    def __init__(self, butler, cameraGeom):
        self.figureButler = butler
        self.cameraGeom   = cameraGeom

    def generateHtml():
        pass

#
###
#

class QaFigure():
    def __init__(self):
        self.butler = None
        self.fig = pylab.figure()

    def makeFigure(self):
        # override
        pass

    def saveFigure(self, outfile, clear = True):
        self.fig.savefig(outfile)
        if clear:
            self.fig.clf()

class SingleChipFigure(QaFigure):
    def __init__(self):
        QaFigure.__init__(self)

class FpaFigure(QaFigure):
    def __init__(self, cameraGeomPaf):
        QaFigure.__init__(self)

        self.cameraGeomPolicy            = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        self.camera                      = cameraGeomUtils.makeCamera(self.cameraGeomPolicy)
        self.rectangles, self.boundaries = self.cameraToRectangles(self.camera)

        # To be filled in by child class
        self.values = {}
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            self.values[rlabel] = {}
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                self.values[rlabel][clabel] = 0.0


    def fillValues(self):
        # override
        pass

    def makeFigure(self, title,
                   DPI = 100., size = (2000, 2000), borderPix = 100,
                   boundaryColors = 'r', doLabel = False):
        self.fig.set_size_inches(size[0] / DPI, size[1] / DPI)
        
        sp     = self.fig.gca()
        values = []  # needs to be synchronized with self.rectangles
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                values.append(self.values[rlabel][clabel][0])

        p = PatchCollection(self.rectangles)
        p.set_array(num.array(values))
        cb = pylab.colorbar(p)
        sp.add_collection(p)

        for b in self.boundaries:
            pylab.plot(b[0], b[1], '%s-' % (boundaryColors), lw=3)

        if doLabel:
            for r in self.rectangles:
                label = r.get_label()
                bbox  = r.get_bbox()
                xplot = 0.5 * (bbox.x0 + bbox.x1)
                yplot = bbox.y1 - size[1]//2
                sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 9, weight = 'bold')

        sp.set_title("Zeropoint %s" % (title), fontsize = 30, weight = 'bold')
        sp.set_xlabel("Focal Plane X", fontsize = 20, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 20, weight = 'bold')

        xmin = +1e10
        ymin = +1e10
        xmax = -1e10
        ymax = -1e10
        for r in self.rectangles:
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

    # Helper function to turn Camera class into Matplotlib objects
    def cameraToRectangles(self, camera):
        rectangles = []
        boundaries = []
        for r in self.camera:
            raft = cameraGeom.cast_Raft(r)
        
            # NOTE: all ccd coords are w.r.t. the *center* of the raft, not its LLC
            rxc     = raft.getCenterPixel().getX()
            ryc     = raft.getCenterPixel().getY()
    
            xmin = +1e10
            ymin = +1e10
            xmax = -1e10
            ymax = -1e10
            for c in raft:
                ccd   = cameraGeom.cast_Ccd(c)
                label = ccd.getId().getName()
        
                cxc     = ccd.getCenterPixel().getX()
                cyc     = ccd.getCenterPixel().getY()
                cbbox   = ccd.getAllPixels(True)
                cwidth  = cbbox.getX1() - cbbox.getX0()
                cheight = cbbox.getY1() - cbbox.getY0()
        
                cx0     = rxc + cxc - cwidth/2
                cy0     = ryc + cyc - cheight/2
                cx1     = cx0 + cwidth
                cy1     = cy0 + cheight
    
                if cx0 < xmin:
                    xmin = cx0
                if cx1 > xmax:
                    xmax = cx1
                if cy0 < ymin:
                    ymin = cy0
                if cy1 > ymax:
                    ymax = cy1
        
                crectangle = Rectangle((cx0, cy0), cwidth, cheight, fill = False, label = label)
                rectangles.append(crectangle)
    
            boundaries.append(((xmin, xmin), (ymin, ymax)))
            boundaries.append(((xmin, xmax), (ymin, ymin)))
            boundaries.append(((xmin, xmax), (ymax, ymax)))
            boundaries.append(((xmax, xmax), (ymin, ymax)))
    
        return rectangles, boundaries
        
#
###
#

class ZeropointFigure(FpaFigure):
    def __init__(self, cameraGeomPaf, database):
        FpaFigure.__init__(self, cameraGeomPaf)
        self.dbInterface = LsstSimDbInterface(database)


    def fillValues(self, visitId, filter):
        filterId = self.dbInterface.filterMap[filter]
        
        self.sql  = 'select rm.raftName, cm.ccdName, sce.fluxMag0'
        self.sql += ' from Science_Ccd_Exposure as sce,'
        self.sql += ' RaftMap as rm,'
        self.sql += ' CcdMap as cm'
        self.sql += ' where sce.raft = rm.raftNum '
        self.sql += ' and sce.ccd = cm.ccdNum'
        self.sql += ' and sce.filterId = %d' % (filterId)
        self.sql += ' and sce.visit = %s' % (visitId)
        
        results  = self.dbInterface.execute(self.sql)
        rafts    = num.array([x[0] for x in results])
        ccds     = num.array([x[1] for x in results])
        fluxMag0 = num.array([x[2] for x in results])
        magMag0  = 2.5 * num.log10(fluxMag0)

        for r in self.camera:
            raft     = cameraGeom.cast_Raft(r)
            rlabel   = raft.getId().getName()
            
            raftId   = '%02d'  % (raft.getId().getSerial())
            raftId2  = '%s,%s' % (raftId[0], raftId[1])

            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                
                ccdId  = str(ccd.getId().getSerial())[-2:]
                ccdId2 = '%s,%s' % (ccdId[0], ccdId[1])

                idx = num.where((rafts == raftId2) & (ccds == ccdId2))
                assert len(idx[0]) == 1
                self.values[rlabel][clabel] = magMag0[idx[0]]

                    
class PhotometricRmsFigure(QaFigure):
    def __init__(self, chip = None):
        QaFigure.__init__(self)
        
        if chip == None:
            self.butler = None
        else:
            self.butler = Butler(chip)
            
