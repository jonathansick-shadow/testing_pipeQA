from .DatabaseQuery import LsstSimDbInterface, DatabaseIdentity

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
from lsst.pex.logging import Trace

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

def IQR(data):
    data  = num.sort(data)
    d25 = data[int(0.25 * len(data))]
    d75 = data[int(0.75 * len(data))]
    return 0.741 * (d75 - d25)

class QaFigure():
    def __init__(self):
        self.butler = None
        self.fig = pylab.figure()

    def binDistrib(self, x, y, dy, binSizeX = 0.5, minPts = 2):
        bx  = []
        bs  = []
        by  = []
        bdy = []
        
        bins  = num.arange(min(x) + binSizeX, max(x), binSizeX)
        for i in range(len(bins) - 1):
            idx = num.where((x >= bins[i]) & (x < bins[i+1]))
            if len(idx[0]) < minPts:
                continue

            # find median x-values in this bin
            xmed = num.median(x[idx])
            bx.append(xmed)

            # find median y-values in this bin
            ymed = num.median(y[idx])
            by.append(ymed)

            # find width about mean y using IQR
            bs.append(IQR(y[idx]))

            # find median dy-values in this bin
            dymed = num.median(dy[idx])
            bdy.append(dymed)

        return num.array(bx), num.array(by), num.array(bs), num.array(bdy)


    def plotSparseContour(self, sp, x, y, binSizeX, binSizeY, minCont = 50, nCont = 7):
        idx   = num.isfinite(x)
        x     = x[idx]
        y     = y[idx]
        
        xmin  = x.min()
        xmax  = x.max()
        ymin  = y.min()
        ymax  = y.max()
        xbins = num.arange(xmin, xmax, binSizeX)
        ybins = num.arange(ymin, ymax, binSizeY)

        cdata = num.zeros((len(ybins), len(xbins)))
        for i in range(len(x)):
            xidx = (x[i] - xmin) // binSizeX
            yidx = (y[i] - ymin) // binSizeY
            cdata[yidx][xidx] += 1
    
        if cdata.max() < minCont:
            minCont = 1
            
        cs    = num.arange(minCont, cdata.max(), (cdata.max() - minCont) // nCont).astype(num.int)
        c     = sp.contour(cdata, cs, origin='lower', linewidths=1, extent=(xmin,xmax,ymin,ymax))
        outer = c.collections[0]._paths
    
        xp = []
        yp = []
        for i in range(len(x)):
            found = [o.contains_point((x[i], y[i])) for o in outer]
            if not (True in found):
                xp.append(x[i])
                yp.append(y[i])
        sp.plot(xp, yp, 'r.', ms = 1)
    
    def makeFigure(self):
        # override
        pass

    def saveFigure(self, outfile, clear = True):
        Trace("lsst.testing.pipeQA.QaFigure", 1, "Saving %s" % (outfile))
        self.fig.savefig(outfile)
        if clear:
            self.fig.clf()

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

        sigZpt = IQR(values)
        
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
                sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 8, weight = 'bold')

        sp.set_title(r"Zeropoint %s: $\sigma = %.3f$ mag" % (title, sigZpt), fontsize = 30, weight = 'bold')
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
        dbId = DatabaseIdentity(database)
        self.dbInterface = LsstSimDbInterface(dbId)


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
    def __init__(self, database, filter):
        QaFigure.__init__(self)
        self.database = database
        self.filter = filter
        
        dbId = DatabaseIdentity(self.database)
        self.dbInterface = LsstSimDbInterface(dbId)
        self.makeFigure(self.filter)

    def makeFigure(self, filter):
        filterId = self.dbInterface.filterMap[filter]
        
        sql  = 'select s.objectId, '                                                         # 0: objectId
        sql += ' sro.%sMag,' % (filter)                                                      # 1: catMag
        sql += ' dnToAbMag(s.apFlux, sce.fluxMag0),'                                         # 2: apMag
        sql += ' dnToAbMagSigma(s.apFlux, s.apFluxErr, sce.fluxMag0, sce.fluxMag0Sigma),'    # 3: apMagErr
        sql += ' dnToAbMag(s.psfFlux, sce.fluxMag0),'                                        # 4: psfMag
        sql += ' dnToAbMagSigma(s.psfFlux, s.psfFluxErr, sce.fluxMag0, sce.fluxMag0Sigma)'   # 5: psfMagErr
        sql += ' from Source as s, Science_Ccd_Exposure as sce,'
        sql += ' RefObjMatch as rom, SimRefObject as sro'
        sql += ' where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
        sql += ' and (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
        sql += ' and (sro.isStar = 1) and (sro.isVar = 0)'
        sql += ' and (s.filterId = %d) and ((s.flagForDetection & 0xa01) = 0)' % (filterId)
        sql += ' and s.objectID is not NULL'        
        sql += ' order by s.objectID'
        
        results = self.dbInterface.execute(sql)
        allPhot = self.parsePhot(results)
        self.plotPhot(allPhot)
        
    def plotPhot(self, allPhot, minmag = 16, maxmag = 23.5, yrange = 1.0,
                 sigmaOffset = 0.5, sigmaRange = 2):
        keys = allPhot.keys()

        allcatMags     = [allPhot[k]['catMag'] for k in keys]
        alldApMags     = [allPhot[k]['apMag'] for k in keys]
        alldApMagErrs  = [allPhot[k]['apMagErr'] for k in keys]
        alldPsfMags    = [allPhot[k]['psfMag'] for k in keys]
        alldPsfMagErrs = [allPhot[k]['psfMagErr'] for k in keys]

        allcatMags     = num.array([item for sublist in allcatMags for item in sublist])
        alldApMags     = num.array([item for sublist in alldApMags for item in sublist])
        alldApMagErrs  = num.array([item for sublist in alldApMagErrs for item in sublist])
        alldPsfMags    = num.array([item for sublist in alldPsfMags for item in sublist])
        alldPsfMagErrs = num.array([item for sublist in alldPsfMagErrs for item in sublist])

        alldApMags     = allcatMags - alldApMags
        alldPsfMags    = allcatMags - alldPsfMags

        # Get RMS around bright end
        idx = num.where( (allcatMags >= (minmag+sigmaOffset)) &
                         (allcatMags <= (minmag+sigmaOffset+sigmaRange)) )
        sigmeanAp  = IQR(alldApMags[idx])
        sigmeanPsf = IQR(alldPsfMags[idx])
        
        binnedAp       = self.binDistrib(allcatMags, alldApMags, alldApMagErrs)
        binnedPsf      = self.binDistrib(allcatMags, alldPsfMags, alldPsfMagErrs)

        pylab.subplots_adjust(hspace=0.0)
        sp1 = self.fig.add_subplot(211)
        sp2 = self.fig.add_subplot(212, sharex = sp1)

        self.plotSparseContour(sp1, allcatMags, alldApMags, binSizeX = 0.25, binSizeY = 0.1)
        self.plotSparseContour(sp2, allcatMags, alldPsfMags, binSizeX = 0.25, binSizeY = 0.1)

        sp1.plot(binnedAp[0], binnedAp[1]+binnedAp[2], 'bv', alpha = 0.5)
        sp1.plot(binnedAp[0], binnedAp[1]-binnedAp[2], 'b^', alpha = 0.5)

        sp2.plot(binnedPsf[0], binnedPsf[1]+binnedPsf[2], 'bv', alpha = 0.5)
        sp2.plot(binnedPsf[0], binnedPsf[1]-binnedPsf[2], 'b^', alpha = 0.5)
        
        sp1.text(0.1, 0.9, r"$\sigma_{\rm %.1f - %.1f} = %.3f$" % ((minmag+sigmaOffset),
                                                                   (minmag+sigmaOffset+sigmaRange),
                                                                   sigmeanAp), transform = sp1.transAxes)
        sp2.text(0.1, 0.9, r"$\sigma_{\rm %.1f - %.1f} = %.3f$" % ((minmag+sigmaOffset),
                                                                   (minmag+sigmaOffset+sigmaRange),
                                                                   sigmeanPsf), transform = sp2.transAxes)
        
        pylab.setp(sp1.get_xticklabels(), visible=False)
        sp1.set_title("Photometric RMS: %s %s" % (self.database, self.filter),
                      fontsize = 12, weight = 'bold')
        sp1.set_ylabel(r"$\Delta\ $Ap (Cat - Meas)", fontsize = 12, weight = 'bold')
        sp2.set_ylabel(r"$\Delta\ $Psf (Cat - Meas)", fontsize = 12, weight = 'bold')
        sp2.set_xlabel("Catalog Mag", fontsize = 12, weight = 'bold')
        sp2.set_xlim(minmag, maxmag)

        sp1.set_ylim(num.median(binnedAp[1]) - yrange/2, num.median(binnedAp[1]) + yrange/2)
        sp2.set_ylim(num.median(binnedPsf[1]) - yrange/2, num.median(binnedPsf[1]) + yrange/2)
        
    def parsePhot(self, results, minPts = 2):
        allPhot    = {}
        for row in results:
            objectId   = row[0]
            if objectId == None:
                # un-associated orphan
                continue
    
            catMag    = row[1]
            apMag     = row[2]
            apMagErr  = row[3]
            psfMag    = row[4]
            psfMagErr = row[5]
            
            if not allPhot.has_key(objectId):
                allPhot[objectId] = {}
                allPhot[objectId]['catMag']    = []
                allPhot[objectId]['psfMag']    = []
                allPhot[objectId]['psfMagErr'] = []
                allPhot[objectId]['apMag']     = []
                allPhot[objectId]['apMagErr']  = []
            
            allPhot[objectId]['catMag'].append(catMag)
            allPhot[objectId]['apMag'].append(apMag)
            allPhot[objectId]['apMagErr'].append(apMagErr)
            allPhot[objectId]['psfMag'].append(psfMag)
            allPhot[objectId]['psfMagErr'].append(psfMagErr)
    
        for key in allPhot.keys():
            if len(allPhot[key]) < minPts:
                del allPhot[key]
                
            allPhot[key]['catMag']    = num.array(allPhot[key]['catMag'])
            allPhot[key]['apMag']     = num.array(allPhot[key]['apMag'])
            allPhot[key]['apMagErr']  = num.array(allPhot[key]['apMagErr'])
            allPhot[key]['psfMag']    = num.array(allPhot[key]['psfMag'])
            allPhot[key]['psfMagErr'] = num.array(allPhot[key]['psfMagErr'])
    
        return allPhot
