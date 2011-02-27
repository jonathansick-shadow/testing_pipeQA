from .DatabaseQuery import LsstSimDbInterface, DatabaseIdentity

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
from lsst.pex.logging import Trace

import numpy as num

import pylab
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse

class HtmlFormatter:
    def __init__(self, butler, cameraGeom):
        self.figureButler = butler
        self.cameraGeom   = cameraGeom

    def generateHtml():
        pass

#
###
#

def pointInsidePolygon(x,y,poly):
    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x,p1y = p2x,p2y
    return inside
    

def sigIQR(data):
    data  = num.sort(data)
    d25 = data[int(0.25 * len(data))]
    d75 = data[int(0.75 * len(data))]
    return 0.741 * (d75 - d25)

#
###
#

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

            # find width about mean y using sigIQR
            bs.append(sigIQR(y[idx]))

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

        sigZpt = sigIQR(values)
        
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

class ZeropointFpaFigure(FpaFigure):
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
        sigmeanAp  = sigIQR(alldApMags[idx])
        sigmeanPsf = sigIQR(alldPsfMags[idx])
        
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

class ZeropointFitFigure(QaFigure):
    def __init__(self, database, visitId, filterName, raftName, ccdName):
        QaFigure.__init__(self)
        dbId = DatabaseIdentity(database)
        self.dbInterface = LsstSimDbInterface(dbId)

        self.database   = database
        self.visitId    = visitId
        self.raftName   = raftName
        self.ccdName    = ccdName
        self.filterName = filterName

        self.getData()
    
    def getData(self):
        axis = self.fig.gca()
        
        # select all reference stars within this field
        # first, get field limits
        scesql  = 'select scienceCcdExposureId,fluxMag0,llcRa,ulcRa,llcDecl,ulcDecl,'
        scesql += ' lrcRa,urcRa,lrcDecl,urcDecl'
        scesql += ' from Science_Ccd_Exposure'
        scesql += ' where visit = %d' % (self.visitId)
        scesql += ' and raftName = "%s"' % (self.raftName)
        scesql += ' and ccdName = "%s"' % (self.ccdName)
        scesql += ' and filterName = "%s"' % (self.filterName)
        sceresults  = self.dbInterface.execute(scesql)
        if len(sceresults) != 1:
            # throw exception or something
            return
        sceId,fmag0,llcRa,ulcRa,llcDecl,ulcDecl,lrcRa,urcRa,lrcDecl,urcDecl = sceresults[0]
        polygon = ( (llcRa,llcDecl),
                    (ulcRa,ulcDecl),
                    (urcRa,urcDecl),
                    (lrcRa,lrcDecl),
                    (llcRa,llcDecl) )
        zpt = -2.5 * num.log10(fmag0)
        
        # select all simRefObjects in this field
        # really need a stored procedure to do this; will fail near poles and around ra=0
        # do some basic filtering
        srosql  = 'select refObjectId,isStar,ra,decl from SimRefObject'
        srosql += ' where (ra >= %f)' % (min( min(llcRa, ulcRa), max(urcRa, lrcRa) ))
        srosql += ' and (ra <= %f)'   % (max( min(llcRa, ulcRa), max(urcRa, lrcRa) ))
        srosql += ' and (decl >= %f)' % (min( min(llcDecl, ulcDecl), max(urcDecl, lrcDecl) ))
        srosql += ' and (decl <= %f)' % (max( min(llcDecl, ulcDecl), max(urcDecl, lrcDecl) ))
        sroresults  = self.dbInterface.execute(srosql)
        refAll = {'s': [], 'g': []}
        for result in sroresults:
            oid, isStar, ra, decl = result
            if pointInsidePolygon(ra, decl, polygon):
                if isStar: refAll['s'].append(oid)
                else: refAll['g'].append(oid)

        # select all matched galaxies
        mrefGsql  = 'select sro.%sMag,' % (self.filterName)
        mrefGsql += ' s.psfFlux, s.psfFluxErr'
        mrefGsql += ' from SimRefObject as sro, RefObjMatch as rom, Source as s'
        mrefGsql += ' where (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
        mrefGsql += ' and (s.scienceCcdExposureId = %d)' % (sceId)
        mrefGsql += ' and (s.objectID is not NULL)'
        mrefGsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['g'])))
        mrefGresults  = self.dbInterface.execute(mrefGsql)
        mrefGmag  = num.array([x[0] for x in mrefGresults])
        mimgGflu  = num.array([x[1] for x in mrefGresults])
        mimgGferr = num.array([x[2] for x in mrefGresults])
        mimgGmag  = -2.5 * num.log10(mimgGflu)
        mimgGmerr =  2.5 / num.log(10.0) * mimgGferr / mimgGflu
        axis.plot(mimgGmag, mrefGmag, '.', color='g', mfc='g', mec='g',
                  alpha=0.5, zorder=10, label = 'Matched Galaxies')
        for i in range(len(mrefGmag)):
            a = Ellipse(xy=num.array([mimgGmag[i], mrefGmag[i]]),
                        width=mimgGmerr[i]/2., height=mimgGmerr[i]/2.,
                        alpha=0.5, fill=True, ec='g', fc='g', zorder=10)
            axis.add_artist(a)            
        
        # select all matched stars
        mrefSsql  = 'select sro.%sMag,' % (self.filterName)
        mrefSsql += ' s.psfFlux, s.psfFluxErr'
        mrefSsql += ' from SimRefObject as sro, RefObjMatch as rom, Source as s'
        mrefSsql += ' where (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
        mrefSsql += ' and (s.scienceCcdExposureId = %d)' % (sceId)
        mrefSsql += ' and (s.objectID is not NULL)'
        mrefSsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['s'])))
        mrefSresults  = self.dbInterface.execute(mrefSsql)
        mrefSmag  = num.array([x[0] for x in mrefSresults])
        mimgSflu  = num.array([x[1] for x in mrefSresults])
        mimgSferr = num.array([x[2] for x in mrefSresults])
        mimgSmag  = -2.5 * num.log10(mimgSflu)
        mimgSmerr =  2.5 / num.log(10.0) * mimgSferr / mimgSflu
        axis.plot(mimgSmag, mrefSmag, '.', color='b', mfc='b', mec='b',
                  alpha=0.5, zorder=12, label = 'Matched Stars')
        for i in range(len(mrefSmag)):
            a = Ellipse(xy=num.array([mimgSmag[i], mrefSmag[i]]),
                        width=mimgSmerr[i]/2., height=mimgSmerr[i]/2.,
                        alpha=0.5, fill=True, ec='b', fc='b', zorder=12)
            axis.add_artist(a)            

        ####

        # Umatched reference catalog objects
        urefsql  = 'select sro.%sMag from SimRefObject as sro, RefObjMatch as rom' % (self.filterName)
        urefsql += ' where (sro.refObjectId = rom.refObjectId)'
        urefsql += ' and (rom.objectId is NULL)'
        urefsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['g'] + refAll['s'])))
        urefresults = self.dbInterface.execute(urefsql)
        urefmag     = num.array([x[0] for x in urefresults])

        # Unmatched detections
        uimgsql  = 'select psfFlux from Source '
        uimgsql += ' where (scienceCcdExposureId = %d)' % (sceId)
        uimgsql += ' and (objectId is NULL)'
        uimgresults  = self.dbInterface.execute(uimgsql)
        uimgmag      = -2.5 * num.log10( num.array([x[0] for x in uimgresults]) )

        # Plot zpt
        xmin, xmax, ymin, ymax = axis.axis()
        xzpt = num.array((xmin, xmax))
        pzpt = axis.plot(xzpt, xzpt - zpt, 'b--', label = 'Zeropoint')

        # Plot ticks

        # Red tick marks show unmatched img sources
        dy       = (ymax - ymin) * 0.05
        y1       = num.ones_like(uimgmag) * ymax
        print 'A', len(y1)
        uimgplot = axis.plot(num.vstack((uimgmag, uimgmag)), num.vstack((y1, y1-dy)), 'r-',
                             alpha=0.5)
        uimgplot[0].set_label('Unmatched Sources')
        
        # Blue tick marks show matched img sources
        yG      = num.ones_like(mimgGmag) * ymax
        print 'B', len(yG)
        miGplot = axis.plot(num.vstack((mimgGmag, mimgGmag)), num.vstack((yG-(0.25*dy), yG-(1.25*dy))),
                         'b-', alpha=0.5)
        yS      = num.ones_like(mimgSmag) * ymax
        print 'C', len(yS)
        miSplot = axis.plot(num.vstack((mimgSmag, mimgSmag)), num.vstack((yS-(0.25*dy), yS-(1.25*dy))),
                         'b-', alpha=0.5)
        miSplot[0].set_label('Matched Sources')

        # Red ticks for unmatched ref sources
        dx = (xmax - xmin) * 0.05
        x1 = num.ones_like(urefmag) * xmax
        print 'D', len(x1)
        #axis.plot(num.vstack((x1, x1-dx)), num.vstack((urefmag, urefmag)), 'r-',
        #          alpha=0.5, label = '_nolegend_')

        # Blue ticks for matched ref sources
        xG      = num.ones_like(mrefGmag) * xmax
        print 'E', len(xG)
        #mrGplot = axis.plot(num.vstack((x1-(0.25*dx), x1-(1.25*dx))), num.vstack((mrefGmag, mrefGmag)),
        #                    'b-', alpha=0.5)
        xS      = num.ones_like(mrefSmag) * xmax
        print 'F', len(xS)
        #mrSplot = axis.plot(num.vstack((x1-(0.25*dx), x1-(1.25*dx))), num.vstack((mrefSmag, mrefSmag)),
        #                    'b-', alpha=0.5)

        # Switch min<->max
        axis.axis((xmax, xmin, ymax, ymin))

        axis.legend(loc = 9, numpoints=1, prop=FontProperties(size='small'))
        axis.set_xlabel('Image instrumental mag')
        axis.set_ylabel('Reference catalog: %s band (mag)' % (self.filterName))
        axis.set_title('%s v%s r%s s%s' %
                       (self.database, self.visitId, self.raftName, self.ccdName),
                       fontsize = 12)
