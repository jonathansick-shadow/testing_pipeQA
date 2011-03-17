from .DatabaseQuery import LsstSimDbInterface, DatabaseIdentity

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
from lsst.pex.logging import Trace

import os
import numpy as num

import pylab
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse

class HtmlFormatter:
    def __init__(self):
        pass

    def generateFileName(self, prefix, raft, ccd):
        return '%s_%s_%s.png' % (prefix, raft, ccd)

    def generateIndex(self, rootdir, width = 75, height = 75):
        files = os.listdir(rootdir)
        print '<h1>Photometric RMS</h1>'
        for file in files:
            if file.startswith('photRms'):
                print '<a href="%s"><img width="%i" height="%i" border="0" src="%s"></a>' % (file, width,
                                                                                             height, file)
        print '<h1>Zeropoint Across Focal Plane</h1>\n'
        for file in files:
            if file.startswith('zptFPA'):
                print '<a href="%s"><img width="%i" height="%i" border="0" src="%s"></a>' % (file, width,
                                                                                             height, file)
        print '<h1>Zeropoint Fit by Chip</h1>'
        for file in files:
            if file.startswith('zptFit') and file.endswith('.html'):
                print '<a href="%s">%s</a>' % (file, file)
            
    
    def generateHtml(self, buff, prefix, width = 75, height = 75):
        buff.write('<html><body><table>\n')

        for rr in [['',    '4,1', '4,2', '4,3',  '',],
                   ['3,0', '3,1', '3,2', '3,3', '3,4',],
                   ['2,0', '2,1', '2,2', '2,3', '2,4',],
                   ['1,0', '1,1', '1,2', '1,3', '1,4',],
                   ['',    '0,1', '0,2', '0,3', '',]]:
            buff.write('  <tr>\n')
            for r in rr:
                if len(r):
                    buff.write('    <th>raft %s</th>\n' % (r))
                else:
                    buff.write('    <th />\n')

            buff.write('  </tr>\n')
            buff.write('  <tr>\n')
            for r in rr:
                if not len(r):
                    buff.write('    <td />\n')
                    continue

                buff.write('    <td><table>\n')
                for cc in [['2,0', '2,1', '2,2'], 
                           ['1,0', '1,1', '1,2'], 
                           ['0,0', '0,1', '0,2']]:
                    for c in cc:
                        imgname = self.generateFileName(prefix, r, c)
                        buff.write('        <td><a href="%s"><img width="%i" height="%i" border="0" src="%s"></a></td>\n' % (imgname, width, height, imgname))
                    buff.write('      </tr>\n')
                buff.write('    </table></td>\n')
            buff.write('  </tr>\n')
        buff.write('</table>\n')

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
    

def sigIQR(data, min = None, max = None):
    data  = num.sort(data)
    if (min != None) and (max != None):
        idx   = num.where( (data > min) & (data < max) )
        data  = data[idx]

    if len(data) == 0:
        return 0.0
    
    d25 = data[int(0.25 * len(data))]
    d75 = data[int(0.75 * len(data))]
    return 0.741 * (d75 - d25)

#
###
#

class QaFigure():
    def __init__(self):
        self.fig      = pylab.figure()
        self.data     = {}
        self.dataType = {}

    def reset(self):
        self.data       = {}
        self.fig.clf()

    def validate(self):
        dkeys = self.data.keys()
        rkeys = self.dataType.keys()
        for key1 in rkeys:
            if not key1 in dkeys:
                return 0

            if type(self.dataType[key1]) == type({}):
                for key2 in self.dataType[key1].keys():
                    if not key2 in self.data[key1].keys():
                        return 0
                    if not type(self.data[key1][key2]) == self.dataType[key1][key2]:
                        return 0
            else:
                if not type(self.data[key1]) == self.dataType[key1]:
                    return 0
                
        return 1

    def retrieveData(self):
        # override
        pass
    
    def makeFigure(self):
        # override
        pass

    def saveFigure(self, outfile, clear = True):
        Trace("lsst.testing.pipeQA.QaFigure", 1, "Saving %s" % (outfile))
        self.fig.savefig(outfile)
        if clear:
            self.fig.clf()

    #
    ################# Helper plotting functions
    #
        
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

            
class FpaFigure(QaFigure):
    def __init__(self, cameraGeomPaf):
        QaFigure.__init__(self)

        self.cameraGeomPolicy            = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        self.camera                      = cameraGeomUtils.makeCamera(self.cameraGeomPolicy)
        self.rectangles, self.boundaries = self.cameraToRectangles(self.camera)

        # To be filled in by child class
        self.data = {}
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            self.data[rlabel] = {}
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                self.data[rlabel][clabel] = 0.0

    def reset(self):
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            self.data[rlabel] = {}
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                self.data[rlabel][clabel] = 0.0
        
    def validate(self):
        # Since we establish the structure of data in __init__, it
        # should always be valid.  Unless someone mucks with self.data
        # directly...
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            if not self.data.has_key(rlabel):
                return 0
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                if not self.data[rlabel].has_key(clabel):
                    return 0
        return 1

    def makeFigure(self, 
                   DPI = 100., size = (2000, 2000), borderPix = 100,
                   boundaryColors = 'r', doLabel = False):

        if not self.validate():
            Trace("lsst.testing.pipeQA.FpaFigure", 1, "Invalid Data")
            return None

        self.fig.set_size_inches(size[0] / DPI, size[1] / DPI)
        
        sp     = self.fig.gca()
        values = []  # needs to be synchronized with self.rectangles
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                values.append(self.data[rlabel][clabel])

        sigVal = sigIQR(values, min = 0, max = 99.99)
        
        p = PatchCollection(self.rectangles)
        p.set_array(num.array(values))
        cb = self.fig.colorbar(p)
        sp.add_collection(p)

        for b in self.boundaries:
            sp.plot(b[0], b[1], '%s-' % (boundaryColors), lw=3)

        if doLabel:
            for r in self.rectangles:
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
    def __init__(self, cameraGeomPaf):
        FpaFigure.__init__(self, cameraGeomPaf)
        
        # set on retrieve; reset on reset()
        self.database = None
        self.visitId  = None
        self.filter   = None

    def reset(self):
        FpaFigure.reset(self)
        self.database = None
        self.visitId  = None
        self.filter   = None

    def retrieveData(self, database, visitId, defZpt = 29.5):
        self.reset()
        self.database = database
        self.visitId  = visitId
        
        dbId        = DatabaseIdentity(database)
        dbInterface = LsstSimDbInterface(dbId)
        
        sql  = 'select rm.raftName, cm.ccdName, sce.fluxMag0, sce.filterName'
        sql += ' from Science_Ccd_Exposure as sce,'
        sql += ' RaftMap as rm,'
        sql += ' CcdMap as cm'
        sql += ' where sce.raft = rm.raftNum '
        sql += ' and sce.ccd = cm.ccdNum'
        sql += ' and sce.visit = %s' % (visitId)
        results  = dbInterface.execute(sql)
        if len(results) == 0:
            return None

        rafts    = num.array([x[0] for x in results])
        ccds     = num.array([x[1] for x in results])
        fluxMag0 = num.array([x[2] for x in results])
        magMag0  = 2.5 * num.log10(fluxMag0)
        self.filter = results[0][-1]

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
                if len(idx[0]) == 0:
                    self.data[rlabel][clabel] = defZpt
                elif len(idx[0]) == 1:
                    self.data[rlabel][clabel] = magMag0[idx[0]][0]
                else:
                    print "ERROR; TRACK ME DOWN"

    def makeFigure(self, doLabel = False):
        FpaFigure.makeFigure(self, doLabel)
        sp     = self.fig.gca()
        sp.set_title(r"Zeropoint %d %s" % (self.visitId, self.filter),
                     fontsize = 30, weight = 'bold')
        
class PhotometricRmsFigure(QaFigure):
    def __init__(self):
        QaFigure.__init__(self)
        self.data     = {}
        self.dataType = {
            "PhotByObject" : {}
            }

        # set on retrieve; reset on reset()
        self.database = None
        self.filter   = None

    def reset(self):
        self.data     = {}
        self.database = None
        self.filter   = None
        
    def retrieveData(self, database, filter, minPts = 2):
        self.reset()
        self.database = database
        self.filter   = filter

        dbId        = DatabaseIdentity(database)
        dbInterface = LsstSimDbInterface(dbId)

        filterId = dbInterface.filterMap[filter]
        sql  = 'select s.objectId, '                                                         # 0: objectId
        sql += ' sro.%sMag,' % (filter)                                                      # 1: catMag
        sql += ' dnToAbMag(s.apFlux, sce.fluxMag0),'                                         # 2: apMag
        sql += ' dnToAbMagSigma(s.apFlux, s.apFluxSigma, sce.fluxMag0, sce.fluxMag0Sigma),'  # 3: apMagSigma
        sql += ' dnToAbMag(s.psfFlux, sce.fluxMag0),'                                        # 4: psfMag
        sql += ' dnToAbMagSigma(s.psfFlux, s.psfFluxSigma, sce.fluxMag0, sce.fluxMag0Sigma)' # 5: psfMagSigma
        sql += ' from Source as s, Science_Ccd_Exposure as sce,'
        sql += ' RefObjMatch as rom, SimRefObject as sro'
        sql += ' where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
        sql += ' and (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
        sql += ' and (sro.isStar = 1) and (sro.isVar = 0)'
        sql += ' and (s.filterId = %d) and ((s.flagForDetection & 0xa01) = 0)' % (filterId)
        sql += ' and s.objectID is not NULL'        
        sql += ' order by s.objectID'
        print sql
        results = dbInterface.execute(sql)

        photByObject    = {}
        for row in results:
            objectId   = row[0]
            if objectId == None:
                # un-associated orphan
                continue
    
            catMag    = row[1]
            apMag     = row[2]
            apMagSig  = row[3]
            psfMag    = row[4]
            psfMagSig = row[5]
            
            if not photByObject.has_key(objectId):
                photByObject[objectId] = {}
                photByObject[objectId]['catMag']    = []
                photByObject[objectId]['psfMag']    = []
                photByObject[objectId]['psfMagSig'] = []
                photByObject[objectId]['apMag']     = []
                photByObject[objectId]['apMagSig']  = []

            if (apMag == None) or (apMagSig == None) or (psfMag == None) or (psfMagSig == None):
                continue
            
            photByObject[objectId]['catMag'].append(catMag)
            photByObject[objectId]['apMag'].append(apMag)
            photByObject[objectId]['apMagSig'].append(apMagSig)
            photByObject[objectId]['psfMag'].append(psfMag)
            photByObject[objectId]['psfMagSig'].append(psfMagSig)
    
        for key in photByObject.keys():
            if len(photByObject[key]) < minPts:
                del photByObject[key]
                
            photByObject[key]['catMag']    = num.array(photByObject[key]['catMag'])
            photByObject[key]['apMag']     = num.array(photByObject[key]['apMag'])
            photByObject[key]['apMagSig']  = num.array(photByObject[key]['apMagSig'])
            photByObject[key]['psfMag']    = num.array(photByObject[key]['psfMag'])
            photByObject[key]['psfMagSig'] = num.array(photByObject[key]['psfMagSig'])
    
        self.data["PhotByObject"] = photByObject

        
    def makeFigure(self, minmag = 16, maxmag = 23.5, yrange = 1.0,
                   sigmaOffset = 0.5, sigmaRange = 2):
        
        if not self.validate():
            Trace("lsst.testing.pipeQA.ZeropointFitFigure", 1, "Invalid Data")
            return None

        photByObject = self.data["PhotByObject"]
        keys = photByObject.keys()

        if len(keys) == 0:
            Trace("lsst.testing.pipeQA.ZeropointFitFigure", 1, "No data for figure")
            return None
            

        allcatMags     = [photByObject[k]['catMag'] for k in keys]
        alldApMags     = [photByObject[k]['apMag'] for k in keys]
        alldApMagSigs  = [photByObject[k]['apMagSig'] for k in keys]
        alldPsfMags    = [photByObject[k]['psfMag'] for k in keys]
        alldPsfMagSigs = [photByObject[k]['psfMagSig'] for k in keys]

        allcatMags     = num.array([item for sublist in allcatMags for item in sublist])
        alldApMags     = num.array([item for sublist in alldApMags for item in sublist])
        alldApMagSigs  = num.array([item for sublist in alldApMagSigs for item in sublist])
        alldPsfMags    = num.array([item for sublist in alldPsfMags for item in sublist])
        alldPsfMagSigs = num.array([item for sublist in alldPsfMagSigs for item in sublist])

        alldApMags     = allcatMags - alldApMags
        alldPsfMags    = allcatMags - alldPsfMags

        # Get RMS around bright end
        idx = num.where( (allcatMags >= (minmag+sigmaOffset)) &
                         (allcatMags <= (minmag+sigmaOffset+sigmaRange)) )
        sigmeanAp  = sigIQR(alldApMags[idx])
        sigmeanPsf = sigIQR(alldPsfMags[idx])
        
        binnedAp       = self.binDistrib(allcatMags, alldApMags, alldApMagSigs)
        binnedPsf      = self.binDistrib(allcatMags, alldPsfMags, alldPsfMagSigs)

        self.fig.subplots_adjust(hspace=0.0)
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

        
class ZeropointFitFigure(QaFigure):
    def __init__(self):
        QaFigure.__init__(self)
        self.data     = {}
        self.dataType = {
            "Zeropoint"         : num.float64,
            "MatchedGalaxies"   : {"Refmag": num.ndarray,
                                   "Imgmag": num.ndarray,
                                   "Imgerr": num.ndarray},
            "MatchedStars"      : {"Refmag": num.ndarray,
                                   "Imgmag": num.ndarray,
                                   "Imgerr": num.ndarray},
            "UnmatchedReference": num.ndarray,
            "UnmatchedImage"    : num.ndarray
            }

        # set on retrieve; reset on reset()
        self.database   = None
        self.visitId    = None
        self.filterName = None
        self.raftName   = None
        self.ccdName    = None
        
    def reset(self):
        self.data       = {}
        self.database   = None
        self.visitId    = None
        self.filterName = None
        self.raftName   = None
        self.ccdName    = None

    def retrieveData(self, database, visitId, filterName, raftName, ccdName):
        self.reset()
        self.database   = database
        self.visitId    = visitId
        self.filterName = filterName
        self.raftName   = raftName
        self.ccdName    = ccdName

        dbId        = DatabaseIdentity(database)
        dbInterface = LsstSimDbInterface(dbId)

        # Select all reference stars within this field
        # First, get field limits and zeropoint
        scesql  = 'select scienceCcdExposureId,fluxMag0,llcRa,ulcRa,llcDecl,ulcDecl,'
        scesql += ' lrcRa,urcRa,lrcDecl,urcDecl'
        scesql += ' from Science_Ccd_Exposure'
        scesql += ' where visit = %d' % (visitId)
        scesql += ' and raftName = "%s"' % (raftName)
        scesql += ' and ccdName = "%s"' % (ccdName)
        scesql += ' and filterName = "%s"' % (filterName)
        sceresults  = dbInterface.execute(scesql)
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
        self.data["Zeropoint"] = zpt

        # Select all simRefObjects in this field
        # Really need a stored procedure to do this; will fail near poles and around ra=0
        # Do some basic filtering for now
        srosql  = 'select refObjectId,isStar,ra,decl from SimRefObject'
        srosql += ' where (ra >= %f)' % (min( min(llcRa, ulcRa), max(urcRa, lrcRa) ))
        srosql += ' and (ra <= %f)'   % (max( min(llcRa, ulcRa), max(urcRa, lrcRa) ))
        srosql += ' and (decl >= %f)' % (min( min(llcDecl, ulcDecl), max(urcDecl, lrcDecl) ))
        srosql += ' and (decl <= %f)' % (max( min(llcDecl, ulcDecl), max(urcDecl, lrcDecl) ))
        sroresults  = dbInterface.execute(srosql)
        refAll = {'s': [], 'g': []}
        for result in sroresults:
            oid, isStar, ra, decl = result
            if pointInsidePolygon(ra, decl, polygon):
                if isStar: refAll['s'].append(oid)
                else: refAll['g'].append(oid)

        # Select all matched galaxies
        mrefGsql  = 'select sro.%sMag,' % (filterName)
        mrefGsql += ' s.psfFlux, s.psfFluxSigma'
        mrefGsql += ' from SimRefObject as sro, RefObjMatch as rom, Source as s'
        mrefGsql += ' where (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
        mrefGsql += ' and (s.scienceCcdExposureId = %d)' % (sceId)
        mrefGsql += ' and (s.objectID is not NULL)'
        mrefGsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['g'])))
        mrefGresults  = dbInterface.execute(mrefGsql)
        mrefGmag  = num.array([x[0] for x in mrefGresults])
        mimgGflu  = num.array([x[1] for x in mrefGresults])
        mimgGferr = num.array([x[2] for x in mrefGresults])
        mimgGmag  = -2.5 * num.log10(mimgGflu)
        mimgGmerr =  2.5 / num.log(10.0) * mimgGferr / mimgGflu
        self.data["MatchedGalaxies"] = {"Refmag": mrefGmag,
                                        "Imgmag": mimgGmag,
                                        "Imgerr": mimgGmerr}
        
        # Select all matched stars
        mrefSsql  = 'select sro.%sMag,' % (filterName)
        mrefSsql += ' s.psfFlux, s.psfFluxSigma'
        mrefSsql += ' from SimRefObject as sro, RefObjMatch as rom, Source as s'
        mrefSsql += ' where (s.objectId = rom.objectId) and (rom.refObjectId = sro.refObjectId)'
        mrefSsql += ' and (s.scienceCcdExposureId = %d)' % (sceId)
        mrefSsql += ' and (s.objectID is not NULL)'
        mrefSsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['s'])))
        mrefSresults  = dbInterface.execute(mrefSsql)
        mrefSmag  = num.array([x[0] for x in mrefSresults])
        mimgSflu  = num.array([x[1] for x in mrefSresults])
        mimgSferr = num.array([x[2] for x in mrefSresults])
        mimgSmag  = -2.5 * num.log10(mimgSflu)
        mimgSmerr =  2.5 / num.log(10.0) * mimgSferr / mimgSflu
        self.data["MatchedStars"] = {"Refmag": mrefSmag,
                                     "Imgmag": mimgSmag,
                                     "Imgerr": mimgSmerr}

        ####

        # Umatched reference catalog objects
        urefsql  = 'select sro.%sMag from SimRefObject as sro, RefObjMatch as rom' % (filterName)
        urefsql += ' where (sro.refObjectId = rom.refObjectId)'
        urefsql += ' and (rom.objectId is NULL)'
        urefsql += ' and (sro.refObjectId in (%s))' % (','.join(map(str, refAll['g'] + refAll['s'])))
        urefresults = dbInterface.execute(urefsql)
        urefmag     = num.array([x[0] for x in urefresults])
        self.data["UnmatchedReference"] = urefmag
            

        # Unmatched detections
        uimgsql  = 'select psfFlux from Source '
        uimgsql += ' where (scienceCcdExposureId = %d)' % (sceId)
        uimgsql += ' and (objectId is NULL)'
        uimgresults  = dbInterface.execute(uimgsql)
        uimgmag      = -2.5 * num.log10( num.array([x[0] for x in uimgresults]) )
        self.data["UnmatchedImage"] = uimgmag


    def makeFigure(self):
        if not self.validate():
            Trace("lsst.testing.pipeQA.ZeropointFitFigure", 1, "Invalid Data")
            return None

        legLines  = []
        legLabels = []
        
        axis = self.fig.add_axes([0.225, 0.225, 0.675, 0.550])

        # Plot all matched galaxies
        mrefGmag  = self.data["MatchedGalaxies"]["Refmag"]
        mimgGmag  = self.data["MatchedGalaxies"]["Imgmag"]
        mimgGmerr = self.data["MatchedGalaxies"]["Imgerr"]
        mimgGplot = axis.plot(mimgGmag, mrefGmag, '.', color='g', mfc='g', mec='g',
                              alpha=0.5, zorder=10, label = 'Matched Galaxies')
        for i in range(len(mrefGmag)):
            a = Ellipse(xy=num.array([mimgGmag[i], mrefGmag[i]]),
                        width=mimgGmerr[i]/2., height=mimgGmerr[i]/2.,
                        alpha=0.5, fill=True, ec='g', fc='g', zorder=10)
            axis.add_artist(a)            
        legLines.append(mimgGplot[0])
        legLabels.append("Matched Galaxies")
        
        # Plot all matched stars
        mrefSmag  = self.data["MatchedStars"]["Refmag"]
        mimgSmag  = self.data["MatchedStars"]["Imgmag"]
        mimgSmerr = self.data["MatchedStars"]["Imgerr"]
        mimgSplot = axis.plot(mimgSmag, mrefSmag, '.', color='b', mfc='b', mec='b',
                              alpha=0.5, zorder=12, label = 'Matched Stars')
        for i in range(len(mrefSmag)):
            a = Ellipse(xy=num.array([mimgSmag[i], mrefSmag[i]]),
                        width=mimgSmerr[i]/2., height=mimgSmerr[i]/2.,
                        alpha=0.5, fill=True, ec='b', fc='b', zorder=12)
            axis.add_artist(a)            
        legLines.append(mimgSplot[0])
        legLabels.append("Matched Stars")

        if len(mrefGmag) == 0 and len(mrefSmag) == 0:
            xmin, xmax, ymin, ymax = -15, -8, 16, 28
        else:
            xmin, xmax, ymin, ymax = axis.axis()

        # Plot zpt
        xzpt = num.array((xmin, xmax))
        pzpt = axis.plot(xzpt, xzpt - self.data["Zeropoint"], 'b--', label = 'Zeropoint')
        legLines.append(pzpt)
        legLabels.append("Zeropoint")

        # Unmatched objects
        urefmag     = self.data["UnmatchedReference"]
        uimgmag     = self.data["UnmatchedImage"]

        # Unmatched & matched reference objects
        ax2  = self.fig.add_axes([0.1,   0.225, 0.125, 0.550], sharey=axis)
        nu, bu, pu = ax2.hist(urefmag, bins=num.arange(ymin, ymax, 0.25),
                              orientation='horizontal', log = True, color = 'r', alpha = 0.5, zorder = 1)
        if len(mrefGmag) > 0 and len(mrefSmag) > 0:
            ax2.hist(num.concatenate((mrefGmag,mrefSmag)), bins=num.arange(ymin, ymax, 0.25),
                     orientation='horizontal', log = True, color = 'b', alpha = 0.5, zorder = 2)
        ax2.set_xlabel('N', fontsize = 10)
        ax2.set_ylabel('Reference catalog: %s band (mag)' % (self.filterName), fontsize = 10)
        legLines.append(pu[0])
        legLabels.append("Unmatched Sources")

        # Unmatched & matched stellar objects
        ax3  = self.fig.add_axes([0.225, 0.1,   0.675, 0.125], sharex=axis)
        ax3.get_yaxis().set_ticks_position('right')
        ax3.get_yaxis().set_label_position('right')
        if len(mimgGmag) > 0 and len(mimgSmag) > 0:
            nm, bm, pm = ax3.hist(num.concatenate((mimgGmag,mimgSmag)), bins=num.arange(xmin, xmax, 0.25),
                                  log = True, color = 'b', alpha = 0.5, zorder = 2)
            legLines.append(pm[0])
            legLabels.append("Matched Sources")
        ax3.hist(uimgmag, bins=num.arange(xmin, xmax, 0.25),
                 log = True, color = 'r', alpha = 0.5, zorder = 1)
        ax3.set_xlabel('Image instrumental mag', fontsize = 10)
        ax3.set_ylabel('N', rotation = 180, fontsize = 10)

        # Mag - Zpt
        ax4  = self.fig.add_axes([0.225, 0.775, 0.675, 0.125], sharex=axis)
        ax4.errorbar(mimgSmag, (mimgSmag - self.data["Zeropoint"]) - mrefSmag, yerr = mimgSmerr, fmt = 'bo',
                     ms = 2, alpha = 0.5, capsize = 0, elinewidth = 0.5)
        ax4.errorbar(mimgGmag, (mimgGmag - self.data["Zeropoint"]) - mrefGmag, yerr = mimgGmerr, fmt = 'go',
                     ms = 2, alpha = 0.5, capsize = 0, elinewidth = 0.5)
        ax4.get_yaxis().set_ticks_position('right')
        ax4.get_yaxis().set_label_position('right')
        ax4.set_ylabel('Cal-Ref', fontsize = 10, rotation = 270)
        ax4.axhline(y = 0, c='k', linestyle='--', alpha = 0.25)

        # Cleaning up figure
        pylab.setp(axis.get_xticklabels()+axis.get_yticklabels(), visible=False)
        pylab.setp(ax2.get_xticklabels()+ax2.get_yticklabels(), fontsize = 8)
        pylab.setp(ax3.get_xticklabels()+ax3.get_yticklabels(), fontsize = 8)
        pylab.setp(ax4.get_xticklabels(), visible=False)
        pylab.setp(ax4.get_yticklabels(), fontsize = 8)

        ax2.set_xlim(1, 999)
        ax3.set_ylim(1, 999)
        ax4.set_ylim(-0.99, 0.99)
        axis.axis((xmax, xmin, ymax, ymin))

        self.fig.legend(legLines, legLabels,
                        numpoints=1, prop=FontProperties(size='small'), loc = 'center right')
        self.fig.suptitle('%s v%s r%s s%s' %
                          (self.database, self.visitId, self.raftName, self.ccdName),
                          fontsize = 12)
        return
