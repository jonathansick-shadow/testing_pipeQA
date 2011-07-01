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
from matplotlib.collections import LineCollection


class PsfShapeQaAnalysis(qaAna.QaAnalysis):

    def __init__(self, ellipMax, fwhmMax, **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.limitsEllip = [0.0, ellipMax]
        self.limitsFwhm = [0.0, fwhmMax]

        self.description = """
         For each CCD, the ellipticity of stars used in the Psf model are
         plotted as a function of position in the focal plane.  The summary FPA
         figures show the median vector (offset and angle) of this ellipticity
         for each chip, as well as the effective FWHM in arcsec for the final
         Psf model.
        """

    def free(self):
        del self.theta
        del self.ellip
        del self.y
        del self.x
        del self.filter
        del self.detector
        del self.ssDict

        del self.ellipMedians
        del self.thetaMedians
        del self.fwhm

        del self.calexpDict
        
    def test(self, data, dataId):
        
        # get data
        self.ssDict        = data.getSourceSetBySensor(dataId)
        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        self.calexpDict    = data.getCalexpBySensor(dataId)

        # create containers for data in the focal plane
        self.x     = raftCcdData.RaftCcdVector(self.detector)
        self.y     = raftCcdData.RaftCcdVector(self.detector)
        self.ellip = raftCcdData.RaftCcdVector(self.detector)
        self.theta = raftCcdData.RaftCcdVector(self.detector)

        # compute values of interest
        filter = None
        for key, ss in self.ssDict.items():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()

            filter = self.filter[key].getName()

            qaAnaUtil.isStar(ss)
            
            for s in ss:
                ixx = s.getIxx()
                iyy = s.getIyy()
                ixy = s.getIxy()

                tmp = 0.25*(ixx-iyy)**2 + ixy**2
                if tmp < 0:
                    continue

                a2 = 0.5*(ixx+iyy) + numpy.sqrt(tmp)
                b2 = 0.5*(ixx+iyy) - numpy.sqrt(tmp)

                if b2/a2 < 0:
                    continue
                
                ellip = 1.0 - numpy.sqrt(b2/a2)
                theta = 0.5*numpy.arctan2(2.0*ixy, ixx-iyy)

                # vectors have no direction, so default to pointing in +ve 'y'
                # - failing to do this caused a stats bug when alignment is near pi/2
                #   both +/- pi/2 arise but are essentially the same, ... and the mean is near zero
                if theta < 0.0:
                    theta += numpy.pi
                    
                #print ixx, iyy, ixy, a2, b2, ellip, theta
                isStar = s.getFlagForDetection() & measAlg.Flags.STAR
                
                if numpy.isfinite(ellip) and numpy.isfinite(theta) and isStar:
                    self.ellip.append(raft, ccd, ellip)
                    self.theta.append(raft, ccd, theta)
                    self.x.append(raft, ccd, s.getXAstrom())
                    self.y.append(raft, ccd, s.getYAstrom())
                
        # create a testset and add values
        testSet = self.getTestSet(data, dataId)

        # gets the stats for each sensor and put the values in the raftccd container
        self.ellipMedians = raftCcdData.RaftCcdData(self.detector)
        self.thetaMedians = raftCcdData.RaftCcdData(self.detector)
        for raft, ccd in self.ellip.raftCcdKeys():
            ellip = self.ellip.get(raft, ccd)
            theta = self.theta.get(raft, ccd)

            if len(ellip) > 0:
                stat = afwMath.makeStatistics(ellip, afwMath.NPOINT | afwMath.MEDIAN)
                ellipMed = stat.getValue(afwMath.MEDIAN)
                stat = afwMath.makeStatistics(theta, afwMath.NPOINT | afwMath.MEDIAN)
                thetaMed = stat.getValue(afwMath.MEDIAN)
                n      = stat.getValue(afwMath.NPOINT)
            else:
                ellipMed = -1.0
                thetaMed = 0.0
                n = 0

            # add a test for acceptible psf ellipticity
            self.ellipMedians.set(raft, ccd, ellipMed)
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "median psf ellipticity "
            comment = "median psf ellipticity (nstar=%d)" % (n)
            testSet.addTest( testCode.Test(label, ellipMed, self.limitsEllip, comment, areaLabel=areaLabel) )

            # stash the angles.  We'll use them to make figures in plot()
            self.thetaMedians.set(raft, ccd, thetaMed)
            

        # And the Fwhm
        self.fwhm  = raftCcdData.RaftCcdData(self.detector)
        for key, item in self.calexpDict.items():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()

            self.fwhm.set(raft, ccd, item['fwhm'])
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "psf fwhm (arcsec) "
            comment = "psf fwhm (arcsec)"
            testSet.addTest( testCode.Test(label, item['fwhm'], self.limitsFwhm, comment, areaLabel=areaLabel) )
                          
    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)

        vLen = 1000.0  # for e=1.0

        # fpa figures
        ellipBase = "medPsfEllip"
        ellipData, ellipMap = testSet.unpickle(ellipBase, default=[None, None])
        ellipFig = qaFig.VectorFpaQaFigure(data.cameraInfo, data=ellipData, map=ellipMap)

        fwhmBase = "psfFwhm"
        fwhmData, fwhmMap = testSet.unpickle(fwhmBase, default=[None, None])
        fwhmFig = qaFig.FpaQaFigure(data.cameraInfo, data=fwhmData, map=fwhmMap)

        fwhmMin =  1e10
        fwhmMax = -1e10
        for raft, ccdDict in ellipFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.ellipMedians.get(raft, ccd) is None:
                    ellipFig.data[raft][ccd] = [self.thetaMedians.get(raft, ccd),
                                                10*vLen*self.ellipMedians.get(raft, ccd),
                                                self.ellipMedians.get(raft, ccd)]
                    ellipFig.map[raft][ccd] = "ell/theta=%.3f/%.0f" % (self.ellipMedians.get(raft, ccd),
                                                                       (180/numpy.pi)*self.thetaMedians.get(raft, ccd))
                if not self.fwhm.get(raft, ccd) is None:
                    fwhm = self.fwhm.get(raft, ccd)
                    if fwhm > fwhmMax:
                        fwhmMax = fwhm
                    if fwhm < fwhmMin:
                        fwhmMin = fwhm
                    fwhmFig.data[raft][ccd] = fwhm
                    fwhmFig.map[raft][ccd] = "fwhm=%.2f asec" % (fwhm)
                    
                
        ellipFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=self.limitsEllip,
                            title="Median PSF Ellipticity", failLimits=self.limitsEllip)
        testSet.addFigure(ellipFig, ellipBase+".png", "Median PSF Ellipticity", navMap=True)
        testSet.pickle(ellipBase, [ellipFig.data, ellipFig.map])

        blue = '#0000ff'
        red = '#ff0000'

        if fwhmMin < 1e10:
            vlimMin = numpy.max([self.limitsFwhm[0], fwhmMin])
        else:
            vlimMin = self.limitsFwhm[0]
        if fwhmMax > -1e10:
            vlimMax = numpy.min([self.limitsFwhm[1], fwhmMax])
        else:
            vlimMax = self.limitsFwhm[1]
            
        fwhmFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=[vlimMin, vlimMax],
                           title="PSF FWHM (arcsec)", cmapOver=red, failLimits=self.limitsFwhm,
                           cmapUnder=blue)
        testSet.addFigure(fwhmFig, fwhmBase + ".png", "FWHM of Psf (arcsec)", navMap=True)
        testSet.pickle(fwhmBase, [fwhmFig.data, fwhmFig.map])

                        
        #
        figsize = (4.0, 4.0)
        
        #xlim = [0, 25.0]
        #ylim = [0, 0.4]

        conv = colors.ColorConverter()
        black = conv.to_rgb('k')

        i = 0
        xmin, xmax = 1.0e99, -1.0e99
        for raft, ccd in self.ellip.raftCcdKeys():
            eLen = vLen*self.ellip.get(raft, ccd)
            
            t = self.theta.get(raft, ccd)
            dx = eLen*numpy.cos(t)
            dy = eLen*numpy.sin(t)
            x = self.x.get(raft, ccd) - 0.5*dx
            y = self.y.get(raft, ccd) - 0.5*dy

            if len(x) == 0:
                x = numpy.array([0.0])
                y = numpy.array([0.0])
                dx = numpy.array([0.0])
                dy = numpy.array([0.0])
                
            xmax, ymax = x.max(), y.max()
            xlim = [0, 1024*int(xmax/1024.0 + 0.5)]
            ylim = [0, 1024*int(ymax/1024.0 + 0.5)]
            
            print "plotting ", ccd
            
            fig = qaFig.QaFigure(size=figsize)
            fig.fig.subplots_adjust(left=0.15)
            ax = fig.fig.add_subplot(111)

            q = ax.quiver(x, y, dx, dy, color='k', scale=4.0*vLen, angles='xy', pivot='middle',
                          headlength=1, headwidth=1)
            ax.quiverkey(q, 0.9, -0.1, 0.5*vLen, "e=0.5", coordinates='axes',
                         fontproperties={'size':"small"}, labelpos='E')
            
            ax.set_title("PSF ellipticity")
            ax.set_xlabel("x [pixels]") #, position=(0.4, 0.0))
            
            ax.set_ylabel("y [pixels]")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            for tic in ax.get_xticklabels() + ax.get_yticklabels():
                tic.set_size("x-small")

            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            testSet.addFigure(fig, "psfEllip.png",
                              "PSF ellipticity (e=1 shown with length %.0f pix))"%(vLen),
                              areaLabel=areaLabel)
            
            i += 1

