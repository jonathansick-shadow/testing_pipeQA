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


class PsfEllipticityQaAnalysis(qaAna.QaAnalysis):

    def __init__(self):
        qaAna.QaAnalysis.__init__(self)


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

    def test(self, data, dataId):
        
        # get data
        self.ssDict        = data.getSourceSetBySensor(dataId)
        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)

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

                a2 = 0.5*(ixx+iyy) + numpy.sqrt(0.25*(ixx-iyy)**2 + ixy**2)
                b2 = 0.5*(ixx+iyy) - numpy.sqrt(0.25*(ixx-iyy)**2 + ixy**2)
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

        self.limits = [0.0, 0.3]
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
                ellipMed = 99.0
                thetaMed = 0.0
                n = 0

            # add a test for acceptible psf ellipticity
            self.ellipMedians.set(raft, ccd, ellipMed)
            areaLabel = re.sub("\s+", "_", ccd)
            label = "median psf ellipticity "+areaLabel
            comment = "median psf ellipticity (nstar=%d)" % (n)
            testSet.addTest( testCode.Test(label, ellipMed, self.limits, comment, areaLabel=areaLabel) )

            # stash the angles.  We'll use them to make figures in plot()
            self.thetaMedians.set(raft, ccd, thetaMed)
            

    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)

        vLen = 1000.0  # for e=1.0

        figFmt = "png"

        # fpa figure
        ellipFig = qaFig.VectorFpaQaFigure(data.cameraInfo.camera)
        for raft, ccdDict in ellipFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.ellipMedians.get(raft, ccd) is None:
                    ellipFig.data[raft][ccd] = [self.thetaMedians.get(raft, ccd),
                                                10*vLen*self.ellipMedians.get(raft, ccd),
                                                self.ellipMedians.get(raft, ccd)]
                    ellipFig.map[raft][ccd] = "ell/theta=%.3f/%.0f" % (self.ellipMedians.get(raft, ccd),
                                                                       (180/numpy.pi)*self.thetaMedians.get(raft, ccd))
                
        ellipFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=self.limits,
                            title="Median PSF Ellipticity", failLimits=self.limits)
        testSet.addFigure(ellipFig, "medPsfEllip."+figFmt, "Median PSF Ellipticity",
                          navMap=True)

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

            xy1 = zip(x, y)
            xy2 = zip(x+dx, y+dy)
            lines = zip(xy1, xy2)
            p = LineCollection(lines, colors=black*len(lines))
            ax.add_collection(p)
            
            ax.set_title("PSF ellipticity")
            ax.set_xlabel("x [pixels]")
            ax.set_ylabel("y [pixels]")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            for tic in ax.get_xticklabels() + ax.get_yticklabels():
                tic.set_size("x-small")

            areaLabel = re.sub("\s+", "_", ccd)
            testSet.addFigure(fig, "psfEllip.png",
                              "PSF ellipticity (e=1 shown with length %.0f pix))"%(vLen),
                              areaLabel=areaLabel)
            
            i += 1

