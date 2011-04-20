import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils

import numpy

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse



def binDistrib(x, y, dy, binSizeX = 0.5, minPts = 2):
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


def plotSparseContour(sp, x, y, binSizeX, binSizeY, minCont = 500, nCont = 7):
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

    #if cdata.max() < minCont:
    #    minCont = 1

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



def cameraToRectangles(camera):
    rectangles = {}
    centers = {}
    raftBoundaries = []
    ccdBoundaries = []
    for r in camera:
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
	    orient  = ccd.getOrientation()
	    nQuart  = ccd.getOrientation().getNQuarter()
	    yaw     = orient.getYaw()

	    cbbox   = ccd.getAllPixels(True)
	    cwidth  = cbbox.getX1() - cbbox.getX0()
	    cheight = cbbox.getY1() - cbbox.getY0()
	    if abs(yaw - numpy.pi/2.0) < 1.0e-3:  # nQuart == 1 or nQuart == 3:
		ctmp = cwidth
		cwidth = cheight
		cheight = ctmp

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
	    rectangles[label] = crectangle
	    centers[label] = (rxc+cxc, ryc+cyc)

	    ccdBoundaries.append(((cx0, cx0), (cy0, cy1)))
	    ccdBoundaries.append(((cx0, cx1), (cy0, cy0)))
	    ccdBoundaries.append(((cx0, cx1), (cy1, cy1)))
	    ccdBoundaries.append(((cx1, cx1), (cy0, cy1)))
	    

	raftBoundaries.append(((xmin, xmin), (ymin, ymax)))
	raftBoundaries.append(((xmin, xmax), (ymin, ymin)))
	raftBoundaries.append(((xmin, xmax), (ymax, ymax)))
	raftBoundaries.append(((xmax, xmax), (ymin, ymax)))

    return centers, rectangles, raftBoundaries, ccdBoundaries

