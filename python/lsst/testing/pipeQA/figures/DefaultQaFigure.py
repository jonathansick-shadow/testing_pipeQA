import os
import sys
import re
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

###############################################
#
# This is a minimal QaFigure class, the real one lives in displayQA
#
# We should be able to run the pipeQA tests without the display setup,
# and we should be able to use displayQa with no dependencies on the pipe.
# However, to make that possible, I'd need a third package containing TestCode and QaFigure
# as both products use those classes.  The cheap solution is to provide watered-down versions
# in the pipeQA, and override them with displayQA versions, if those are setup.
#
###############################################


class QaFigure(object):
    """A wrapper for a matplotlib figure, and a Baseclass for more complicated QA figures. """

    def __init__(self, size=(4.0, 4.0), dpi=100):  # (512, 512), DPI=100):
        """
        @param size  Figure size in inches
        @param dpi   Dots per inch to use.
        """

        self.fig = figure.Figure(figsize=size)
        self.fig.set_dpi(dpi)
        self.canvas = FigCanvas(self.fig)
        self.map = {}
        #self.fig.set_size_inches(size[0] / DPI, size[1] / DPI)
        self.mapAreas = []
        self.mapTransformed = True

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
        """Save figure."""
        self.fig.savefig(path, dpi=self.fig.get_dpi(), **kwargs)

    def savemap(self, path):
        """Save internal map area data to .map file. """
        pass

    def getTransformedMap(self):
        """Take plot coordinates for map areas and convert to figure coordinates."""
        return []

    def addMapArea(self, label, area, info, axes=None):
        """Add a map area to this figure.

        @param label    an areaLabel to associate the area to other QaFigures or Tests.
        @param area     [x0, y0, x1, y1]  llc, urc corners of the area in plot coordinates.
        @param info     string to show no mouseover (no whitespace)
        @param axes     axes to use for map transformation
        """
        pass

    def getMapInfo(self):
        return self.mapAreas


