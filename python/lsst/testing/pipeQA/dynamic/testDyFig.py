import sys, os, re

import numpy

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse
from matplotlib import cm
from matplotlib import colors

def plot(data, dataId):

    n = 100
    x = 2.0*3.142*numpy.arange(n)/n
    y = numpy.sin(x)

    fig = figure.Figure()
    canvas = FigCanvas(fig)

    ax = fig.add_subplot(111)
    ax.plot(x, y)

    fig.savefig("testDyFig.png")

    
