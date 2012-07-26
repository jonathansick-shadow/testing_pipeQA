#!/usr/bin/env python

import sys
import numpy
import matplotlib

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import colors
import matplotlib.font_manager as fm
from  matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle


import QaPlotUtils as qaPlotUtil

def plot(data):
    pass

if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, n = qaPlotUtil.unshelveGlob(filename)
    if n > 1:
        data['gridVectors'] = True
    fig = plot(data)
    fig.savefig(filename)
