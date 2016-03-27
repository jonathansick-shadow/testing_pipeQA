#!/usr/bin/env python

import sys
import numpy as num
import matplotlib

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import colors
import matplotlib.font_manager as fm
from matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle


import QaPlotUtils as qaPlotUtil


def plot(data):

    dmags = data['dmags']
    radii = data['radii']
    ids = data['ids']
    med, std = data['offsetStats']
    magTypes = data['magTypes']
    summary = data['summary']

    ymin, ymax = -0.5, 0.5
    if len(dmags) > 0:
        ymin = num.max([dmags.min(), -0.5])
        ymax = num.min([dmags.max(), 0.5])
    ylim = [ymin, ymax]

    if len(dmags) == 0:
        dmags = num.array([0.0])
        radii = num.array([0.0])
        ids = num.array([0])

    figsize = (4.0, 4.0)
    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)

    if not summary:

        sp1 = fig.add_subplot(111)
        sp1.plot(radii, dmags, 'ro', ms=2.0)
        sp1.set_ylim(ylim)

        ddmag = 0.001
        drad = 0.01 * (max(radii) - min(radii))
        if False:
            for i in range(len(dmags)):
                info = "nolink:sourceId=%s" % (ids[i])
                area = (radii[i]-drad, dmags[i]-ddmag, radii[i]+drad, dmags[i]+ddmag)
                fig.addMapArea("no_label_info", area, info, axes=sp1)

        sp1.axhline(y=0, c = 'k', linestyle = ':', alpha = 0.25)
        sp1.axhline(y=med, c = 'b', linestyle = '-')
        sp1.axhspan(ymin = med-std, ymax = med+std, fc = 'b', alpha = 0.15)
        sp1x2 = sp1.twinx()
        ylab = sp1x2.set_ylabel('Delta magnitude (%s-%s)' % tuple(magTypes), fontsize=10)
        ylab.set_rotation(-90)
        sp1.set_xlabel('Dist from focal plane center (pixels)', fontsize=10)
        qaPlotUtil.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize = 8)
        qaPlotUtil.qaSetp(sp1x2.get_xticklabels()+sp1x2.get_yticklabels(), visible=False)

    else:
        ymin = num.max([dmags.min(), -0.5])
        ymax = num.min([dmags.max(), 0.5])
        ylim = [ymin, ymax]
        if ymin == ymax:
            ylim = [ymin - 0.1, ymax + 0.1]

        sp1 = fig.add_subplot(111)
        sp1.plot(radii, dmags, 'ro', ms=2, alpha = 0.5)

        sp1.set_ylim(ylim)

        sp1.axhline(y=0, c = 'k', linestyle = ':', alpha = 0.25)
        sp1x2 = sp1.twinx()
        ylab = sp1x2.set_ylabel('Delta magnitude (%s-%s)' % tuple(magTypes), fontsize=10)
        ylab.set_rotation(-90)
        sp1.set_xlabel('Dist from focal plane center (pixels)', fontsize=10)
        qaPlotUtil.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize = 8)
        qaPlotUtil.qaSetp(sp1x2.get_xticklabels()+sp1x2.get_yticklabels(), visible=False)

    return fig


if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename)
    if isSummary:
        data['summary'] = True
    fig = plot(data)
    fig.savefig(filename)
