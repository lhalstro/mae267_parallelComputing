"""RESIDUAL PLOTTER
Logan Halstrom
MAE 267
Created:  11/05/15

Description:  plot residual history
"""

import sys
sys.path.append('/Users/Logan/lib/python')
from lplot import *
import numpy as np
import matplotlib.pyplot as plt
import os


markers = ['x', 'o', '.']
colors = ['blue', 'green', 'red']


def main():


    # READ DATA
    filename = 'res_hist.dat'
    folders = ['finalSingle', 'finalMulti']
    savedir = 'Results'
    savetype = '.png'

    sing = {}
    path = '{}/{}/{}'.format(savedir, folders[0], filename)
    sing['iter'], sing['res'] = np.loadtxt(path, skiprows=1, unpack=True)

    mult = {}
    path = '{}/{}/{}'.format(savedir, folders[1], filename)
    mult['iter'], mult['res'] = np.loadtxt(path, skiprows=1, unpack=True)


    title = 'Residual History of Single and Multi-Block Solvers'
    _, ax = PlotStart(title, 'ITER', 'RESID', horzy='vertical', figsize='tex')

    # ax.axis('equal')
    ax.grid(True)
    xmin = 1500
    ax.set_ylim([0, mult['res'][xmin] * 1.25])
    ax.set_xlim([50, 80000])


    spacing = 100
    ax.plot(sing['iter'], sing['res'], color='red', label='Single Block',
                    linewidth=line, marker='x', markersize=mark, markevery=spacing)
    ax.plot(mult['iter'], mult['res'], color='blue', label='Multi Block',
                    linewidth=line, marker='.', markersize=mark, markevery=spacing)


    PlotLegend(ax, loc='best')

    savename = '{}/ResHist{}'.format(savedir, savetype)
    SavePlot(savename)
    plt.show()





if __name__ == "__main__":
    main()
