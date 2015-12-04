"""PARALLEL SPEEDUP PLOTTER
Logan Halstrom
MAE 267
Created:  12/03/15

Description:  plot parallel speedup
"""

import sys
sys.path.append('/Users/Logan/lib/python')
from lplot import *
import numpy as np
import matplotlib.pyplot as plt
import os

from move import CaseDir, savedir

def ReadWalltime(casedir):
    """Read CPU wall time for given run
    """
    #OPEN PERFORMANCE FILE
    perffile = 'SolnInfo.dat'
    file = open('{}/{}'.format(casedir, perffile), 'r')
    for i, line in enumerate(file):

        #READ WALLTIME
        if i == 3:
            walltime = float( line.split()[0] )
    return walltime

def PlotResid(ax, iter, res, label, color='blue', marker='.', spacing=100):
    """Plot residual convergence history on log-scale y axis

    """
    ax.semilogy(iter, res, label=label,
                        color=color, linewidth=line,
                        marker=marker, markersize=mark, markevery=spacing)





#INPUTS
from plotResid import markers, colors, savetype

def main(nprocs, nx, N, M, name=''):
    """
    name --> text to include in file name after ResHist
    """


    # READ DATA

    #serial walltime
    case = CaseDir(1, nx, N, M, savedir)
    ts = ReadWalltime(case)

    #Parallel walltimes
    tps = []
    #Parallel speedups
    Sps = []
    #Parallel Computational Efficiencies
    Eps = []
    for nproc in nprocs:
        #File path to load
        case = CaseDir(nproc, nx, N, M, savedir)
        #PARALLEL WALLTIME
        tps.append( ReadWalltime(case) )
        #PARALLEL SPEEDUP
        Sps.append( ts / tps[-1] )
        #PARALLEL COMPUTATIONAL EFFICIENCY
        Eps.append( Sps[-1] / nproc )

    #CALC SPEEDUP AND EFFICIENCY

    for nproc, tp, Sp, Ep in zip(nprocs, tps, Sps, Eps):
        print('{} Processors:'.format(nproc))
        print('walltime: {}, Speedup: {}, Efficiency: {}'.format(tp, Sp, Ep))


    #SPEEDUP PLOT
    title = 'Parallel Speedup\n({}x{} Blocks, {}x{} Mesh)'.format(N, M, nx, nx)
    _, ax = PlotStart(title, 'NPROCS', '$S_P$', horzy='horizontal', figsize='tex')


    # ax.axis('equal')
    ax.grid(True)
    ax.set_xlim( [0, max(nprocs)] )
    ax.set_ylim( [0, max(nprocs)] )
    allprocs = np.linspace(0, max(nprocs), max(nprocs)+1 )
    plt.xticks( allprocs )
    plt.yticks( allprocs )

    # for c, clr, mkr in zip(cases, colors, markers):

    #PLOT IDEAL SPEEDUP
    ax.plot(allprocs, allprocs, label='Ideal',
                    color='black', linewidth=line,
                    marker='x', markersize=mark*1.5)
    #PLOT ACTUAL SPEEDUP
    ax.plot(list([1]) + list(nprocs), list([1]) + list(Sps), label='Parallel',
                    color='blue', linewidth=line,
                    marker='^', markersize=mark*1.5)

    PlotLegend(ax, loc='best')

    text = 'Speedup' + name
    savename = '{}/{}{}'.format(savedir, text, savetype)
    SavePlot(savename)
    plt.show()

    #EFFICIENCY PLOT
    title = 'Parallel Computational Efficiency\n({}x{} Blocks, {}x{} Mesh)'.format(N, M, nx, nx)
    _, ax = PlotStart(title, 'NPROCS', '$E_P$', horzy='horizontal', figsize='tex')

    # ax.axis('equal')
    ax.grid(True)
    ax.set_xlim( [0, max(nprocs)] )
    allprocs = np.linspace(0, max(nprocs), max(nprocs)+1 )
    plt.xticks( allprocs )


    # #PLOT IDEAL EFFICIENCY
    ax.plot(allprocs, np.ones(len(allprocs)), label='Ideal',
                    color='black', linewidth=line)
    #PLOT ACTUAL EFFICIENCY
    ax.plot(nprocs, Eps, label='Parallel',
                    color='blue', linewidth=line,
                    marker='^', markersize=mark*1.5)

    PlotLegend(ax, loc='best')

    text = 'Efficiency' + name
    savename = '{}/{}{}'.format(savedir, text, savetype)
    SavePlot(savename)
    plt.show()



if __name__ == "__main__":

    # NPROCS = [4, 8]
    # NX = 501
    # N = 10
    # M = 10

    NPROCS = [2, 4, 6, 8]
    NX = 501
    N = 10
    M = 10

    main(NPROCS, NX, N, M)
