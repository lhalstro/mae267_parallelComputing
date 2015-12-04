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

from move import CaseDir, savedir


def PlotResid(ax, iter, res, label, color='blue', marker='.', spacing=100):
    """Plot residual convergence history on log-scale y axis

    """
    ax.semilogy(iter, res, label=label,
                        color=color, linewidth=line*.75,
                        marker=marker, markersize=mark, markevery=spacing)



#INPUTS


markers = ['o',   '^',     'x',   '.'    'sq']
colors = ['blue', 'green', 'red', 'cyan', 'magenta']
#point spacing
spacing = 2000
# savedir = 'Results'
savetype = '.pdf'

def main(nprocs, nx, Ns, Ms, name=''):
    """
    name --> text to include in file name after ResHist
    """


    # READ DATA
    filename = 'res_hist.dat'

    #ADD SERIAL, MULTI-BLOCK CASE
    nprocs = list([1 ]) + list(nprocs)
    Ns     = list([10]) + list(Ns)
    Ms     = list([10]) + list(Ms)

    #ADD SERIAL, SINGLE-BLOCK CASE
    nprocs = list([1]) + list(nprocs)
    Ns     = list([1]) + list(Ns)
    Ms     = list([1]) + list(Ms)

    #LOCATE CASE DIRECTORIES
    folders = []
    cases = []
    for nproc, N, M in zip(nprocs, Ns, Ms):
        #File path to load
        case = CaseDir(nproc, nx, N, M, savedir)
        path = '{}/{}'.format(case, filename)
        #Add empty dict to cases
        cases.append( {} )
        #Populate dict with data
        cases[-1]['itr'], cases[-1]['res'] = np.loadtxt(path, skiprows=1, unpack=True)
        #make name for case
        cases[-1]['label'] = 'NP={}, NxM=({}x{})'.format(nproc, N, M)

    #Special labels for simple cases
    for c, l in zip(cases[:2],['Serial, Single-Block', 'Serial, Multi (10x10)'] ):
        c['label'] = l


    #RESIDUAL PLOT

    title = 'Residual Convergence History\n({}x{} Mesh)'.format(nx, nx)
    _, ax = PlotStart(title, 'ITER', 'RESID', horzy='vertical', figsize='tex')

    # ax.axis('equal')
    ax.set_yscale("log", nonposy='clip')
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.grid(True)
    # xmin = 1500
    # ax.set_ylim([0, mult['res'][xmin] * 1.25])
    # ax.set_xlim([50, 80000])

    for c, clr, mkr in zip(cases, colors, markers):
        PlotResid(ax, c['itr'], c['res'], c['label'], clr, mkr, spacing)

    PlotLegend(ax, loc='best')

    text = 'ResHist' + name
    savename = '{}/{}{}'.format(savedir, text, savetype)
    SavePlot(savename)
    plt.show()





if __name__ == "__main__":

    NPROCS = [4]
    NX = 501
    NS = [10] * len(NPROCS)
    MS = [10] * len(NPROCS)

    main(NPROCS, NX, NS, MS)
