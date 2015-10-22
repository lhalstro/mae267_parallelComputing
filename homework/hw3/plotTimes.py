"""PLOT CPU TIMES
Logan Halstrom
21 OCTOBER 2015

DESCRIPTION:
"""

import sys
sys.path.append('/Users/Logan/lib/python')
from lplot import *
from lutil import TexTable
import os
import matplotlib.pyplot as plt
import numpy as np

def ReadData(path):
    """Read cpu times and results for given results file"""
    ifile = open(path, 'r')
    for i, line in enumerate(ifile):
        if i == 4:
            #READ NUMBER OF PROCESSORS
            nproc = int(line)
        if i == 6:
            pi = float(line)
            #READ PI
        if i == 8:
            #READ CPU WALLTIME
            walltime = float(line)
        if i == 10:
            #READ WALLTIME PER PROCESSOR
            walltimeper = float(line)
    return nproc, pi, walltime, walltimeper

def DictData(inttype, N, ncpus):
    """Read data and store in dictionary"""
    #INITIALIZE DICT OF EMPTY LISTS WITH GIVEN KEYS
    data = {k: [] for k in keys}
    for np in ncpus:
        filename = '{}{}_i{}_np{}.dat'.format(savedir, inttype, N, np)
        out = ReadData(filename)
        #MAKE LISTS FOR ALL RESULTS VARIABLES
        for o, key in zip(out, keys):
            data[key].append(o)
    return data

savedir = 'Results/'
savetype = '.pdf'
#DICTIONARY KEYS (# of procs, pi result, walltime, walltime per cpu)
keys = ['np', 'pi', 't', 'tper']
mark=10

def main():

    #TYPE OF INTEGRATION
        #'simp' or 'trap'
    inttype = 'simp'
    #NUMBER OF SUB-INTERVALS
        #integer
    N = 256
    #LIST OF NUMBER OF CPUS TO PLOT
    ncpus = [1, 2, 4, 8]

    #READ DATA
        #dictionary in a dictionary
    data = {}
    #simpson's rule results
    data['simp'] = DictData('simp', N, ncpus)
    #simpson's with lots of procs
    data['simplong'] = DictData('simp', N, [1, 2, 4, 8, 16, 32, 64, 128])
    #trapezoid rule results
    data['trap'] = DictData('trap', N, ncpus)
    for d in data.values():
        #error in calculating pi
        d['err'] = np.subtract(d['pi'], np.pi)

    #TABULATE DATA
    A = np.zeros((len(data['simp']['np']), 4))
    for i, n in enumerate(data['simp']['np']):
        for j, k in enumerate(['t', 'tper', 'pi', 'err']):
            A[i,j] = data['simp'][k][i]
    print(A)
    rows = [str(s) for s in data['simp']['np']]
    columns = ['# of\\\\CPUs', 'Wall Time (s)', 'Wall Time\\\\Per Processor',
                    '$\pi$ Result', '$\pi$ Error']
    TexTable('{}results_simp.tex'.format(savedir), A, rows, columns, 4,
                'results', 'Timing and Solver Results')


    #PLOT WALL TIME VS # OF CPUS
    for key in ['simp', 'simplong']:
        title = 'Wall Times for Solving $\pi$ with Simpson''s Rule'
        xlbl = '$n_{CPU}$'
        ylbl = 'Total Wall Time (s)'
        _, ax = PlotStart(title, xlbl, ylbl)
        ax.plot(data[key]['np'], data[key]['t'], color='blue',
                    linestyle='-', linewidth=line, marker='.', markersize=mark)
        ax.grid()
        SavePlot('{}walltimes_{}{}'.format(savedir, key, savetype))

    #PLOT WALL TIME PER CPU VS # OF CPUS
    title = 'Wall Time Per Processor for Solving $\pi$ with Simpson''s Rule'
    xlbl = '$n_{CPU}$'
    ylbl = 'Wall Time Per CPU(s)'
    _, ax = PlotStart(title, xlbl, ylbl)
    ax.plot(data['simplong']['np'], data['simplong']['tper'], color='blue',
                linestyle='-', linewidth=line, marker='.', markersize=mark)
    ax.grid()
    SavePlot('{}walltimesPer_{}{}'.format(savedir, 'simp', savetype))

    #PLOT WALL TIME VS # OF CPUS BOTH METHODS
    title = 'Wall Times for Solving $\pi$'
    xlbl = '$n_{CPU}$'
    ylbl = 'Total Wall Time (s)'
    _, ax = PlotStart(title, xlbl, ylbl)
    ax.plot(data['simp']['np'], data['simp']['t'], color='blue',
                linestyle='-', linewidth=line, marker='.', markersize=mark,
                label='Simpson''s')
    ax.plot(data['trap']['np'], data['trap']['t'], color='green',
                linestyle='-', linewidth=line, marker='x', markersize=mark,
                label='Trapezoid')
    ax.grid()
    leg = ax.legend(loc='best', fancybox=True, framealpha=0.5)
    SavePlot('{}walltimes_both{}'.format(savedir, savetype))

    #PLOT PI RESULTS
    title = 'Error In $\pi$ (Simpson''s Rule)'
    xlbl = '$n_{CPU}$'
    ylbl = 'Error (\pi_{Solver} - $\pi$)'
    _, ax = PlotStart(title, xlbl, ylbl)
    # ax.plot(data['np'], data['t'], color=color,
    #             linestyle='-', linewidth=line, marker=marker, markersize=mark,
    #             label='{} $\delta^*$'.format(self.info['name']),)
    # ax.plot([ncpus[0], ncpus[-1]], [np.pi, np.pi], color='black',
    #             linewidth=line, label='$\pi$')
    ax.plot(data['simplong']['np'], data['simplong']['err'], color='blue',
                linestyle='-', linewidth=line, marker='.', markersize=mark,
                label='Simpson''s')
    # ax.plot(data['trap']['np'], data['trap']['pi'], color='green',
    #             linestyle='-', linewidth=line, marker='x', markersize=mark,
    #             label='Trapezoid')
    ax.grid()
    # plt.ylim([min(data['simp']['pi']), max(data['simp']['pi'])])
    leg = ax.legend(loc='best', fancybox=True, framealpha=0.5)
    SavePlot('{}piErr_simp{}'.format(savedir, savetype))



if __name__ == "__main__":




    main()
