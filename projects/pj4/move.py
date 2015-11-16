#! /usr/local/bin/python3
#! /usr/bin/python3
"""FILE MOVER
Logan Halstrom
MAE 267
Created:  11/14/15

Description:  Move results files from top directory to results directory
for specific case.
(outputing files directoy to case directory in fortran was too nit-picky)

NOTE:  Update the sha-bang (top line) depending on running
linux (usr/bin/python3) or MacOS (usr/local/bin/python3)
"""

import os
import subprocess
import glob

def cmd(command):
    """Execute a shell command.
    Execute multiple commands by separating with semicolon+space: '; '
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    #print proc_stdout
    proc_stdout = process.communicate()[0].strip()
    return proc_stdout

def MakeOutputDir(savedir):
    """make results output directory if it does not already exist.
    instring --> directory path from script containing folder
    """
    #split individual directories
    splitstring = savedir.split('/')
    prestring = ''
    for string in splitstring:
        prestring += string + '/'
        try:
            os.mkdir(prestring)
        except Exception:
            pass

def GetCaseInfo(configfile='config.in'):
    """Read input file to determine info about specific case."""
    #OPEN CONFIG FILE
    config = open(configfile, 'r')
    for i, line in enumerate(config):

        #READ NUMBER OF GRIDPOINTS
        if i == 3:
            nx = int(line)
        #READ NUMBER OF BLOCKS IN J-DIRECTION
        if i == 5:
            M = int(line)
        #READ NUMBER OF BLOCKS IN I-DIRECTION
        if i == 7:
            N = int(line)
        #GET NUMBER OF PROCESSORS FROM PROC CONFIG FILES
        #list of config file names for each proc
        configs = glob.glob('p*.config')
        nproc = len(configs)
    return nx, M, N, nproc, configs

def main():

    #GET CASE INFO
    nx, M, N, nproc, configs = GetCaseInfo()
    #MAKE CASE DIRECTORY (i.e. np4_nx101_n5_m4)
    outdir = 'Results/np{}_nx{}_n{}_m{}'.format(nproc, nx, N, M)
    MakeOutputDir(outdir)
    #MOVE FILES
    filelist = [
                    # 'blockconfig.dat',
                    'a.out', 'SolnInfo.dat', 'res_hist.dat',
                    'grid.xyz', 'grid_form.xyz', 'T.dat', 'T_form.dat'
               ]
    #add config files to move
    filelist += configs
    for f in filelist:
        cmd( 'mv {} {}/.'.format(f, outdir) )

    print( 'Case files moved to {}'.format(outdir) )

if __name__ == "__main__":
    main()
