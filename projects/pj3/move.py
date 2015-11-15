"""FILE MOVER
Logan Halstrom
MAE 267
Created:  11/14/15

Description:  Move results files from top directory to results directory
for specific case.
(outputing files directoy to case directory in fortran was too nit-picky)
"""

import sys
sys.path.append('/Users/Logan/lib/python')
import lutil as lu
import numpy as np
import os





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
    return nx, M, N



def main():

    nx, M, N = GetCaseInfo()
    print(nx, M, N)


if __name__ == "__main__":
    main()
