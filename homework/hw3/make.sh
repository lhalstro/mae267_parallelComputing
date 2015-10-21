#!/bin/bash
echo "Compiling calcpip..."
mpif90 -o calcpip -O3 simp.f calcpip.f
#mpif90 -o calcpip -O3 trap.f calcpip.f

