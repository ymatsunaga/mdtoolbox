#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_36.in \
 -c ../eq2/run_36.rst \
 -o run_36.out \
 -r run_36.rst \
 -x run_36.nc

