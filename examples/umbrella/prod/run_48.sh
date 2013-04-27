#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_48.in \
 -c ../eq2/run_48.rst \
 -o run_48.out \
 -r run_48.rst \
 -x run_48.nc

