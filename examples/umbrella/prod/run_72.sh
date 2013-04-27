#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_72.in \
 -c ../eq2/run_72.rst \
 -o run_72.out \
 -r run_72.rst \
 -x run_72.nc

