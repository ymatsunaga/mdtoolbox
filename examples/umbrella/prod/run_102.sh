#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_102.in \
 -c ../eq2/run_102.rst \
 -o run_102.out \
 -r run_102.rst \
 -x run_102.nc

