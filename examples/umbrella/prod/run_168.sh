#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_168.in \
 -c ../eq2/run_168.rst \
 -o run_168.out \
 -r run_168.rst \
 -x run_168.nc

