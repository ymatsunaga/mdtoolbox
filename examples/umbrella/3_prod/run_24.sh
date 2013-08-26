#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_24.in \
 -c ../2_eq2/run_24.rst \
 -o run_24.out \
 -r run_24.rst \
 -x run_24.nc

