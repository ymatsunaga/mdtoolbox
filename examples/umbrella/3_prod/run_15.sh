#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_15.in \
 -c ../2_eq2/run_15.rst \
 -o run_15.out \
 -r run_15.rst \
 -x run_15.nc

