#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_78.in \
 -c ../2_eq2/run_78.rst \
 -o run_78.out \
 -r run_78.rst \
 -x run_78.nc

