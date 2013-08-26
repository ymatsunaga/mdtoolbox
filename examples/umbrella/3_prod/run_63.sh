#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_63.in \
 -c ../2_eq2/run_63.rst \
 -o run_63.out \
 -r run_63.rst \
 -x run_63.nc

