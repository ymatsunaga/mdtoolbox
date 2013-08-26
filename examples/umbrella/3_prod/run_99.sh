#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_99.in \
 -c ../2_eq2/run_99.rst \
 -o run_99.out \
 -r run_99.rst \
 -x run_99.nc

