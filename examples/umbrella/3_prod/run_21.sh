#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_21.in \
 -c ../2_eq2/run_21.rst \
 -o run_21.out \
 -r run_21.rst \
 -x run_21.nc

