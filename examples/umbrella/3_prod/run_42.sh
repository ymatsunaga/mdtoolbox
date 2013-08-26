#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_42.in \
 -c ../2_eq2/run_42.rst \
 -o run_42.out \
 -r run_42.rst \
 -x run_42.nc

