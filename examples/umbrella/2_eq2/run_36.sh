#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_36.in \
 -c ../1_eq1/run.rst \
 -o run_36.out \
 -r run_36.rst \
 -x run_36.nc

