#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_90.in \
 -c ../1_eq1/run.rst \
 -o run_90.out \
 -r run_90.rst \
 -x run_90.nc

