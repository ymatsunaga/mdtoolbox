#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_12.in \
 -c ../1_eq1/run.rst \
 -o run_12.out \
 -r run_12.rst \
 -x run_12.nc

