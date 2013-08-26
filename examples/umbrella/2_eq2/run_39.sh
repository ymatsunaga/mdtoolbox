#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_39.in \
 -c ../1_eq1/run.rst \
 -o run_39.out \
 -r run_39.rst \
 -x run_39.nc

