#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_78.in \
 -c ../1_eq1/run.rst \
 -o run_78.out \
 -r run_78.rst \
 -x run_78.nc

