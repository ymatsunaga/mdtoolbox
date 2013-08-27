#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_150.in \
 -c ../1_eq1/run.rst \
 -o run_150.out \
 -r run_150.rst \
 -x run_150.nc
