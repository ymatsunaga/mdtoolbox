#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_96.in \
 -c ../1_eq1/run.rst \
 -o run_96.out \
 -r run_96.rst \
 -x run_96.nc

