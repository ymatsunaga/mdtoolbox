#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_6.in \
 -c ../1_eq1/run.rst \
 -o run_6.out \
 -r run_6.rst \
 -x run_6.nc

