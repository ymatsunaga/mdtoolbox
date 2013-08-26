#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_45.in \
 -c ../1_eq1/run.rst \
 -o run_45.out \
 -r run_45.rst \
 -x run_45.nc

