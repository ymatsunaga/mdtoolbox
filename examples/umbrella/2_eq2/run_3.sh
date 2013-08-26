#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_3.in \
 -c ../1_eq1/run.rst \
 -o run_3.out \
 -r run_3.rst \
 -x run_3.nc

