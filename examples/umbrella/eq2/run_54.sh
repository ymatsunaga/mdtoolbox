#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_54.in \
 -c ../eq1/run.rst \
 -o run_54.out \
 -r run_54.rst \
 -x run_54.nc

