#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_165.in \
 -c ../eq1/run.rst \
 -o run_165.out \
 -r run_165.rst \
 -x run_165.nc

