#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_57.in \
 -c ../eq1/run.rst \
 -o run_57.out \
 -r run_57.rst \
 -x run_57.nc

