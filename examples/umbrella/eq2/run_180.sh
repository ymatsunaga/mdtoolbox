#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_180.in \
 -c ../eq1/run.rst \
 -o run_180.out \
 -r run_180.rst \
 -x run_180.nc

