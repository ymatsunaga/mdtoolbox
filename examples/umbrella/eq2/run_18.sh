#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_18.in \
 -c ../eq1/run.rst \
 -o run_18.out \
 -r run_18.rst \
 -x run_18.nc

