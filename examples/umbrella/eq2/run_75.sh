#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_75.in \
 -c ../eq1/run.rst \
 -o run_75.out \
 -r run_75.rst \
 -x run_75.nc

