#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_114.in \
 -c ../eq1/run.rst \
 -o run_114.out \
 -r run_114.rst \
 -x run_114.nc

