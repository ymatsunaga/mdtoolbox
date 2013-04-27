#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_84.in \
 -c ../eq1/run.rst \
 -o run_84.out \
 -r run_84.rst \
 -x run_84.nc

