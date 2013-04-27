#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_84.in \
 -c ../eq2/run_84.rst \
 -o run_84.out \
 -r run_84.rst \
 -x run_84.nc

