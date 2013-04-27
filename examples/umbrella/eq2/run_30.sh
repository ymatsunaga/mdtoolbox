#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_30.in \
 -c ../eq1/run.rst \
 -o run_30.out \
 -r run_30.rst \
 -x run_30.nc

