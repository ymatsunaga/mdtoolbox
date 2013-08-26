#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_117.in \
 -c ../2_eq2/run_117.rst \
 -o run_117.out \
 -r run_117.rst \
 -x run_117.nc

