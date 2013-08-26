#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_174.in \
 -c ../2_eq2/run_174.rst \
 -o run_174.out \
 -r run_174.rst \
 -x run_174.nc

