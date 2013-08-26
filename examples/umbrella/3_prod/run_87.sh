#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_87.in \
 -c ../2_eq2/run_87.rst \
 -o run_87.out \
 -r run_87.rst \
 -x run_87.nc

