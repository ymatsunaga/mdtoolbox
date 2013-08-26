#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_147.in \
 -c ../2_eq2/run_147.rst \
 -o run_147.out \
 -r run_147.rst \
 -x run_147.nc

