#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_138.in \
 -c ../2_eq2/run_138.rst \
 -o run_138.out \
 -r run_138.rst \
 -x run_138.nc

