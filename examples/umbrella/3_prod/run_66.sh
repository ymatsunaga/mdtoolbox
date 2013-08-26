#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_66.in \
 -c ../2_eq2/run_66.rst \
 -o run_66.out \
 -r run_66.rst \
 -x run_66.nc

