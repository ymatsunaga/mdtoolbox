#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_171.in \
 -c ../2_eq2/run_171.rst \
 -o run_171.out \
 -r run_171.rst \
 -x run_171.nc

