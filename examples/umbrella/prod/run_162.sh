#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_162.in \
 -c ../eq2/run_162.rst \
 -o run_162.out \
 -r run_162.rst \
 -x run_162.nc

