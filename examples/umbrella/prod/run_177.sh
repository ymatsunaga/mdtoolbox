#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_177.in \
 -c ../eq2/run_177.rst \
 -o run_177.out \
 -r run_177.rst \
 -x run_177.nc

