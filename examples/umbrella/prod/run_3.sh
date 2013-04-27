#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_3.in \
 -c ../eq2/run_3.rst \
 -o run_3.out \
 -r run_3.rst \
 -x run_3.nc

