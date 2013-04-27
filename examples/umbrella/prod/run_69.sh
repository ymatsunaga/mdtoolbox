#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_69.in \
 -c ../eq2/run_69.rst \
 -o run_69.out \
 -r run_69.rst \
 -x run_69.nc

