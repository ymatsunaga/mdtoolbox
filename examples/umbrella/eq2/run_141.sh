#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_141.in \
 -c ../eq1/run.rst \
 -o run_141.out \
 -r run_141.rst \
 -x run_141.nc

