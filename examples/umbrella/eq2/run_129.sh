#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_129.in \
 -c ../eq1/run.rst \
 -o run_129.out \
 -r run_129.rst \
 -x run_129.nc

