#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_60.in \
 -c ../eq1/run.rst \
 -o run_60.out \
 -r run_60.rst \
 -x run_60.nc

