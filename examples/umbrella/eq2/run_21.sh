#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_21.in \
 -c ../eq1/run.rst \
 -o run_21.out \
 -r run_21.rst \
 -x run_21.nc

