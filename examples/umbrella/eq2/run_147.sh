#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_147.in \
 -c ../eq1/run.rst \
 -o run_147.out \
 -r run_147.rst \
 -x run_147.nc

