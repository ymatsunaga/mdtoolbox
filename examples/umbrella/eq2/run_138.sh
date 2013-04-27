#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_138.in \
 -c ../eq1/run.rst \
 -o run_138.out \
 -r run_138.rst \
 -x run_138.nc

