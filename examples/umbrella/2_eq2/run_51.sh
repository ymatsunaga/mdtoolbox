#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_51.in \
 -c ../1_eq1/run.rst \
 -o run_51.out \
 -r run_51.rst \
 -x run_51.nc

