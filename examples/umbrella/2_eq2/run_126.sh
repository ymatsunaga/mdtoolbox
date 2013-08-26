#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_126.in \
 -c ../1_eq1/run.rst \
 -o run_126.out \
 -r run_126.rst \
 -x run_126.nc

