#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_66.in \
 -c ../eq1/run.rst \
 -o run_66.out \
 -r run_66.rst \
 -x run_66.nc

