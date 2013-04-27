#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_111.in \
 -c ../eq1/run.rst \
 -o run_111.out \
 -r run_111.rst \
 -x run_111.nc

