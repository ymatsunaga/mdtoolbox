#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_132.in \
 -c ../eq1/run.rst \
 -o run_132.out \
 -r run_132.rst \
 -x run_132.nc

