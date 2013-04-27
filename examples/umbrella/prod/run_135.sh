#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_135.in \
 -c ../eq2/run_135.rst \
 -o run_135.out \
 -r run_135.rst \
 -x run_135.nc

