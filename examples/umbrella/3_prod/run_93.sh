#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_93.in \
 -c ../2_eq2/run_93.rst \
 -o run_93.out \
 -r run_93.rst \
 -x run_93.nc

