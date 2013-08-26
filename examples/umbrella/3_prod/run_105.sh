#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_105.in \
 -c ../2_eq2/run_105.rst \
 -o run_105.out \
 -r run_105.rst \
 -x run_105.nc

