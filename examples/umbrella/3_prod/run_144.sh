#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_144.in \
 -c ../2_eq2/run_144.rst \
 -o run_144.out \
 -r run_144.rst \
 -x run_144.nc

