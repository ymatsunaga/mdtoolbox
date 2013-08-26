#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run_144.in \
 -c ../1_eq1/run.rst \
 -o run_144.out \
 -r run_144.rst \
 -x run_144.nc

