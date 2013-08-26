#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run.in \
 -c ../0_init/ala.crd \
 -o run.out \
 -r run.rst \
 -x run.nc

