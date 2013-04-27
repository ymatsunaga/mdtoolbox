#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p parm \
 -i run.in \
 -c ../init/ala.crd \
 -o run.out \
 -r run.rst \
 -x run.nc

