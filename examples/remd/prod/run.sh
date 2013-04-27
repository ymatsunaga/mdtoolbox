#!/bin/bash

NPROC=16
mpirun -np $NPROC sander.MPI -ng 8 -groupfile groupfile

