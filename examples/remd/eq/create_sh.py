#!/usr/bin/python
# coding: utf-8

import os

############# define functions
def print_lines(f, value):
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("NPROC=2\n")
    f.write("mpirun -np $NPROC sander.MPI -O \\\n")
    f.write(" -p parm \\\n")
    f.write(" -i run.in.%03d \\\n" % (value))
    f.write(" -c ../init/ala.crd \\\n")
    f.write(" -o run.out.%03d \\\n" % (value))
    f.write(" -r run.rst.%03d \\\n" % (value))
    f.write(" -x run.nc.%03d\n" % (value))
    f.write("\n")

############# main
temperature = [300.00, 331.00, 364.48, 401.69, 443.05, 489.01, 540.04, 596.76]

for i in range(len(temperature)):
    filename = "run.sh.%03d" % (i+1)
    print "writing %s..." % (filename)
    f = open(filename, 'w')
    print_lines(f, i+1)
    f.close()
    os.chmod(filename, 0755)


