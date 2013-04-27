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
    f.write(" -i run_%d.in \\\n" % (value))
    f.write(" -c ../eq1/run.rst \\\n")
    f.write(" -o run_%d.out \\\n" % (value))
    f.write(" -r run_%d.rst \\\n" % (value))
    f.write(" -x run_%d.nc\n" % (value))
    f.write("\n")

############# main
torsion = range(0, 181, 3)

for i in torsion:
    filename = "run_%d.sh" % (i)
    print "writing %s..." % (filename)
    f = open(filename, 'w')
    print_lines(f, i)
    f.close()
    os.chmod(filename, 0755)


