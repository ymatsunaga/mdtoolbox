#!/usr/bin/python
# coding: utf-8

import random

############# define functions
def print_lines(f, value):
    ig = random.randint(0,1000000)
    f.write("production with restraint\n")
    f.write(" &cntrl\n")
    f.write("   ig=%d, \n" % (ig))
    f.write("   irest=1, ntx=5,\n")
    f.write("   igb=8, gbsa=1,\n")
    f.write("   cut=9999.0, rgbmax=9998.0,\n")
    f.write("   ntc=2, ntf=1, tol=0.000001,\n")
    f.write("   ntt=3, gamma_ln=2.0, temp0=300.0,\n")
    f.write("   ntb=0, nscm=10000,\n")
    f.write("   ioutfm=1,\n")
    f.write("   nstlim=500000, dt=0.002,\n")
    f.write("   ntpr=5000, ntwx=5000, ntwv=0, ntwr=500000,\n")
    f.write("   nmropt=1,\n")
    f.write(" /\n")
    f.write(" &wt\n")
    f.write("  type='DUMPFREQ', istep1=5000,\n")
    f.write(" /\n")
    f.write(" &wt\n")
    f.write("  type='END',\n")
    f.write(" /\n")
    f.write("DISANG=run_%d.disang\n" % (value))
    f.write("DUMPAVE=run_%d.dat\n" % (value))
    f.write("\n")

############# main
torsion = range(0, 181, 3)

for i in torsion:
    filename = "run_%d.in" % (i)
    print "writing %s..." % (filename)
    f = open(filename, 'w')
    print_lines(f, i)
    f.close()

