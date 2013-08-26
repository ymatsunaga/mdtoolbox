#!/usr/bin/python
# coding: utf-8

import random

############# define functions
def print_lines(f, value):
    ig = random.randint(0, 1000000)
    f.write("production\n")
    f.write(" &cntrl\n")
    f.write("   numexchg=3000,\n")
    f.write("   ig=%d,\n" % (ig))
    f.write("   irest=1, ntx=5,\n")
    f.write("   igb=8, gbsa=1,\n")
    f.write("   cut=9999.0, rgbmax=9998.0,\n")
    f.write("   ntc=2, ntf=1, tol=0.000001,\n")
    f.write("   ntt=3, gamma_ln=2.0, temp0=%f,\n" % (value))
    f.write("   ntb=0, nscm=10000,\n")
    f.write("   ioutfm=1,\n")
    f.write("   nstlim=1000, dt=0.001,\n")
    f.write("   ntpr=1000, ntwx=1000, ntwv=0, ntwr=1000,\n")
    f.write("   nmropt=1,\n")
    f.write(" /\n")
    f.write(" &wt\n")
    f.write("   type='END'\n")
    f.write(" /\n")
    f.write("DISANG=chirality.disang\n")
    f.write("\n")

############# main
temperature = [300.00, 331.00, 364.48, 401.69, 443.05, 489.01, 540.04, 596.76]

for i in range(len(temperature)):
    filename = "run.in.%03d" % (i+1)
    print "writing %s..." % (filename)
    f = open(filename, 'w')
    print_lines(f, temperature[i])
    f.close()

