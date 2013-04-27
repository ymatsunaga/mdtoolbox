#!/usr/bin/python
# coding: utf-8

############# define functions
def print_lines(f, value):
    f.write("harmonic restraint changing spring constant\n")
    f.write(" &rst\n")
    f.write("   iat=9,15,17,19,\n")
    f.write("   r0=%f, r0a=%f, k0=0.01, k0a=200.0,\n" % (value, value))
    f.write("   ifvari=1, nstep1=0, nstep2=250000,\n")
    f.write(" /\n")
    f.write(" &rst\n")
    f.write("   iat=9,15,17,19,\n")
    f.write("   r0=%f, r0a=%f, k0=200.0, k0a=200.0,\n" % (value, value))
    f.write("   ifvari=1, nstep1=250001, nstep2=500000,\n")
    f.write(" /\n")
    f.write("\n")

############# main
torsion = range(0, 181, 3)

for i in torsion:
    filename = "run_%d.disang" % (i)
    print "writing %s..." % (filename)
    f = open(filename, 'w')
    print_lines(f, i)
    f.close()

