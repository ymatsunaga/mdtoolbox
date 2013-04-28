#!/usr/bin/python
# coding: utf-8

############# define functions
def print_lines(f, value):
    f.write("harmonic restraint fixed spring constant\n")
    f.write(" &rst\n")
    f.write("   iat=9,15,17,19,\n")
    f.write("   r0=%f, k0=200.0,\n" % (value))
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

