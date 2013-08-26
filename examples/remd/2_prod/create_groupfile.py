#!/usr/bin/python
# coding: utf-8

############# define functions
def print_line(f, value):
    s = ""
    s = s + "-O -remlog rem.log -remtype rem.type -rem 1 "
    s = s + "-p parm "
    s = s + "-i run.in.%03d " % (value)
    s = s + "-c ../1_eq/run.rst.%03d " % (value)
    s = s + "-o run.out.%03d " % (value)
    s = s + "-r run.rst.%03d " % (value)
    s = s + "-x run.nc.%03d\n" % (value)
    f.write(s)

############# main
temperature = [300.00, 331.00, 364.48, 401.69, 443.05, 489.01, 540.04, 596.76]

filename = "groupfile"
print "writing %s..." % (filename)
f = open(filename, 'w')
for i in range(len(temperature)):
    print_line(f, i+1)
f.close()

