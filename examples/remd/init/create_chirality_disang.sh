#!/bin/bash -x

ambpdb -p ala.parm <ala.crd >ala.ambpdb
makeCHIR_RST ala.ambpdb chirality.disang


