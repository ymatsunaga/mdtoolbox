#!/bin/bash
#
# Exctract a trajectory with constant temperature from replica trajectories
# Usage: 
# $ ./extract_temperature_trj.sh 2>&1 | tee extract_temperature_trj.log
#
# temperature = [300.00, 331.00, 364.48, 401.69, 443.05, 489.01, 540.04, 596.76]

ptraj parm << EOF
trajin ../prod/run.nc.001 remdtraj remdtrajtemp 300.00
trajout run.300.00K.nc netcdf

go
EOF

