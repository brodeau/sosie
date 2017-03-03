#!/bin/bash

# File containing field to be drowned:
ftd=data/zoom_SA_sst_reyn_025.nc

# Output file:
fd=`basename ${ftd} | sed -e s/'.nc'/'_drowned.nc'/g`

echo "../bin/mask_drown_field.x -D -i ${ftd} -v sst -m 0 -p -1 -o ${fd}"
../bin/mask_drown_field.x -D -i ${ftd} -v sst -m 0 -p -1 -o ${fd}

echo "Check file ${fd} the drowned version of file ${ftd} !"; echo

