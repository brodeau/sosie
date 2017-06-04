#!/bin/bash

../bin/sosie.x -f namelist.example4_x
../bin/sosie.x -f namelist.example4_y

../bin/corr_vect.x -x u10 -y v10 -f namelist.example4_x -m data/mesh_mask_ORCA1_light.nc

rm -f uraw*.nc* vraw*.nc*

echo "Check files u10_1x1-deg-ORCA1_vct.nc4 and v10_1x1-deg-ORCA1_vct.nc4 !"; echo
