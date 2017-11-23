#!/bin/bash

../bin/sosie.x -f namelist.example4_x
../bin/sosie.x -f namelist.example4_y


# On U,V grids:
../bin/corr_vect.x -G U -f namelist.example4 -m data/mesh_mask_ORCA1v2_light.nc4


#rm -f uraw*.nc* vraw*.nc*


