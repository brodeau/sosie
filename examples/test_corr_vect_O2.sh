#!/bin/bash

../bin/sosie.x -f namelist.example4_O2_x
../bin/sosie.x -f namelist.example4_O2_y


# On U,V grids:
../bin/corr_vect.x -G U -f namelist.example4_O2 -m data/mesh_mask_ORCA2.nc4


#rm -f uraw*.nc* vraw*.nc*


