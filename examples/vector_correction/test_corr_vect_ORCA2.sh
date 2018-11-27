#!/bin/bash

# This example does the interpolation of a surface 10m wind vector field from a
# regular grid (ECMWF) to the ORCA2 tri-polar NEMO grid, including the rotation of
# the vector to account for the distortion of the ORCA grid in the North.
#
# *** u10, v10 regular => ORCA2:
#   => u10 on grid_T and v10 on grid_T
#   => u10 on grid_U and v10 on grid_V
#
#
#
#
# On T grid:
../../bin/sosie.x -f namelist.example4_O2t_x
../../bin/sosie.x -f namelist.example4_O2t_y

# Correcting on T grid:
../../bin/corr_vect.x -G T -f namelist.example4_O2t -m ../data/mesh_mask_ORCA2.nc4
# Correcting on U,V grid:
../../bin/corr_vect.x -G U -f namelist.example4_O2t -m ../data/mesh_mask_ORCA2.nc4

#rm -f uraw*.nc* vraw*.nc*



