#!/bin/bash

# This example does the interpolation of a surface 10m wind vector field from a
# regular grid (ECMWF) to the ORCA1 tri-polar NEMO grid, including the rotation of
# the vector to account for the distortion of the ORCA grid in the North.
#
# *** u10, v10 regular => ORCA1:
#   => u10 on grid_T and v10 on grid_T
#   => u10 on grid_U and v10 on grid_V
#
# On T grid:
../../bin/sosie3.x -f namelist.example4_O1t_x
../../bin/sosie3.x -f namelist.example4_O1t_y

# Correcting on T grid:
../../bin/corr_vect.x -G T -f namelist.example4_O1t -m ../data/mesh_mask_ORCA1v2_light.nc4

# Correcting on U,V grid:
../../bin/corr_vect.x -G U -f namelist.example4_O1t -m ../data/mesh_mask_ORCA1v2_light.nc4

rm -f uraw*.nc* vraw*.nc*
