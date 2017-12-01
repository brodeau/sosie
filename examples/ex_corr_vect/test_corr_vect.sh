#!/bin/bash

# On U,V grids:
../../bin/sosie.x -f namelist.example4_O2uv_x
../../bin/sosie.x -f namelist.example4_O2uv_y
../../bin/corr_vect.x -G U -f namelist.example4_O2uv -m ../data/mesh_mask_ORCA2.nc4

# On T grid:
../../bin/sosie.x -f namelist.example4_O2t_x
../../bin/sosie.x -f namelist.example4_O2t_y
../../bin/corr_vect.x -G T -f namelist.example4_O2t -m ../data/mesh_mask_ORCA2.nc4


#rm -f uraw*.nc* vraw*.nc*


