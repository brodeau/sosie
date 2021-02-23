#!/bin/bash

# This example does the interpolation of a surface 10m wind vector field from a
# regular grid (ECMWF) to the ORCA1 tri-polar NEMO grid, including the rotation of
# the vector to account for the distortion of the ORCA grid in the North.
#
# *** u10, v10 regular => ORCA1:
#   => u10 on grid_T and v10 on grid_T
#   => u10 on grid_U and v10 on grid_V


fu="u10_1x1-deg-ORCA1_gridT_grid_T.nc"
fv="v10_1x1-deg-ORCA1_gridT_grid_T.nc"


if [ ! -f ${fu} ] || [ ! -f ${fv} ]; then
    echo "Run script 'test_corr_vect_ORCA1.sh' prior to this one so that these files exist:"
    echo "   ${fu} , ${fv} "
    echo
    exit
fi

CMD="../../bin/corr_vect.x -I -i ${fu} ${fv} -v u10 v10 -t time_counter -m ../data/mesh_mask_ORCA1v2_light.nc4"
echo; echo ${CMD}; ${CMD}; echo
