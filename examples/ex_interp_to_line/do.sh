#!/bin/bash

F_IN="../data/sss_ORCA1_example.nc"
V_IN="sos"
M_IN="../data/mesh_mask_ORCA1v2_light.nc4"

echo "../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}"
../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}
echo
