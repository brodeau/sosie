#!/bin/bash

#../../bin/interp_to_line.x -i ../data/ssh_ORCA025.nc4 -v ssh  -m ../data/mesh_mask_ORCA025_light.nc4 \
#                           -p ephem_calval_june2015_sph_short.txt -a
##                           -p ephem_calval_june2015_sph.txt

# ORCA1:
../../bin/interp_to_line.x -i ../data/sss_ORCA1_example.nc -v sos  -m ../data/mesh_mask_ORCA1v2_light.nc4 \
                           -p ephem_calval_june2015_sph.txt -a



exit

#F_IN="../data/sss_ORCA1_example.nc"
#V_IN="sos"
#M_IN="../data/mesh_mask_ORCA1v2_light.nc4"

#echo "../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}"
#../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}
#echo




