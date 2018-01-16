#!/bin/bash
#../../bin/interp_to_ephem.x -i ../data/ssh_ORCA025.nc4 -v ssh  -m ../data/mesh_mask_ORCA025_light.nc4 \
#                           -p ephem_calval_june2015_sph_short.txt -a
##                           -p ephem_calval_june2015_sph.txt

# ORCA1:
#../../bin/interp_to_ephem.x -i ../data/sss_ORCA1_example.nc -v sos  -m ../data/mesh_mask_ORCA1v2_light.nc4 \
#                           -p ephem_calval_june2015_sph.txt -f 600.,3600.
#exit

../../bin/interp_to_ephem.x -i ../data/sss_ORCA1_example.nc -v sos  -m ../data/mesh_mask_ORCA1v2_light.nc4 \
    -p /home/users/brodeau/sat/saral_sentinel2_max_cls/SARAL-Altika/dt_global_alg_sla_vxxc_20170401_20180104_tr_41496-42490.nc4 \
    -f 600.,3600. -n
#-p /home/users/brodeau/sat/saral_sentinel2_max_cls/SARAL-Altika/dt_global_alg_sla_vxxc_20170401_20180104_tr_42490-44404.nc4 \
#    -p /home/users/brodeau/sat/saral_sentinel2_max_cls/SARAL-Altika/dt_global_alg_sla_vxxc_20170501_20180104_tr.nc4 \


exit

#F_IN="../data/sss_ORCA1_example.nc"
#V_IN="sos"
#M_IN="../data/mesh_mask_ORCA1v2_light.nc4"

#echo "../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}"
#../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}
#echo



echo "../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}"
../../bin/interp_to_line.x -i ${F_IN} -v ${V_IN} -x nav_lon -y nav_lat -t time_counter -m ${M_IN}
echo



