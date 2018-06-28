#!/bin/bash

F_E="../data/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4"
#F_E="${HOME}/DATA/SOSIE/dt_global_s3a_sla_vxxc_20170402_Sentinel3.nc4"

#F_E="../dt_short.nc"

F_M="../data/ssh_ORCA1_20170101_20171231_grid_T.nc4"

CMD="../../bin/interp_to_ground_track.x \
    -i ${F_M} -v ssh \
    -m ../data/mesh_mask_ORCA1v2_light.nc4 \
    -p ${F_E} -n adt_unfiltered -S"


echo
echo ${CMD}
echo
${CMD}
    

