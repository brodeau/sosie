#!/bin/bash

#F_E="../data/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4"
#F_E="${HOME}/DATA/SOSIE/dt_global_s3a_sla_vxxc_20170402_Sentinel3.nc4"
F_E="/data/SOSIE/dt_global_s3a_sla_vxxc_20170402_Sentinel3.nc4"

F_M="${HOME}/DATA/SOSIE/NATL60-CJM165_y2017m04d02.1h_SSH.nc4"

../../bin/interp_to_ground_track.x \
    -i ${F_M} -v sossheig \
    -m /data/gcm_setup/NATL60/NATL60-I/mesh_mask_2D_NATL60.nc \
    -p ${F_E} -n adt_unfiltered
#    -p ${F_E} -n sla_unfiltered
