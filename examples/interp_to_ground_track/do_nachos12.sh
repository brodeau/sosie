#!/bin/bash

F_E="../data/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4"

F_M="${HOME}/DATA/SOSIE/NACHOS12.L75-MAA13_y2017m04d02.1h_SSH.nc4"

../../bin/interp_to_ground_track.x \
    -i ${F_M} -v sossheig \
    -m ${HOME}/DATA/SOSIE/mesh_mask_1st_level_NACHOS12.nc4 \
    -p ${F_E} -n adt_unfiltered -S
