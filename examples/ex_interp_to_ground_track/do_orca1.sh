#!/bin/bash

F_E="../data/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4"

F_M="../data/ssh_ORCA1_20170101_20171231_grid_T.nc4"

../../bin/interp_to_ground_track.x -i ${F_M} -v ssh \
                                   -m ../data/mesh_mask_ORCA1v2_light.nc4 \
                                   -p ${F_E} -n sla_unfiltered
