#!/bin/bash

PATH_NATL="${HOME}/Dropbox/tmp/NATL60"

F_E="../data/dt_global_al_phy_vxxc_l3_20130401_20170110_tr.nc4"
F_M="${PATH_NATL}/NATL60-CJM165_y2013m04d01.1h_SSH_jt_21-23_time_20h-22h.nc4"


../../bin/interp_to_ground_track.x \
    -i ${F_M} -v sossheig \
    -m ${PATH_NATL}/mesh_t_minimum_NATL60.nc4 \
    -p ${F_E} -n sla_unfiltered
