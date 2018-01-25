#!/bin/bash

#PATH_NATL="${HOME}/DATA/NEMO/NATL60_data" ; #meolkerg
PATH_NATL="${HOME}/Dropbox/tmp/NATL60"

../../bin/interp_to_ground_track.x -i ${PATH_NATL}/NATL60-CJM165_y2013m04d01.1h_SSH_jt_21-23_time_20h-22h.nc4 -v sossheig \
                            -m ${PATH_NATL}/mesh_t_minimum_NATL60.nc4 \
                            -p ../data/dt_global_al_phy_vxxc_l3_20130401_20170110_tr.nc4 -n

