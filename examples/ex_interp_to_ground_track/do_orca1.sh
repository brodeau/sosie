#!/bin/bash


F_E="../data/dt_global_al_phy_vxxc_l3_20130401_20170110_tr.nc4"

../../bin/interp_to_ground_track.x -i ../data/sss_ORCA1_example.nc -v sos \
                                   -m ../data/mesh_mask_ORCA1v2_light.nc4 \
                                   -p ${F_E} -n sla_unfiltered
