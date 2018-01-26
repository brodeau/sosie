#!/bin/bash


F_E="../data/dt_global_al_phy_vxxc_l3_20130401_20170110_tr.nc4"
#F_E="/home/users/brodeau/DATA/sat/duacs-rep-global-al-phy-unfiltered-l3-v3/2013/dt_global_al_phy_vxxc_l3_20130401_20170110_tr.nc"


#F_E="/home/users/brodeau/DATA/sat/duacs-rep-global-al-phy-unfiltered-l3-v3/monthly/dt_global_al_phy_vxxc_l3_201304_SMALL_SLICE.nc4"
#F_E="/home/users/brodeau/DATA/sat/duacs-rep-global-al-phy-unfiltered-l3-v3/2013/dt_global_al_phy_vxxc_l3_20130404_20170110_tr.nc"



../../bin/interp_to_ground_track.x -i ../data/sss_ORCA1_example.nc -v sos \
                                   -m ../data/mesh_mask_ORCA1v2_light.nc4 \
                                   -p ${F_E} -n
