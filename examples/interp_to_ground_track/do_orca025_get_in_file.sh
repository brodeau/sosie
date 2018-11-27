#!/bin/bash

# we read everything in the model file, no mesh-mask provided!

F_E="../data/dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4"

F_M="../data/ssh_ORCA025.nc4"

../../bin/interp_to_ground_track.x \
    -i ${F_M} -v ssh -x 'nav_lon' -y 'nav_lat' \
    -p ${F_E} -n adt_unfiltered -S
