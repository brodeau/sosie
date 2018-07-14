#!/bin/bash

FMM="/MEDIA/data/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
FIN="/MEDIA/data/eNATL4/sossheig_snap_eNATL4.nc"

#for cv in "sossheig" "sosaline" "sosstsst"; do
for cv in "sossheig" ; do

    ./bin/nemo_coarsener_2d.x -m ${FMM} -i ${FIN} -v ${cv} -o ${cv}_test.nc

done
