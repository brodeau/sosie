#!/bin/bash


if [ "`hostname`" = "merlat" ]; then
    FMM="/MEDIA/data/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    FIN="/MEDIA/data/eNATL4/eNATL4-BLB400_1h_20080101_20080110_gridT_test.nc4"
elif [ "`hostname`" = "login3" ]; then
    FMM="/data/gcm_setup/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    FIN="/data/gcm_output/NEMO/eNATL4/eNATL4-BLB400_1h_20080101_20080110_gridT_merg.nc"
else
    echo " Host unknown!"; exit
fi


./bin/nemo_coarsener.x -m ${FMM} -i ${FIN} -v sossheig -o 2d_test.nc

#for cv in "sossheig" "sosaline" "sosstsst"; do
#for cv in "sossheig" ; do
    #./bin/nemo_coarsener_2d.x -m ${FMM} -i ${FIN} -v ${cv} -o ${cv}_test.nc
#done
