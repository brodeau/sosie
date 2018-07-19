#!/bin/bash


if [ "`hostname`" = "merlat" ]; then
    FMM="/MEDIA/data/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    #FIN="/MEDIA/data/eNATL4/eNATL4-BLB400_1h_20080101_20080110_gridT_test.nc4"
    FIN="/MEDIA/data/eNATL4/eNATL4-BLB400_5d_20080111_20080120_gridT.nc"
elif [ "`hostname`" = "login3" ]; then
    FMM="/data/gcm_setup/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    FIN="/data/gcm_output/NEMO/eNATL4/eNATL4-BLB400_1h_20080101_20080110_gridT_merg.nc"
elif [ "`hostname`" = "luitel" ]; then
    FMM="/data/gcm_setup/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    #FIN="/data/gcm_output/NEMO/eNATL4/eNATL4-BLB400_1h_20080101_20080110_gridT_merg.nc"
    FIN="/data/gcm_output/NEMO/eNATL4/eNATL4-BLB400_5d_20080111_20080120_gridT.nc"
else
    echo " Host unknown!"; exit
fi


./bin/nemo_coarsener.x -m ${FMM} -i ${FIN} -P T -o test.nc
