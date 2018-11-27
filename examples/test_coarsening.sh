#!/bin/bash


if [ "`hostname`" = "merlat" ]; then
    SDIR="/MEDIA/data/eNATL4"
    FMM="${SDIR}/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    FT_2D="${SDIR}/00000721-00001440/eNATL4-BLB400_1h_20080111_20080120_gridT.nc"
    FT_3D="${SDIR}/00000721-00001440/eNATL4-BLB400_5d_20080111_20080120_gridT.nc"
    FU_2D="${SDIR}/00000721-00001440/eNATL4-BLB400_1h_20080111_20080120_gridU.nc"
    FV_2D="${SDIR}/00000721-00001440/eNATL4-BLB400_1h_20080111_20080120_gridV.nc"
    
elif [ "`hostname`" = "login3" ]; then
    FMM="/data/gcm_setup/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    FT_2D="/data/gcm_output/NEMO/eNATL4/eNATL4-BLB400_1h_20080101_20080110_gridT_merg.nc"


elif [ "`hostname`" = "luitel" ]; then
    SDIR="/data/gcm_output/NEMO/eNATL4/eNATL4-BLB400-S"
    FMM="/data/gcm_setup/eNATL4/eNATL4-I/mesh_mask_eNATL4_3.6.nc"
    FT_2D="${SDIR}/00002161-00002880/eNATL4-BLB400_1h_20080131_20080209_gridT.nc"
    FT_3D="${SDIR}/00002161-00002880/eNATL4-BLB400_5d_20080131_20080209_gridT.nc"
    FU_2D="${SDIR}/00002161-00002880/eNATL4-BLB400_1h_20080131_20080209_gridU.nc"
    FV_2D="${SDIR}/00002161-00002880/eNATL4-BLB400_1h_20080131_20080209_gridV.nc"
    FU_3D="${SDIR}/00002161-00002880/eNATL4-BLB400_5d_20080131_20080209_gridU.nc"
    FV_3D="${SDIR}/00002161-00002880/eNATL4-BLB400_5d_20080131_20080209_gridV.nc"
    
else
    echo " Host unknown!"; exit
fi


../bin/nemo_coarsener.x -m ${FMM} -i ${FT_3D} -P T -o test_T_3d.tmp -M
../bin/nemo_coarsener.x -m ${FMM} -i ${FT_2D} -P T -o test_T_2d.tmp

../bin/nemo_coarsener.x -m ${FMM} -i ${FV_2D} -P V -o test_V_2d.tmp
../bin/nemo_coarsener.x -m ${FMM} -i ${FU_2D} -P U -o test_U_2d.tmp

../bin/nemo_coarsener.x -m ${FMM} -i ${FU_3D} -P U -o test_U_3d.tmp
../bin/nemo_coarsener.x -m ${FMM} -i ${FV_3D} -P V -o test_V_3d.tmp
