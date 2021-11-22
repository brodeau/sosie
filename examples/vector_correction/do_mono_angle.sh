#!/bin/bash

CONF="eORCA12"
NEMOv="4.0.6"
NZ="75"

DIRI="/mnt/meom/workdir/brodeau/${CONF}/${CONF}.L${NZ}-I"

./bin/create_angle_file.x -m ${DIRI}/domain_cfg_L${NZ}_${NEMOv}_${CONF}.nc


fout="cos_sin_angles_${CONF}_L${NZ}_${NEMOv}.nc"

#fmm='/mnt/meom/workdir/brodeau/NANUK1/NANUK1.L${NZ}-I/domain_cfg_L${NZ}_${NEMOv}_NANUK1.nc'

rm -f ${fout}

rsync -avP cost_angle.nc ${fout}

for cv in cosu cosv cosf sint sinu sinv sinf; do

    ncks -A -v ${cv} ${cv}_angle.nc -o ${fout}

done


    

