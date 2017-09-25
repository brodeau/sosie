#!/bin/bash

BRKD_DIR="${HOME}/DEV/barakuda"


for ff in \
        bw_400x225-0.25degx0.25deg_akima \
        bw_400x225-0.25degx0.25deg_bilin \
        bw_1600x900-1degx1deg_akima \
        bw_1600x900-1degx1deg_smooth+akima \
        bw_1600x900-1degx1deg_bilin \
        bw_1600x900-1degx1deg_smooth+bilin; do

    python ${BRKD_DIR}/python/exec/netcdf_to_image_bw.py ${ff}.nc4 bw bmp

done
