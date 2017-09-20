#!/bin/bash

BRKD_DIR="${HOME}/DEV/barakuda"

python ${BRKD_DIR}/python/exec/netcdf_to_image_bw.py bw_1600x900-1degx1deg_bilin.nc4 bw bmp
python ${BRKD_DIR}/python/exec/netcdf_to_image_bw.py bw_1600x900-1degx1deg_smooth+bilin.nc4 bw bmp

