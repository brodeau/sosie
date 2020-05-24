#!/bin/bash

F_MOD="../temp_EN.4.2.0_360x173L42-ORCA1L75_january_clim_1990-2010.nc"

F_HYS="../data/coord_section_OVIDE2012.nc"

HERE=`pwd`


rm -f *.nc *.dat

echo

    if [ ! -f ${F_MOD} ]; then
        echo "Need to generate `basename ${F_MOD}` !"
        cd ../
        ../bin/sosie3.x -f namelist.example2
        cd ${HERE}
        echo
    fi

../../bin/interp_to_hydro_section.x -i  ${F_MOD} -v temp \
                                -p ${F_HYS} -n temp_OVIDE -S


./plot.gp
