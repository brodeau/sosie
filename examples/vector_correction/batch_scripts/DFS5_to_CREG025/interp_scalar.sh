#!/bin/bash

#PBS -N itrpsc
#PBS -q mpi_1
#PBS -l select=1:ncpus=28:mem=60G
#PBS -l walltime=01:59:00
#PBS -o out_itrpsc.out
#PBS -e err_itrpsc.err
#PBS -m n

DIR_I="/home/datawork-lops-drakkarcom/DATA-REFERENCE/DFS5.2_RD/ALL"
DIR_O="/home1/scratch/ctalandi/FORCING/DFS5.2"

FMM_I="/home3/datawork/lbrodeau/ECMWF/mask_N128_corrected_LB.nc4"
FMM_O="/home1/scratch/ctalandi/FORCING/CREG025.L75_mesh_mask_clean.nc"


LIST_VAR="t2 q2 radlw radsw precip snow"

#Y1=1958
Y1=1986
Y2=$((Y1+28-1))

#Y1=2014
#Y2=2015



cd /home3/scratch/lbrodeau/SOSIE/




for vv in ${LIST_VAR}; do

    yy=${Y1}

    while [ ${yy} -le ${Y2} ]; do

        fnml="${vv}_${yy}.tmp"

        cat > ${fnml} <<EOF
!! -------------------
!! Namelist for SOSIE
!! -------------------
!!--------------------------------------------------------------------------
!!
&ndom_src
csource    = 'DFS5.2'
ivect      = 0
l_reg_src  = .true.
cf_src     = '${DIR_I}/drowned_${vv}_DFS5.2_y${yy}.nc'
cv_src     = '${vv}'
cv_t_src   = 'time'
cf_x_src   = '${DIR_I}/drowned_${vv}_DFS5.2_y${yy}.nc'
cv_lon_src = 'lon0'
cv_lat_src = 'lat0'
cf_lsm_src = '${FMM_I}'
cv_lsm_src = 'lsm'
ewper_src  = 0
/
!!
!!
&ndom_trg
ctarget    = 'CREG025'
l_reg_trg  = .false.
cf_x_trg   = '${FMM_O}'
cv_lon_trg = 'nav_lon'
cv_lat_trg = 'nav_lat'
cf_lsm_trg = '${FMM_O}'
cv_lsm_trg = 'tmask'
ewper_trg  = -1
/
!!
!!
&ninterp
cmethod     = 'akima'
!!
idrown      = 70,50
l_save_drwn = .false.
ismooth     = 0
jt1         = 0
jt2         = 0
!!jt1         = 1
!!jt2         = 16
jplev       = 1
vmax        =  1.E6
vmin        = -1.E6
ismooth_out = 0
ibx_xtra_sm = 0, 0,0, 0,0  ! Extra-smoothing on a given rectangular region: ibx_xtra_sm = ntimes, i1,j1, i2,j2
/
!!
&noutput
cv_out    = '${vv}'
cu_out    = ''
cln_out   = ''
cv_t_out  = 'time_counter'
cd_out    = '${DIR_O}'
cextra    = 'y${yy}'
lmout     = .true.
rmiss_val = -9999.
lct       = .false.
t0        = 0.
t_stp     = 0.
/

EOF

        echo "Launching: sosie3.x -f ${fnml} &"
        sosie3.x -f ${fnml} &
        echo

        yy=$((yy + 1))
    done

    wait

done
