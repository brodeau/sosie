#!/bin/bash

#PBS -N itrpW
#PBS -q mpi_1
#PBS -l select=1:ncpus=28:mem=120G
#PBS -l walltime=02:59:00
#PBS -o out_itrpW.out
#PBS -e err_itrpW.err
#PBS -m n

DIR_I="/home/datawork-lops-drakkarcom/DATA-REFERENCE/DFS5.2_RD/ALL"
DIR_O="/home1/scratch/ctalandi/FORCING/DFS5.2"

FMM_I="/home3/datawork/lbrodeau/ECMWF/mask_N128_corrected_LB.nc4"
FMM_O="/home1/scratch/ctalandi/FORCING/CREG025.L75_mesh_mask_clean.nc"

LIST_VAR="u10 v10"


#Y1=1958
Y1=1986
Y2=$((Y1+28-1))

#Y1=2014
#Y2=2015




cd /home3/scratch/lbrodeau/SOSIE/

ivect=1

for vv in ${LIST_VAR}; do

    echo; echo; echo
    echo " #### Doing ${vv} ($ivect)"
    echo

    yy=${Y1}

    while [ ${yy} -le ${Y2} ]; do


        if [ ${ivect} -eq 1 ]; then
            fnml="w10_${yy}.tmp_x"
        fi

        if [ ${ivect} -eq 2 ]; then
            fnml="w10_${yy}.tmp_y"
        fi

        cat > ${fnml} <<EOF
!! -------------------
!! Namelist for SOSIE
!! -------------------
!!--------------------------------------------------------------------------
!!
&ndom_src
csource    = 'DFS5.2'
ivect      = ${ivect}
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
        if [ ! -f ${DIR_O}/uraw_DFS5.2-CREG025_y${yy}.nc ] || [ ! -f ${DIR_O}/vraw_DFS5.2-CREG025_y${yy}.nc ]; then
            echo; echo; echo; echo
            echo "Launching: sosie3.x -f ${fnml} &"; echo
            sosie3.x -f ${fnml} &
            echo; echo
        fi
        
        yy=$((yy + 1))
    done

    wait

    ivect=$((ivect+1))

done



yy=${Y1}
while [ ${yy} -le ${Y2} ]; do
    
    fnml0="w10_${yy}.tmp"
    
    echo; echo; echo
    echo "Launching: corr_vect.x -G T -f ${fnml0} -m ${FMM_O} &"
    corr_vect.x -G T -f ${fnml0} -m ${FMM_O} &
    echo; echo

    yy=$((yy + 1))
done
wait

rm -f uraw_*.nc vraw_*.nc
