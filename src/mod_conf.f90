MODULE MOD_CONF

   USE io_ezcdf, ONLY: nbatt_max, var_attr

   IMPLICIT NONE

   !CHARACTER(len=16), PARAMETER :: csosie_version = 'svn_trunk' !: SOSIE version
   CHARACTER(len=16), PARAMETER :: csosie_version = '3.0 beta' !: SOSIE version

   INTEGER, PARAMETER :: wpl = 4        !: local working precision

   LOGICAL, SAVE :: l_first_call_interp_routine, &
      &             l_glob_lon_wize, l_glob_lat_wize
   

   LOGICAL, DIMENSION(2) :: l_2d_grid_yet_regular = (/ .FALSE. , .FALSE. /) !: in case input (i=1) and/or output (i=2) domains are
   !!                                      defined with 2D longitude and latitude arrays but their grid is actually regular...

   CHARACTER(len=400) :: cf_nml_sosie = 'namelist' !: namelist for sosie.x

   TYPE, PUBLIC :: grid_type
      INTEGER          :: ifld_nord
      CHARACTER(len=1) :: cgrd_type
   END TYPE grid_type

   TYPE(grid_type), SAVE  :: gt_orca_in, gt_orca_out   ! Source, target, grid is not an ORCA grid (0), ORCA2 (4), ORCA1 (6), to be completed!
   INTEGER, SAVE          :: i_orca_in, i_orca_out   ! Source, target, grid is not an ORCA grid (0), ORCA2 (4), ORCA1 (6), to be completed!
   CHARACTER(len=1), SAVE :: c_orca_in, c_orca_out   ! Source, target, grid is not an ORCA grid (0), ORCA2 (4), ORCA1 (6), to be completed!

   
   INTEGER :: &
      &   Ntr, &                      !: time dimmension == number of time records to interpolate!
      &   ni_in, nj_in, nk_in,   &    !: dimension of input arrays
      &   Ntr0, nlev, j_start, j_stop, &
      &   ni_out, nj_out, nk_out, &      !: dimension of output arrays
      
      &   nlon_inc_in,  &  !: wether input longitude increases with i
      &   nlat_inc_in,  &  !: wether input latitude  increases with j
      &   nlon_inc_out, &  !:  //    output  //
      &   nlat_inc_out, &  !:  //    output  //
      
      &   jj_ex_top, jj_ex_btm, &
      &   i_chg_lon

   REAL(8) ::  &
      &   lon_min_1, lon_max_1, & !lolo
      &   max_lon_in,  min_lon_in, &
      &   max_lat_in,  min_lat_in, &
      &   max_lat_out, min_lat_out

   LOGICAL    :: &
      &   l3d, l_int_3d, &  !: wether input variables is 3D
      &   ltime             !: wether input field contain a time variable





   !!       NAMELIST PARAMETERS
   !!       ###################

   !! Logicals
   !! --------
   LOGICAL :: &
      &    lregin,  & ! whether source grid is regular or not
      &    lregout, & ! whether target grid is regular or not
      &    ldrown,  & ! whether to use the DROWN procedure to extrapolate sea values on earth
      &    lmout,   & ! masking output or not
      &    lct        ! time control


   ! &    ldrown_out ! drowning out put (only if lmout=FALSE and cv_lsm_out /= ''

   !! Integers
   !! --------
   INTEGER  :: &
      &    jt1=0, jt2=0, ismooth=0

   !! Name of files or directories:
   !! -----------------------------
   CHARACTER(LEN=400) :: &
      &    cf_in = '',      &
      &    cf_x_in = '',    &
      &    cf_z_in = '',    &
      &    cf_coor_in = '', &
      &    cf_z_out = '',   &
      &    cf_lsm_in = '',  &
      &    cf_x_out = '',   &
      &    cf_lsm_out = '', &
      &    cf_out = '',     &
      &    cd_out = '',     &
      &    csource = '',    &
      &    ctarget = '',    &
      &    cextra = '',     &
      &    cf_short = '',   &
      &    cf_misc = '',    &
      &    cpat = ''


   !! Name of variables :
   !! -------------------
   CHARACTER(LEN=80) ::  &
      &    cv_t_in = '',    &
      &    cv_lon_in = '',  &
      &    cv_lat_in = '',  &
      &    cv_z_in = '',    &
      &    cv_z_out = '',   &
      &    cv_z_out_name = '',   &
      &    cv_lsm_in = '',  &
      &    cv_in = '',      &
      &    cv_lon_out = 'lon', &
      &    cv_lat_out = 'lat', &
      &    cv_lsm_out = '', &
      &    ca_missval

   !! S-coordinates specific
   !! -----------------------
   CHARACTER(LEN=80) ::  &
      &    cv_bathy_in = '',&
      &    cv_bathy_out= '',&
      &    cf_bathy_in = '',&
      &    cf_bathy_out= '',&
      &    ctype_z_in='z', &
      &    ctype_z_out='z'  

   TYPE scoord_params
            integer :: Vtransform
            integer :: Vstretching
            integer :: Nlevels
            real(wpl) :: theta_s
            real(wpl) :: theta_b
            real(wpl) :: Tcline
            real(wpl) :: hmin
   END TYPE scoord_params

   TYPE (scoord_params) :: ssig_in, ssig_out

   REAL(8),  DIMENSION(:),    ALLOCATABLE ::   &
       &   Cs_rho, Sc_rho
  
   !! Netcdf output strings :
   !! -----------------------
   CHARACTER(LEN=5)  ::  cmethod
   CHARACTER(LEN=80) ::  &
      &    cv_t_out = '',   &
      &    cv_out = '',     &
      &    cu_out
   !!
   CHARACTER(LEN=400) :: cln_out

   INTEGER           :: &
      &     ivect,    &
      &     jplev,    &
      &     ewper,    &
      &     ewper_out

   REAL(4)       ::   &
      &     rmaskvalue,  &
      &     t0, t_stp,   &
      &     vmax,        &
      &     vmin

   INTEGER, PARAMETER :: numnam = 13


   !! Attribute info for some variables LOLO: move to io_ezcdf instead????
   !! ---------------------------------
   !! Time vector:
   INTEGER                              :: nb_att_t, nb_att_lon, nb_att_lat, nb_att_z, nb_att_F
   TYPE(var_attr), DIMENSION(nbatt_max) :: vatt_info_t, vatt_info_lon, vatt_info_lat, &
      &                                    vatt_info_z, vatt_info_F


   !! Data arrays
   !! -----------

   REAL(8),  DIMENSION(:),    ALLOCATABLE ::   &
      &   vt0, vt          !: time arrays

   !! Arrays on source grid :
   REAL(8),  DIMENSION(:,:,:),    ALLOCATABLE ::  &
      &   depth_in,   &
      &   depth_out

  REAL(4),  DIMENSION(:,:,:),    ALLOCATABLE ::  &    
       &   depth_in_tmp

  REAL(wpl),  DIMENSION(:,:),    ALLOCATABLE ::  &    
       &   bathy_in,  &
       &   bathy_out

   REAL(wpl),  DIMENSION(:,:),  ALLOCATABLE ::  &
      &   data_in,    &   !: data array on source grid
      &   data_in_b

   REAL(8),  DIMENSION(:,:),  ALLOCATABLE ::  &
      &   lon_in,   &
      &   lat_in

   REAL(wpl),  DIMENSION(:,:,:),  ALLOCATABLE ::  &
      &      data3d_in,  &
      &      data3d_out, &
      &      data3d_tmp  ! horizontal target resol. + vertical source resol.

   INTEGER(1),   DIMENSION(:,:),  ALLOCATABLE :: IGNORE !: point of target domain to disregard (IGNORE==0) 
   
   INTEGER(2),   DIMENSION(:,:,:),  ALLOCATABLE ::   &
      &   mask_in, mask_in_b, &  !: land-sea mask on input grid
      &   mask_out               !: land-sea mask on output grid

   !! Arrays on target grid :
   REAL(wpl) , DIMENSION(:,:),  ALLOCATABLE ::             &
      &   data_out     ! data array interpolated

   REAL(8) , DIMENSION(:,:),  ALLOCATABLE ::             &
      &   lon_out,   & ! longitude array on target grid
      &   lon_out_b, & !  //   backup
      &   lat_out      ! latitude array on target grid

END MODULE MOD_CONF
