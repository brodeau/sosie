MODULE MOD_CONF

   USE io_ezcdf, ONLY: nbatt_max, var_attr

   IMPLICIT NONE

   CHARACTER(len=16), PARAMETER :: csosie_version = '3.0' !: SOSIE version

   INTEGER, PARAMETER :: wpl = 4        !: local working precision

   LOGICAL, SAVE :: l_first_call_interp_routine, &
      &             l_drown_src, & ! DROWN source field
      &             l_glob_lon_wize, l_glob_lat_wize, &
      &             l_identical_levels=.FALSE.  ! if true and 3D interpolation => means that source and target vertical grids are identical!


   LOGICAL, DIMENSION(2) :: l_2d_grid_yet_regular = (/ .FALSE. , .FALSE. /) !: in case input (i=1) and/or output (i=2) domains are
   !!                                      defined with 2D longitude and latitude arrays but their grid is actually regular...

   CHARACTER(len=400) :: cf_nml_sosie = 'namelist' !: namelist for sosie.x

   TYPE, PUBLIC :: grid_type
      INTEGER          :: ifld_nord
      CHARACTER(len=1) :: cgrd_type
   END TYPE grid_type

   TYPE(grid_type), SAVE  :: gt_orca_src, gt_orca_trg   ! Source, target, grid is not an ORCA grid (0), ORCA2 (4), ORCA1 (6), to be completed!
   INTEGER, SAVE          :: i_orca_src, i_orca_trg   ! Source, target, grid is not an ORCA grid (0), ORCA2 (4), ORCA1 (6), to be completed!
   CHARACTER(len=1), SAVE :: c_orca_src, c_orca_trg   ! Source, target, grid is not an ORCA grid (0), ORCA2 (4), ORCA1 (6), to be completed!


   INTEGER :: &
      &   Ntr, &                      !: time dimmension == number of time records to interpolate!
      &   ni_src, nj_src, nk_src,   &    !: dimension of input arrays
      &   Ntr0, nlev, j_start, j_stop, &
      &   ni_trg, nj_trg, nk_trg, &      !: dimension of output arrays
      
      &   nlon_icr_src,  &  !: wether input longitude increases with i
      &   nlat_icr_src,  &  !: wether input latitude  increases with j
      &   nlon_icr_trg, &  !:  //    output  //
      &   nlat_icr_trg, &  !:  //    output  //
      
      &   jj_ex_top, jj_ex_btm, &
      &   i_chg_lon

   REAL(8) ::  &
      &   lon_min_1, lon_max_1, & !lolo
      &   max_lon_src,  min_lon_src, &
      &   max_lat_src,  min_lat_src, &
      &   max_lat_trg, min_lat_trg

   LOGICAL    :: &
      &   l3d, l_itrp_3d, &  !: wether input variables is 3D
      &   ltime             !: wether input field contain a time variable





   !!       NAMELIST PARAMETERS
   !!       ###################

   !! Logicals
   !! --------
   LOGICAL :: &
      &    l_reg_src=.TRUE.,    & ! whether source grid is regular or not
      &    l_reg_trg=.FALSE.,   & ! whether target grid is regular or not
      &    l_save_drwn=.FALSE., & ! save the input field drowned (on the source grid)
      &    lmout=.FALSE.,       & ! masking output or not
      &    lct=.FALSE.            ! time control

   !! Integers
   !! --------
   INTEGER  :: &
      &    jt1=0,       &
      &    jt2=0,       &
      &    ismooth=0,   &
      &    ismooth_out=0


   INTEGER, DIMENSION(0:4) :: ibx_xtra_sm = (/ 0, 0,0, 0,0 /)  ! box for exta-smoothing => (/ntimes, i1,j1, i2,j2/)


   !! Name of files or directories:
   !! -----------------------------
   CHARACTER(LEN=400) ::    &
      &    cf_src = '',     &
      &    cf_x_src = '',   &
      &    cf_z_src = '',   &
      &    cf_z_trg = '',   &
      &    cf_lsm_src = '', &
      &    cf_x_trg = '',   &
      &    cf_lsm_trg = '', &
      &    cf_out = '',     &
      &    cd_out = '',     &
      &    cv_z_out = 'depth', &
      &    csource = '',    &
      &    ctarget = '',    &
      &    cextra = '',     &
      &    cf_short = '',   &
      &    cf_misc = '',    &
      &    cpat = ''


   !! Name of variables :
   !! -------------------
   CHARACTER(LEN=80) ::  &
      &    cv_t_src = '',    &
      &    cv_lon_src = '',  &
      &    cv_lat_src = '',  &
      &    cv_z_src = '',    &
      &    cv_z_trg = '',    &
      &    cv_lsm_src = '',  &
      &    cv_src = '',      &
      &    cv_lon_trg = 'lon', &
      &    cv_lat_trg = 'lat', &
      &    cv_lsm_trg = '', &
      &    ca_missval


   TYPE idrown_info
      INTEGER :: np_penetr
      INTEGER :: nt_smooth
      LOGICAL :: l_msk_chg  ! mask is changing with time
   END TYPE idrown_info
   TYPE(idrown_info) ::  idrown    ! number of pixels to propagate sea-values onto land (DROWN), by default, no DROWN = 0!
   !                               ! and how many times to smooth the drowned area


   !! S-coordinates specific
   !! -----------------------
   CHARACTER(LEN=80) ::  &
      &    cv_bathy_src = '',&
      &    cv_bathy_trg= '',&
      &    cf_bathy_src = '',&
      &    cf_bathy_trg= '',&
      &    ctype_z_src='z', &
      &    ctype_z_trg='z'

   TYPE scoord_params
      integer :: Vtransform
      integer :: Vstretching
      integer :: Nlevels
      real(wpl) :: theta_s
      real(wpl) :: theta_b
      real(wpl) :: Tcline
      real(wpl) :: hmin
   END TYPE scoord_params

   TYPE (scoord_params) :: ssig_src, ssig_trg

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
      &     ewper_src,    &
      &     ewper_trg

   REAL(4)       ::   &
      &     rmiss_val,  &
      &     t0, t_stp,   &
      &     vmax,        &
      &     vmin

   INTEGER, PARAMETER :: numnam = 13


   !! Attribute info for some variables LOLO: move to io_ezcdf instead????
   !! ---------------------------------
   !! Time vector:
   INTEGER     :: nb_att_t, &
      &           nb_att_lon_src, nb_att_lat_src, nb_att_lon_trg, nb_att_lat_trg, nb_att_z_trg, nb_att_z_src, nb_att_F
   TYPE(var_attr), DIMENSION(nbatt_max) :: vatt_info_t, &
      &                                    vatt_info_lon_src, vatt_info_lat_src, &
      &                                    vatt_info_lon_trg, vatt_info_lat_trg, &
      &                                    vatt_info_z_trg, vatt_info_z_src, vatt_info_F


   !! Data arrays
   !! -----------

   REAL(8),  DIMENSION(:),    ALLOCATABLE ::   &
      &   vt0, vt          !: time arrays

   !! Arrays on source grid :
   REAL(4), DIMENSION(:,:,:),    ALLOCATABLE ::  &
      &   depth_src,   &
      &   depth_trg,  &
      &   depth_src_trgt2d

   REAL(wpl),  DIMENSION(:,:),    ALLOCATABLE ::  &
      &   bathy_src,  &
      &   bathy_trg

   REAL(wpl),  DIMENSION(:,:),  ALLOCATABLE ::  &
      &   data_src,    &   !: data array on source grid
      &   data_src_b

   REAL(8),  DIMENSION(:,:),  ALLOCATABLE ::  &
      &   lon_src,   &
      &   lat_src

   REAL(wpl),  DIMENSION(:,:,:),  ALLOCATABLE ::  &
      &      data3d_src, &
      &      data3d_trg, &
      &      data3d_tmp, &  ! horizontal target resol. + vertical source resol.
      &    data_src_drowned

   INTEGER(1),   DIMENSION(:,:),  ALLOCATABLE :: IGNORE !: point of target domain to disregard (IGNORE==0)
   INTEGER(1),   DIMENSION(:,:,:),  ALLOCATABLE ::   &
      &   mask_src, mask_src_b, &  !: land-sea mask on input grid
      &   mask_trg               !: land-sea mask on output grid

   !! Arrays on target grid :
   REAL(wpl) , DIMENSION(:,:),  ALLOCATABLE ::             &
      &   data_trg     ! data array interpolated

   REAL(8) , DIMENSION(:,:),  ALLOCATABLE ::             &
      &   lon_trg,   & ! longitude array on target grid
      &   lon_trg_b, & !  //   backup
      &   lat_trg      ! latitude array on target grid

END MODULE MOD_CONF
