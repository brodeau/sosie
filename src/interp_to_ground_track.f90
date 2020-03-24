PROGRAM INTERP_TO_GROUND_TRACK

   USE io_ezcdf
   USE mod_conf
   USE mod_bilin_2d
   USE mod_manip
   USE mod_drown

   !!========================================================================
   !! Purpose :
   !! ---------
   !!
   !! Author :   Laurent Brodeau, 2018
   !! --------
   !!
   !!========================================================================

   IMPLICIT NONE

   !! ************************ Configurable part ****************************
   !!
   LOGICAL, PARAMETER :: &
      &   l_debug_SARAL = .FALSE., &
      &   l_debug = .FALSE., &
      &   l_debug_mapping = .FALSE., &
      &   l_drown_in = .FALSE., & ! Not needed since we ignore points that are less than 1 point away from land... drown the field to avoid spurious values right at the coast!
      &   l_akima = .TRUE., &
      &   l_bilin = .FALSE.
   !!
   REAL(8), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER :: Nt0, Nti, Ntf, iP, jP, iquadran
   !!
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: xlon_gt_0, xlat_gt_0, xlon_gt_i, xlat_gt_i, xlon_gt_f, xlat_gt_f, xdum_r8
   !!


   !! Coupe stuff:
   REAL(8), DIMENSION(:),   ALLOCATABLE :: Ftrack_mod, Ftrack_mod_np, Ftrack_obs, rcycle_obs, vdistance
   REAL(8), DIMENSION(:),   ALLOCATABLE :: vtf, vt_mod, vt_obs, F_gt_0, F_gt_f, rcycle
   REAL(4), DIMENSION(:,:), ALLOCATABLE :: RES_2D_MOD, RES_2D_OBS

   REAL(8),    DIMENSION(:,:,:), ALLOCATABLE :: RAB       !: alpha, beta
   INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: IMETRICS  !: iP, jP, iquadran at each point
   INTEGER(2), DIMENSION(:,:),   ALLOCATABLE :: IPB       !: ID of problem


   !! Grid, default name :
   CHARACTER(len=80) :: &
      &    cv_mod = '', &
      &    cv_obs = '', &
      &    cv_t   = 'time_counter',  &
      &    cv_mt  = 'tmask',         &
      &    cv_z   = 'deptht',        &
      &    cv_lon = 'glamt',         & ! input grid longitude name, T-points
      &    cv_lat = 'gphit'            ! input grid latitude name,  T-points

   CHARACTER(len=256)  :: cr, cunit, cnm_fill
   CHARACTER(len=512)  :: cdum, cconf
   !!
   !!
   !!******************** End of conf for user ********************************
   !!
   !!               ** don't change anything below **
   !!
   LOGICAL ::  &
      &     l_get_mask_metrics_from_meshmask = .FALSE., &
      &     l_mask_input_data = .FALSE., &
      &   l_write_nc_show_track = .FALSE., &
      &     l_exist   = .FALSE., &
      &     l_use_anomaly = .FALSE., &  ! => will transform a SSH into a SLA (SSH - MEAN(SSH))
      &     l_loc1, l_loc2, &
      &     l_obs_ok, lfillval
   !!
   !!
   CHARACTER(len=400)  :: &
      &    cf_obs = 'ground_track.nc', &
      &    cf_mod = '', &
      &    cf_msk = 'mask.nc', &
      &    cf_mm  = 'mesh_mask.nc', &
      &    cf_mpg = ''
   !!
   INTEGER      :: &
      &    jarg,   &
      &    idot, ip1, jp1, im1, jm1, &
      &    i0, j0,  &
      &    ni, nj, nk=0, Ntm=0, Ntdum, &
      &    ni1, nj1, ni2, nj2, &
      &    iargc, id_f1, id_v1
   !!
   !!
   INTEGER :: ji_min, ji_max, jj_min, jj_max, nib, njb

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xdum_r4, show_obs, xvar, xvar1, xvar2, xmean
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: xlont, xlatt
   !!
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJidx, JIidx    ! debug
   !!
   INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: imask, imask_ignr
   INTEGER(1), DIMENSION(:),   ALLOCATABLE :: Fmask, icycle
   !!
   INTEGER :: jt, it, jt0, jtf, jt_s, jtm_1, jtm_2, jtm_1_o, jtm_2_o
   !!
   REAL(8) :: rt, t_min_e, t_max_e, t_min_m, t_max_m, &
      &       alpha, beta, t_min, t_max
   !!
   CHARACTER(LEN=2), DIMENSION(12), PARAMETER :: &
      &            clist_opt = (/ '-h','-v','-x','-y','-z','-t','-i','-p','-n','-m','-S','-M' /)

   REAL(8) :: lon_min_2, lon_max_2, lat_min, lat_max, r_obs
   REAL(4) :: zdt, zdst, rrr, rfillval_mod
   TYPE(t_unit_t0) :: tut_epoch, tut_obs, tut_mod

   INTEGER :: it1, it2

   CHARACTER(80), PARAMETER :: cunit_time_trg = 'seconds since 1970-01-01 00:00:00'

   OPEN(UNIT=6, FORM='FORMATTED', RECL=512)



   !! Epoch is our reference time unit, it is "seconds since 1970-01-01 00:00:00" which translates into:
   tut_epoch%unit   = 's'
   tut_epoch%year   = 1970
   tut_epoch%month  = 1
   tut_epoch%day    = 1
   tut_epoch%hour   = 0
   tut_epoch%minute = 0
   tut_epoch%second = 0

   PRINT *, ''






   !PRINT *, 'Distance Paris - New-York =', DISTANCE(2.35_8, 360._8-74._8, 48.83_8, 40.69_8)


   !! Getting string arguments :
   !! --------------------------
   
   l_get_mask_metrics_from_meshmask = .FALSE.
   jarg = 0

   DO WHILE ( jarg < iargc() )

      jarg = jarg + 1
      CALL getarg(jarg,cr)

      SELECT CASE (trim(cr))

      CASE('-h')
         call usage()

      CASE('-i')
         CALL GET_MY_ARG('input file', cf_mod)

      CASE('-v')
         CALL GET_MY_ARG('model input variable', cv_mod)

      CASE('-x')
         CALL GET_MY_ARG('longitude', cv_lon)

      CASE('-y')
         CALL GET_MY_ARG('latitude', cv_lat)

      CASE('-z')
         CALL GET_MY_ARG('depth', cv_z)

      CASE('-t')
         CALL GET_MY_ARG('time', cv_t)

      CASE('-p')
         CALL GET_MY_ARG('orbit track file', cf_obs)

      CASE('-m')
         l_get_mask_metrics_from_meshmask = .TRUE.
         CALL GET_MY_ARG('mesh_mask file', cf_mm)

      CASE('-S')
         l_write_nc_show_track = .TRUE.

      CASE('-n')
         CALL GET_MY_ARG('ground track input variable', cv_obs)

      CASE('-M')
         l_mask_input_data = .TRUE.
         CALL GET_MY_ARG('masking file', cf_msk)


      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
         CALL usage()

      END SELECT

   END DO

   IF ( (trim(cv_mod) == '').OR.(trim(cf_mod) == '') ) THEN
      PRINT *, ''
      PRINT *, 'You must at least specify input file (-i) and input variable (-v)!!!'
      CALL usage()
   END IF
   
   IF ( TRIM(cv_obs) == '' ) THEN
      PRINT *, ''
      PRINT *, 'You must specify the name of which variable to look at in orbit file ! => -n <name> !!!'
      CALL usage()
   END IF
   
   

   PRINT *, ''
   PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   PRINT *, ''

   PRINT *, ' * Input file = ', trim(cf_mod)
   PRINT *, '   => associated variable names = ', trim(cv_mod)
   PRINT *, '   => associated longitude/latitude/time = ', TRIM(cv_lon), ', ', TRIM(cv_lat), ', ', TRIM(cv_t)
   IF (l_get_mask_metrics_from_meshmask) PRINT *, '   => mesh_mask file = ', TRIM(cf_mm)


   PRINT *, ''

   !! Name of config: lulu
   idot = SCAN(cf_mod, '/', back=.TRUE.)
   cdum = cf_mod(idot+1:)
   idot = SCAN(cdum, '.', back=.TRUE.)
   cconf = cdum(:idot-1)

   idot = SCAN(cf_obs, '/', back=.TRUE.)
   cdum = cf_obs(idot+1:)
   idot = SCAN(cdum, '.', back=.TRUE.)
   cconf = TRIM(cconf)//'__to__'//cdum(:idot-1)


   PRINT *, ' *** CONFIG: cconf ='//TRIM(cconf)



   IF (.NOT.l_get_mask_metrics_from_meshmask) cf_mm = cf_mod


   !! testing longitude and latitude
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CALL DIMS(cf_mm, cv_lon, ni1, nj1, nk, Ntdum)
   CALL DIMS(cf_mm, cv_lat, ni2, nj2, nk, Ntdum)

   IF ( (nj1==-1).AND.(nj2==-1) ) THEN
      ni = ni1 ; nj = ni2
      PRINT *, 'Grid is 1D: ni, nj =', ni, nj
      l_reg_src = .TRUE.
   ELSE
      IF ( (ni1==ni2).AND.(nj1==nj2) ) THEN
         ni = ni1 ; nj = nj1
         PRINT *, 'Grid is 2D: ni, nj =', ni, nj
         l_reg_src = .FALSE.
      ELSE
         PRINT *, 'ERROR: problem with grid!' ; STOP
      END IF
   END IF

   PRINT *, ''
   

   !! testing variable dimensions
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CALL DIMS(cf_mod, cv_mod, ni1, nj1, nk, Ntm)

   IF ( (ni1/=ni).AND.(nj1/=nj) ) THEN
      PRINT *, 'ERROR: dimension of ',trim(cv_mod), 'does not agree with lon/lat' ; STOP
   END IF

   IF ( nk < 1 ) nk = 1

   IF ( Ntm < 1 ) THEN
      PRINT *, 'ERROR: ',trim(cv_mod),' must have at least a time record!' ; STOP
   END IF


   PRINT *, 'Dimension for ',trim(cv_mod),':'
   PRINT *, '   => ni =', ni ;   PRINT *, '   => nj =', nj
   PRINT *, '   => nk =', nk ;   PRINT *, '   => Ntm =', Ntm
   PRINT *, ''


   PRINT *, ''
   PRINT *, ' *** Allocating ni x nj arrays...'
   ALLOCATE ( xlont(ni,nj), xlatt(ni,nj), xdum_r4(ni,nj), &
      &       xvar(ni,nj), xvar1(ni,nj), xvar2(ni,nj),    &
      &       imask(ni,nj), vt_mod(Ntm) )
   IF ( l_debug ) THEN
      ALLOCATE ( RES_2D_MOD(ni,nj), RES_2D_OBS(ni,nj) )
      RES_2D_MOD(:,:) = 0.
      RES_2D_OBS(:,:) = 0.
   END IF
   PRINT *, ' *** Done!'; PRINT *, ''


   !! In case we mask input data with field 'mask' found into file 'cf_msk' (option: "-M")
   IF ( l_mask_input_data ) THEN
      CALL DIMS(cf_msk, 'mask', ni1, nj1, nk, Ntdum)
      IF ( (ni1/=ni).OR.(nj1/=nj) ) THEN
         PRINT *, 'ERROR: shape of mask not consistent with your domain!', ni1,ni, nj1,nj ; STOP
      END IF
      ALLOCATE ( imask_ignr(ni,nj) )
      PRINT *, 'Reading mask "mask" into file '//TRIM(cf_msk)//' !'
      CALL GETVAR_2D   (i0, j0, cf_msk,  'mask', 0, 0, 0, xdum_r4)
      imask_ignr(:,:) = INT(xdum_r4, 1) ; i0=0 ; j0=0
   END IF


   IF ( l_reg_src ) THEN
      PRINT *, 'Regular case not supported yet! Priority to ORCA grids...'
      STOP
   END IF




   !! The first important stage is to compare time slices in the OGCM 2D input field
   !! w.r.t the one in the track orbit file!
   !! Since we are dealing with satellite data, useing the UNIX "epoch" time seems
   !! appropriate:
   !!
   !! The Unix epoch (or Unix time or POSIX time or Unix timestamp) is the
   !! number of seconds that have elapsed since January 1, 1970 (midnight
   !! UTC/GMT), not counting leap seconds (in ISO 8601: 1970-01-01T00:00:00Z).
   !!
   !! As such, step #1 is to convert the time vector in both files to epoch time
   !! If there is no overlapping period of time between the two file, then no
   !! need to go further...
   !!
   CALL GET_VAR_INFO(cf_mod, cv_t, cunit, cdum)
   tut_mod  = GET_TIME_UNIT_T0(TRIM(cunit))
   PRINT *, ' *** Unit and reference time in model file:'
   PRINT *, tut_mod

   CALL GET_VAR_INFO(cf_obs, 'time', cunit, cdum)
   tut_obs  = GET_TIME_UNIT_T0(TRIM(cunit))
   PRINT *, ' *** Unit and reference time in track file:'
   PRINT *, tut_obs
   PRINT *, ''





   !! Getting coordinates
   !! ~~~~~~~~~~~~~~~~~~~

   !!IF ( nk > 1 ) CALL GETVAR_1D(cf_mod, cv_z, vdepth(:,1))


   !! Reading it in input file:
   CALL GETVAR_1D(cf_mod, cv_t, vt_mod)


   !! Getting longitude, latitude and mask in mesh_mask file:
   ! Longitude array:
   CALL GETVAR_2D (i0, j0, cf_mm,  cv_lon, 0, 0, 0, xlont(:,:))  ; i0=0 ; j0=0
   !xlont(:,:) = xdum_r4(:,:) ; i0=0 ; j0=0
   !!


   !! Min an max lon:
   lon_min_1 = MINVAL(xlont)
   lon_max_1 = MAXVAL(xlont)
   PRINT *, ' *** Minimum longitude on source domain before : ', lon_min_1
   PRINT *, ' *** Maximum longitude on source domain before : ', lon_max_1
   !!
   WHERE ( xdum_r4 < 0. ) xlont = xlont + 360.0_8
   !!
   lon_min_2 = MINVAL(xlont)
   lon_max_2 = MAXVAL(xlont)
   PRINT *, ' *** Minimum longitude on source domain: ', lon_min_2
   PRINT *, ' *** Maximum longitude on source domain: ', lon_max_2

   ! lolo: disgusting!:
   !l_loc1 = (lon_min_1 <  0.).AND.(lon_min_1 > -170.).AND.(lon_max_1 >  0. ).AND.(lon_min_1 <  170.)
   l_loc1 = ((lon_min_1 <  0.).AND.(lon_min_1 > -170.)).OR.((lon_max_1 >  0. ).AND.(lon_min_1 <  170.))
   l_loc2 = (lon_min_2 >= 0.).AND.(lon_min_2 <   2.5).AND.(lon_max_2 >357.5).AND.(lon_max_2 <= 360.)
   IF (.NOT. l_loc1) THEN
      IF ( l_loc2 ) THEN
         l_glob_lon_wize = .TRUE.
         PRINT *, 'Looks like global setup (longitude-wise at least...)'
      ELSE
         PRINT *, 'ERROR: cannot find if regional or global source domain...'; STOP
      END IF
   ELSE
      PRINT *, 'Looks like regional setup (longitude-wise at least...)'
      l_glob_lon_wize = .FALSE.
      !!
      WRITE(*,'("  => going to disregard points of target domain with lon < ",f7.2," and lon > ",f7.2)'), lon_min_1,lon_max_1
   END IF
   PRINT *, ''




   ! Latitude array:
   CALL GETVAR_2D   (i0, j0, cf_mm,  cv_lat, 0, 0, 0, xlatt(:,:)) ; i0=0 ; j0=0
   !xlatt(:,:) = xdum_r4(:,:) ; i0=0 ; j0=0

   !! Min an max lat:
   lat_min = MINVAL(xlatt)
   lat_max = MAXVAL(xlatt)
   PRINT *, ' *** Minimum latitude on source domain : ', lat_min
   PRINT *, ' *** Maximum latitude on source domain : ', lat_max


   l_glob_lat_wize = .TRUE.
   IF ( lat_max < 88. ) THEN
      l_glob_lat_wize =.FALSE.
      WRITE(*,'("  => going to disregard points of target domain with lat < ",f7.2," and lat > ",f7.2)'), lat_min,lat_max
   END IF
   PRINT *, ''



   !! 2D land-sea mask on model grid:
   IF (l_get_mask_metrics_from_meshmask) THEN
      CALL GETMASK_2D(cf_mm, cv_mt, imask, jlev=1)
   ELSE
      !! getting mask from _FillValue on first field of file:
      CALL CHECK_4_MISS(cf_mod, cv_mod, lfillval, rfillval_mod, cnm_fill)
      CALL GETVAR_2D(i0, j0, cf_mod, cv_mod, Ntm, 0, 1, xdum_r4, jt1=1, jt2=1)
      imask(:,:) = 1
      WHERE ( xdum_r4 == rfillval_mod ) imask = 0
      i0=0 ; j0=0
   END IF
   
   IF ( l_mask_input_data ) THEN
      ! taking into consideration the forced masked region 'imask_ignr'
      WHERE ( imask_ignr(:,:) == 0 ) imask(:,:) = 0
      DEALLOCATE ( imask_ignr )
   END IF
   
   !CALL DUMP_FIELD(REAL(imask), 'mask_in.nc', 'lsm') !, xlont, xlatt, 'nav_lon', 'nav_lat', rfill=-9999.)









   !! Reading along-track from file:

   INQUIRE(FILE=TRIM(cf_obs), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: please provide the file containing definition of orbit track'; STOP
   END IF

   CALL DIMS(cf_obs, 'time', Nt0, nj1, nk, ni1)
   PRINT *, ' *** Nb. time records in NetCDF track file:', Nt0
   ALLOCATE ( xlon_gt_0(1,Nt0), xlat_gt_0(1,Nt0), vt_obs(Nt0), F_gt_0(Nt0), rcycle(Nt0) )
  
   CALL GETVAR_1D(cf_obs, 'time', vt_obs)
   CALL GETVAR_1D(cf_obs, 'longitude', xlon_gt_0(1,:))
   CALL GETVAR_1D(cf_obs, 'latitude',  xlat_gt_0(1,:))
   
   CALL GETVAR_1D(cf_obs, 'cycle',     rcycle)

   CALL GETVAR_1D(cf_obs, cv_obs,      F_gt_0)

   PRINT *, 'Done!'; PRINT *, ''

   !IF ( l_debug_SARAL ) THEN
   !   OPEN(16, FILE='debug_saral.txt', FORM='FORMATTED', RECL=256, STATUS='unknown')
   !   WRITE(16,*) '#     Fucked-up points! '
   !   WRITE(16,*) '# time rec. in file  | time (d since 1950)  |   dt (s)     | dl (km)        | speed (km/s)'
   !   DO jt = 2, Nt0
   !      zdt = (vt_obs(jt) - vt_obs(jt-1))*3600.*24. ! dt
   !      zdst = DISTANCE( xlon_gt_0(1,jt), xlon_gt_0(1,jt-1), xlat_gt_0(1,jt), xlat_gt_0(1,jt-1) )! dl
   !      rrr = zdst/zdt
   !      IF (rrr > 8.) WRITE(16,*) jt, vt_obs(jt), zdt, zdst, zdst/zdt
   !   END DO
   !   CLOSE(16)
   !   STOP
   !END IF





   nib = ni ; njb = nj ; ji_min=1 ; ji_max=ni ; jj_min=1 ; jj_max=nj

   PRINT *, ''
   PRINT *, ' Time vector in track file:'
   CALL to_epoch_time_vect( tut_obs, vt_obs, l_dt_below_sec=.true. )
   PRINT *, '' ;  PRINT *, ''

   !! Defaults:
   Nti = Nt0
   it1  = 1
   it2  = Nt0

   IF ( .NOT. l_debug_mapping ) THEN

      PRINT *, ' Time vector in model file:'
      CALL to_epoch_time_vect( tut_mod, vt_mod, l_dt_below_sec=.FALSE. )
      !PRINT *, vt_mod(:)
      PRINT *, ''

      IF ( l_use_anomaly ) THEN
         PRINT *, ''
         PRINT *, ' *** Computing mean of field '//TRIM(cv_mod)//' for considering anomaly later...'
         ALLOCATE ( xmean(ni,nj) )
         xmean = 0.
         DO jt = 1, Ntm
            CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jt, xdum_r4)
            xmean = xmean + xdum_r4/REAL(Ntm,4)
         END DO
         id_f1=0 ;  id_v1=0
         !WHERE ( imask == 0 ) xmean = -9999.
         CALL DUMP_FIELD(xmean, 'mean_'//TRIM(cv_mod)//'.nc', cv_mod, xlont, xlatt, 'nav_lon', 'nav_lat', rfill=-9999.)
         !STOP 'lolo'
      END IF

      t_min_e = MINVAL(vt_obs)
      t_max_e = MAXVAL(vt_obs)
      t_min_m = MINVAL(vt_mod)
      t_max_m = MAXVAL(vt_mod)

      PRINT *, ''
      PRINT *, ' *** Max min time for track:', t_min_e, t_max_e
      PRINT *, ' *** Max min time for model:', t_min_m, t_max_m
      PRINT *, ''

      IF ( (t_min_m >= t_max_e).OR.(t_min_e >= t_max_m).OR.(t_max_m <= t_min_e).OR.(t_max_e <= t_min_m) ) THEN
         PRINT *, ' No time overlap for Model and Track file! '
         STOP
      END IF

      t_min = MAX(t_min_e, t_min_m)
      t_max = MIN(t_max_e, t_max_m)
      PRINT *, ' *** Time overlap for Model and Track file:', t_min, t_max


      !! Finding when we can start and stop when scanning the track file:
      !! it1, it2
      DO it1 = 1, Nt0-1
         IF ( (vt_obs(it1) <= t_min).AND.(vt_obs(it1+1) > t_min) ) EXIT
      END DO
      DO it2 = it1, Nt0-1
         IF ( (vt_obs(it2) <= t_max).AND.(vt_obs(it2+1) > t_max) ) EXIT
      END DO

      Nti = it2 - it1 + 1

      PRINT *, ' it1, it2 =',it1, it2
      PRINT *, Nti, '  out of ', Nt0
      PRINT *, ' => ', vt_obs(it1), vt_obs(it2)
      PRINT *, ''
   END IF ! IF ( .NOT. l_debug_mapping )

   ALLOCATE ( IGNORE(1,Nti), xlon_gt_i(1,Nti), xlat_gt_i(1,Nti) )

   IGNORE(:,:) = 1 !lolo

   !! Main time loop is on time vector in track file!


   cf_mpg = 'MAPPING__'//TRIM(cconf)//'.nc'


   !!
   xlon_gt_i(:,:) = xlon_gt_0(:,it1:it2)
   xlat_gt_i(:,:) = xlat_gt_0(:,it1:it2)

   DEALLOCATE ( xlon_gt_0, xlat_gt_0 )

   IF ( .NOT. l_glob_lon_wize ) THEN
      ALLOCATE ( xdum_r8(1,Nti) )
      xdum_r8 = degE_to_degWE(xlon_gt_i)
      
      !xdum_r8 = SIGN(1._8, 180._8-xlon_gt_i)*MIN(xlon_gt_i,ABS(xlon_gt_i-360._8)) ! like xlon_gt_i but between -180 and +180 !
      WHERE ( xdum_r8 < lon_min_1 ) IGNORE=0
      WHERE ( xdum_r8 > lon_max_1 ) IGNORE=0
      DEALLOCATE ( xdum_r8 )
   END IF

   IF ( .NOT. l_glob_lat_wize ) THEN
      WHERE ( xlat_gt_i < lat_min ) IGNORE=0
      WHERE ( xlat_gt_i > lat_max ) IGNORE=0
   END IF


   !! We are going to shorten our 1D input arrays, only keeping values included
   !! in target domain (i.e. where IGNORE==1):

   PRINT *, ''
   PRINT *, ' Intially we had Nt0 ', Nt0, ' points'
   PRINT *, ' - excluding non relevant time led to Nti', Nti, ' points', SIZE(IGNORE(1,:),1)

   !Ntf = SUM(INT4(IGNORE)) !lolo wtf Gfortran ???
   Ntf = SUM(INT(IGNORE))
   PRINT *, ' - and in the end we only retain Ntf ', Ntf , ' points!'
   PRINT *, ''



   !! Allocate arrays on the final retained size
   ALLOCATE ( IMETRICS(1,Ntf,3), RAB(1,Ntf,2), IPB(1,Ntf), xlon_gt_f(1,Ntf), xlat_gt_f(1,Ntf), vtf(Ntf), F_gt_f(Ntf), icycle(Ntf) )
   ALLOCATE ( Ftrack_mod(Ntf), Ftrack_obs(Ntf), Fmask(Ntf), Ftrack_mod_np(Ntf), rcycle_obs(Ntf), vdistance(Ntf) )

   xlon_gt_f(1,:) = SHRINK_VECTOR(xlon_gt_i(1,:),  IGNORE(1,:), Ntf)
   xlat_gt_f(1,:) = SHRINK_VECTOR(xlat_gt_i(1,:),  IGNORE(1,:), Ntf)
   vtf(:)         = SHRINK_VECTOR(rcycle(it1:it2), IGNORE(1,:), Ntf)
   icycle = INT2(vtf)
   vtf(:)         = SHRINK_VECTOR(vt_obs(it1:it2), IGNORE(1,:), Ntf) !
   F_gt_f(:)      = SHRINK_VECTOR(F_gt_0(it1:it2), IGNORE(1,:), Ntf)

   ! 
   DEALLOCATE ( xlon_gt_i , xlat_gt_i , vt_obs , F_gt_0 , rcycle )


   
   IF ( l_debug_SARAL ) THEN
      OPEN(16, FILE='debug_saral.txt', FORM='FORMATTED', RECL=256, STATUS='unknown')
      !WRITE(16,*) '#     Fucked-up points! '
      WRITE(16,*) '# time rec. in file  | time (d since 1950)  |   dt (s)     | dl (km)        | speed (km/s)'
      DO jt = 2, Ntf
         it  = it1-1+jt
         zdt = (vtf(it) - vtf(it-1))*3600.*24. ! dt
         zdst = DISTANCE( xlon_gt_f(1,jt), xlon_gt_f(1,jt-1), xlat_gt_f(1,jt), xlat_gt_f(1,jt-1) )! dl
         rrr = zdst/zdt
         !IF (rrr > 8.) WRITE(16,*) jt, vt_obs(jt), zdt, zdst, zdst/zdt
         WRITE(16,*) it, vtf(it), zdt, zdst, zdst/zdt
      END DO
      CLOSE(16)
      STOP 'You are in SARAL debug mode (l_debug_SARAL=.TRUE.), we stop here! An check file "debug_saral.txt"!!!'
   END IF

   




   INQUIRE(FILE=trim(cf_mpg), EXIST=l_exist ) !
   IF ( .NOT. l_exist ) THEN
      PRINT *, ' *** Creating mapping file...' !
      CALL MAPPING_BL(-1, xlont, xlatt, xlon_gt_f, xlat_gt_f, cf_mpg ) !,  mask_domain_trg=IGNORE) don't need ignore, points have been removed!
      PRINT *, ' *** Done!'; PRINT *, ''
   ELSE
      PRINT *, ' *** File "',trim(cf_mpg),'" found in current directory, using it!'
      PRINT *, ''
   END IF

   CALL RD_MAPPING_AB(cf_mpg, IMETRICS, RAB, IPB)
   PRINT *, ''; PRINT *, ' *** Mapping and weights read into "',trim(cf_mpg),'"'; PRINT *, ''

   ALLOCATE (JIidx(1,Ntf) , JJidx(1,Ntf) )
   JIidx(1,:) = IMETRICS(1,:,1)
   JJidx(1,:) = IMETRICS(1,:,2)



   !! Showing iy in file mask_+_nearest_points.nc:
   IF ( l_write_nc_show_track ) THEN
      !! Finding and storing the nearest points of NEMO grid to track points:
      !CALL FIND_NEAREST_POINT(xlon_gt_0, xlat_gt_0, xlont, xlatt,  JIidx, JJidx)
      ALLOCATE ( show_obs(nib,njb) )
      show_obs(:,:) = -9999.
      DO jtf = 1, Ntf
         IF ( (JIidx(1,jtf)>0).AND.(JJidx(1,jtf)>0) )  show_obs(JIidx(1,jtf), JJidx(1,jtf)) = REAL(jtf,4)
      END DO
      WHERE (imask == 0) show_obs = -100.
      CALL DUMP_FIELD(REAL(show_obs(:,:),4), 'mask_+_nearest_points__'//TRIM(cconf)//'.nc', 'mask', xlont, xlatt, cv_lon, cv_lat, rfill=-9999.)
      !lolo:
      !CALL DUMP_FIELD(REAL(xlont(:,:),4), 'lon_360.nc', 'lon')
      !show_obs = SIGN(1.,180.-xlont)*MIN(xlont,ABS(xlont-360.))
      !CALL DUMP_FIELD(REAL(show_obs(:,:),4), 'lon_-180-180.nc', 'lon')
      !WHERE ( (show_obs > 10.).OR.(show_obs < -90.) ) show_obs = -800.
      !CALL DUMP_FIELD(REAL(show_obs(:,:),4), 'lon_masked.nc', 'lon')
      !STOP 'interp_to_ground_obs.f90'
      !lolo.

      DEALLOCATE ( show_obs )
   END IF

   IF ( l_debug_mapping ) STOP 'l_debug_mapping'



   !STOP 'mapping done!'



   Ftrack_mod_np(:) = -9999.
   Ftrack_obs(:)    = -9999.
   Ftrack_mod(:)    = -9999.
   Fmask(:)         = 0
   rcycle_obs(:)    = -9999.

   jtm_1_o = -100
   jtm_2_o = -100
   jt_s    = 1

   DO jtf = 1, Ntf

      rt = vtf(jtf)

      IF ( (rt >= t_min_m).AND.(rt < t_max_m) ) THEN

         !! Two surrounding time records in model file => jtm_1 & jtm_2
         DO jt=jt_s, Ntm-1
            IF ( (rt >= vt_mod(jt)).AND.(rt < vt_mod(jt+1)) ) EXIT
         END DO
         !!
         jtm_1 = jt
         jtm_2 = jt+1
         IF (jtf==1) jt0 = jtm_1 ! Saving the actual first useful time step of the model!

         PRINT *, ' * Track time =>', rt, '/ model jtm_1,jtm_2 =', INT2(jtm_1), INT2(jtm_2)

         !! If first time we have these jtm_1 & jtm_2, getting the two surrounding fields:
         IF ( (jtm_1>jtm_1_o).AND.(jtm_2>jtm_2_o) ) THEN
            IF ( (jtm_1_o == -100).OR.(jtm_1 > jtm_2_o) ) THEN
               PRINT *, ' *** Reading field '//TRIM(cv_mod)//' in '//TRIM(cf_mod)
               PRINT *, '    => at jtm_1=', jtm_1, '  (starting from jt1=',jt0,')'
               CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jtm_1, xvar1, jt1=jt0)
               IF ( l_use_anomaly ) xvar1 = xvar1 - xmean
               IF ( l_drown_in ) CALL DROWN(-1, xvar1, imask, nb_inc=5)
            ELSE
               xvar1(:,:) = xvar2(:,:)
            END IF
            PRINT *, ' *** Reading field '//TRIM(cv_mod)//' in '//TRIM(cf_mod)
            PRINT *, '    => at jtm_2=', jtm_2, '  (starting from jt1=',jt0,')'
            CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jtm_2, xvar2, jt1=jt0)
            IF ( l_use_anomaly ) xvar2 = xvar2 - xmean
            IF ( l_drown_in ) CALL DROWN(-1, xvar2, imask, nb_inc=5)

            xdum_r4 = (xvar2 - xvar1) / (vt_mod(jtm_2) - vt_mod(jtm_1)) ! xdum_r4 is the slope here !!!

            PRINT *, ''
         END IF

         !! Linear interpolation of field at time rt:
         xvar(:,:) = xvar1(:,:) + xdum_r4(:,:)*(rt - vt_mod(jtm_1))

         !! Performing bilinear interpolation:
         iP       = IMETRICS(1,jtf,1)
         jP       = IMETRICS(1,jtf,2)
         iquadran = IMETRICS(1,jtf,3)
         alpha    = RAB(1,jtf,1)
         beta     = RAB(1,jtf,2)

         !LOLO: IF ( (iP/=INT(rflg)).AND.(jP/=INT(rflg)) ) THEN
         IF ( (iP>0).AND.(jP>0) ) THEN
            IF ( imask(iP,jP)==1 ) THEN
               r_obs    = F_gt_f(jtf)
               l_obs_ok = ( r_obs > -20.).AND.( r_obs < 20.)
               IF ( l_obs_ok ) THEN
                  !! Ignore points that are just 1 point away from land:
                  ip1 = MIN(iP+1,ni) ; jp1 = MIN(jP+1,nj)
                  im1 = MAX(iP-1,1)  ; jm1 = MAX(jP-1,1)
                  idot = imask(ip1,jP) + imask(ip1,jp1) + imask(iP,jp1) + imask(im1,jp1) &
                     & + imask(im1,jP) + imask(im1,jm1) + imask(iP,jm1) + imask(ip1,jm1)
                  !! => idot == 8 if in that case...
                  !!
                  IF (idot==8) THEN
                     !! Model, nearest point:
                     Ftrack_mod_np(jtf) =  xvar(JIidx(1,jtf),JJidx(1,jtf)) ! NEAREST POINT interpolation
                     !! Model, 2D bilinear interpolation:
                     Ftrack_mod(jtf) = REAL( INTERP_BL(-1, iP, jP, iquadran, alpha, beta, xvar) , 8)
                     !! Observations as on their original point:
                     Ftrack_obs(jtf) = r_obs
                     IF ( l_debug ) THEN
                        !! On the model grid for info:
                        RES_2D_MOD(iP,jP) = Ftrack_mod(jtf)
                        RES_2D_OBS(iP,jP) = Ftrack_obs(jtf)
                     END IF
                     !!
                     Fmask(jtf) = 1 ! That was a valid point!
                     !!
                     rcycle_obs(jtf) = REAL( icycle(jtf), 8 )
                     !!
                  END IF
               END IF
            END IF
         END IF

         jtm_1_o = jtm_1
         jtm_2_o = jtm_2
         jt_s    = jtm_1 ! so we do not rescan from begining...

      END IF

   END DO


   !! Vector distance (in km)
   vdistance(:) = 0.
   DO jt = 2, Ntf
      IF ( (Fmask(jt)==1).AND.(Fmask(jt-1)==1) ) vdistance(jt) = vdistance(jt-1) + DISTANCE( xlon_gt_f(1,jt), xlon_gt_f(1,jt-1), xlat_gt_f(1,jt), xlat_gt_f(1,jt-1) )
   END DO

   WHERE ( Fmask == 0 )
      Ftrack_mod    = -9999.
      Ftrack_mod_np = -9999.
      Ftrack_obs    = -9999.
      rcycle_obs    = -9999.
      vdistance     = -9999.
   END WHERE

   WHERE ( Ftrack_mod < -9990. ) Ftrack_mod = -9999.


   PRINT *, ''
   !WRITE(cf_out, '("track_",a,"_",a,".nc")') TRIM(cv_mod), TRIM(cf_obs)
   cf_out = 'result__'//TRIM(cconf)//'.nc'
   PRINT *, ' * Output file = ', trim(cf_out)
   PRINT *, ''

   CALL PT_SERIES(vtf(:), REAL(Ftrack_mod,4), cf_out, 'time', cv_mod, 'm', 'Model data, bi-linear interpolation', -9999., &
      &           ct_unit=TRIM(cunit_time_trg), &
      &           vdt2=REAL(Ftrack_mod_np,4),cv_dt2=TRIM(cv_mod)//'_np',cln2='Model data, nearest-point interpolation', &
      &           vdt3=REAL(Ftrack_obs,4),   cv_dt3=cv_obs,             cln3='Original data as in track file...',   &
      &           vdt4=REAL(xlon_gt_f(1,:),4), cv_dt4='longitude',        cln4='Longitude (as in track file)',  &
      &           vdt5=REAL(xlat_gt_f(1,:),4), cv_dt5='latitude',         cln5='Latitude (as in track file)' ,  &
      &           vdt6=REAL(Fmask,4),          cv_dt6='mask',             cln6='Mask', &
      &           vdt7=REAL(rcycle_obs,4),     cv_dt7='cycle',            cln7='cycle', &
      &           vdt8=REAL(vdistance,4),      cv_dt8='distance',         cln8='Distance (in km) from first point of segment' )

   IF ( l_debug ) THEN
      WHERE ( imask == 0 )
         RES_2D_MOD = -9999.
         RES_2D_OBS = -9999.
      END WHERE
      CALL DUMP_FIELD(RES_2D_MOD, 'RES_2D_MOD__'//TRIM(cconf)//'.nc', cv_mod, xlont, xlatt, 'nav_lon', 'nav_lat', rfill=-9999.)
      CALL DUMP_FIELD(RES_2D_OBS, 'RES_2D_OBS__'//TRIM(cconf)//'.nc', cv_obs, xlont, xlatt, 'nav_lon', 'nav_lat', rfill=-9999.)
   END IF

   !IF ( l_debug ) DEALLOCATE ( JIidx, JJidx )
   !DEALLOCATE ( F_gt_0 )
   !DEALLOCATE ( Ftrack_mod, Ftrack_mod_np, Ftrack_obs )
   !DEALLOCATE ( xlont, xlatt, xvar, xvar1, xvar2, xdum_r4, imask, xdum_r8 )


   PRINT *, ''
   PRINT *, 'Written!'
   PRINT *, ' => check:'
   PRINT *, TRIM(cf_out)
   PRINT *, ''

   CLOSE(6)

CONTAINS






   SUBROUTINE GET_MY_ARG(cname, cvalue)
      CHARACTER(len=*), INTENT(in)    :: cname
      CHARACTER(len=*), INTENT(inout) :: cvalue
      !!
      IF ( jarg + 1 > iargc() ) THEN
         PRINT *, 'ERROR: Missing ',trim(cname),' name!' ; call usage()
      ELSE
         jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         IF ( ANY(clist_opt == trim(cr)) ) THEN
            PRINT *, 'ERROR: Missing',trim(cname),' name!'; call usage()
         ELSE
            cvalue = trim(cr)
         END IF
      END IF
   END SUBROUTINE GET_MY_ARG




   SUBROUTINE usage()
      !!
      !OPEN(UNIT=6, FORM='FORMATTED', RECL=512)
      !!
      WRITE(6,*) ''
      WRITE(6,*) '   List of command line options:'
      WRITE(6,*) '   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(6,*) ''
      WRITE(6,*) ' -i <input_file.nc>   => INPUTE FILE'
      WRITE(6,*) ''
      WRITE(6,*) ' -v  <name>           => Specify variable name in input file'
      WRITE(6,*) ''
      WRITE(6,*) ' -p  <track_file>     => Specify name of NetCDF file containing orbit tack'
      WRITE(6,*) ''
      WRITE(6,*) ' -n  <name>           => name of variable of interest in orbit tack file'
      WRITE(6,*) ''
      !!
      WRITE(6,*) ''
      WRITE(6,*) '    Optional:'
      WRITE(6,*)  ''
      WRITE(6,*) ' -x  <name>           => Specify longitude name in input file (default: lon)'
      WRITE(6,*) ''
      WRITE(6,*) ' -y  <name>           => Specify latitude  name in input file  (default: lat)'
      WRITE(6,*) ''
      WRITE(6,*) ' -z  <name>           => Specify depth name in input file (default: depth)'
      WRITE(6,*) ''
      WRITE(6,*) ' -t  <name>           => Specify time name in input file (default: time)'
      WRITE(6,*) ''
      WRITE(6,*) ' -m  <mesh_mask_file> => Specify mesh_mask file to be used (default: mesh_mask.nc)'
      WRITE(6,*) ''
      WRITE(6,*) ' -S                => dump boxes on 2D output field "mask_+_nearest_points.nc" '
      WRITE(6,*) ''
      WRITE(6,*) ' -M <masking_file> => ignore regions of input field where field "mask"==0 in "masking_file"'
      WRITE(6,*) ''
      WRITE(6,*) ' -h                   => Show this message'
      WRITE(6,*) ''
      !!
      CLOSE(6)
      STOP
      !!
   END SUBROUTINE usage


END PROGRAM INTERP_TO_GROUND_TRACK
