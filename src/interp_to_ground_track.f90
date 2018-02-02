PROGRAM INTERP_TO_GROUND_TRACK

   USE io_ezcdf
   USE mod_conf
   USE mod_bilin_2d
   USE mod_manip
   USE mod_drown

   !!========================================================================
   !! Purpose :
   !!
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
      &   l_debug = .TRUE., &
      &   l_debug_mapping = .FALSE., &
      &   l_drown_in = .FALSE., & ! Not needed since we ignore points that are less than 1 point away from land... drown the field to avoid spurious values right at the coast!
      &   l_akima = .TRUE., &
      &   l_bilin = .FALSE.
   !!
   LOGICAL :: &
      &      l_orbit_file_is_nc    = .FALSE.
   !!
   REAL(8), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER :: Nte, Nten, io, idx, iP, jP, iquadran
   !!
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xgt, Ygt, Fgt, xlon_gt, xlat_gt, xdum_r8
   !!
   !! Coupe stuff:
   REAL(8), DIMENSION(:), ALLOCATABLE :: Ftrack_mod, Ftrack_mod_np, Ftrack_obs

   REAL(8), DIMENSION(:),     ALLOCATABLE :: vte, vt_mod, vt_obs   ! in seconds

   REAL(4), DIMENSION(:,:),   ALLOCATABLE :: RES_2D_MOD, RES_2D_OBS

   REAL(8),    DIMENSION(:,:,:), ALLOCATABLE :: RAB       !: alpha, beta
   INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: IMETRICS  !: iP, jP, iquadran at each point
   INTEGER,    DIMENSION(:,:),   ALLOCATABLE :: IPB       !: ID of problem

   !! Grid, default name :
   CHARACTER(len=80) :: &
      &    cv_mod, &
      &    cv_obs, &
      &    cv_t   = 'time_counter',  &
      &    cv_mt  = 'tmask',         &
      &    cv_z   = 'deptht',        &
      &    cv_lon = 'glamt',         & ! input grid longitude name, T-points
      &    cv_lat = 'gphit'            ! input grid latitude name,  T-points

   CHARACTER(len=256)  :: cr, cunit
   CHARACTER(len=512)  :: cdum, cconf
   !!
   !!
   !!******************** End of conf for user ********************************
   !!
   !!               ** don't change anything below **
   !!
   LOGICAL ::  &
      &     l_exist   = .FALSE., &
      &     l_use_anomaly = .FALSE., &  ! => will transform a SSH into a SLA (SSH - MEAN(SSH))
      &     l_loc1, l_loc2, &
      &     l_obs_ok
   !!
   !!
   CHARACTER(len=400)  :: &
      &    cf_obs   = 'ground_track.nc', &
      &    cf_mod, &
      &    cf_mm='mesh_mask.nc', &
      &    cf_mapping, &
      &    cs_force_tv_m='', &
      &    cs_force_tv_e=''
   !!
   INTEGER      :: &
      &    jarg,   &
      &    idot, ip1, jp1, im1, jm1, &
      &    i0, j0,  &
      &    ni, nj, nk=0, Ntm=0, &
      &    ni1, nj1, ni2, nj2, &
      &    iargc, id_f1, id_v1
   !!
   !!
   INTEGER :: ji_min, ji_max, jj_min, jj_max, nib, njb

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xvar, xvar1, xvar2, xmean

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xdum_r4, show_obs
   REAL(8), DIMENSION(:,:), ALLOCATABLE ::    &
      &    xlont, xlatt
   !!
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJidx, JIidx    ! debug
   !!
   INTEGER(2), DIMENSION(:,:), ALLOCATABLE :: mask
   INTEGER(2), DIMENSION(:),   ALLOCATABLE :: Fmask
   !!
   INTEGER :: jt, jte, jt_s, jtm_1, jtm_2, jtm_1_o, jtm_2_o
   !!
   REAL(8) :: rt, rt0, rdt, &
      &       t_min_e, t_max_e, t_min_m, t_max_m, &
      &       alpha, beta, t_min, t_max
   !!
   CHARACTER(LEN=2), DIMENSION(12), PARAMETER :: &
      &            clist_opt = (/ '-h','-v','-x','-y','-z','-t','-i','-p','-n','-m','-f','-g' /)

   REAL(8) :: lon_min_2, lon_max_2, lat_min, lat_max, r_obs

   TYPE(t_unit_t0) :: tut_epoch, tut_obs, tut_mod

   INTEGER :: it1, it2

   CHARACTER(80), PARAMETER :: cunit_time_out = 'seconds since 1970-01-01 00:00:00'

   !! Epoch is our reference time unit, it is "seconds since 1970-01-01 00:00:00" which translates into:
   tut_epoch%unit   = 's'
   tut_epoch%year   = 1970
   tut_epoch%month  = 1
   tut_epoch%day    = 1
   tut_epoch%hour   = 0
   tut_epoch%minute = 0
   tut_epoch%second = 0

   PRINT *, ''


   !! Getting string arguments :
   !! --------------------------

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
         CALL GET_MY_ARG('mesh_mask file', cf_mm)

      CASE('-f')
         CALL GET_MY_ARG('forced time vector construction for model', cs_force_tv_m)

      CASE('-g')
         CALL GET_MY_ARG('forced time vector construction for track', cs_force_tv_e)

      CASE('-n')
         l_orbit_file_is_nc = .TRUE.
         CALL GET_MY_ARG('ground track input variable', cv_obs)

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

   PRINT *, ''
   PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   PRINT *, ''

   PRINT *, ' * Input file = ', trim(cf_mod)
   PRINT *, '   => associated variable names = ', trim(cv_mod)
   PRINT *, '   => associated longitude/latitude/time = ', trim(cv_lon), ', ', trim(cv_lat), ', ', trim(cv_t)
   PRINT *, '   => mesh_mask file = ', trim(cf_mm)


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


   !! testing longitude and latitude
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CALL DIMS(cf_mm, cv_lon, ni1, nj1, nk, Ntm)
   CALL DIMS(cf_mm, cv_lat, ni2, nj2, nk, Ntm)

   IF ( (nj1==-1).AND.(nj2==-1) ) THEN
      ni = ni1 ; nj = ni2
      PRINT *, 'Grid is 1D: ni, nj =', ni, nj
      lregin = .TRUE.
   ELSE
      IF ( (ni1==ni2).AND.(nj1==nj2) ) THEN
         ni = ni1 ; nj = nj1
         PRINT *, 'Grid is 2D: ni, nj =', ni, nj
         lregin = .FALSE.
      ELSE
         PRINT *, 'ERROR: problem with grid!' ; STOP
      END IF
   END IF

   ALLOCATE ( xlont(ni,nj), xlatt(ni,nj), xdum_r4(ni,nj) )
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



   ALLOCATE ( xvar(ni,nj), xvar1(ni,nj), xvar2(ni,nj), mask(ni,nj), vt_mod(Ntm) )

   ALLOCATE ( RES_2D_MOD(ni,nj), RES_2D_OBS(ni,nj) )
   RES_2D_MOD(:,:) = 0.
   RES_2D_OBS(:,:) = 0.


   IF ( lregin ) THEN
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

   IF ( l_orbit_file_is_nc ) THEN
      CALL GET_VAR_INFO(cf_obs, 'time', cunit, cdum)
      tut_obs  = GET_TIME_UNIT_T0(TRIM(cunit))
      PRINT *, ' *** Unit and reference time in track file:'
      PRINT *, tut_obs
   END IF
   PRINT *, ''





   !! Getting coordinates
   !! ~~~~~~~~~~~~~~~~~~~

   !!IF ( nk > 1 ) CALL GETVAR_1D(cf_mod, cv_z, vdepth(:,1))


   IF ( TRIM(cs_force_tv_m) /= '' ) THEN
      !! Building new time vector!
      idx = SCAN(TRIM(cs_force_tv_m),',')
      cdum = cs_force_tv_m(1:idx-1)
      READ(cdum,'(f)') rt0
      cdum = cs_force_tv_m(idx+1:)
      READ(cdum,'(f)') rdt
      PRINT *, ' *** MODEL: OVERIDING time vector with t0 and dt =', REAL(rt0,4), REAL(rdt,4)
      DO jt=1, Ntm
         vt_mod(jt) = rt0 + REAL(jt-1)*rdt
      END DO
   ELSE
      !! Reading it in input file:
      CALL GETVAR_1D(cf_mod, cv_t, vt_mod)
   END IF


   !! Getting longitude, latitude and mask in mesh_mask file:
   ! Longitude array:
   CALL GETVAR_2D   (i0, j0, cf_mm, cv_lon, 0, 0, 0, xdum_r4)
   xlont(:,:) = xdum_r4(:,:) ; i0=0 ; j0=0
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
   l_loc1 = (lon_min_1 <  0.).AND.(lon_min_1 > -170.).AND.(lon_max_1 >  0. ).AND.(lon_min_1 <  170.)
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
   CALL GETVAR_2D   (i0, j0, cf_mm, cv_lat, 0, 0, 0, xdum_r4)
   xlatt(:,:) = xdum_r4(:,:) ; i0=0 ; j0=0

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



   !! 3D LSM
   CALL GETMASK_2D(cf_mm, cv_mt, mask, jlev=1)











   !! Reading along-track from file:

   INQUIRE(FILE=TRIM(cf_obs), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: please provide the file containing definition of orbit track'; STOP
   END IF

   IF ( .NOT. l_orbit_file_is_nc ) THEN
      !! Getting number of lines:
      Nte = -1 ; io = 0
      OPEN (UNIT=13, FILE=TRIM(cf_obs))
      DO WHILE (io==0)
         READ(13,*,iostat=io)
         Nte = Nte + 1
      END DO
      PRINT*, Nte, ' points in '//TRIM(cf_obs)//'...'
      ALLOCATE ( Xgt(1,Nte), Ygt(1,Nte), vt_obs(Nte), Fgt(1,Nte) )
      !!
      REWIND(13)
      DO jte = 1, Nte
         READ(13,*) vt_obs(jte), Xgt(1,jte), Ygt(1,jte)
      END DO
      CLOSE(13)

   ELSE
      PRINT *, ''
      PRINT *, 'NetCDF orbit track!'
      CALL DIMS(cf_obs, 'time', Nte, nj1, nk, ni1)
      PRINT *, ' *** Nb. time records in NetCDF track file:', Nte
      ALLOCATE ( Xgt(1,Nte), Ygt(1,Nte), vt_obs(Nte), Fgt(1,Nte))
      CALL GETVAR_1D(cf_obs, 'time', vt_obs)
      CALL GETVAR_1D(cf_obs, 'longitude', Xgt(1,:))
      CALL GETVAR_1D(cf_obs, 'latitude',  Ygt(1,:))
      CALL GETVAR_1D(cf_obs, cv_obs,  Fgt(1,:))
      PRINT *, 'Done!'; PRINT *, ''
   END IF


   IF ( TRIM(cs_force_tv_e) /= '' ) THEN
      !! Building new time vector!
      idx = SCAN(TRIM(cs_force_tv_e),',')
      cdum = cs_force_tv_e(1:idx-1)
      READ(cdum,'(f)') rt0
      cdum = cs_force_tv_e(idx+1:)
      READ(cdum,'(f)') rdt
      PRINT *, ' *** TRACK: OVERIDING time vector with t0 and dt =', REAL(rt0,4), REAL(rdt,4)
      DO jt=1, Nte
         vt_obs(jt) = rt0 + REAL(jt-1)*rdt
         !PRINT *, ' vt_obs(jt)= ', vt_obs(jt)
      END DO
   END IF




   nib = ni ; njb = nj ; ji_min=1 ; ji_max=ni ; jj_min=1 ; jj_max=nj








   PRINT *, ''
   PRINT *, ' Time vector in track file:'
   CALL to_epoch_time_vect( tut_obs, vt_obs, l_dt_below_sec=.true. )
   PRINT *, '' ;  PRINT *, ''

   !! Defaults:
   Nten = Nte
   it1  = 1
   it2  = Nte

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
         !WHERE ( mask == 0 ) xmean = -9999.
         CALL PRTMASK(xmean, 'mean_'//TRIM(cv_mod)//'.nc', cv_mod, xlont, xlatt, 'nav_lon', 'nav_lat', rfill=-9999.)
         !STOP'lolo'
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


      !! Findin when we can start and stop when scanning the track file:
      !! it1, it2
      DO it1 = 1, Nte-1
         IF ( (vt_obs(it1) <= t_min).AND.(vt_obs(it1+1) > t_min) ) EXIT
      END DO
      DO it2 = it1, Nte-1
         IF ( (vt_obs(it2) <= t_max).AND.(vt_obs(it2+1) > t_max) ) EXIT
      END DO

      Nten = it2 - it1 + 1

      PRINT *, ' it1, it2 =',it1, it2
      PRINT *, Nten, '  out of ', Nte
      PRINT *, ' => ', vt_obs(it1), vt_obs(it2)
      PRINT *, ''
   END IF ! IF ( .NOT. l_debug_mapping )

   ALLOCATE ( IMETRICS(1,Nten,3), RAB(1,Nten,2), IPB(1,Nten), IGNORE(1,Nten), xlon_gt(1,Nten), xlat_gt(1,Nten) )

   IGNORE(:,:) = 1 !lolo

   !! Main time loop is on time vector in track file!


   cf_mapping = 'MAPPING__'//TRIM(cconf)//'.nc'


   !!
   xlon_gt(:,:) = Xgt(:,it1:it2)
   xlat_gt(:,:) = Ygt(:,it1:it2)

   DEALLOCATE ( Xgt, Ygt )

   IF ( .NOT. l_glob_lon_wize ) THEN
      ALLOCATE ( xdum_r8(1,Nten) )
      xdum_r8 = SIGN(1.,180.-xlon_gt)*MIN(xlon_gt,ABS(xlon_gt-360.)) ! like xlon_gt but between -180 and +180 !
      WHERE ( xdum_r8 < lon_min_1 ) IGNORE=0
      WHERE ( xdum_r8 > lon_max_1 ) IGNORE=0
      DEALLOCATE ( xdum_r8 )
   END IF

   IF ( .NOT. l_glob_lat_wize ) THEN
      WHERE ( xlat_gt < lat_min ) IGNORE=0
      WHERE ( xlat_gt > lat_max ) IGNORE=0
   END IF

   INQUIRE(FILE=trim(cf_mapping), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, ' *** Creating mapping file...'
      CALL MAPPING_BL(-1, xlont, xlatt, xlon_gt, xlat_gt, cf_mapping,  mask_domain_out=IGNORE)
      PRINT *, ' *** Done!'; PRINT *, ''
   ELSE
      PRINT *, ' *** File "',trim(cf_mapping),'" found in current directory, using it!'
      PRINT *, ''
   END IF

   CALL RD_MAPPING_AB(cf_mapping, IMETRICS, RAB, IPB)
   PRINT *, ''; PRINT *, ' *** Mapping and weights read into "',trim(cf_mapping),'"'; PRINT *, ''

   ALLOCATE (JIidx(1,Nten) , JJidx(1,Nten) )
   JIidx(1,:) = IMETRICS(1,:,1)
   JJidx(1,:) = IMETRICS(1,:,2)



   !! Showing iy in file mask_+_nearest_points.nc:
   IF ( l_debug ) THEN
      !! Finding and storing the nearest points of NEMO grid to track points:
      !CALL FIND_NEAREST_POINT(Xgt, Ygt, xlont, xlatt,  JIidx, JJidx)
      ALLOCATE ( show_obs(nib,njb) )
      show_obs(:,:) = 0.
      DO jte = 1, Nten
         IF ( (JIidx(1,jte)>0).AND.(JJidx(1,jte)>0) )  show_obs(JIidx(1,jte), JJidx(1,jte)) = REAL(jte,4)
      END DO
      WHERE (mask == 0) show_obs = -9999.
      CALL PRTMASK(REAL(show_obs(:,:),4), 'mask_+_nearest_points__'//TRIM(cconf)//'.nc', 'mask', xlont, xlatt, 'lon0', 'lat0', rfill=-9999.)
      !lolo:
      !CALL PRTMASK(REAL(xlont(:,:),4), 'lon_360.nc', 'lon')
      !show_obs = SIGN(1.,180.-xlont)*MIN(xlont,ABS(xlont-360.))
      !CALL PRTMASK(REAL(show_obs(:,:),4), 'lon_-180-180.nc', 'lon')
      !WHERE ( (show_obs > 10.).OR.(show_obs < -90.) ) show_obs = -800.
      !CALL PRTMASK(REAL(show_obs(:,:),4), 'lon_masked.nc', 'lon')
      !STOP 'interp_to_ground_obs.f90'
      !lolo.

      DEALLOCATE ( show_obs )
   END IF

   IF ( l_debug_mapping ) STOP'l_debug_mapping'



   !STOP 'mapping done!'

   ALLOCATE ( vte(Nten), Ftrack_mod(Nten), Ftrack_obs(Nten), Fmask(Nten), Ftrack_mod_np(Nten) )


   vte(:) = vt_obs(it1:it2)

   Ftrack_mod_np(:) = -9999.
   Ftrack_obs(:)    = -9999.
   Ftrack_mod(:)    = -9999.
   Fmask(:)         = 0

   jt_s = 1 ; ! time step model!

   jtm_1_o = -100
   jtm_2_o = -100

   DO jte = 1, Nten

      rt = vte(jte)

      IF ( (rt >= t_min_m).AND.(rt < t_max_m) ) THEN

         !! Two surrounding time records in model file => jtm_1 & jtm_2
         DO jt=jt_s, Ntm-1
            IF ( (rt >= vt_mod(jt)).AND.(rt < vt_mod(jt+1)) ) EXIT
         END DO
         !!
         jtm_1 = jt
         jtm_2 = jt+1
         IF (jte==1) jt_s = jtm_1 ! Saving the actual first useful time step of the model!

         PRINT *, 'Treating track time =>', rt, '     model jtm_1 =', jtm_1

         !! If first time we have these jtm_1 & jtm_2, getting the two surrounding fields:
         IF ( (jtm_1>jtm_1_o).AND.(jtm_2>jtm_2_o) ) THEN
            IF ( jtm_1_o == -100 ) THEN
               PRINT *, ' *** Reading field '//TRIM(cv_mod)//' in '//TRIM(cf_mod)
               PRINT *, '    => at jtm_1=', jtm_1, '  (starting from jt1=',jt_s,')'
               CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jtm_1, xvar1, jt1=jt_s)
               IF ( l_use_anomaly ) xvar1 = xvar1 - xmean
               IF ( l_drown_in ) CALL DROWN(-1, xvar1, mask, nb_inc=5, nb_smooth=2)
            ELSE
               xvar1(:,:) = xvar2(:,:)
            END IF
            PRINT *, ' *** Reading field '//TRIM(cv_mod)//' in '//TRIM(cf_mod)
            PRINT *, '    => at jtm_2=', jtm_2, '  (starting from jt1=',jt_s,')'
            CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jtm_2, xvar2, jt1=jt_s)
            IF ( l_use_anomaly ) xvar2 = xvar2 - xmean
            IF ( l_drown_in ) CALL DROWN(-1, xvar2, mask, nb_inc=5, nb_smooth=2)

            xdum_r4 = (xvar2 - xvar1) / (vt_mod(jtm_2) - vt_mod(jtm_1)) ! xdum_r4 is the slope here !!!

         END IF

         !! Linear interpolation of field at time rt:
         xvar(:,:) = xvar1(:,:) + xdum_r4(:,:)*(rt - vt_mod(jtm_1))

         !! Performing bilinear interpolation:
         iP       = IMETRICS(1,jte,1)
         jP       = IMETRICS(1,jte,2)
         iquadran = IMETRICS(1,jte,3)
         alpha    = RAB(1,jte,1)
         beta     = RAB(1,jte,2)





         IF ( mask(iP,jP)==1 ) THEN

            r_obs    = Fgt(1,it1+jte-1)
            l_obs_ok = ( r_obs > -20.).AND.( r_obs < 20.)
            
            IF ( l_obs_ok .AND. (mask(iP,jP)==1).AND.(iP/=INT(rflg)).AND.(jP/=INT(rflg)) ) THEN
               
               !! Ignore points that are just 1 point away from land:
               ip1 = MIN(iP+1,ni) ; jp1 = MIN(jP+1,nj)
               im1 = MAX(iP-1,1)  ; jm1 = MAX(jP-1,1)
               idot = mask(ip1,jP) + mask(ip1,jp1) + mask(iP,jp1) + mask(im1,jp1) &
                  & + mask(im1,jP) + mask(im1,jm1) + mask(iP,jm1) + mask(ip1,jm1)
               !! => idot == 8 if in that case...
               IF (idot==8) THEN
                  !! Model, nearest point:
                  Ftrack_mod_np(jte) =  xvar(JIidx(1,jte),JJidx(1,jte)) ! NEAREST POINT interpolation
                  !! Model, 2D bilinear interpolation:
                  Ftrack_mod(jte) = INTERP_BL(-1, iP, jP, iquadran, alpha, beta, REAL(xvar,8))
                  !! Observations as on their original point:
                  Ftrack_obs(jte) = r_obs
                  !! On the model grid for info:
                  RES_2D_MOD(iP,jP) = Ftrack_mod(jte)
                  RES_2D_OBS(iP,jP) = Ftrack_obs(jte)
                  !!
                  Fmask(jte) = 1 ! That was a valid point!
                  !!
               END IF
            END IF
         END IF
         
         !IF ( (JIidx(1,jte)>0).AND.(JJidx(1,jte)>0) ) &
         !   &

         jtm_1_o = jtm_1
         jtm_2_o = jtm_2
         jt_s    = jtm_1 ! so we do not rescan from begining...

      END IF

   END DO

   WHERE ( Fmask == 0 )
      Ftrack_mod    = -9999.
      Ftrack_mod_np = -9999.
      Ftrack_obs    = -9999.
   END WHERE


   PRINT *, ''
   !WRITE(cf_out, '("track_",a,"_",a,".nc")') TRIM(cv_mod), TRIM(cf_obs)
   cf_out = 'result__'//TRIM(cconf)//'.nc4'
   PRINT *, ' * Output file = ', trim(cf_out)
   PRINT *, ''

   CALL PT_SERIES(vte(:), REAL(Ftrack_mod,4), cf_out, 'time', cv_mod, 'm', 'Model data, bi-linear interpolation', -9999., &
      &           ct_unit=TRIM(cunit_time_out), &
      &           vdt2=REAL(Ftrack_mod_np,4),cv_dt2=TRIM(cv_mod)//'_np',cln2='Model data, nearest-point interpolation', &
      &           vdt3=REAL(Ftrack_obs,4),   cv_dt3=cv_obs,             cln3='Original data as in track file...',   &
      &           vdt4=REAL(xlon_gt(1,:),4), cv_dt4='longitude',        cln4='Longitude (as in track file)',  &
      &           vdt5=REAL(xlat_gt(1,:),4), cv_dt5='latitude',         cln5='Latitude (as in track file)' ,  &
      &           vdt6=REAL(Fmask,4),        cv_dt6='mask',             cln6='Mask', &
      &           vdt7=REAL(IGNORE(1,:),4),  cv_dt7='ignore_out',       cln7='Ignore mask on target track (ignored where ignore_out==0)')

   WHERE ( mask == 0 )
      RES_2D_MOD = -9999.
      RES_2D_OBS = -9999.
   END WHERE

   CALL PRTMASK(RES_2D_MOD, 'RES_2D_MOD__'//TRIM(cconf)//'.nc', cv_mod, xlont, xlatt, 'nav_lon', 'nav_lat', rfill=-9999.)
   CALL PRTMASK(RES_2D_OBS, 'RES_2D_OBS__'//TRIM(cconf)//'.nc', cv_obs, xlont, xlatt, 'nav_lon', 'nav_lat', rfill=-9999.)



   !IF ( l_debug ) DEALLOCATE ( JIidx, JJidx )
   !DEALLOCATE ( Fgt )
   !DEALLOCATE ( Ftrack_mod, Ftrack_mod_np, Ftrack_obs )
   !DEALLOCATE ( xlont, xlatt, xvar, xvar1, xvar2, xdum_r4, mask, xdum_r8 )


   PRINT *, ''
   PRINT *, 'Written!'
   PRINT *, ' => check:'
   PRINT *, TRIM(cf_out)
   PRINT *, ''


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


END PROGRAM INTERP_TO_GROUND_TRACK




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
   WRITE(6,*) ' -p  <track_file>     => Specify name of file containing orbit tack (columns: time, lon, lat)'
   WRITE(6,*) ''
   WRITE(6,*) ' -n  <name>           => file containing orbit track is in NetCDF, and this is the name of var'
   WRITE(6,*) '                         (default is columns in ASCII file <time> <lon> <lat>'
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
   WRITE(6,*) ' -f  <t0,dt>          => overide time vector in input NEMO file with one of same length'
   WRITE(6,*) '                         based on t0 and dt (in seconds!) (ex: " ... -f 0.,3600.")'
   WRITE(6,*) ''
   WRITE(6,*) ' -g  <t0,dt>          => overide time vector in track file with one of same length'
   WRITE(6,*) '                         based on t0 and dt (in seconds!) (ex: " ... -f 0.,3600.")'
   WRITE(6,*) ''
   WRITE(6,*) ' -h                   => Show this message'
   WRITE(6,*) ''
   !!
   !CLOSE(6)
   STOP
   !!
END SUBROUTINE usage
!!
