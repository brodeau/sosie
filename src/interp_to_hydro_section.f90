PROGRAM INTERP_TO_HYDRO_SECTION

   USE io_ezcdf
   USE mod_conf
   USE mod_bilin_2d
   USE mod_manip
   USE mod_drown

   !!========================================================================
   !! Purpose :
   !! ---------
   !!
   !! Author :   Laurent Brodeau, 2020
   !! --------
   !!
   !!
   !!  Requirements for the hydrographic Section netCDF file to be fed in:
   !!
   !!  - the record is named "record" and is "UNLIMITED" (1 record == 1 vertical profile)
   !!
   !!  - coordinates: "longitude" and "latitude" variables at each record => "latitude(record)", "longitude(record)"
   !!
   !!  - since it's a section and not a "track" a depth dimension is needed
   !!    => depth dimension is "z"
   !!    => depth variable is a 2D array called "depth" => "depth(record,z)"
   !!
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
   REAL(8), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER :: Nrec, jr, jt_m, nz, iP, jP, iquadran, jk
   !!
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: xlon_o, xlat_o, xdum_r8
   !!


   !! Coupe stuff:
   REAL(8), DIMENSION(:),   ALLOCATABLE :: Fhs_m, Fhs_m_np, Fhs_o, vdist
   REAL(8), DIMENSION(:),   ALLOCATABLE :: vzo, vfo, vrec, vt_mod
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: z_o, F_o
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
      &    cv_z   = 'depth',         &
      &    cv_lon = 'nav_lon',         & ! input grid longitude name, T-points
      &    cv_lat = 'nav_lat'            ! input grid latitude name,  T-points

   CHARACTER(len=256)  :: cr, cnm_fill
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
      &    l_write_nc_show_route = .FALSE., &
      &     l_exist   = .FALSE., &
      &     l_use_anomaly = .FALSE., &  ! => will transform a SSH into a SLA (SSH - MEAN(SSH))
      &     l_loc1, l_loc2, &
      &     l_z_o_upd, &
      &     l_obs_ok, lfillval
   !!
   !!
   CHARACTER(len=400)  :: &
      &    cf_obs = 'hydro_section.nc', &
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

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xdum_r4, show_obs, xf_m, xmean
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: xlon_m, xlat_m
   REAL(8), DIMENSION(:),   ALLOCATABLE :: vdepth_mod
   !!
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJidx, JIidx    ! debug
   !!
   INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: imask, imask_ignr
   INTEGER(1), DIMENSION(:),   ALLOCATABLE :: Fmask
   !!
   INTEGER :: jt, jt0, jt_s
   !!
   REAL(8) :: alpha, beta
   !!
   CHARACTER(LEN=2), DIMENSION(12), PARAMETER :: &
      &            clist_opt = (/ '-h','-v','-x','-y','-z','-t','-i','-p','-n','-m','-S','-M' /)

   REAL(8) :: lon_min_2, lon_max_2, lat_min, lat_max, r_obs
   REAL(4) :: rfillval_mod

   CHARACTER(80), PARAMETER :: cunit_time_trg = 'seconds since 1970-01-01 00:00:00'

   OPEN(UNIT=6, FORM='FORMATTED', RECL=512)



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
         CALL GET_MY_ARG('input depth', cv_z)

      CASE('-t')
         CALL GET_MY_ARG('time', cv_t)

      CASE('-p')
         CALL GET_MY_ARG('orbit track file', cf_obs)

      CASE('-m')
         l_get_mask_metrics_from_meshmask = .TRUE.
         CALL GET_MY_ARG('mesh_mask file', cf_mm)

      CASE('-S')
         l_write_nc_show_route = .TRUE.

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

   !! Name of config:
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
   IF ( nk < 2 )  THEN
      PRINT *, 'ERROR: ',trim(cv_mod),' must have a depth dimension of at least 2!' ; STOP
   END IF
   IF ( Ntm < 1 ) THEN
      PRINT *, 'ERROR: ',trim(cv_mod),' must have at least a time record!' ; STOP
   END IF
   PRINT *, 'Dimension for "'//TRIM(cv_mod)//'" into "'//TRIM(cf_mod)//'":'
   PRINT *, '   => ni =', ni ;   PRINT *, '   => nj =', nj
   PRINT *, '   => nk =', nk ;   PRINT *, '   => Ntm =', Ntm
   PRINT *, ''


   PRINT *, ''
   PRINT *, ' *** Allocating ni x nj arrays...'
   ALLOCATE ( xlon_m(ni,nj), xlat_m(ni,nj), xdum_r4(ni,nj), &
      &       xf_m(ni,nj),  &
      &       vdepth_mod(nk), imask(ni,nj), vt_mod(Ntm) )
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






   !! Getting coordinates and depth vector of model data
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   ! Longitude array:
   CALL GETVAR_2D (i0, j0, cf_mm,  cv_lon, 0, 0, 0, xlon_m(:,:))  ; i0=0 ; j0=0

   !! Min an max lon:
   lon_min_1 = MINVAL(xlon_m)
   lon_max_1 = MAXVAL(xlon_m)
   PRINT *, ' *** Minimum longitude on source domain before : ', lon_min_1
   PRINT *, ' *** Maximum longitude on source domain before : ', lon_max_1
   !!
   WHERE ( xlon_m < 0. ) xlon_m = xlon_m + 360.0_8
   !!
   lon_min_2 = MINVAL(xlon_m)
   lon_max_2 = MAXVAL(xlon_m)
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
   CALL GETVAR_2D   (i0, j0, cf_mm,  cv_lat, 0, 0, 0, xlat_m(:,:)) ; i0=0 ; j0=0

   !! Min an max lat:
   lat_min = MINVAL(xlat_m)
   lat_max = MAXVAL(xlat_m)
   PRINT *, ' *** Minimum latitude on source domain : ', lat_min
   PRINT *, ' *** Maximum latitude on source domain : ', lat_max

   l_glob_lat_wize = .TRUE.
   IF ( lat_max < 88. ) THEN
      l_glob_lat_wize =.FALSE.
      WRITE(*,'("  => going to disregard points of target domain with lat < ",f7.2," and lat > ",f7.2)'), lat_min,lat_max
   END IF
   PRINT *, ''


   !! Depth array
   CALL GETVAR_1D(cf_mm, cv_z, vdepth_mod(:))

   PRINT *, ''
   PRINT *, ' *** Model depths:'
   PRINT *, vdepth_mod
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

   IF(l_debug) CALL DUMP_FIELD(REAL(imask), 'mask_surface_in.nc', 'lsm') !, xlon_m, xlat_m, 'nav_lon', 'nav_lat', rfill=-9999.)









   !! Time to read the hydrographic section file:

   INQUIRE(FILE=TRIM(cf_obs), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: please provide the file containing hydrographic section'; STOP
   END IF

   !! Geting dimmensions:
   CALL DIMS(cf_obs, cv_obs, nz, ni1, nj1, Nrec)
   IF( (ni1 /= -1).OR.(nj1 /= -1) ) THEN
      PRINT *, 'ERROR: variable ',TRIM(cv_obs), 'must have only 2 dimmensions [record,z] !' ; STOP
   END IF
   PRINT *, ' *** Nb. of records in hydrographic section:', Nrec
   PRINT *, ' *** Nb. of depths  in hydrographic section:', nz
   PRINT *, ''


   ALLOCATE ( xlon_o(1,Nrec), xlat_o(1,Nrec), vzo(nz), vfo(nz), z_o(nz,Nrec), F_o(nz,Nrec) )

   CALL GETVAR_1D(cf_obs, 'longitude', xlon_o(1,:))
   CALL GETVAR_1D(cf_obs, 'latitude',  xlat_o(1,:))

   
   !DO jr=1, Nrec
   !   
   !END DO

   PRINT *, 'Done!'; PRINT *, ''



   nib = ni ; njb = nj ; ji_min=1 ; ji_max=ni ; jj_min=1 ; jj_max=nj


   !IF ( .NOT. l_debug_mapping ) THEN
   !IF ( l_use_anomaly ) THEN
   !   PRINT *, ''
   !   PRINT *, ' *** Computing mean of field '//TRIM(cv_mod)//' for considering anomaly later...'
   !   ALLOCATE ( xmean(ni,nj) )
   !   xmean = 0.
   !   DO jt = 1, Ntm
   !      CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jt, xdum_r4)
   !      xmean = xmean + xdum_r4/REAL(Ntm,4)
   !   END DO
   !   id_f1=0 ;  id_v1=0
   !   !WHERE ( imask == 0 ) xmean = -9999.
   !   CALL DUMP_FIELD(xmean, 'mean_'//TRIM(cv_mod)//'.nc', cv_mod, xlon_m, xlat_m, 'nav_lon', 'nav_lat', rfill=-9999.)
   !   !STOP 'lolo'
   !END IF
   !END IF ! IF ( .NOT. l_debug_mapping )


   !! Main time loop is on time vector in track file!


   cf_mpg = 'MAPPING__'//TRIM(cconf)//'.nc'


   !PRINT *, '  lon =  20. =>', to_degE(20._8)
   !PRINT *, '  lon = -20. =>', to_degE(-20._8)
   !STOP   

   PRINT *, 'LOLO: xlon_o =', xlon_o
   PRINT *, ''
   xlon_o = to_degE(xlon_o) ! to degrees East
   PRINT *, 'LOLO: xlon_o =', xlon_o

   
   !! Allocate arrays on the final retained size
   ALLOCATE ( IMETRICS(1,Nrec,3), RAB(1,Nrec,2), IPB(1,Nrec),  vrec(Nrec) )
   ALLOCATE ( Fhs_m(Nrec), Fhs_o(Nrec), Fmask(Nrec), Fhs_m_np(Nrec), vdist(Nrec) )


   INQUIRE(FILE=trim(cf_mpg), EXIST=l_exist ) !
   IF ( .NOT. l_exist ) THEN
      PRINT *, ' *** Creating mapping file...' !
      CALL MAPPING_BL(-1, xlon_m, xlat_m, xlon_o, xlat_o, cf_mpg ) !,  mask_domain_trg=IGNORE) don't need ignore, points have been removed!
      PRINT *, ' *** Done!'; PRINT *, ''
   ELSE
      PRINT *, ' *** File "',trim(cf_mpg),'" found in current directory, using it!'
      PRINT *, ''
   END IF

   CALL RD_MAPPING_AB(cf_mpg, IMETRICS, RAB, IPB)
   PRINT *, ''; PRINT *, ' *** Mapping and weights read into "',trim(cf_mpg),'"'; PRINT *, ''

   ALLOCATE (JIidx(1,Nrec) , JJidx(1,Nrec) )
   JIidx(1,:) = IMETRICS(1,:,1)
   JJidx(1,:) = IMETRICS(1,:,2)


   !! Showing route in file mask_+_nearest_points.nc:
   IF ( l_write_nc_show_route ) THEN
      !! Finding and storing the nearest points of NEMO grid to track points:
      !CALL FIND_NEAREST_POINT(xlon_o, xlat_o, xlon_m, xlat_m,  JIidx, JJidx)
      ALLOCATE ( show_obs(nib,njb) )
      show_obs(:,:) = -9999.
      DO jr = 1, Nrec
         IF ( (JIidx(1,jr)>0).AND.(JJidx(1,jr)>0) )  show_obs(JIidx(1,jr), JJidx(1,jr)) = REAL(jr,4)
      END DO
      WHERE (imask == 0) show_obs = -100.
      CALL DUMP_FIELD(REAL(show_obs(:,:),4), 'mask_+_nearest_points__'//TRIM(cconf)//'.nc', 'mask', xlon_m, xlat_m, cv_lon, cv_lat, rfill=-9999.)
      !lolo:
      !CALL DUMP_FIELD(REAL(xlon_m(:,:),4), 'lon_360.nc', 'lon')
      !show_obs = SIGN(1.,180.-xlon_m)*MIN(xlon_m,ABS(xlon_m-360.))
      !CALL DUMP_FIELD(REAL(show_obs(:,:),4), 'lon_-180-180.nc', 'lon')
      !WHERE ( (show_obs > 10.).OR.(show_obs < -90.) ) show_obs = -800.
      !CALL DUMP_FIELD(REAL(show_obs(:,:),4), 'lon_masked.nc', 'lon')
      !STOP 'interp_to_ground_obs.f90'
      !lolo.

      DEALLOCATE ( show_obs )
   END IF

   IF ( l_debug_mapping ) STOP 'l_debug_mapping'




   !STOP 'mapping done!'



   Fhs_m_np(:) = -9999.
   Fhs_o(:)    = -9999.
   Fhs_m(:)    = -9999.
   Fmask(:)         = 0

   jt_s    = 1



   ! Reading model field:
   jt_m = 1   
   CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jt_m, xf_m) !#lolo: whole domain! => improve

   
   DO jr = 1, Nrec

      !! Reading the input vertical grid at record jr into F_o:
      CALL GETVAR_1D(cf_obs, 'depth', vzo, jrec=jr)
      vzo = ABS(vzo)      ! we want depth to be positive:
                  
      !! Reading the input vertical profiles at record jr into F_o:
      CALL GETVAR_1D(cf_obs, cv_obs, vfo, jrec=jr)
      !! => technically useless, only for extra info in result file and to test if was a valid profile...

      PRINT *, vfo
      l_z_o_upd = ( vzo(1) > vzo(nz) )
      
      IF (l_z_o_upd) THEN
         PRINT *, 'Depth and data are upside-down into input file...'
         z_o(:,jr) = vzo(nz:1:-1)
         F_o(:,jr) = vfo(nz:1:-1)
         PRINT *, ''
      END IF
      
      !PRINT *, 'DEPTHS at record jr =',jr
      !PRINT *, z_o(:,jr)      
      !PRINT *, 'PROFILE at record jr =',jr
      !PRINT *, F_o(:,jr)
      !STOP'LALA'
      
      !! If first time we have these jtm_1 & jtm_2, getting the two surrounding fields:
      !IF ( (jtm_1>jtm_1_o).AND.(jtm_2>jtm_2_o) ) THEN
      !   IF ( (jtm_1_o == -100).OR.(jtm_1 > jtm_2_o) ) THEN
      !      PRINT *, ' *** Reading field '//TRIM(cv_mod)//' in '//TRIM(cf_mod)
      !      PRINT *, '    => at jtm_1=', jtm_1, '  (starting from jt1=',jt0,')'!
      !
      !      IF ( l_use_anomaly ) xf_m1 = xf_m1 - xmean
      !      IF ( l_drown_in ) CALL DROWN(-1, xf_m1, imask, nb_inc=5)
      !   ELSE
      !      xf_m1(:,:) = xf_m2(:,:)
      !   END IF
      !   PRINT *, ' *** Reading field '//TRIM(cv_mod)//' in '//TRIM(cf_mod)
      !   PRINT *, '    => at jtm_2=', jtm_2, '  (starting from jt1=',jt0,')'
      !   CALL GETVAR_2D(id_f1, id_v1, cf_mod, cv_mod, Ntm, 0, jtm_2, xf_m2, jt1=jt0)
      !   IF ( l_use_anomaly ) xf_m2 = xf_m2 - xmean
      !   IF ( l_drown_in ) CALL DROWN(-1, xf_m2, imask, nb_inc=5)
      !
      !   xdum_r4 = (xf_m2 - xf_m1) / (vt_mod(jtm_2) - vt_mod(jtm_1)) ! xdum_r4 is the slope here !!!
      !
      !   PRINT *, ''
      !END IF

      !! Linear interpolation of field at time rt:
      !xf_m(:,:) = xf_m1(:,:) + xdum_r4(:,:)*(rt - vt_mod(jtm_1))

      !! Performing bilinear interpolation:
      iP       = IMETRICS(1,jr,1)
      jP       = IMETRICS(1,jr,2)
      iquadran = IMETRICS(1,jr,3)
      alpha    = RAB(1,jr,1)
      beta     = RAB(1,jr,2)



      jk = 1   ! #lolo only first level for now
      
      !LOLO: IF ( (iP/=INT(rflg)).AND.(jP/=INT(rflg)) ) THEN
      IF ( (iP>0).AND.(jP>0) ) THEN
         IF ( imask(iP,jP)==1 ) THEN
            r_obs    = F_o(jk,jr)
            l_obs_ok = ( r_obs > -20.).AND.( r_obs < 50.)
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
                  Fhs_m_np(jr) =  xf_m(JIidx(1,jr),JJidx(1,jr)) ! NEAREST POINT interpolation
                  !! Model, 2D bilinear interpolation:
                  Fhs_m(jr) = REAL( INTERP_BL(-1, iP, jP, iquadran, alpha, beta, xf_m) , 8)
                  !! Observations as on their original point:
                  Fhs_o(jr) = r_obs
                  IF ( l_debug ) THEN
                     !! On the model grid for info:
                     RES_2D_MOD(iP,jP) = Fhs_m(jr)
                     RES_2D_OBS(iP,jP) = Fhs_o(jr)
                  END IF
                  !!
                  Fmask(jr) = 1 ! That was a valid point!
                  !!
                  !rcycle_obs(jr) = REAL( icycle(jr), 8 )
                  !!
               END IF
            END IF
         END IF
      END IF


   END DO



   
   !! Vector distance (in km)
   vdist(:) = 0.
   DO jt = 2, Nrec
      IF ( (Fmask(jt)==1).AND.(Fmask(jt-1)==1) ) vdist(jt) = vdist(jt-1) + DISTANCE( xlon_o(1,jt), xlon_o(1,jt-1), xlat_o(1,jt), xlat_o(1,jt-1) )
   END DO

   WHERE ( Fmask == 0 )
      Fhs_m    = -9999.
      Fhs_m_np = -9999.
      Fhs_o    = -9999.
      vdist     = -9999.
   END WHERE

   WHERE ( Fhs_m < -9990. ) Fhs_m = -9999.


   PRINT *, 'LOLO Fhs_m = ', Fhs_m
   

   PRINT *, ''
   !WRITE(cf_out, '("track_",a,"_",a,".nc")') TRIM(cv_mod), TRIM(cf_obs)
   cf_out = 'result__'//TRIM(cconf)//'.nc'
   PRINT *, ' * Output file = ', trim(cf_out)
   PRINT *, ''


   vrec(:) = 0. ! => forces no variable record into PT_SERIES
   
   CALL PT_SERIES(vrec(:), REAL(Fhs_m,4), cf_out, 'record', cv_mod, 'm', 'Model data, bi-linear interpolation', -9999.,      &
      &           ct_unit=TRIM(cunit_time_trg),                                                                              &
      &           vdt02=REAL(Fhs_m_np,4),    cv_dt02=TRIM(cv_mod)//'_np',cln02='Model data, nearest-point interpolation',    &
      &           vdt03=REAL(Fhs_o,4),       cv_dt03=cv_obs,             cln03='Original data as in track file...',          &
      &           vdt04=REAL(xlon_o(1,:),4), cv_dt04='longitude',        cln04='Longitude (as in track file)',               &
      &           vdt05=REAL(xlat_o(1,:),4), cv_dt05='latitude',         cln05='Latitude (as in track file)' ,               &
      &           vdt06=REAL(Fmask,4),       cv_dt06='mask',             cln06='Mask' ,                                      &
      &           vdt07=REAL(vdist,4),       cv_dt07='distance',         cln07='Distance (in km) from first point of segment' )
   

   !#lulu
   STOP'LOLO0'



   
   IF ( l_debug ) THEN
      WHERE ( imask == 0 )
         RES_2D_MOD = -9999.
         RES_2D_OBS = -9999.
      END WHERE
      CALL DUMP_FIELD(RES_2D_MOD, 'RES_2D_M__'//TRIM(cconf)//'.nc', cv_mod, xlon_m, xlat_m, 'nav_lon', 'nav_lat', rfill=-9999.)
      CALL DUMP_FIELD(RES_2D_OBS, 'RES_2D_OBS__'//TRIM(cconf)//'.nc', cv_obs, xlon_m, xlat_m, 'nav_lon', 'nav_lat', rfill=-9999.)
   END IF

   !IF ( l_debug ) DEALLOCATE ( JIidx, JJidx )
   DEALLOCATE ( F_o )
   !DEALLOCATE ( Fhs_m, Fhs_m_np, Fhs_o )
   !DEALLOCATE ( xlon_m, xlat_m, xf_m, xf_m1, xf_m2, xdum_r4, imask, xdum_r8 )


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
      WRITE(6,*) ' -v  <name>           => Specify variable name of interest in input file'
      WRITE(6,*) ''
      WRITE(6,*) ' -p  <track_file>     => Specify name of NetCDF file containing hydrographic section'
      WRITE(6,*) ''
      WRITE(6,*) ' -n  <name>           => name of variable of interest in hydrographic section'
      WRITE(6,*) ''
      !!
      WRITE(6,*) ''
      WRITE(6,*) '    Optional:'
      WRITE(6,*)  ''
      WRITE(6,*) ' -x  <name>           => Specify longitude name in input file (default: '//TRIM(cv_lon)//')'
      WRITE(6,*) ''
      WRITE(6,*) ' -y  <name>           => Specify latitude  name in input file  (default: '//TRIM(cv_lat)//')'
      WRITE(6,*) ''
      WRITE(6,*) ' -z  <name>           => Specify depth name in input file (default: '//TRIM(cv_z)//')'
      WRITE(6,*) ''
      WRITE(6,*) ' -t  <name>           => Specify time name in input file (default: '//TRIM(cv_t)//')'
      WRITE(6,*) ''
      WRITE(6,*) ' -m  <mesh_mask_file> => Specify mesh_mask file to be used (default: '//TRIM(cf_mm)//')'
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


END PROGRAM INTERP_TO_HYDRO_SECTION
