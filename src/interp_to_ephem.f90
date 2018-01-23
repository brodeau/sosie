PROGRAM INTERP_TO_EPHEM

   !USE datetime_module, ONLY:datetime ! => https://github.com/wavebitscientific/datetime-fortran/releases

   USE io_ezcdf
   USE mod_conf
   USE mod_drown
   USE mod_akima_2d
   USE mod_bilin_2d
   USE mod_manip ! debug

   !!========================================================================
   !! Purpose :
   !!
   !! ---------
   !!
   !! Author :   Laurent Brodeau
   !! --------
   !!
   !!========================================================================

   IMPLICIT NONE

   !! ************************ Configurable part ****************************
   !!
   LOGICAL, PARAMETER :: &
      &   l_debug = .TRUE., &
      &   l_debug_mapping = .FALSE., &
      &   l_akima = .true., &
      &   l_bilin = .false.
   !!
   LOGICAL :: &
      &      l_orbit_file_is_nc    = .FALSE.
   !!
   REAL(8), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER :: Nte, Nten, io, idx, iP, jP, iquadran
   !!
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xtar, Ytar, Ztar
   REAL(4), DIMENSION(:,:), ALLOCATABLE :: Zpt4
   !!
   !! Coupe stuff:
   REAL(8), DIMENSION(:), ALLOCATABLE :: Ftrack, Fmask, Ftrack_np, Ftrack_ephem
   REAL(4), DIMENSION(:),   ALLOCATABLE :: xcmask
   REAL(8), DIMENSION(:,:),   ALLOCATABLE :: vposition

   REAL(8), DIMENSION(:,:),   ALLOCATABLE :: vdepth
   REAL(8), DIMENSION(:),     ALLOCATABLE :: vt_model, vt_ephem   ! in seconds

   REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: RAB       !: alpha, beta
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IMETRICS  !: iP, jP, iquadran at each point
   INTEGER, DIMENSION(:,:),   ALLOCATABLE :: IPB       !: ID of problem


   !! Grid, default name :
   CHARACTER(len=80) :: &
      &    cv_t   = 'time_counter',  &
      &    cv_mt  = 'tmask',         &
      &    cv_z   = 'deptht',        &
      &    cv_lon = 'glamt',       & ! input grid longitude name, T-points
      &    cv_lat = 'gphit'          ! input grid latitude name,  T-points

   CHARACTER(len=256)  :: cr, cunit, ctrack, cdum
   !!
   !!
   !!******************** End of conf for user ********************************
   !!
   !!               ** don't change anything below **
   !!
   LOGICAL ::  &
      &     l_exist   = .FALSE.
   !!
   !!
   CHARACTER(len=400)  :: &
      &    cf_track   = 'track.dat', &
      &    cf_mm='mesh_mask.nc', &
      &    cs_force_tv_m='', &
      &    cs_force_tv_e=''
   !!
   INTEGER      :: &
      &    jarg,   &
      &    i0, j0, ifo, ivo,   &
      &    ni, nj, nk=0, Ntm=0, &
      &    ni1, nj1, ni2, nj2, &
      &    iargc, id_f1, id_v1
   !!
   !!
   INTEGER :: imin, imax, jmin, jmax, isav, jsav, ji_min, ji_max, jj_min, jj_max, nib, njb

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xvar, xvar1, xvar2, xslp

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xdum2d, show_track
   REAL(8), DIMENSION(:,:), ALLOCATABLE ::    &
      &    xlont, xlatt
   !!
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJidx, JIidx    ! debug
   !!
   INTEGER(2), DIMENSION(:,:), ALLOCATABLE :: mask
   !!
   INTEGER :: jt, jte, jt_s, jtm_1, jtm_2, jtm_1_o, jtm_2_o
   !!
   REAL(8) :: rA, rB, dlon, dlat, dang, lon_min, lon_max, lat_min, lat_max, rt, rt0, rdt, &
      &       t_min_e, t_max_e, t_min_m, t_max_m, &
      &       alpha, beta, t_min, t_max
   !!
   CHARACTER(LEN=2), DIMENSION(12), PARAMETER :: &
      &            clist_opt = (/ '-h','-v','-x','-y','-z','-t','-i','-p','-n','-m','-f','-g' /)


   INTEGER(4) :: itime
   TYPE(t_unit_t0) :: tut_epoch, tut_ephem, tut_model

   INTEGER :: it1, it2

   !TYPE(datetime) :: dt_r0

   !dt_r0 = datetime(1977,4,19,0,0,0)

   !PRINT *, ' T epoch for 1977-04-19 =>', dt_r0%secondsSinceEpoch()




   !! Epoch is our reference time unit, it is "seconds since 1970-01-01 00:00:00" which translates into:
   tut_epoch%unit   = 's'
   tut_epoch%year   = 1970
   tut_epoch%month  = 1
   tut_epoch%day    = 1
   tut_epoch%hour   = 0
   tut_epoch%minute = 0
   tut_epoch%second = 0


   PRINT *, ''


   !tut_ephem  = GET_TIME_UNIT_T0('days since 1950-01-01 00:00:00')
   !PRINT *, ' tut_ephem =', tut_ephem
   !tut_model  = GET_TIME_UNIT_T0('seconds since 1950-01-01 00:00:00')
   !PRINT *, ' tut_model =', tut_model
   !STOP


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
         CALL GET_MY_ARG('input file', cf_in)

      CASE('-v')
         CALL GET_MY_ARG('input variable', cv_in)

      CASE('-x')
         CALL GET_MY_ARG('longitude', cv_lon)

      CASE('-y')
         CALL GET_MY_ARG('latitude', cv_lat)

      CASE('-z')
         CALL GET_MY_ARG('depth', cv_z)

      CASE('-t')
         CALL GET_MY_ARG('time', cv_t)

      CASE('-p')
         CALL GET_MY_ARG('orbit ephem track file', cf_track)

      CASE('-m')
         CALL GET_MY_ARG('mesh_mask file', cf_mm)

      CASE('-f')
         CALL GET_MY_ARG('forced time vector construction for model', cs_force_tv_m)

      CASE('-g')
         CALL GET_MY_ARG('forced time vector construction for ephem', cs_force_tv_e)

      CASE('-n')
         l_orbit_file_is_nc = .TRUE.

      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
         CALL usage()

      END SELECT

   END DO

   IF ( (trim(cv_in) == '').OR.(trim(cf_in) == '') ) THEN
      PRINT *, ''
      PRINT *, 'You must at least specify input file (-i) and input variable (-v)!!!'
      CALL usage()
   END IF

   PRINT *, ''
   PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   PRINT *, ''

   PRINT *, ' * Input file = ', trim(cf_in)
   PRINT *, '   => associated variable names = ', trim(cv_in)
   PRINT *, '   => associated longitude/latitude/time = ', trim(cv_lon), ', ', trim(cv_lat), ', ', trim(cv_t)
   PRINT *, '   => mesh_mask file = ', trim(cf_mm)
   PRINT *, ' * Output file = ', trim(cf_out)

   PRINT *, ''


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

   ALLOCATE ( xlont(ni,nj), xlatt(ni,nj), xdum2d(ni,nj) )
   PRINT *, ''



   !! testing variable dimensions
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~
   CALL DIMS(cf_in, cv_in, ni1, nj1, nk, Ntm)

   IF ( (ni1/=ni).AND.(nj1/=nj) ) THEN
      PRINT *, 'ERROR: dimension of ',trim(cv_in), 'does not agree with lon/lat' ; STOP
   END IF

   IF ( nk < 1 ) nk = 1

   IF ( Ntm < 1 ) THEN
      PRINT *, 'ERROR: ',trim(cv_in),' must have at least a time record!' ; STOP
   END IF


   PRINT *, 'Dimension for ',trim(cv_in),':'
   PRINT *, '   => ni =', ni ;   PRINT *, '   => nj =', nj
   PRINT *, '   => nk =', nk ;   PRINT *, '   => Ntm =', Ntm
   PRINT *, ''

   ALLOCATE ( xvar(ni,nj), xvar1(ni,nj), xvar2(ni,nj), xslp(ni,nj), mask(ni,nj), vdepth(nk,1), vt_model(Ntm) )

   IF ( lregin ) THEN
      PRINT *, 'Regular case not supported yet! Priority to ORCA grids...'
      STOP
   END IF




   !! The first important stage is to compare time slices in the OGCM 2D input field
   !! w.r.t the one in the ephem orbit file!
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
   CALL GET_VAR_INFO(cf_in, cv_t, cunit, cdum)
   tut_model  = GET_TIME_UNIT_T0(TRIM(cunit))
   PRINT *, ' *** Unit and reference time in model file:'
   PRINT *, tut_model

   IF ( l_orbit_file_is_nc ) THEN
      CALL GET_VAR_INFO(cf_track, 'time', cunit, cdum)
      tut_ephem  = GET_TIME_UNIT_T0(TRIM(cunit))
      PRINT *, ' *** Unit and reference time in ephem file:'
      PRINT *, tut_ephem
   END IF
   PRINT *, ''





   !! Getting coordinates
   !! ~~~~~~~~~~~~~~~~~~~

   IF ( nk > 1 ) CALL GETVAR_1D(cf_in, cv_z, vdepth(:,1))


   IF ( TRIM(cs_force_tv_m) /= '' ) THEN
      !! Building new time vector!
      idx = SCAN(TRIM(cs_force_tv_m),',')
      cdum = cs_force_tv_m(1:idx-1)
      READ(cdum,'(f)') rt0
      cdum = cs_force_tv_m(idx+1:)
      READ(cdum,'(f)') rdt
      PRINT *, ' *** MODEL: OVERIDING time vector with t0 and dt =', REAL(rt0,4), REAL(rdt,4)
      DO jt=1, Ntm
         vt_model(jt) = rt0 + REAL(jt-1)*rdt
      END DO
   ELSE
      !! Reading it in input file:
      CALL GETVAR_1D(cf_in, cv_t, vt_model)
   END IF

   IF ( l_debug ) THEN
      PRINT *, ''
      PRINT *, 'Time vector in NEMO input file is (s), (h), (d):'
      DO jt=1, Ntm
         PRINT *, vt_model(jt), vt_model(jt)/3600., vt_model(jt)/(3600.*24.)
      END DO
      PRINT *, ''
      PRINT *, ''
   END IF



   !! Getting longitude, latitude and mask in mesh_mask file:
   ! Longitude array:
   CALL GETVAR_2D   (i0, j0, cf_mm, cv_lon, 0, 0, 0, xdum2d)
   xlont(:,:) = xdum2d(:,:) ; i0=0 ; j0=0
   WHERE ( xdum2d < 0. ) xlont = xlont + 360.0_8
   ! Latitude array:
   CALL GETVAR_2D   (i0, j0, cf_mm, cv_lat, 0, 0, 0, xdum2d)
   xlatt(:,:) = xdum2d(:,:) ; i0=0 ; j0=0
   !! 3D LSM
   CALL GETMASK_2D(cf_mm, cv_mt, mask, jlev=1)



   !! Reading along-track from file:

   INQUIRE(FILE=TRIM(cf_track), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: please provide the file containing definition of orbit ephem track'; STOP
   END IF

   IF ( .NOT. l_orbit_file_is_nc ) THEN
      !! Getting number of lines:
      Nte = -1 ; io = 0
      OPEN (UNIT=13, FILE=TRIM(cf_track))
      DO WHILE (io==0)
         READ(13,*,iostat=io)
         Nte = Nte + 1
      END DO
      PRINT*, Nte, ' points in '//TRIM(cf_track)//'...'
      ALLOCATE ( Xtar(1,Nte), Ytar(1,Nte), vt_ephem(Nte), Ztar(1,Nte) )
      !!
      REWIND(13)
      DO jte = 1, Nte
         READ(13,*) vt_ephem(jte), Xtar(1,jte), Ytar(1,jte)
      END DO
      CLOSE(13)

   ELSE
      PRINT *, ''
      PRINT *, 'NetCDF orbit ephem!'
      CALL DIMS(cf_track, 'time', Nte, nj1, nk, ni1)
      PRINT *, ' *** Nb. time records in NetCDF ephem file:', Nte
      ALLOCATE ( Xtar(1,Nte), Ytar(1,Nte), vt_ephem(Nte), Ztar(1,Nte))
      CALL GETVAR_1D(cf_track, 'time', vt_ephem)
      CALL GETVAR_1D(cf_track, 'longitude', Xtar(1,:))
      CALL GETVAR_1D(cf_track, 'latitude',  Ytar(1,:))
      CALL GETVAR_1D(cf_track, 'sla_unfiltered',  Ztar(1,:))
      PRINT *, 'Done!'; PRINT *, ''
   END IF


   IF ( TRIM(cs_force_tv_e) /= '' ) THEN
      !! Building new time vector!
      idx = SCAN(TRIM(cs_force_tv_e),',')
      cdum = cs_force_tv_e(1:idx-1)
      READ(cdum,'(f)') rt0
      cdum = cs_force_tv_e(idx+1:)
      READ(cdum,'(f)') rdt
      PRINT *, ' *** EPHEM: OVERIDING time vector with t0 and dt =', REAL(rt0,4), REAL(rdt,4)
      DO jt=1, Nte
         vt_ephem(jt) = rt0 + REAL(jt-1)*rdt
         !PRINT *, ' vt_ephem(jt)= ', vt_ephem(jt)
      END DO
   END IF




   nib = ni ; njb = nj ; ji_min=1 ; ji_max=ni ; jj_min=1 ; jj_max=nj

   ALLOCATE ( Ftrack(Nte), Ftrack_ephem(Nte), Fmask(Nte), xcmask(Nte), vposition(Nte,1) )


   !ctrack = TRIM(cf_track(1:LEN(cf_track)-4))
   WRITE(cf_out, '("track_",a,"_",a,".nc")') TRIM(cv_in), TRIM(cf_track)






   !PRINT *, ''
   !PRINT *, 'First time record for model:', vt_model(1)
   !itime = to_epoch_time_scalar( tut_model, vt_model(1) )
   !PRINT *, '     ==> in epoch time =>',  itime
   !PRINT *, ''
   !PRINT *, 'Last time record for model:', vt_model(Ntm)
   !itime = to_epoch_time_scalar( tut_model, vt_model(Ntm) )
   !PRINT *, '     ==> in epoch time =>',  itime

   !PRINT *, '' ; PRINT *, ''

   !PRINT *, 'First time record for ephem:', vt_ephem(1)
   !itime = to_epoch_time_scalar( tut_ephem, vt_ephem(1), dt=0.1_8 )
   !PRINT *, '     ==> in epoch time =>',  itime
   !PRINT *, ''
   !PRINT *, 'Last time record for ephem:', vt_ephem(Nte)
   !itime = to_epoch_time_scalar( tut_ephem, vt_ephem(Nte), dt=0.1_8 )
   !PRINT *, '     ==> in epoch time =>',  itime

   !PRINT *, ''

   !!
   !! Converting time vectors to epoch:
   !CALL time_vector_to_epoch_time( tut_ephem, vt_ephem )
   !CALL time_vector_to_epoch_time( tut_model, vt_model )


   !CALL to_epoch_time_vect( tut_model, vt_model )
   !PRINT *, vt_model(:)
   !PRINT *, ''


   PRINT *, ''
   PRINT *, ' Time vector in ephem file:'
   CALL to_epoch_time_vect( tut_ephem, vt_ephem, l_dt_below_sec=.true. )
   !PRINT *, vt_ephem(:)
   PRINT *, ''
   PRINT *, ''


   IF ( .NOT. l_debug_mapping ) THEN
      PRINT *, ' Time vector in model file:'
      CALL to_epoch_time_vect( tut_model, vt_model, l_dt_below_sec=.FALSE. )
      !PRINT *, vt_model(:)
      PRINT *, ''




      t_min_e = MINVAL(vt_ephem)
      t_max_e = MAXVAL(vt_ephem)
      t_min_m = MINVAL(vt_model)
      t_max_m = MAXVAL(vt_model)

      PRINT *, ''
      PRINT *, ' *** Max min time for ephem:', t_min_e, t_max_e
      PRINT *, ' *** Max min time for model:', t_min_m, t_max_m
      PRINT *, ''

      IF ( (t_min_m >= t_max_e).OR.(t_min_e >= t_max_m).OR.(t_max_m <= t_min_e).OR.(t_max_e <= t_min_m) ) THEN
         PRINT *, ' No time overlap for Model and Ephem file! '
         STOP
      END IF

      t_min = MAX(t_min_e, t_min_m)
      t_max = MIN(t_max_e, t_max_m)
      PRINT *, ' *** Time overlap for Model and Ephem file:', NINT(t_min), NINT(t_max)


      !! Findin when we can start and stop when scanning the ephem file:
      !! it1, it2
      DO it1 = 1, Nte-1
         IF ( (vt_ephem(it1) <= t_min).AND.(vt_ephem(it1+1) > t_min) ) EXIT
      END DO
      it1 = it1 + 1
      DO it2 = it1, Nte-1
         IF ( (vt_ephem(it2) <= t_max).AND.(vt_ephem(it2+1) > t_max) ) EXIT
      END DO

      Nten = it2 - it1 + 1

      PRINT *, ' it1, it2 =',it1, it2
      PRINT *, Nten, '  out of ', Nte
      PRINT *, ' => ', vt_ephem(it1), vt_ephem(it2)
      PRINT *, ''

   END IF ! IF ( .NOT. l_debug_mapping )

   !ALLOCATE ( Xtar(1,Nten), Ytar(1,Nten), vt_e(Nten), Ztar(1,Nten) )

   ALLOCATE ( IMETRICS(1,Nten,3), RAB(1,Nten,2), IPB(1,Nten) )



   IF ( l_debug_mapping ) THEN
      it1  = 1
      it2  = Nte
      Nten = Nte
   END IF


   !! Main time loop is on time vector in ephem file!


   !ALLOCATE ( Xpt(1,1), Ypt(1,1), Zpt4(1,1) )

   INQUIRE(FILE='mapping.nc', EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, ' *** Creating mapping file...'
      CALL MAPPING_BL(-1, xlont, xlatt, Xtar(:,it1:it2), Ytar(:,it1:it2), 'mapping.nc')
      PRINT *, ' *** Done!'; PRINT *, ''
   ELSE
      PRINT *, ' *** File "mapping.nc" found in current directory, using it!'
      PRINT *, ''
   END IF

   CALL RD_MAPPING_AB('mapping.nc', IMETRICS, RAB, IPB)
   PRINT *, ''; PRINT *, ' *** Mapping and weights read into "mapping.nc"'; PRINT *, ''


   !lulu
   !! Showing iy in file mask_+_nearest_points.nc:
   IF ( l_debug ) THEN
      ALLOCATE (JIidx(1,Nten) , JJidx(1,Nten) , Ftrack_np(Nten) )
      !! Finding and storing the nearest points of NEMO grid to ephem points:
      !CALL FIND_NEAREST_POINT(Xtar, Ytar, xlont, xlatt,  JIidx, JJidx)
      JIidx(1,:) = IMETRICS(1,:,1)
      JJidx(1,:) = IMETRICS(1,:,2)
      ALLOCATE ( show_track(nib,njb) )
      show_track(:,:) = 0.
      DO jte = 1, Nten
         show_track(JIidx(1,jte), JJidx(1,jte)) = REAL(jte,4)
      END DO
      WHERE (mask == 0) show_track = -9999.
      CALL PRTMASK(REAL(show_track(:,:),4), 'mask_+_nearest_points.nc', 'mask', xlont, xlatt, 'lon0', 'lat0', rfill=-9999.)
      DEALLOCATE ( show_track )
   END IF

   IF ( l_debug_mapping ) STOP'l_debug_mapping'

   STOP'#lolo'


   !STOP 'mapping done!'


   IF ( l_debug ) Ftrack_np(:) = -9999.
   Ftrack_ephem(:) = -9999.
   Ftrack(:) = -9999.
   Fmask(:) = -9999.

   jt_s = 1 ; ! time step model!

   jtm_1_o = -100
   jtm_2_o = -100

   DO jte = 1, Nte
      !!
      rt = vt_ephem(jte)
      PRINT *, 'Treating ephem time =>', rt
      !!
      IF ( (rt >= t_min_m).AND.(rt < t_max_m) ) THEN
         !!
         !! Two surrounding time records in model file => jtm_1 & jtm_2
         DO jt=jt_s, Ntm-1
            IF ( (rt >= vt_model(jt)).AND.(rt < vt_model(jt+1)) ) EXIT
         END DO
         !!
         jtm_1 = jt
         jtm_2 = jt+1
         PRINT *, ' rt, vt_model(jtm_1), vt_model(jtm_2) =>', rt, vt_model(jtm_1), vt_model(jtm_2)
         !!
         !! If first time we have these jtm_1 & jtm_2, getting the two surrounding fields:
         IF ( (jtm_1>jtm_1_o).AND.(jtm_2>jtm_2_o) ) THEN
            IF ( jtm_1_o == -100 ) THEN
               PRINT *, 'Reading field '//TRIM(cv_in)//' in '//TRIM(cf_in)//' at jtm_1=', jtm_1
               CALL GETVAR_2D(id_f1, id_v1, cf_in, cv_in, Ntm, 0, jtm_1, xvar1)
            ELSE
               PRINT *, 'Getting field '//TRIM(cv_in)//' at jtm_1=', jtm_1,' from previous jtm_2 !'
               xvar1(:,:) = xvar2(:,:)
            END IF
            PRINT *, 'Reading field '//TRIM(cv_in)//' in '//TRIM(cf_in)//' at jtm_2=', jtm_2
            CALL GETVAR_2D(id_f1, id_v1, cf_in, cv_in, Ntm, 0, jtm_2, xvar2)
            xslp = (xvar2 - xvar1) / (vt_model(jtm_2) - vt_model(jtm_1)) ! slope...

         END IF

         !! Linear interpolation of field at time rt:
         xvar(:,:) = xvar1(:,:) + xslp(:,:)*(rt - vt_model(jtm_1))

         !! Performing bilinear interpolation:
         iP       = IMETRICS(1,jte,1)
         jP       = IMETRICS(1,jte,2)
         iquadran = IMETRICS(1,jte,3)

         alpha    = RAB(1,jte,1)
         beta     = RAB(1,jte,2)

         IF ( (iP == INT(rflg)).OR.(jP == INT(rflg)) ) THEN
            Ftrack(jte) = -9999. ; ! masking
            Ftrack_ephem(jte) = -9999. ; ! masking
            Fmask(jte) = -9999. ; ! masking
         ELSE
            !! INTERPOLATION !
            Ftrack(jte) = INTERP_BL(-1, iP, jP, iquadran, alpha, beta, REAL(xvar,8))
            Fmask(jte)  = INTERP_BL(-1, iP, jP, iquadran, alpha, beta, REAL(mask,8))
            !!
            Ftrack_ephem(jte) = Ztar(1,jte) ! Input ephem data
            !!
         END IF

         IF ( l_debug ) Ftrack_np(jte) =  xvar(JIidx(1,jte),JJidx(1,jte)) ! NEAREST POINT interpolation

         jtm_1_o = jtm_1
         jtm_2_o = jtm_2
         jt_s    = jtm_1 ! so we do not rescan from begining...

      END IF

   END DO


   !WHERE ( IPB(1,:) > 0 ) Ftrack = -9999.
   WHERE ( Ftrack > 1.E9 ) Ftrack = -9999.
   WHERE ( Fmask < 1.    ) Ftrack = -9999.
   WHERE ( Ftrack_ephem > 1.E9 ) Ftrack_ephem = -9999.

   CALL PT_SERIES(vt_ephem, REAL(Ftrack,4), 'result.nc', 'time', cv_in, 'boo', 'ta mere', -9999.)
   CALL PT_SERIES(vt_ephem, REAL(Fmask,4), 'result_mask.nc', 'time', 'lsm', 'boo', 'ta mere', -9999.)

   CALL PT_SERIES(vt_ephem, REAL(Ftrack_ephem,4), 'data_ephem.nc', 'time', cv_in, 'boo', 'ta mere', -9999.)
   CALL PT_SERIES(vt_ephem, REAL(Ytar(1,:),4), 'lat_ephem.nc', 'time', 'latitude', 'boo', 'ta mere', -9999.)
   CALL PT_SERIES(vt_ephem, REAL(Xtar(1,:),4), 'lon_ephem.nc', 'time', 'longitude', 'boo', 'ta mere', -9999.)

   IF ( l_debug ) THEN
      DEALLOCATE ( JIidx, JJidx )
      WHERE ( Ftrack_np > 1.E9 ) Ftrack_np = -9999.
      WHERE ( Fmask < 1.    ) Ftrack_np = -9999.
      CALL PT_SERIES(vt_ephem, REAL(Ftrack_np,4), 'result_np.nc', 'time', cv_in, 'boo', 'ta mere', -9999.)
   END IF

   !DO jte = 1, Nte
   !  PRINT *, ''

   STOP 'LOLO: stop for now...'

   DEALLOCATE ( Xtar, Ytar, Ztar )
   DEALLOCATE ( Ftrack, xcmask, vposition )
   DEALLOCATE ( xlont, xlatt, xvar, xvar1, xvar2, xslp, mask ) !, xtmp4 )
   !lolo




CONTAINS



   SUBROUTINE BUILD_TRACK()

      IF ( imax < imin ) THEN ! switching points,
         isav = imin ; jsav = jmin
         imin = imax ; jmin = jmax
         imax = isav ; jmax = jsav
      END IF

      !! Creating the local box with 4 extrapoints:
      !!   => lolo: should mind limits (1,1,ni,nj) !!!!
      ji_min = MIN(imin,imax) - 2
      jj_min = MIN(jmin,jmax) - 2
      ji_max = MAX(imin,imax) + 2
      jj_max = MAX(jmin,jmax) + 2
      PRINT *, 'ji_min, ji_max, jj_min, jj_max:'; PRINT *, ji_min, ji_max, jj_min, jj_max
      nib = ji_max - ji_min + 1
      njb = jj_max - jj_min + 1
      PRINT *, 'nib, njb =', nib, njb

      !lolo: change with actual min and max (use min_val ..)
      lon_max = MAX(xlont(ji_max,jj_min),xlont(ji_max,jj_max))
      lon_min = MIN(xlont(ji_min,jj_min),xlont(ji_min,jj_max))
      lat_max = MAX(xlatt(ji_max,jj_max),xlatt(ji_min,jj_max))
      lat_min = MIN(xlatt(ji_max,jj_min),xlatt(ji_min,jj_min))

      PRINT *, ''
      dlon = lon_max - lon_min ; PRINT *, 'long. range =', dlon
      dlat = lat_max - lat_min ; PRINT *, 'latg. range =', dlat
      dang = SQRT(dlon*dlon + dlat*dlat) ; PRINT *, 'Ang. range =', dang

      Nte = INT(dang/res) + 1
      IF ( MOD(Nte,2) == 0 ) Nte = Nte - 1 ! we want odd integer...
      PRINT *, 'Number of points to create on segment:', Nte ; PRINT *, ''

      ALLOCATE ( Xtar(1,Nte), Ytar(1,Nte), Ztar(1,Nte) )
      ALLOCATE ( xcmask(Nte), vposition(Nte,1) )

      IF ( ABS(dlon) < 1.E-12 ) THEN
         PRINT *, 'ERROR: Section seems to be vertical!'; STOP
      END IF

      rA = (xlatt(imax,jmax) - xlatt(imin,jmin))/ dlon
      rB = xlatt(imin,jmin) - rA*xlont(imin,jmin)

      PRINT *, 'rA, rB = ', rA, rB
      PRINT *, 'Lat1 =', rA*xlont(imin,jmin) + rB
      PRINT *, 'Lat2 =', rA*xlont(imax,jmax) + rB

      dlon = dlon / (Nte-1) ;  ; PRINT *, 'dlon =', dlon

      DO jte = 1, Nte
         Xtar(1,jte) = xlont(imin,jmin) + (jte-1)*dlon
         Ytar(1,jte) = rA*Xtar(1,jte) + rB
      END DO

      ! Only positive longitudes:
      WHERE( Xtar(1,:) < 0. ) Xtar(1,:) = Xtar(1,:) + 360.

      ! BAD! lolo should find good metrics...
      vposition(:,1) = Xtar(1,:)

      ALLOCATE( xlont(nib,njb), xlatt(nib,njb), mask(nib,njb)) !, xtmp4(nib,njb) )

      xlont = xlont(ji_min:ji_max,jj_min:jj_max)
      xlatt = xlatt(ji_min:ji_max,jj_min:jj_max)

      mask = mask(ji_min:ji_max,jj_min:jj_max)

   END SUBROUTINE BUILD_TRACK










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


END PROGRAM INTERP_TO_EPHEM




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
   WRITE(6,*) ' -n                   => file containing orbit ephem is in NetCDF'
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
   WRITE(6,*) ' -g  <t0,dt>          => overide time vector in ephem file with one of same length'
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
