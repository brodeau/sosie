PROGRAM INTERP_TO_EPHEM

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
      &   l_akima = .true., &
      &   l_bilin = .false.
   !!
   LOGICAL :: &
      &      l_file_is_ascii = .FALSE., &
      &      l_file_is_nc    = .FALSE.
   !!
   REAL(8), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER :: nval, io, idx, iP, jP, iquadran
   !!
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xtar, Ytar, Xpt, Ypt
   REAL(4), DIMENSION(:,:), ALLOCATABLE :: Ztar4, Zpt4
   !!
   !! Coupe stuff:
   REAL(8), DIMENSION(:), ALLOCATABLE :: Ftrack, Fmask, Ftrack_np
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

   CHARACTER(len=256)  :: cr, ctrack, cdum
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
      &    cs_force_tv=''
   !!
   INTEGER      :: &
      &    jarg,   &
      &    jd,  &
      &    i0, j0, ifo, ivo,   &
      &    ni, nj, nk=0, Ntm=0, &
      &    ni1, nj1, ni2, nj2, &
      &    iargc, id_f1, id_v1
   !!
   !!
   INTEGER :: imin, imax, jmin, jmax, isav, jsav, ji_min, ji_max, jj_min, jj_max, nib, njb

   REAL(8), DIMENSION(:,:), ALLOCATABLE :: xvar, xvar1, xvar2, xslp

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xdum2d
   REAL(8), DIMENSION(:,:), ALLOCATABLE ::    &
      &    xlont, xlatt
   !!
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJidx, JIidx, mask_show_track    ! debug
   !!
   INTEGER(2), DIMENSION(:,:), ALLOCATABLE :: mask
   !!
   INTEGER :: jt, jte, jl, jt_s, jtm_1, jtm_2, jtm_1_o, jtm_2_o
   !!
   REAL(8) :: rA, rB, dlon, dlat, dang, lon_min, lon_max, lat_min, lat_max, rt, rt0, rdt, &
      &       t_min_e, t_max_e, t_min_m, t_max_m, &
      &       alpha, beta
   !!
   CHARACTER(LEN=2), DIMENSION(12), PARAMETER :: &
      &            clist_opt = (/ '-h','-v','-x','-y','-z','-t','-i','-p','-a','-n','-m','-f' /)

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
         CALL GET_MY_ARG('forced time vector construction', cs_force_tv)

      CASE('-a')
         l_file_is_ascii = .true.

      CASE('-n')
         l_file_is_nc    = .true.

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


   IF ( ((.NOT. l_file_is_ascii).AND.(.NOT. l_file_is_nc)).OR.(l_file_is_ascii .AND.  l_file_is_nc) ) THEN
      PRINT *, ''
      PRINT *, 'When using an orbit, you must specify whether input file is in ASCII (-a) or netcdf (-n)!'
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
   CALL DIMS(cf_in, 'nav_lon', ni1, nj1, nk, Ntm)
   CALL DIMS(cf_in, 'nav_lat', ni2, nj2, nk, Ntm)

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





   !! Getting coordinates
   !! ~~~~~~~~~~~~~~~~~~~

   IF ( nk > 1 ) CALL GETVAR_1D(cf_in, cv_z, vdepth(:,1))


   IF ( TRIM(cs_force_tv) /= '' ) THEN
      !! Building new time vector!
      idx = SCAN(TRIM(cs_force_tv),',')
      cdum = cs_force_tv(1:idx-1)
      READ(cdum,'(f)') rt0
      cdum = cs_force_tv(idx+1:)
      READ(cdum,'(f)') rdt
      PRINT *, ' *** OVERIDING time vector with t0 and dt =', REAL(rt0,4), REAL(rdt,4)
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

   IF ( l_file_is_ascii ) THEN
      !! Getting number of lines:
      nval = -1 ; io = 0
      OPEN (UNIT=13, FILE=TRIM(cf_track))
      DO WHILE (io==0)
         READ(13,*,iostat=io)
         nval = nval + 1
      END DO
      PRINT*, nval, ' points in '//TRIM(cf_track)//'...'
      ALLOCATE ( Xtar(1,nval), Ytar(1,nval), vt_ephem(nval) )
      !!
      REWIND(13)
      DO jl = 1, nval
         READ(13,*) vt_ephem(jl), Xtar(1,jl), Ytar(1,jl)
      END DO
      CLOSE(13)

   ELSE
      PRINT *, 'Netcdf orbit ephem case not supported yet!!!'
      STOP
   END IF



   nib = ni ; njb = nj ; ji_min=1 ; ji_max=ni ; jj_min=1 ; jj_max=nj

   ALLOCATE ( Ztar4(1,nval), Ftrack(nval), Fmask(nval), xcmask(nval), vposition(nval,1) )

   ALLOCATE ( IMETRICS(1,nval,3), RAB(1,nval,2), IPB(1,nval) )

   !ctrack = TRIM(cf_track(1:LEN(cf_track)-4))
   WRITE(cf_out, '("track_",a,"_",a,".nc")') TRIM(cv_in), TRIM(cf_track)





   !! Showing iy in file mask_+_nearest_points.nc:
   IF ( l_debug ) THEN
      ALLOCATE (JIidx(1,nval) , JJidx(1,nval) , Ftrack_np(nval) )
      !! Finding and storing the nearest points of NEMO grid to ephem points:
      CALL FIND_NEAREST_POINT(Xtar, Ytar, xlont, xlatt,  JIidx, JJidx)
      ALLOCATE ( mask_show_track(nib,njb) )
      mask_show_track(:,:) = mask(:,:)
      DO jd = 1, nval
         mask_show_track(JIidx(1,jd), JJidx(1,jd)) = -5
      END DO
      CALL PRTMASK(REAL(mask_show_track(:,:),4), 'mask_+_nearest_points.nc', 'mask', xlont, xlatt, 'lon0', 'lat0')
      DEALLOCATE ( mask_show_track )
   END IF


   !PRINT *, ''
   !PRINT *, ' Time vector in ephem file:'
   !PRINT *, vt_ephem(:)
   PRINT *, ''
   PRINT *, ' Time vector for model file:'
   PRINT *, vt_model(:)
   PRINT *, ''

   t_min_e = MINVAL(vt_ephem)
   t_max_e = MAXVAL(vt_ephem)
   t_min_m = MINVAL(vt_model)
   t_max_m = MAXVAL(vt_model)


   !! Main time loop is on time vector in ephem file!


   ALLOCATE ( Xpt(1,1), Ypt(1,1), Zpt4(1,1) )

   INQUIRE(FILE='mapping.nc', EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, ' *** Creating mapping file...'
      CALL MAPPING_BL(-1, xlont, xlatt, Xtar, Ytar, 'mapping.nc')
      PRINT *, ' *** Done!'; PRINT *, ''
   ELSE
      PRINT *, ' *** File "mapping.nc" found in current directory, using it!'
      PRINT *, ''
   END IF

   CALL RD_MAPPING_AB('mapping.nc', IMETRICS, RAB, IPB)
   PRINT *, ''; PRINT *, ' *** Mapping and weights read into "mapping.nc"'; PRINT *, ''
   !STOP 'mapping done!'


   IF ( l_debug ) Ftrack_np(:) = -9999.
   Ftrack(:) = -9999.
   Fmask(:) = -9999.

   jt_s = 1 ; ! time step model!

   jtm_1_o = -100
   jtm_2_o = -100

   DO jte = 1, nval
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
               CALL GETVAR_2D(id_f1, id_v1, cf_in, cv_in, Ntm, 0, jtm_1, xdum2d(:,:))
               xvar1(:,:) = REAL(xdum2d(:,:),8)
            ELSE
               PRINT *, 'Getting field '//TRIM(cv_in)//' at jtm_1=', jtm_1,' from previous jtm_2 !'
               xvar1(:,:) = xvar2(:,:)
            END IF
            PRINT *, 'Reading field '//TRIM(cv_in)//' in '//TRIM(cf_in)//' at jtm_2=', jtm_2
            CALL GETVAR_2D(id_f1, id_v1, cf_in, cv_in, Ntm, 0, jtm_2, xdum2d(:,:))
            xvar2(:,:) = REAL(xdum2d(:,:),8)
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
            Fmask(jte) = -9999. ; ! masking
         ELSE
            !! INTERPOLATION !
            Ftrack(jte) = INTERP_BL(-1, iP, jP, iquadran, alpha, beta, xvar)
            Fmask(jte)  = INTERP_BL(-1, iP, jP, iquadran, alpha, beta, REAL(mask,8))
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

   CALL PT_SERIES(vt_ephem, REAL(Ftrack,4), 'result.nc', 'time', cv_in, 'boo', 'ta mere', -9999.)
   CALL PT_SERIES(vt_ephem, REAL(Fmask,4), 'result_mask.nc', 'time', 'lsm', 'boo', 'ta mere', -9999.)

   IF ( l_debug ) THEN
      DEALLOCATE ( JIidx, JJidx )
      WHERE ( Ftrack_np > 1.E9 ) Ftrack_np = -9999.
      WHERE ( Fmask < 1.    ) Ftrack_np = -9999.
      CALL PT_SERIES(vt_ephem, REAL(Ftrack_np,4), 'result_np.nc', 'time', cv_in, 'boo', 'ta mere', -9999.)
   END IF

   !DO jte = 1, nval
    !  PRINT *, ''

   STOP 'LOLO: stop for now...'









   !! Filling arrays for "local box" (defined by the section)
   !!
   !DO jt = 1, Ntm
   !   !!
   !   xtmp4(:,:) = REAL(xvar(ji_min:ji_max,jj_min:jj_max,jt), 8)
   !   !!
   !   !! Extrapolating values over continents:
   !   CALL DROWN(-1, xtmp4(:,:), mask(:,:)) ! ORCA
   !   !!
   !   xvar(:,:,jt) = REAL(xtmp4, 8)
   !   !!
   !END DO


   l_first_call_interp_routine = .TRUE.
   ! Interpolating the mask on target section
   IF ( l_akima ) CALL AKIMA_2D(-1, xlont, xlatt, REAL(mask(:,:),4), &
      &                        Xtar, Ytar, Ztar4,  icall=1)
   !!
   IF ( l_bilin) CALL BILIN_2D(-1, xlont, xlatt, REAL(mask(:,:),4), &
      &                        Xtar, Ytar, Ztar4, TRIM(ctrack))
   !!
   xcmask(:) = Ztar4(1,:)

   ifo=0 ; ivo=0

   l_first_call_interp_routine = .TRUE.
   ! Interpolating the Ntm snapshots of field on target section
   DO jt = 1, Ntm
      !!
      IF ( l_akima ) CALL AKIMA_2D(-1, xlont, xlatt, REAL(xvar(:,:),4), &
         &                           Xtar,    Ytar,    Ztar4,  icall=1)
      !!
      IF ( l_bilin ) CALL BILIN_2D(-1, xlont, xlatt, REAL(xvar(:,:),4), &
         &                           Xtar,    Ytar,    Ztar4, trim(ctrack))
      !!
      !xcoupe(:,jt) = Ztar4(1,:)
      !!
      !WHERE( xcmask < 0.25 ) xcoupe(:,jt) = -9999.

      !CALL P2D_T(ifo, ivo, Ntm, jt, vposition, vdepth, vt_model, xcoupe(:,jt), cf_out, &
      !   &       'position', 'profo', cv_t, cv_in, -9999.)

   END DO

   !LOLO: CALL PT_SERIES(vtime, vseries, cf_in, cv_t, cv_in, cunit, cln, vflag, &
   !LOLO: &                 lpack)



   !lulu
   !IF ( l_debug ) THEN
   !ivo=0 ; ifo=0
   !CALL P2D_T(ifo, ivo, 1, 1, vposition, vdepth, vt_model, xcmask(:), 'MASK_'//TRIM(cf_out), &
   !   &       'position', 'profo', cv_t, 'lsm', -9999.)
   !END IF




   PRINT *, 'File created => ', trim(cf_out)

   DEALLOCATE ( Xtar, Ytar, Ztar4 )
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

      nval = INT(dang/res) + 1
      IF ( MOD(nval,2) == 0 ) nval = nval - 1 ! we want odd integer...
      PRINT *, 'Number of points to create on segment:', nval ; PRINT *, ''

      ALLOCATE ( Xtar(1,nval), Ytar(1,nval), Ztar4(1,nval) )
      ALLOCATE ( xcmask(nval), vposition(nval,1) )

      IF ( ABS(dlon) < 1.E-12 ) THEN
         PRINT *, 'ERROR: Section seems to be vertical!'; STOP
      END IF

      rA = (xlatt(imax,jmax) - xlatt(imin,jmin))/ dlon
      rB = xlatt(imin,jmin) - rA*xlont(imin,jmin)

      PRINT *, 'rA, rB = ', rA, rB
      PRINT *, 'Lat1 =', rA*xlont(imin,jmin) + rB
      PRINT *, 'Lat2 =', rA*xlont(imax,jmax) + rB

      dlon = dlon / (nval-1) ;  ; PRINT *, 'dlon =', dlon

      DO jd = 1, nval
         Xtar(1,jd) = xlont(imin,jmin) + (jd-1)*dlon
         Ytar(1,jd) = rA*Xtar(1,jd) + rB
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
   WRITE(6,*) ' -a                   => file containing orbit ephem is in ASCII'
   WRITE(6,*) ''
   WRITE(6,*) ' -n                   => file containing orbit ephem is in NetCDF'
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
   WRITE(6,*) ' -h                   => Show this message'
   WRITE(6,*) ''
   !!
   !CLOSE(6)
   STOP
   !!
END SUBROUTINE usage
!!
