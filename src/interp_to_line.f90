PROGRAM INTERP_TO_LINE

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
   INTEGER :: nval, io, idx
   !!
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xtar, Ytar
   REAL(4), DIMENSION(:,:), ALLOCATABLE :: Ztar4
   !!
   !! Coupe stuff:
   REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: xcoupe
   REAL(4), DIMENSION(:,:),   ALLOCATABLE :: xcmask
   REAL(8), DIMENSION(:,:),   ALLOCATABLE :: vposition

   REAL(8), DIMENSION(:,:),   ALLOCATABLE :: vdepth
   REAL(8), DIMENSION(:),     ALLOCATABLE :: vt_model, vt_ephem   ! in seconds



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
      &    jd, jk, &
      &    i0, j0, ifo, ivo,   &
      &    ni, nj, nk=0, nt=0, &
      &    ni1, nj1, ni2, nj2, &
      &    iargc, id_f1, id_v1
   !!
   !!
   INTEGER :: imin, imax, jmin, jmax, isav, jsav, ji_min, ji_max, jj_min, jj_max, nib, njb

   REAL(4), DIMENSION(:,:,:,:), ALLOCATABLE :: xvar
   REAL(4), DIMENSION(:,:,:),   ALLOCATABLE :: xtmp4_b
   REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: xvar_b

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xdum2d
   REAL(8), DIMENSION(:,:), ALLOCATABLE ::    &
      &    xlont, xlatt, &
      &    xlont_b, xlatt_b
   !!
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJidx, JIidx, mask_show_track    ! debug
   !!
   INTEGER(2), DIMENSION(:,:,:), ALLOCATABLE :: mask, mask_b
   !!
   INTEGER :: jt, jl
   !!
   REAL(8) :: rA, rB, dlon, dlat, dang, lon_min, lon_max, lat_min, lat_max, rt0, rdt
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
   CALL DIMS(cf_in, 'nav_lon', ni1, nj1, nk, nt)
   CALL DIMS(cf_in, 'nav_lat', ni2, nj2, nk, nt)

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
   CALL DIMS(cf_in, cv_in, ni1, nj1, nk, nt)

   IF ( (ni1/=ni).AND.(nj1/=nj) ) THEN
      PRINT *, 'ERROR: dimension of ',trim(cv_in), 'does not agree with lon/lat' ; STOP
   END IF

   IF ( nk < 1 ) nk = 1

   IF ( nt < 1 ) THEN
      PRINT *, 'ERROR: ',trim(cv_in),' must have at least a time record!' ; STOP
   END IF


   PRINT *, 'Dimension for ',trim(cv_in),':'
   PRINT *, '   => ni =', ni ;   PRINT *, '   => nj =', nj
   PRINT *, '   => nk =', nk ;   PRINT *, '   => nt =', nt
   PRINT *, ''

   ALLOCATE ( xvar(ni,nj,nk,nt), mask(ni,nj,nk), vdepth(nk,1), vt_model(nt) )

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
      DO jt=1, nt
         vt_model(jt) = rt0 + REAL(jt-1)*rdt
      END DO
   ELSE
      !! Reading it in input file:
      CALL GETVAR_1D(cf_in, cv_t, vt_model)
   END IF
   
   IF ( l_debug ) THEN
      PRINT *, ''
      PRINT *, 'Time vector in NEMO input file is (s), (h), (d):'
      DO jt=1, nt
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
   CALL GETMASK_3D(cf_mm, cv_mt, mask, jz1=1, jz2=nk)



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
   ALLOCATE( xlont_b(nib,njb), xlatt_b(nib,njb), xvar_b(nib,njb,nk,nt), mask_b(nib,njb,nk), xtmp4_b(nib,njb,nk) )
   xlont_b(:,:) =  xlont(:,:)
   mask_b(:,:,:)  =  mask(:,:,:)
   xlatt_b(:,:) =  xlatt(:,:)

   ALLOCATE ( Ztar4(1,nval), xcoupe(nval,nk,nt), xcmask(nval,nk), vposition(nval,1), &
      &       JIidx(1,nval) , JJidx(1,nval) )

   !ctrack = TRIM(cf_track(1:LEN(cf_track)-4))
   WRITE(cf_out, '("track_",a,"_",a,".nc")') TRIM(cv_in), TRIM(cf_track)




   !! Finding and storing the nearest points of NEMO grid to ephem points:
   
   CALL FIND_NEAREST_POINT(Xtar, Ytar, xlont_b, xlatt_b,  JIidx, JJidx)

   !! Showing iy in file mask_+_nearest_points.nc:
   IF ( l_debug ) THEN
      ALLOCATE ( mask_show_track(nib,njb) )
      mask_show_track(:,:) = mask_b(:,:,1)
      DO jd = 1, nval
         mask_show_track(JIidx(1,jd), JJidx(1,jd)) = -5
      END DO
      CALL DUMP_2D_FIELD(REAL(mask_show_track(:,:),4), 'mask_+_nearest_points.nc', 'mask', xlont_b, xlatt_b, 'lon0', 'lat0')
      DEALLOCATE ( mask_show_track )
   END IF


   PRINT *, ''
   PRINT *, ' Time vector in ephem file:'
   PRINT *, vt_ephem(:)
   PRINT *, ''
   PRINT *, ' Time vector for model file:'
   PRINT *, vt_model(:)
   PRINT *, ''



   STOP 'LOLO: stop for now...'


   !! Filling xvar once for all, !BAD lolo, if few virtual memory on the machine...
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~
   DO jt = 1, nt
      !!
      PRINT *, ' *** Reading record', jt
      CALL GETVAR_3D(id_f1, id_v1, cf_in, cv_in, nt, jt, xvar(:,:,1:nk,jt), jz1=1, jz2=nk)
      !!
      !CALL DUMP_2D_FIELD(xvar(:,:,1,jt), TRIM(cv_in)//'_stage_1.nc', cv_in,   xlont, xlatt, 'lon0', 'lat0')
      !!
   END DO
   !!
   





   !! Filling arrays for "local box" (defined by the section)
   !!
   DO jt = 1, nt
      !!
      xtmp4_b(:,:,:) = REAL(xvar(ji_min:ji_max,jj_min:jj_max,:,jt), 8)
      !!
      !! Extrapolating values over continents:
      DO jk = 1, nk
         CALL DROWN(-1, xtmp4_b(:,:,jk), mask_b(:,:,jk)) ! ORCA
      END DO
      !!
      xvar_b(:,:,:,jt) = REAL(xtmp4_b, 8)
      !!
   END DO


   l_first_call_interp_routine = .TRUE.
   ! Interpolating the mask on target section
   DO jk = 1, nk
      !!
      IF ( l_akima ) CALL AKIMA_2D(-1, xlont_b, xlatt_b, REAL(mask_b(:,:,jk),4), &
         &                        Xtar, Ytar, Ztar4,  icall=1)
      !!
      IF ( l_bilin) CALL BILIN_2D(-1, xlont_b, xlatt_b, REAL(mask_b(:,:,jk),4), &
         &                        Xtar, Ytar, Ztar4, trim(ctrack))
      !!
      xcmask(:,jk) = Ztar4(1,:)
      !!
   END DO


   ifo=0 ; ivo=0

   l_first_call_interp_routine = .TRUE.
   ! Interpolating the nt snapshots of field on target section
   DO jt = 1, nt
      !!
      DO jk = 1, nk
         !!
         IF ( l_akima ) CALL AKIMA_2D(-1, xlont_b, xlatt_b, REAL(xvar_b(:,:,jk,jt),4), &
            &                           Xtar,    Ytar,    Ztar4,  icall=1)
         !!
         IF ( l_bilin ) CALL BILIN_2D(-1, xlont_b, xlatt_b, REAL(xvar_b(:,:,jk,jt),4), &
            &                           Xtar,    Ytar,    Ztar4, trim(ctrack))
         !!
         xcoupe(:,jk,jt) = Ztar4(1,:)
         !!
      END DO
      !!
      WHERE( xcmask < 0.25 ) xcoupe(:,:,jt) = -9999.

      CALL P2D_T(ifo, ivo, nt, jt, vposition, vdepth, vt_model, xcoupe(:,:,jt), cf_out, &
         &       'position', 'profo', cv_t, cv_in, -9999.)

   END DO

   !lulu
   IF ( l_debug ) THEN
      ivo=0 ; ifo=0
      CALL P2D_T(ifo, ivo, 1, 1, vposition, vdepth, vt_model, xcmask(:,:), 'MASK_'//TRIM(cf_out), &
         &       'position', 'profo', cv_t, 'lsm', -9999.)
   END IF


   l_first_call_interp_routine = .TRUE.

   PRINT *, 'File created => ', trim(cf_out)

   DEALLOCATE ( Xtar, Ytar, Ztar4 )
   DEALLOCATE ( xcoupe, xcmask, vposition )
   DEALLOCATE ( xlont_b, xlatt_b, xvar_b, mask_b, xtmp4_b )
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
      ALLOCATE ( xcoupe(nval,nk,nt), xcmask(nval,nk), vposition(nval,1) )

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

      ALLOCATE( xlont_b(nib,njb), xlatt_b(nib,njb), xvar_b(nib,njb,nk,nt), mask_b(nib,njb,nk), xtmp4_b(nib,njb,nk) )

      xlont_b = xlont(ji_min:ji_max,jj_min:jj_max)
      xlatt_b = xlatt(ji_min:ji_max,jj_min:jj_max)

      mask_b = mask(ji_min:ji_max,jj_min:jj_max,:)

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


END PROGRAM INTERP_TO_LINE








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
