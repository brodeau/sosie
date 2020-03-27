PROGRAM IJ_FROM_LON_LAT

   USE io_ezcdf
   USE mod_manip

   IMPLICIT NONE

   !! ************************ Configurable part ****************************
   LOGICAL, PARAMETER :: &
      &   l_debug = .FALSE., &
      &   l_drown_in = .FALSE. ! Not needed since we ignore points that are less than 1 point away from land... drown the field to avoid spurious values right at the coast!
   !!
   REAL(8), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER :: io,  npoints, jl
   !!
   REAL(8), DIMENSION(:,:), ALLOCATABLE :: vpt_lon, vpt_lat
   !!


   !! Coupe stuff:
   !REAL(8), DIMENSION(:), ALLOCATABLE :: Ftrack_mod, Ftrack_mod_np, Ftrack_obs, rcycle_obs, vdistance


   !! Grid, default name :
   CHARACTER(len=80) :: &
      &    cv_mod, &
      &    cv_lon = 'glamt',         & ! input grid longitude name, T-points
      &    cv_lat = 'gphit'            ! input grid latitude name,  T-points

   CHARACTER(len=256)  :: cr
   CHARACTER(len=512)  :: cdir_home, cdir_out, cdir_tmpdir, cdum, cconf
   !!
   !!
   !!******************** End of conf for user ********************************
   !!
   !!               ** don't change anything below **
   !!
   LOGICAL ::  &
      &     l_reg_src, &
      &     l_exist   = .FALSE.

   !!
   !!
   CHARACTER(len=400)  :: &
      &    cf_mod, &
      &    cf_ascii='file_in.txt'
   !!
   CHARACTER(len=512), DIMENSION(:), ALLOCATABLE :: cb_name
   !!
   INTEGER      :: &
      &    jarg, &
      &    idot, &
      &    i0, j0,  &
      &    ni, nj, Ntm=0, nk=0, &
      &    ni1, nj1, ni2, nj2, &
      &    iargc
   !!

   !!


   REAL(4), DIMENSION(:,:), ALLOCATABLE :: xdum_r4
   REAL(8), DIMENSION(:,:), ALLOCATABLE ::    &
      &    xlont, xlatt, xlont_tmp
   !!
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: JJidx, JIidx    ! debug
   !!
   !INTEGER :: jt, jt0, jtf, jt_s, jtm_1, jtm_2, jtm_1_o, jtm_2_o, jb, js
   !!
   
   CHARACTER(LEN=2), DIMENSION(5), PARAMETER :: &
      &            clist_opt = (/ '-h','-i','-x','-y','-p' /)

   REAL(8) :: lon_min_1, lon_max_1, lon_min_2, lon_max_2, lat_min, lat_max

   REAL(8) :: lon_min_trg, lon_max_trg, lat_min_trg, lat_max_trg



   CALL GET_ENVIRONMENT_VARIABLE("HOME", cdir_home)
   CALL GET_ENVIRONMENT_VARIABLE("TMPDIR", cdir_tmpdir)


   cdir_out = TRIM(cdir_tmpdir)//'/EXTRACTED_BOXES' ! where to write data!



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

      CASE('-x')
         CALL GET_MY_ARG('longitude', cv_lon)

      CASE('-y')
         CALL GET_MY_ARG('latitude', cv_lat)

      CASE('-p')
         CALL GET_MY_ARG('list lon,lat ascii file', cf_ascii)

      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
         CALL usage()

      END SELECT

   END DO

   IF ( (trim(cv_mod) == '').OR.(trim(cf_mod) == '') ) THEN
      PRINT *, ''
      PRINT *, 'You must at least specify input file (-i) !!!'
      CALL usage()
   END IF

   PRINT *, ''
   PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   PRINT *, ''

   PRINT *, ' * Input file = ', trim(cf_mod)
   PRINT *, '   => associated variable names = ', trim(cv_mod)
   PRINT *, '   => associated longitude/latitude/time = ', trim(cv_lon), ', ', trim(cv_lat)


   PRINT *, ''

   !! Name of config: lulu
   idot = SCAN(cf_mod, '/', back=.TRUE.)
   cdum = cf_mod(idot+1:)
   idot = SCAN(cdum, '.', back=.TRUE.)
   cconf = cdum(:idot-1)

   PRINT *, ' *** CONFIG: cconf ='//TRIM(cconf) ; PRINT *, ''


   !! testing longitude and latitude
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INQUIRE(FILE=TRIM(cf_mod), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: please provide the file with model grid!'
      call usage()
   END IF

   CALL DIMS(cf_mod, cv_lon, ni1, nj1, nk, Ntm)
   CALL DIMS(cf_mod, cv_lat, ni2, nj2, nk, Ntm)

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

   ALLOCATE ( xlont(ni,nj), xlatt(ni,nj), xdum_r4(ni,nj) )
   PRINT *, ''



   !! testing variable dimensions
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~

   IF ( l_reg_src ) THEN
      PRINT *, 'Regular case not supported yet! Priority to ORCA grids...'
      STOP
   END IF




   ALLOCATE ( xlont_tmp(ni,nj) )

   !! Getting model longitude & latitude:
   ! Longitude array:
   CALL GETVAR_2D(i0, j0, cf_mod, cv_lon, 0, 0, 0, xdum_r4)
   xlont(:,:) = xdum_r4(:,:) ; i0=0 ; j0=0
   !!


   !! Min an max lon:
   lon_min_1 = MINVAL(xlont)
   lon_max_1 = MAXVAL(xlont)
   PRINT *, ' *** Minimum longitude on model grid before : ', lon_min_1
   PRINT *, ' *** Maximum longitude on model grid before : ', lon_max_1
   !!
   xlont_tmp = xlont
   WHERE ( xdum_r4 < 0. ) xlont_tmp = xlont + 360.0_8
   !!
   lon_min_2 = MINVAL(xlont_tmp)
   lon_max_2 = MAXVAL(xlont_tmp)
   PRINT *, ' *** Minimum longitude on model grid: ', lon_min_2
   PRINT *, ' *** Maximum longitude on model grid: ', lon_max_2


   ! Latitude array:
   CALL GETVAR_2D   (i0, j0, cf_mod, cv_lat, 0, 0, 0, xdum_r4)
   xlatt(:,:) = xdum_r4(:,:) ; i0=0 ; j0=0

   !! Min an max lat:
   lat_min = MINVAL(xlatt)
   lat_max = MAXVAL(xlatt)
   PRINT *, ' *** Minimum latitude on model grid : ', lat_min
   PRINT *, ' *** Maximum latitude on model grid : ', lat_max











   !! Reading along-track from file:

   INQUIRE(FILE=TRIM(cf_ascii), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: please provide the file containing lists of lat lon...'
      call usage()
   END IF



   npoints = -1 ; io = 0
   OPEN (UNIT=13, FILE=TRIM(cf_ascii))
   DO WHILE (io==0)
      READ(13,*,iostat=io)
      npoints = npoints + 1
   END DO
   PRINT*, ' *** Found ', npoints, ' points in '//TRIM(cf_ascii)//'...'

   ALLOCATE ( vpt_lon(1,npoints), vpt_lat(1,npoints), cb_name(npoints) )
   !!
   REWIND(13)
   DO jl = 1, npoints
      !READ(13,'(a,f,f)') cb_name(jl), vpt_lon(1,jl), vpt_lat(1,jl)
      READ(13,*) cb_name(jl), vpt_lon(1,jl), vpt_lat(1,jl)
      IF ( vpt_lon(1,jl) < 0. ) vpt_lon(1,jl) = vpt_lon(1,jl) + 360.0      
   END DO
   CLOSE(13)


   PRINT *, ''; PRINT *, ''
   DO jl = 1, npoints
      WRITE(*,'(" *** Point #",i4," -> Box ",a," : ",f9.4,f9.4)') &
         &     jl, TRIM(cb_name(jl)), vpt_lon(1,jl), vpt_lat(1,jl)
   END DO
   PRINT *, ''

   lat_min_trg = MINVAL(vpt_lat)
   lat_max_trg = MAXVAL(vpt_lat)
   lon_min_trg = MINVAL(vpt_lon)
   lon_max_trg = MAXVAL(vpt_lon)


   IF ( (lat_min_trg<lat_min).OR.(lat_max_trg>lat_max)) THEN
      PRINT *, 'ERROR: based on latitude, 1 of your points might be outside of model grid domain!'
      STOP
   END IF

   !IF ( (lon_min_trg<lon_min_1).OR.(lon_max_trg>lon_max_1)) THEN
   !   PRINT *, 'ERROR: based on longitude, 1 of your points might be outside of model grid domain!'
   !   STOP
   !END IF


   ALLOCATE (JIidx(1,npoints) , JJidx(1,npoints) )

   CALL FIND_NEAREST_POINT(vpt_lon, vpt_lat, xlont_tmp, xlatt,  JIidx, JJidx)

   DEALLOCATE ( xlont_tmp )


   
   PRINT *, ''
   PRINT *, '#  NEAREST POINTS ON GRID:'
   
   DO jl = 1, npoints

      IF ( (JJidx(1,jl)<0).OR.(JIidx(1,jl)<0) ) THEN
         PRINT *, 'ERROR #2 : 1 of your points might be outside of model grid domain!'
         PRINT *, ' *** Point ",a," ->  [i,j]:', TRIM(cb_name(jl)), JIidx(1,jl), JJidx(1,jl)
         STOP
      END IF

      IF ( (JJidx(1,jl)>3.E7).OR.(JIidx(1,jl)>3.E7) ) THEN
         PRINT *, 'ERROR #3: 1 of your points might be outside of model grid domain!'
         PRINT *, ' *** Point [i,j]:', JIidx(1,jl), JJidx(1,jl)
         STOP
      END IF

      WRITE(*,'(" *** Point #",i4,", ",a," -> i,j: ",i4.4,",",i4.4," (lon,lat:",f9.4,",",f9.4,")")') &
         &     jl, TRIM(cb_name(jl)), JIidx(1,jl), JJidx(1,jl), vpt_lon(1,jl), vpt_lat(1,jl)
      
   END DO
   PRINT *, ''

   IF ( npoints == 2 ) THEN
      PRINT *, ''
      PRINT *, '  ==> for XIOS-2.0 :'
      WRITE(*,'("     ibegin=",i4,"  jbegin=",i4," ni=",i4," nj=",i4)') JIidx(1,1)-1, JJidx(1,1)-1, JIidx(1,2)-JIidx(1,1), JJidx(1,2)-JJidx(1,1)
      PRINT *, ''
   END IF




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
      WRITE(6,*) ' -i <input_file.nc>   => file containing grid (longitude and latitude) of model'
      WRITE(6,*) ''
      WRITE(6,*) ' -p <ascii_file>      => Specify name of ASCII file containing "name", lon,lat (3 columns)'
      WRITE(6,*) '                         example of a valid line:'
      WRITE(6,*) '                             "Denmark_Strait" -26.49 66.32'
      WRITE(6,*) ''
      WRITE(6,*) '    Optional:'
      WRITE(6,*) ' -h                   => Show this message'
      WRITE(6,*) ''
      WRITE(6,*) ' -x  <name>           => Specify longitude name in input file (default: '//TRIM(cv_lon)//')'
      WRITE(6,*) ''
      WRITE(6,*) ' -y  <name>           => Specify latitude  name in input file  (default: '//TRIM(cv_lon)//')'
      WRITE(6,*) ''
      WRITE(6,*) ''
      !!
      STOP
      !!
   END SUBROUTINE usage




END PROGRAM IJ_FROM_LON_LAT




!!
