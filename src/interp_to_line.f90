PROGRAM INTERP_TO_LINE
  !!
  USE io_ezcdf
  USE mod_drown
  USE mod_akima_2d
  USE mod_bilin_2d
  !!
  !!========================================================================
  !! Purpose :  
  !!            
  !! ---------  
  !!
  !! Author :   Laurent Brodeau
  !! -------- 
  !!
  !!========================================================================
  !!
  IMPLICIT NONE
  !!
  !!
  !! ************************ Configurable part ****************************
  !!
  LOGICAL, PARAMETER :: &
       &   l_akima = .FALSE., &
       &   l_bilin = .TRUE.

  !!
  REAL(8), PARAMETER :: res = 0.1  ! resolution in degree
  !!
  INTEGER :: ndiv
  !!
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xtar, Ytar
  REAL(4), DIMENSION(:,:), ALLOCATABLE :: Ztar4
  !!
  !! Coupe stuff:
  REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: xcoupe
  REAL(4), DIMENSION(:,:),   ALLOCATABLE :: xcmask
  REAL(8), DIMENSION(:,:),   ALLOCATABLE :: vposition

  REAL(8), DIMENSION(:,:),   ALLOCATABLE :: vdepth
  REAL(8), DIMENSION(:),     ALLOCATABLE :: vtime



  !! Grid :
  CHARACTER(len=80) :: &
       &    cv_in  = '',      &
       &    cv_t   = 'time',  &
       &    cv_mt  = 'tmask', &
       &    cv_z   = 'depth', &
       &    cv_lon = 'lon',   &   ! input grid longitude name, T-points
       &    cv_lat = 'lat'       ! input grid latitude name,  T-points

  CHARACTER(len=256)  :: cr, csection
  !!
  !!
  !!******************** End of conf for user ********************************
  !!
  !!               ** don't change anything below **
  !!
  LOGICAL :: lregin, &
       &     lcontinue = .TRUE., &
       &     l_exist   = .FALSE.
  !!
  !!
  CHARACTER(len=400)  :: &
       &    cf_in = '',  cf_out, &
       &    cf_section = 'transportiz.dat', &
       &    cf_mm='mesh_mask.nc'
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
  !! For section:
  INTEGER :: imin, imax, jmin, jmax, isav, jsav, ji_min, ji_max, jj_min, jj_max, nib, njb

  REAL(4), DIMENSION(:,:,:,:), ALLOCATABLE :: xvar
  REAL(4), DIMENSION(:,:,:),   ALLOCATABLE :: xtmp4_b
  REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: xvar_b


  REAL(8), DIMENSION(:,:), ALLOCATABLE ::    &
       &    xlont, xlatt, &
       &   xlont_b, xlatt_b
  !!
  INTEGER(2), DIMENSION(:,:,:), ALLOCATABLE :: mask, mask_b


  !!
  INTEGER :: jt
  !!
  REAL(8) :: rA, rB, dlon, dlat, dang, lon_min, lon_max, lat_min, lat_max
  !!
  CHARACTER(LEN=2), DIMENSION(9), PARAMETER :: &
       &            clist_opt = (/ '-h','-v','-x','-y','-z','-t','-i','-s','-m' /)
  !!
  !!

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

     CASE('-s')
        CALL GET_MY_ARG('section ascii file', cf_section)

     CASE('-m')
        CALL GET_MY_ARG('mesh_mask file', cf_mm)

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
  !!


  !! testing longitude and latitude
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  CALL DIMS(cf_in, cv_lon, ni1, nj1, nk, nt)
  CALL DIMS(cf_in, cv_lat, ni2, nj2, nk, nt)

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

  ALLOCATE ( xlont(ni,nj), xlatt(ni,nj) )
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

  ALLOCATE ( xvar(ni,nj,nk,nt), mask(ni,nj,nk), vdepth(nk,1), vtime(nt) )

  IF ( lregin ) THEN
     PRINT *, 'Regular case not supported yet! Priority to ORCA grids...'
     STOP
  END IF





  !! Getting coordinates
  !! ~~~~~~~~~~~~~~~~~~~
  CALL GETVAR_1D(        cf_in, cv_t, vtime)
  CALL GETVAR_1D(        cf_in, cv_z, vdepth(:,1))
  !!
  CALL GETVAR_2D_R8(i0, j0, cf_in, cv_lon, 1, 1, 1, xlont) ; i0=0 ; j0=0
  CALL GETVAR_2D_R8(i0, j0, cf_in, cv_lat, 1, 1, 1, xlatt)

  WHERE ( xlont < 0. ) xlont = xlont + 360.

  CALL GETMASK_3D(cf_mm, cv_mt, mask)





  !! Filling xvar once for all, !BAD lolo, if few virtual memory on the machine...
  !! ~~~~~~~~~~~~~~~~~~~~~~~~~
  DO jt = 1, nt
     !!
     CALL GETVAR_3D(id_f1, id_v1, cf_in, cv_in, nt, jt, xvar(:,:,1:nk,jt))
     !!
  END DO
  !!








  !! LOOP ALONG SECTIONS:
  !! ~~~~~~~~~~~~~~~~~~~~

  INQUIRE(FILE=trim(cf_section), EXIST=l_exist )
  IF ( .NOT. l_exist ) THEN
     PRINT *, 'ERROR: please provide the file containing definition of sections'; STOP
  END IF

  OPEN(UNIT=13, FILE=trim(cf_section), FORM='FORMATTED', STATUS='old')  ! 11 and 12 busy with bilin module...
  
  DO WHILE ( lcontinue ) 
     
     READ(13,'(a)') csection
     IF (TRIM(csection) == 'EOF' ) THEN
        PRINT *, 'Closing ', trim(cf_section) ; PRINT *, ''
        CLOSE(13)
        lcontinue = .FALSE.
        EXIT
     END IF
     
     IF ( .NOT. lcontinue ) EXIT

     PRINT *, '';  PRINT *, ''
     PRINT *, 'Doing section ', trim(csection)
     
     READ(13,*) imin, imax, jmin, jmax
     PRINT*, '   =>', imin, imax, jmin, jmax
     PRINT*, '   => Point 1:', xlont(imin,jmin), xlatt(imin,jmin)
     PRINT*, '   => Point 2:', xlont(imax,jmax), xlatt(imax,jmax)
     PRINT *, ''

     WRITE(cf_out, '("section_",a,"_",a,".nc")') trim(cv_in), trim(csection)


     !! Creating the line joining the 2 points defining the section
     !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     CALL PREPARE_LINE()


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


     ! Interpolating the mask on target section
     DO jk = 1, nk
        !!
        IF ( l_akima ) CALL AKIMA_2D(-1, xlont_b, xlatt_b, REAL(mask_b(:,:,jk),4), &
             &                        Xtar, Ytar, Ztar4)
        !!
        IF ( l_bilin) CALL BILIN_2D(-1, xlont_b, xlatt_b, REAL(mask_b(:,:,jk),4), &
             &                        Xtar, Ytar, Ztar4, trim(csection))
        !!
        xcmask(:,jk) = Ztar4(1,:)
        !!
     END DO


     ifo=0 ; ivo=0

     ! Interpolating the nt snapshots of field on target section
     DO jt = 1, nt
        !!
        DO jk = 1, nk
           !!
           IF ( l_akima ) CALL AKIMA_2D(-1, xlont_b, xlatt_b, REAL(xvar_b(:,:,jk,jt),4), &
                &                           Xtar,    Ytar,    Ztar4)
           !!
           IF ( l_bilin ) CALL BILIN_2D(-1, xlont_b, xlatt_b, REAL(xvar_b(:,:,jk,jt),4), &
                &                           Xtar,    Ytar,    Ztar4, trim(csection))
           !!
           xcoupe(:,jk,jt) = Ztar4(1,:)
           !!
        END DO
        !!
        WHERE( xcmask < 0.25 ) xcoupe(:,:,jt) = -9999.

        CALL P2D_T(ifo, ivo, nt, jt, vposition, vdepth, vtime, xcoupe(:,:,jt), cf_out, &
             &     'position', 'profo', cv_t, cv_in, 'boo', 'bla bla bla', -9999.)

     END DO

     l_first_call_bl = .TRUE. ! next is another conf (needed by mod_bilin_2d and mod_akima_2d)
     l_first_call_akima = .TRUE. ! next is another conf (needed by mod_bilin_2d and mod_akima_2d)

     PRINT *, 'File created => ', trim(cf_out)



     DEALLOCATE ( Xtar, Ytar, Ztar4 )
     DEALLOCATE ( xcoupe, xcmask, vposition )
     DEALLOCATE ( xlont_b, xlatt_b, xvar_b, mask_b, xtmp4_b )
     !lolo



  END DO  ! loop along section


CONTAINS



  SUBROUTINE PREPARE_LINE()

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

    lon_max = MAX(xlont(ji_max,jj_min),xlont(ji_max,jj_max))
    lon_min = MIN(xlont(ji_min,jj_min),xlont(ji_min,jj_max))
    lat_max = MAX(xlatt(ji_max,jj_max),xlatt(ji_min,jj_max))
    lat_min = MIN(xlatt(ji_max,jj_min),xlatt(ji_min,jj_min))

    PRINT *, ''
    dlon = lon_max - lon_min ; PRINT *, 'long. range =', dlon
    dlat = lat_max - lat_min ; PRINT *, 'latg. range =', dlat
    dang = SQRT(dlon*dlon + dlat*dlat) ; PRINT *, 'Ang. range =', dang

    ndiv = INT(dang/res) + 1
    IF ( MOD(ndiv,2) == 0 ) ndiv = ndiv - 1 ! we want odd integer...
    PRINT *, 'Number of points to create on segment:', ndiv ; PRINT *, ''

    ALLOCATE ( Xtar(1,ndiv), Ytar(1,ndiv), Ztar4(1,ndiv) )
    ALLOCATE ( xcoupe(ndiv,nk,nt), xcmask(ndiv,nk), vposition(ndiv,1) )

    IF ( ABS(dlon) < 1.E-12 ) THEN
       PRINT *, 'ERROR: Section seems to be vertical!'; STOP
    END IF

    rA = (xlatt(imax,jmax) - xlatt(imin,jmin))/ dlon
    rB = xlatt(imin,jmin) - rA*xlont(imin,jmin)

    PRINT *, 'rA, rB = ', rA, rB
    PRINT *, 'Lat1 =', rA*xlont(imin,jmin) + rB
    PRINT *, 'Lat2 =', rA*xlont(imax,jmax) + rB

    dlon = dlon / (ndiv-1) ;  ; PRINT *, 'dlon =', dlon

    DO jd = 1, ndiv
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

  END SUBROUTINE PREPARE_LINE










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







SUBROUTINE angle_dist(lx, ly, xlon_t, xlat_t, xlon_u, xlat_u, cos_t, sin_t)
  !!
  !!
  !!========================================================================
  !!
  !!                     *******  ANGLE_DIST  *******
  !!
  !! Given the coordinates at T-points and U-points, this routine computes 
  !! sinus and cosinus of the angle of the local grid distorsion at T-points
  !!
  !!
  !! INPUT :     - lx     = x dimension                     [integer]
  !! -------     - ly     = y dimension                     [integer]
  !!             - xlon_t = longitude array at T-points     [array(lx,ly) of real]
  !!             - xlat_t = latitude array at T-points      [array(lx,ly) of real]
  !!             - xlon_u = longitude array at U-points     [array(lx,ly) of real]
  !!             - xlat_t = latitude array at U-points      [array(lx,ly) of real]
  !!
  !! OUTPUT :
  !! --------    - cos_t  = cosinus of the distortion angle [array(lx,ly) of real]
  !!             - sin_t  = sininus of the distortion angle [array(lx,ly) of real]
  !!
  !!
  !! Author : Laurent Brodeau (brodeau@hmg.inpg.fr)
  !!          directly inspired from 'geo2ocean.F90' (OPA/Madec)
  !!
  !!========================================================================
  !!
  !!
  IMPLICIT NONE
  !!
  !!
  !! INPUT :
  !! -------
  INTEGER, INTENT(in) :: lx, ly       ! dimension of arrays
  !!
  REAL(8), DIMENSION(lx,ly), INTENT(in) ::    &
       &       xlon_t, xlon_u,   &  ! latitudes  at point T and U
       &       xlat_t, xlat_u       ! longitudes at point T and U 
  !!
  !!
  !! OUTPUT :
  !! --------
  REAL(8), DIMENSION(lx,ly), INTENT(out) :: cos_t, sin_t
  !!
  !!
  !! LOCAL :
  !! -------
  INTEGER :: ji
  !!
  REAL(8), DIMENSION(ly) :: &
       &       zlon, zlat, zxnpt, znnpt, zlan, &
       &       zphh, zxuut, zyuut, zmnput, zynpt
  !!
  REAL(8), PARAMETER :: &
       &       Pi  = 3.141592653,     &
       &       rad = Pi/180.0
  !!
  !!
  !!
  DO ji = 1, lx
     !!
     !! North pole direction & modulous (at T-point) :
     !! ----------------------------------------------
     zlon  = xlon_t(ji,:)
     zlat  = xlat_t(ji,:)
     zxnpt = 0. - 2*COS( rad*zlon )*TAN(Pi/4 - rad*zlat/2)
     zynpt = 0. - 2*SIN( rad*zlon )*TAN(Pi/4 - rad*zlat/2)
     znnpt = zxnpt*zxnpt + zynpt*zynpt
     !!
     !!
        !! "i" direction & modulous (at T-point) :
     !! ---------------------------------------
     !!   ( since we deal with T points we look at U point 
     !!     on a V point we would look at the F point )
     !!
     zlon = xlon_u(ji,:)
     zlat = xlat_u(ji,:)
     !!
     IF ( ji == 1 ) THEN    
        zlan = xlon_u(lx-3,:)          !! periodicity of ORCA grid
        zphh = xlat_u(lx-3,:)          !! with overlap of 2 points
     ELSE
        zlan = xlon_u(ji-1,:)
        zphh = xlat_u(ji-1,:)
     END IF
     !!
     !!
     zxuut  = 2*COS(rad*zlon) &
          &    *TAN(Pi/4 - rad*zlat/2) - 2*COS(rad*zlan)*TAN(Pi/4 - rad*zphh/2)
     zyuut  = 2*SIN(rad*zlon) &
          &    *TAN(Pi/4 - rad*zlat/2) - 2*SIN(rad*zlan)*TAN(Pi/4 - rad*zphh/2)
     zmnput = SQRT(znnpt*(zxuut*zxuut + zyuut*zyuut))
     !!
     WHERE ( zmnput < 1.e-14 ) zmnput = 1.e-14
     !!
     !!
     !! Cosinus and sinus using scalar and vectorial products :
     !! -------------------------------------------------------
     !!   (caution, rotation of 90 degres)
     sin_t(ji,:) =  ( zxnpt*zxuut + zynpt*zyuut ) / zmnput
     cos_t(ji,:) = -( zxnpt*zyuut - zynpt*zxuut ) / zmnput
     !!
  END DO
  !!
  !!
END SUBROUTINE angle_dist
!!
!!
!!
SUBROUTINE usage()
  !!
  OPEN(UNIT=6, FORM='FORMATTED', RECL=512)
  !!
  WRITE(6,*) ''
  WRITE(6,*) '   List of command line options:'
  WRITE(6,*) '   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  WRITE(6,*) ''
  WRITE(6,*) ' -i <input_file.nc>   => INPUTE FILE'
  WRITE(6,*) ''
  WRITE(6,*) ' -v  <name>           => Specify variable name in input file'
  WRITE(6,*) ''
  WRITE(6,*) '    Optional:'
  WRITE(6,*)  ''
  WRITE(6,*) ' -x  <name>           => Specify longitude name in input file (default: lon)'
  WRITE(6,*) ''
  WRITE(6,*) ' -y  <name>           => Specify latitude  name in input file (default: lat)'
  WRITE(6,*) ''
  WRITE(6,*) ' -z  <name>           => Specify depth name in input file (default: depth)'
  WRITE(6,*) ''
  WRITE(6,*) ' -t  <name>           => Specify time name in input file (default: time)'
  WRITE(6,*) ''
  WRITE(6,*) ' -s  <section_file>   => Specify name of ASCII file containing sections (default: transportiz.dat)'
  WRITE(6,*) ''
  WRITE(6,*) ' -m  <mesh_mask_file> => Specify mesh_mask file to be used (default: mesh_mask.nc)'
  WRITE(6,*) ''
  WRITE(6,*) ' -h                   => Show this message'
  WRITE(6,*) ''
  !!
  CLOSE(6)
  STOP
  !!
END SUBROUTINE usage
!!
