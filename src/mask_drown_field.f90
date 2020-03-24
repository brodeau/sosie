PROGRAM mask_drown_field

   !! Brodeau, 2011
   !!
   !! MASK or DROWN a given field with a special (missing) value
   !!         Extrapolate sea values over continents thanks to the DROWN algorithm

   use io_ezcdf
   use mod_drown

   IMPLICIT NONE


   REAL(4), PARAMETER :: rmv = -9999.

   !! Grid :
   CHARACTER(len=80) :: &
      &    cv_in  = '', &
      &    cv_t   = 'time', &
      &    cv_mm  = 'lsm', &
      &    cv_ewp = '0', &
      &    cv_how_far = '200', &
      &    cv_z   = 'depth', &
      &    cv_lon = 'lon',  &
      &    cv_lat = 'lat',  &
      &    cdum   = '',    &
      &    clbound = '',   &
      &    cubound = '',   &
      &    cu     = ''

   CHARACTER(len=400) :: &
      &      cf_in,         &
      &      cf_out = 'fout.nc', &
      &      cf_mm  = 'mask.nc'

   CHARACTER(len=400) :: clnm, cr

   REAL :: rlbound, rubound

   INTEGER :: &
      &     iewper = 0, &
      &     i_how_far = 200, &
      &     jarg, &
      &     ni, nj, nk, nt, jk, &
      &     ni_m, nj_m, nt_m, nk_m, &
      &     ifx, ivx, ify, ivy, ifi, ivi, ifo, ivo

   REAL(8), DIMENSION(:,:),   ALLOCATABLE :: xlon, xlat
   REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: data
   REAL(8), DIMENSION(:),     ALLOCATABLE :: vtime, vdpth, vlon, vlat
   INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: mask

   INTEGER :: jt, iargc
   REAL(4) :: zsf, zao, rmissv

   LOGICAL :: lreg, l3d, l_missv, &
      &   l_mask_f = .FALSE., &
      &   l_drwn_f = .FALSE.
   CHARACTER(LEN=14) :: cmissv
   CHARACTER(LEN=2), DIMENSION(16), PARAMETER :: &
      &  clist_opt = (/ '-h','-v','-x','-y','-z','-t','-i','-m','-q','-p','-g','-o','-M','-D','-l','-u' /)



   l3d = .FALSE.

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
         CALL GET_MY_ARG('input file', cf_in, 1)

      CASE('-v')
         CALL GET_MY_ARG('input variable', cv_in, 1)

      CASE('-x')
         CALL GET_MY_ARG('longitude', cv_lon, 1)

      CASE('-y')
         CALL GET_MY_ARG('latitude', cv_lat, 1)

      CASE('-z')
         CALL GET_MY_ARG('depth', cv_z, 1)

      CASE('-t')
         CALL GET_MY_ARG('time', cv_t, 1)

      CASE('-m')
         CALL GET_MY_ARG('mask file', cf_mm, 1)

      CASE('-q')
         CALL GET_MY_ARG('mask', cv_mm, 1)

      CASE('-p')
         CALL GET_MY_ARG('east west periodicity', cv_ewp, 1)

      CASE('-g')
         CALL GET_MY_ARG('how far into land', cv_how_far, 1)

      CASE('-o')
         CALL GET_MY_ARG('output file', cf_out, 1)

      CASE('-l')
         CALL GET_MY_ARG('lower bound', clbound, 1)

      CASE('-u')
         CALL GET_MY_ARG('upper bound', cubound, 1)

      CASE('-M')
         CALL GET_MY_ARG('', cdum, 0)
         l_mask_f = .TRUE.

      CASE('-D')
         CALL GET_MY_ARG('', cdum, 0)
         l_drwn_f = .TRUE.


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

   IF ( (.NOT. l_mask_f) .AND. (.NOT. l_drwn_f) ) THEN
      PRINT *, ''
      PRINT *, 'Should we mask ( => -M) or drown ( => -D) the field?'; PRINT *, ''
      STOP
   END IF

   IF ( l_mask_f .AND. l_drwn_f ) THEN
      PRINT *, ''
      PRINT *, 'Cannot DROWN and MASK at the same time! => -M OR -D'; PRINT *, ''
      STOP
   END IF




   IF ( trim(clbound) /= '' ) THEN
      READ(clbound, *) rlbound
      PRINT *, ''; PRINT *, 'Lower bound to apply to the field:', rlbound ; PRINT *, ''
   END IF
   IF ( trim(cubound) /= '' ) THEN
      READ(cubound, *) rubound
      PRINT *, ''; PRINT *, 'Upper bound to apply to the field:', rubound ; PRINT *, ''
   END IF



   READ(cv_how_far, '(i3)') i_how_far
   IF ( l_drwn_f ) WRITE (*, '("How far (points) from sea to extrapolate into land ",i3)')  i_how_far


   READ(cv_ewp, '(i2)') iewper
   IF ( l_drwn_f ) WRITE (*, '("East-West periodicity is set to ",i2)')  iewper


   PRINT *, ''
   PRINT *, 'File in      :'  ;   PRINT *, cf_in
   PRINT *, 'Variable in  :'  ;   PRINT *, cv_in
   PRINT *, 'File mask    :'  ;   PRINT *, cf_mm
   PRINT *, ''


   !! Testing dimensions of coordinates to see if regular or not
   !CALL DIMS(cf_in, cv_lon, ni, nj, nt, l0)
   CALL DIMS(cf_in, cv_lon, ni, nj, nk, nt)
   IF ( ni > 0 ) THEN
      IF ( (nj == -1).and.(nt == -1).and.(nk == -1) )  lreg = .TRUE.
      IF ( (nj > 0)  .and.(nt == -1).and.(nk == -1) )  lreg = .FALSE.
      IF ( (nt > 0).or.(nk > 0) ) THEN
         PRINT *, 'Problem found 1! Dimensions of coordinates are suspicious!' ; STOP
      END IF
   ELSE
      PRINT *, 'Problem found 2! Dimensions of coordinates are suspicious!' ; STOP
   END IF


   !CALL DIMS(cf_in, cv_in, ni, nj, nt, l0)
   CALL DIMS(cf_in, cv_in, ni, nj, nk, nt)

   IF ( nk > 0 ) THEN
      PRINT *, 'The variable is 3D+T!'
      l3d = .TRUE.
      !nk = nt
      !nt = l0
      PRINT *, 'Dimension = ', ni, nj, nk, nt
      !!
   ELSE
      !! 2D or 2D + T
      PRINT *, 'Dimension = ', ni,nj,nt
      IF ( nt == -1 )   nt = 1
      nk = 1
      !!
   END IF

   PRINT *, ''

   CALL GET_SF_AO(cf_in, cv_in, zsf, zao)


   ALLOCATE ( vtime(nt) )

   IF ( .NOT. l3d ) nk = 1
   ALLOCATE ( data(ni,nj,nk), mask(ni,nj,nk), vdpth(nk) )



   IF ( nt /= 1 ) THEN  ! there is a time dimension
      CALL GETVAR_1D(cf_in, cv_t, vtime)
   ELSE
      vtime = 1.
   END IF

   IF ( l3d ) CALL GETVAR_1D(cf_in, cv_z, vdpth)

   !! Getting long name and unit :
   CALL GET_VAR_INFO(cf_in, cv_in, cu, clnm)
   PRINT *, 'Unit      = ', trim(cu)
   PRINT *, 'Long name = ', trim(clnm); PRINT *, ''



   IF ( lreg ) THEN
      ALLOCATE ( vlon(ni), vlat(nj), xlon(ni,1), xlat(nj,1) )
      CALL GETVAR_1D(cf_in, cv_lon, vlon)
      CALL GETVAR_1D(cf_in, cv_lat, vlat)
      xlon(:,1) = vlon(:) ; xlat(:,1) = vlat(:)
   ELSE
      ALLOCATE ( xlon(ni,nj), xlat(ni,nj) )
      CALL GETVAR_2D(ifx, ivx, cf_in, cv_lon, 0, 0, 0, xlon)   !lolo => not cool should be able to DOUBLE!
      CALL GETVAR_2D(ify, ivy, cf_in, cv_lat, 0, 0, 0, xlat)   !lolo => not cool should be able to DOUBLE!
   END IF


   PRINT *, ''

   IF ( TRIM(cf_mm) == '0' ) THEN
      CALL CHECK_4_MISS(cf_in, cv_in, l_missv, rmissv, cmissv)
      PRINT *, 'Will extract mask from treated file from missing value "'//TRIM(cmissv)//'"....'
      IF ( .NOT. l_missv ) THEN
         PRINT *, 'PROBLEM: variable ',TRIM(cv_in),' of file ',TRIM(cf_in),' doesnt have a missing-value attribute!'
         STOP
      END IF
      PRINT *, ' => missing value to use is', rmissv
      mask(:,:,:) = 1
   ELSE
      PRINT *, 'Will extract mask "',trim(cv_mm),'" from file ', trim(cf_mm)
      CALL DIMS(cf_mm, cv_mm, ni_m, nj_m, nk_m, nt_m)
      PRINT *, 'Mask dim =', ni_m, nj_m, nk_m, nt_m

      IF ( (ni_m /= ni).OR.(nj_m /= nj) ) THEN
         PRINT *, 'Error: the mask has not the same ni x nj than input variable!'; STOP
      END IF

      IF ( l3d .AND. nk_m == -1 ) THEN
         PRINT *, 'Error: your variable is 3D and your mask is 2D!!!'; STOP
      END IF

      IF ( l3d .AND. nk_m /= nk ) THEN
         PRINT *, 'Error: your mask does not have the same number of levels than the input variable!'
         STOP
      END IF

      IF ( (l3d) .AND. (nk > 1) ) THEN
         CALL GETMASK_3D(cf_mm, cv_mm, mask)
      ELSE
         CALL GETMASK_2D(cf_mm, cv_mm, mask(:,:,1), jlev=1)
      END IF

   END IF

   PRINT *, ''; PRINT *, ''





   DO jt = 1, nt

      IF ( l3d ) THEN
         CALL GETVAR_3D(ifi, ivi, cf_in, cv_in, nt,     jt, data(:,:,:))
      ELSE
         CALL GETVAR_2D(ifi, ivi, cf_in, cv_in, nt,  0, jt, data(:,:,1))
      END IF

      IF ( trim(cf_mm) == '0' ) THEN
         mask(:,:,:) = 1  ! cuz mask can change for each time record...
         WHERE ( data(:,:,:) == rmissv ) mask(:,:,:) = 0
      END IF

      data = data*zsf + zao

      IF ( trim(clbound) /= '' ) THEN
         WHERE ( (mask /= 0.).AND.(data < rlbound) ) data = rlbound
      END IF
      IF ( trim(cubound) /= '' ) THEN
         WHERE ( (mask /= 0.).AND.(data > rubound) ) data = rubound
      END IF




      IF ( l_mask_f ) THEN
         PRINT *, ' *** masking field at time =', jt
         WHERE ( mask == 0 ) data = rmv
      END IF

      IF ( l_drwn_f ) THEN
         PRINT *, ' *** drowning field at time =', jt
         DO jk = 1, nk
            CALL DROWN(iewper, DATA(:,:,jk), mask(:,:,jk), nb_inc=i_how_far)
         END DO
      END IF


      IF ( l3d ) THEN
         CALL P3D_T(ifo, ivo, nt, jt, xlon, xlat, vdpth, vtime, data, &
            &     cf_out, cv_lon, cv_lat, cv_z, cv_t, cv_in, rmv)
      ELSE
         CALL P2D_T(ifo, ivo, nt, jt, xlon, xlat,        vtime, data(:,:,1), &
            &     cf_out, cv_lon, cv_lat,       cv_t, cv_in, rmv)
      END IF


   END DO




CONTAINS

   SUBROUTINE GET_MY_ARG(cname, cvalue, nbarg)

      CHARACTER(len=*), INTENT(in)    :: cname
      CHARACTER(len=*), INTENT(inout) :: cvalue
      INTEGER,          INTENT(in)    :: nbarg !: number of expected character string after a -* option
      !                                        !: so far: 0 or 1

      IF ( (jarg + 1 > iargc()).AND.(nbarg > 0) ) THEN
         PRINT *, 'ERROR: Missing ',trim(cname),' name!' ; call usage()
      ELSE

         CALL getarg(jarg+1,cr) ! reading next argument

         !! We do not want a -* after a -* that expects an argument:
         IF ( (nbarg > 0) .AND. ANY(clist_opt == trim(cr)) ) THEN
            PRINT *, 'ERROR: Missing',trim(cname),' name!'; call usage()

         ELSEIF (nbarg > 0) THEN
            cvalue = trim(cr)
            jarg = jarg + 1 ! as we just read the argument we must advance 1 argument

            !! We do not want a following argument so we want a -* after or nothing
         ELSE IF ( ( (jarg + 1 <= iargc()).AND.(nbarg == 0) ).AND. (.NOT. ANY(clist_opt == trim(cr)) ) ) THEN
            PRINT *, 'ERROR: ',trim(cname),' does not expect an argument!'; call usage()

         END IF
      END IF
   END SUBROUTINE GET_MY_ARG


END PROGRAM mask_drown_field





SUBROUTINE usage()
   !!
   !! OPEN(UNIT=6, FORM='FORMATTED', RECL=512)  !! MB: Pb with gfortran
   !!
   WRITE(6,*) ''
   WRITE(6,*) '   List of command line options:'
   WRITE(6,*) '   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
   WRITE(6,*) ''
   WRITE(6,*) ' -M                   => MASK the field'
   WRITE(6,*) ''
   WRITE(6,*) ' -D                   => DROWN the field'
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
   WRITE(6,*) ' -m  <mask_file>      => Specify mask file to be used (default: mask.nc)'
   WRITE(6,*) '                         or "0" to use missing-value attribute value of input variable'
   WRITE(6,*) ''
   WRITE(6,*) ' -q  <name>           => Specify mask name in mask file (default: lsm)'
   WRITE(6,*) ''
   WRITE(6,*) ' -p  <integer>        => DROWN: east-west periodicity in points (default: 0)'
   WRITE(6,*) '                         * no periodicity       => -1'
   WRITE(6,*) '                         * no overlaping point  =>  0'
   WRITE(6,*) '                         * N  overlaping points =>  N'
   WRITE(6,*) ''
   WRITE(6,*) ' -l  <real>           => to set a lower bound for the field (ex: "-10." not "-10")'
   WRITE(6,*) ''
   WRITE(6,*) ' -u  <real>           => to set an upper bound for the field (ex: "100." not "100")'
   WRITE(6,*) ''
   WRITE(6,*) ' -g  <integer>        => DROWN: how far into land the variable is drowned' 
   WRITE(6,*) ''
   WRITE(6,*) ' -o  <output_file.nc> => Output file (default: fout.nc)'
   WRITE(6,*) ''
   WRITE(6,*) ' -h                   => Show this message'
   WRITE(6,*) ''
   !!
   CLOSE(6)
   STOP
   !!
END SUBROUTINE usage
!!
