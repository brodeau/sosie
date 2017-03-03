MODULE io_ezcdf

   USE netcdf

   !! Netcdf input/output
   !!
   !! Author: Laurent Brodeau, 2010
   !!


   IMPLICIT NONE

   PRIVATE


   !! List of public routines
   !! =======================
   PUBLIC :: dims,      &
      &    get_sf_ao,        &
      &    getvar_1d,        &
      &    getvar_2d,        &
      &    getvar_2d_r8,     & ! Mostly for 2D coordinates
      &    getvar_3d,        &
      &    getmask_2d,       &
      &    getmask_3d,       &
      &    pt_series,        &
      &    p2d_t,            &
      &    p3d_t,            &
      &    check_4_miss,     &
      &    get_var_info,     &
      &    prtmask,          &
      &    p2d_mapping_ab,    &
      &    rd_mapping_ab,     &
      &    phovmoller,       &
      &    who_is_mv
   !!===========================

   CHARACTER(len=80) :: cv_misc

   CHARACTER(len=2) :: cdt  ! '1d' or '2d'

   REAL(8), DIMENSION(3,2) :: vextrema

   INTEGER :: &
      &    nd,    &
      &    id_f, id_v,           &   !: ID for a variable
      &    id_x, id_y, id_z, id_t,     &
      &    id_lo, id_la, &
      &    id_dpt, id_tim

   CHARACTER(len=400)    :: crtn, cu

   CHARACTER(len=8), PARAMETER :: cdum = 'dummy'

   CHARACTER(LEN=400), PARAMETER   ::     &
      &    cabout = 'Created by SOSIE interpolation environement (http://sosie.sourceforge.net/) / Laurent Brodeau, 2015'

   INTEGER :: ji, jj, jk

   INTEGER, PARAMETER :: nmval = 3
   CHARACTER(len=80), DIMENSION(nmval), PARAMETER :: &
      &     c_nm_missing_val = (/ 'FillValue ', '_FillValue', '_Fillvalue' /)


CONTAINS




   SUBROUTINE DIMS(cf_in, cv_in, lx, ly, lz, lt)

      !!-----------------------------------------------------------------------
      !!
      !! Updated August 2012, L. Brodeau
      !!
      !! This routine opens a netcdf file 'cf_in' to check the dimension
      !! of the variable 'cv_in'. It then gives the length of each of the dimension,
      !! if the length returns '-1' that means that the dimension does not exist
      !!
      !! example : if the variable has only 1 dimension, of length 132,
      !!           DIMS will return lx=132, ly=-1, lz=-1, lt=-1
      !!
      !!
      !! INPUT :
      !! -------
      !!          * cf_in       : name of the input file          (character)
      !!          * cv_in       : name of the variable            (character)
      !!
      !! OUTPUT :
      !! --------
      !!          * lx      : first dimension                     (integer)
      !!          * ly      : second dimension  (-1 if none)      (integer)
      !!          * lz      : third dimension   (-1 if none)      (integer)
      !!          * lt      : number of records (-1 if none)      (integer)
      !!
      !!------------------------------------------------------------------------

      CHARACTER(len=*),   INTENT(in)  :: cf_in, cv_in
      INTEGER,            INTENT(out) :: lx, ly, lz, lt

      INTEGER, DIMENSION(:), ALLOCATABLE :: id_dim, nlen
      INTEGER :: jdim, id_unlim_dim


      crtn = 'DIMS'

      lx = -1 ; ly = -1 ; lz = -1 ; lt = -1

      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)

      ! Get ID of unlimited dimension
      CALL sherr(  NF90_INQUIRE(id_f, unlimitedDimId = id_unlim_dim), crtn,cf_in,cv_in)



      !! Getting variable ID:
      CALL sherr( NF90_INQ_VARID(id_f,cv_in,id_v),              crtn,cf_in,cv_in)

      !! nd => number of dimensions for the variable:
      CALL sherr( NF90_INQUIRE_VARIABLE(id_f, id_v, ndims=nd),  crtn,cf_in,cv_in)

      ALLOCATE ( id_dim(nd) , nlen(nd) )

      !! Vector containing the IDs of each dimension (id_dim):
      CALL sherr( NF90_INQUIRE_VARIABLE(id_f, id_v, dimids=id_dim),  crtn,cf_in,cv_in)

      DO jdim = 1, nd
         CALL sherr( NF90_INQUIRE_DIMENSION(id_f, id_dim(jdim), len=nlen(jdim)),   crtn,cf_in,cv_in)
      END DO

      IF ( (nd > 2).AND.(id_unlim_dim < 1) ) THEN
         WRITE(6,*) 'WARNING: DIMS of io_ezcdf.f90'
         WRITE(6,*) '   => variable '//TRIM(cv_in)//' in file:'
         WRITE(6,*) '      '//TRIM(cf_in)
         WRITE(6,*) '      does not have an UNLIMITED dimension !!!'
         WRITE(6,*) '   => if it is supposed to depend on a time-record, this time-record should be'
         WRITE(6,*) '      an unlimited DIMENSION into the netcdf file!!!'
         WRITE(6,*) '   => otherwize time dimension might be confused with a space dimension!!!'
         WRITE(6,*) '   => maybe SOSIE can overcome this, but at your own risks... ;)'
         WRITE(6,*) ''
      END IF

      SELECT CASE(nd)

      CASE(1)
         !! Purely 1D
         lx = nlen(1)

      CASE(2)
         IF ( id_dim(2) == id_unlim_dim ) THEN
            !! 1D with time records
            lx = nlen(1) ; lt = nlen(2)
         ELSE
            !! 2D with no time records
            lx = nlen(1) ; ly = nlen(2)
         END IF

      CASE(3)
         IF ( id_dim(3) == id_unlim_dim ) THEN
            !! 2D with time records
            lx = nlen(1) ; ly = nlen(2) ; lt = nlen(3)
         ELSE
            !! 3D with no time records
            lx = nlen(1) ; ly = nlen(2) ; lz = nlen(3)
         END IF

      CASE(4)
         IF ( id_unlim_dim < 1 ) THEN
            WRITE(6,*) 'ERROR: file ',trim(cf_in),' doesnt have an unlimited dimension (time record)!'
         END IF

         lx = nlen(1) ; ly = nlen(2)
         IF ( id_dim(3) == id_unlim_dim ) THEN
            lz = nlen(4) ; lt = nlen(3)   ! time record (unlimited dim) comes as 3rd dim and lz as 4th
         ELSE
            lz = nlen(3) ; lt = nlen(4)   ! time record (unlimited dim) comes as last dim and lz as 3rd
         END IF

      CASE DEFAULT
         CALL print_err(crtn, 'the dimension is not realistic')

      END SELECT

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE DIMS










   SUBROUTINE GETVAR_1D(cf_in, cv_in, X)

      !!-----------------------------------------------------------------------
      !! This routine extract a variable 1D from a netcdf file
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file             (character l=100)
      !!          * cv_in      : name of the variable               (character l=20)
      !!
      !! OUTPUT :
      !! --------
      !!          * X         : 1D array contening the variable   (double)
      !!
      !!------------------------------------------------------------------------

      CHARACTER(len=*),       INTENT(in)  :: cf_in, cv_in
      REAL(8), DIMENSION (:), INTENT(out) ::  X

      crtn = 'GETVAR_1D'

      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)

      CALL sherr( NF90_INQ_VARID(id_f, trim(cv_in), id_v),  crtn,cf_in,cv_in)
      CALL sherr( NF90_GET_VAR(id_f, id_v, X),              crtn,cf_in,cv_in)

      CALL sherr( NF90_CLOSE(id_f),                         crtn,cf_in,cv_in)

   END SUBROUTINE GETVAR_1D




   SUBROUTINE GETVAR_2D(idx_f, idx_v, cf_in, cv_in, lt, kz, kt, X, jt1, jt2, lz)

      !!-----------------------------------------------------------------------------
      !! This routine extract a 2D field from a netcdf file
      !! at a given time
      !!
      !! INPUT :
      !! -------
      !!          * idx_f    : ID of current file                  (integer)
      !!          * idx_v    : ID of current variable              (integer)
      !!          * cf_in      : name of the input file              (character)
      !!          * cv_in      : name of the variable                (character)
      !!          * lt        : time dimension of the variable     (integer)
      !!          * kz        : level to extract                    (integer)
      !!                      0 => input file does not have levels ( 1 would work anyway...)
      !!
      !!          * kt        : time snapshot to extract            (integer)
      !!                      0 => input file does not have a time snapshot
      !!                           (= old GETVAR_2D_NOTIME)
      !!
      !!
      !! OUTPUT :
      !! --------
      !!          * X         : 2D array contening the variable     (real)
      !!
      !! OPTIONAL INPUT :
      !! ----------------
      !!          * jt1, jt2  : first and last time snapshot to extract
      !!          *       lz  : number of levels to know when they are all read
      !!                        so we can close the file
      !!
      !!------------------------------------------------------------------------

      INTEGER,                   INTENT(inout) :: idx_f, idx_v
      CHARACTER(len=*),          INTENT(in)    :: cf_in, cv_in
      INTEGER,                   INTENT(in)    :: lt, kz, kt
      REAL(4),  DIMENSION (:,:), INTENT(out)   :: X
      INTEGER,       OPTIONAL,   INTENT(in)    :: jt1, jt2, lz

      INTEGER :: &
         & lx, &    ! x dimension of the variable        (integer)
         & ly       ! y dimension of the variable        (integer)

      INTEGER :: n1, n2, n3, n4, jlev, its, ite, kz_stop = 0

      LOGICAL :: l_okay

      crtn = 'GETVAR_2D'

      lx = size(X,1)
      ly = size(X,2)

      jlev = kz ! so we can modify jlev without affecting kz...

      CALL DIMS(cf_in, cv_in, n1, n2, n3, n4)

      IF ( (lx /= n1).OR.(ly /= n2) ) CALL print_err(crtn, ' PROBLEM #1 => '//TRIM(cv_in)//' in '//TRIM(cf_in))
      !PRINT *, 'lt =', lt
      !PRINT *, ' n1, n2, n3, n4 =>', n1, n2, n3, n4
      !lolo IF ( (lt > 0).AND.(n4 < lt)  ) CALL print_err(crtn, ' PROBLEM #2  => '//TRIM(cv_in)//' in '//TRIM(cf_in))

      IF ( present(jt1).AND.present(jt2) ) THEN
         its = jt1 ; ite = jt2
      ELSE
         its = 1   ; ite = lt
      END IF

      IF ( present(lz) ) kz_stop = lz

      IF ( (kt == its).OR.(kt == 0) ) THEN   ! Opening file and defining variable :
         CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, idx_f),  crtn,cf_in,cv_in)
         CALL sherr( NF90_INQ_VARID(idx_f, cv_in, idx_v),  crtn,cf_in,cv_in)
      END IF


      l_okay = .FALSE.
      DO WHILE ( .NOT. l_okay )

         IF ( jlev == 0 ) THEN    ! No levels

            IF ( kt == 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1/), count=(/lx,ly/)), &
                  &      crtn,cf_in,cv_in)
            ELSEIF ( kt > 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,kt/), count=(/lx,ly,1/)), &
                  &      crtn,cf_in,cv_in)
            END IF
            l_okay = .TRUE. ! we can exit the WHILE loop...

         ELSEIF ( jlev > 0 ) THEN

            !! User possibly specified jlev = 1 and there is not an existing level dimension:
            IF ( n3 == -1 ) THEN
               PRINT *, ' *** warning: ',trim(crtn),' => there is actually no levels for ', trim(cv_in),' in ',trim(cf_in)
               PRINT *, '              => fixing it...'
               jlev = 0 ! => should be treated at next "while" loop...

            ELSE
               IF ( jlev >  n3 ) CALL print_err(crtn, ' you want extract a level greater than max value')
               IF ( kt == 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev/), count=(/lx,ly,1/)), crtn,cf_in,cv_in)
               ELSEIF ( kt > 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev,kt/), count=(/lx,ly,1,1/)), crtn,cf_in,cv_in)
               END IF
               l_okay = .TRUE.  ! we can exit the WHILE loop...
            END IF

         END IF

      END DO

      ! Closing when needed:
      IF ( (( kt == ite ).OR.( kt == 0 )).AND.(jlev == kz_stop) )  THEN
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
         idx_f = 0 ; idx_v = 0
      END IF

   END SUBROUTINE GETVAR_2D



   SUBROUTINE GETVAR_2D_R8(idx_f, idx_v, cf_in, cv_in, lt, kz, kt, X, jt1, jt2, lz)

      !!-----------------------------------------------------------------------------
      !! This routine extract a 2D field from a netcdf file
      !! at a given time
      !!
      !! INPUT :
      !! -------
      !!          * idx_f    : ID of current file                  (integer)
      !!          * idx_v    : ID of current variable              (integer)
      !!          * cf_in      : name of the input file              (character)
      !!          * cv_in      : name of the variable                (character)
      !!          * lt        : time dimension of the variable     (integer)
      !!          * kz        : level to extract                    (integer)
      !!                      0 => input file does not have levels ( 1 would work anyway...)
      !!
      !!          * kt        : time snapshot to extract            (integer)
      !!                      0 => input file does not have a time snapshot
      !!                           (= old GETVAR_2D_NOTIME)
      !!
      !!
      !! OUTPUT :
      !! --------
      !!          * X         : 2D array contening the variable     (double)
      !!
      !! OPTIONAL INPUT :
      !! ----------------
      !!          * jt1, jt2  : first and last time snapshot to extract
      !!          *       lz  : number of levels to know when they are all read
      !!                        so we can close the file
      !!
      !!------------------------------------------------------------------------

      INTEGER,                   INTENT(inout) :: idx_f, idx_v
      CHARACTER(len=*),          INTENT(in)    :: cf_in, cv_in
      INTEGER,                   INTENT(in)    :: lt, kz, kt
      REAL(8),  DIMENSION (:,:), INTENT(out)   :: X
      INTEGER,       OPTIONAL,   INTENT(in)    :: jt1, jt2, lz

      INTEGER :: &
         & lx, &    ! x dimension of the variable        (integer)
         & ly       ! y dimension of the variable        (integer)

      INTEGER :: n1, n2, n3, n4, jlev, its, ite, kz_stop = 0

      LOGICAL :: l_okay

      crtn = 'GETVAR_2D_R8'

      lx = size(X,1)
      ly = size(X,2)

      jlev = kz ! so we can modify jlev without affecting kz...

      CALL DIMS(cf_in, cv_in, n1, n2, n3, n4)

      IF ( (lx /= n1).OR.(ly /= n2) ) CALL print_err(crtn, ' PROBLEM #1 => '//trim(cv_in)//' in '//trim(cf_in))
      !PRINT *, ' n4, lt =>', n4, lt
      IF ( (lt > 0).AND.(n4 /= lt)  ) CALL print_err(crtn, ' PROBLEM #2  => '//trim(cv_in)//' in '//trim(cf_in))

      IF ( present(jt1).AND.present(jt2) ) THEN
         its = jt1 ; ite = jt2
      ELSE
         its = 1   ; ite = lt
      END IF

      IF ( present(lz) ) kz_stop = lz

      IF ( (kt == its).OR.(kt == 0) ) THEN   ! Opening file and defining variable :
         CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, idx_f),  crtn,cf_in,cv_in)
         CALL sherr( NF90_INQ_VARID(idx_f, cv_in, idx_v),  crtn,cf_in,cv_in)
      END IF


      l_okay = .FALSE.
      DO WHILE ( .NOT. l_okay )

         IF ( jlev == 0 ) THEN    ! No levels

            IF ( kt == 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1/), count=(/lx,ly/)), &
                  &      crtn,cf_in,cv_in)
            ELSEIF ( kt > 0 ) THEN
               CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,kt/), count=(/lx,ly,1/)), &
                  &      crtn,cf_in,cv_in)
            END IF
            l_okay = .TRUE. ! we can exit the WHILE loop...

         ELSEIF ( jlev > 0 ) THEN

            !! User possibly specified jlev = 1 and there is not an existing level dimension:
            IF ( n3 == -1 ) THEN
               PRINT *, ' *** warning: ',trim(crtn),' => there is actually no levels for ', trim(cv_in),' in ',trim(cf_in)
               PRINT *, '              => fixing it...'
               jlev = 0 ! => should be treated at next "while" loop...

            ELSE
               IF ( jlev >  n3 ) CALL print_err(crtn, ' you want extract a level greater than max value')
               IF ( kt == 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev/), count=(/lx,ly,1/)), crtn,cf_in,cv_in)
               ELSEIF ( kt > 0 ) THEN
                  CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,jlev,kt/), count=(/lx,ly,1,1/)), crtn,cf_in,cv_in)
               END IF
               l_okay = .TRUE.  ! we can exit the WHILE loop...
            END IF

         END IF

      END DO

      ! Closing when needed:
      IF ( (( kt == ite ).OR.( kt == 0 )).AND.(jlev == kz_stop) )  THEN
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
         idx_f = 0 ; idx_v = 0
      END IF

   END SUBROUTINE GETVAR_2D_R8













   SUBROUTINE GETVAR_3D(idx_f, idx_v, cf_in, cv_in, lt, kt, X, jt1, jt2)

      !!------------------------------------------------------------------
      !! This routine extract a 3D field from a netcdf file
      !! at a given time
      !!
      !! INPUT :
      !! -------
      !!          * idx_f    : ID of current file                  (integer)
      !!          * idx_v    : ID of current variable              (integer)
      !!          * cf_in      : name of the input file              (character)
      !!          * cv_in      : name of the variable                (character)
      !!          * lt        : time dimension of the variable     (integer)
      !!
      !!          * kt        : time snapshot to extract            (integer)
      !!                      0 => input file does not have a time snapshot
      !!                           (= old GETVAR_2D_NOTIME)
      !!
      !! OUTPUT :
      !! --------
      !!          * X         : 3D array contening the variable     (real)
      !!
      !! OPTIONAL INPUT :
      !! ----------------
      !!          * jt1, jt2  : first and last time snapshot to extract
      !!
      !!------------------------------------------------------------------------

      INTEGER,                    INTENT(inout) :: idx_f, idx_v
      CHARACTER(len=*),           INTENT(in)    :: cf_in, cv_in
      INTEGER,                    INTENT(in)    :: lt, kt
      REAL(4), DIMENSION (:,:,:), INTENT(out)   :: X
      INTEGER,       OPTIONAL,    INTENT(in)    :: jt1, jt2

      INTEGER ::  &
         & lx,  &        ! x dimension of the variable        (integer)
         & ly,  &        ! y dimension of the variable        (integer)
         & lz            ! z dimension of the variable        (integer)

      INTEGER :: n1, n2, n3, n4, its, ite

      crtn = 'GETVAR_3D'

      lx = size(X,1)
      ly = size(X,2)
      lz = size(X,3)

      CALL DIMS(cf_in, cv_in, n1, n2, n3, n4)

      IF ( (lx /= n1).OR.(ly /= n2) ) CALL print_err(crtn, ' PROBLEM #1 => '//TRIM(cv_in)//' in '//TRIM(cf_in))

      IF ( PRESENT(jt1).AND.PRESENT(jt2) ) THEN
         its = jt1 ; ite = jt2
      ELSE
         its = 1   ; ite = lt
      END IF

      !IF ( kt <= 1 ) THEN   ! Opening file and defining variable :
      IF ( (kt == its).OR.(kt == 0) ) THEN   ! Opening file and defining variable :
         CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE,  idx_f),  crtn,cf_in,cv_in)
         CALL sherr( NF90_INQ_VARID(idx_f, cv_in, idx_v),  crtn,cf_in,cv_in)
      END IF

      IF ( kt == 0 ) THEN
         CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,1/), count=(/lx,ly,lz/)), &
            &      crtn,cf_in,cv_in)
      ELSEIF ( kt > 0 ) THEN
         CALL sherr( NF90_GET_VAR(idx_f, idx_v, X, start=(/1,1,1,kt/), count=(/lx,ly,lz,1/)), &
            &      crtn,cf_in,cv_in)
      END IF

      IF ( ( kt == ite ).OR.( kt == 0 ) )  THEN
         !IF (( kt == lt ).OR.( kt == 0 ))  THEN
         CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)
         idx_f = 0 ; idx_v = 0
      END IF

   END SUBROUTINE GETVAR_3D




   SUBROUTINE GETMASK_2D(cf_in, cv_in, IX, jlev)

      !!-----------------------------------------------------------------------
      !!  Get mask (variable 'cv_in') from a netcdf file.
      !! - mask is stored in integer array IX
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file          (character)
      !!          * cv_in      : name of mask variable           (character)
      !!
      !!          * jlev      : level to get (0 if no levels) |OPTIONAL|     (integer)
      !!
      !! OUTPUT :
      !! --------
      !!          * IX        :  2D array contening mask       (integer)
      !!
      !!------------------------------------------------------------------------

      CHARACTER(len=*),        INTENT(in)  :: cf_in, cv_in
      INTEGER, OPTIONAL,       INTENT(in)  :: jlev
      INTEGER(2), DIMENSION(:,:), INTENT(out) :: IX

      INTEGER :: &
         & lx, &    ! x dimension of the mask
         & ly       ! y dimension of the mask

      INTEGER :: nx, ny, nk, nt, icz

      crtn = 'GETMASK_2D'

      lx = size(IX,1)
      ly = size(IX,2)

      icz = 1 ! getting mask at level 1 for default

      IF ( present(jlev) ) THEN
         IF ( jlev > 0 ) THEN
            icz = jlev ; WRITE(6,*) 'Getting mask at level', icz
         ELSE
            CALL print_err(crtn, 'you cannot specify a level jlev <= 0')
         END IF
      END IF

      CALL DIMS(cf_in, cv_in, nx, ny, nk, nt)

      IF ( (nx /= lx).OR.(ny /= ly) ) CALL print_err(crtn, 'data and mask file dont have same horizontal dimensions')



      !!    Opening MASK netcdf file
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, trim(cv_in), id_v),  crtn,cf_in,cv_in)


      IF ( nk > 0 ) THEN

         !! Mask is 3D
         !! ~~~~~~~~~~

         IF ( .NOT. present(jlev) ) THEN
            WRITE(6,*) trim(crtn),': WARNING => mask is 3D, should specify a level to extract, defaulting to 1st level!'
         END IF

         !!  3D+T
         IF ( nt > 0 ) THEN
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,icz,1/), count=(/nx,ny,1,1/)),  &
               &      crtn,cf_in,cv_in)
         ELSE
            !! 3D
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,icz/), count=(/nx,ny,1/)), &
               &      crtn,cf_in,cv_in)
         END IF

      ELSE

         !! Mask is 2D
         !! ~~~~~~~~~~
         IF ( present(jlev) ) THEN
            IF (jlev > 1) CALL print_err(crtn, 'you want mask at a given level (jlev > 1) but mask is 2D!!!')
         END IF


         IF ( nt > 0 ) THEN
            !!  2D+T
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,1/), count=(/nx,ny,1/)), crtn,cf_in,cv_in)
         ELSE
            !! 2D
            CALL sherr( NF90_GET_VAR(id_f, id_v, IX),  crtn,cf_in,cv_in)
         END IF


      END IF

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE GETMASK_2D





   SUBROUTINE GETMASK_3D(cf_in, cv_in, IX)

      !!-----------------------------------------------------------------------
      !!  Get mask (variable 'cv_in') from a netcdf file.
      !! - mask is stored in integer array IX
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file          (character)
      !!          * cv_in      : name of mask variable           (character)
      !!
      !! OUTPUT :
      !! --------
      !!          * IX        :  3D array contening mask       (integer)
      !!
      !!  Author :            Laurent Brodeau
      !!  --------
      !!------------------------------------------------------------------------

      CHARACTER(len=*),          INTENT(in)  :: cf_in, cv_in
      INTEGER(2), DIMENSION(:,:,:), INTENT(out) :: IX

      INTEGER :: &
         & lx, &    ! x dimension of the mask
         & ly, &    ! y dimension of the mask
         & lz       ! z dimension of the mask

      INTEGER :: nx, ny, nk, nt

      crtn = 'GETMASK_3D'

      lx = size(IX,1)
      ly = size(IX,2)
      lz = size(IX,3)


      CALL DIMS(cf_in, cv_in, nx, ny, nk, nt)

      IF ( nk < 1 ) THEN
         WRITE(6,*) 'mask 3D file => ', trim(cf_in)
         WRITE(6,*) 'mask 3D name => ', trim(cv_in)
         CALL print_err(crtn, 'mask is not 3D')
      END IF

      IF ( (nx /= lx).OR.(ny /= ly).OR.(nk /= lz) ) THEN
         !&   CALL print_err(crtn, 'data and mask file dont have same dimensions => '\\)
         PRINT *, ''
         WRITE(6,*) 'ERROR in ',TRIM(crtn),' (io_ezcdf.f90): '
         WRITE(6,*) 'data and mask file dont have same dimensions'
         WRITE(6,*) '  => nx, ny, nk =', nx, ny, nk
         WRITE(6,*) '  => lx, ly, lz =', lx, ly, lz
         PRINT *, ''
         STOP
      END IF

      !!    Opening MASK netcdf file
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),     crtn,cf_in,cv_in)
      CALL sherr( NF90_INQ_VARID(id_f, trim(cv_in), id_v),  crtn,cf_in,cv_in)

      IF ( nt > 0 ) THEN
         !! 3D+T
         CALL sherr( NF90_GET_VAR(id_f, id_v, IX, start=(/1,1,1,1/), count=(/nx,ny,nk,1/)), &
            &      crtn,cf_in,cv_in)

      ELSE
         !! 3D
         CALL sherr( NF90_GET_VAR(id_f, id_v, IX),  crtn,cf_in,cv_in)
      END IF

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE GETMASK_3D




   SUBROUTINE PT_SERIES(vtime, vseries, cf_in, cv_t, cv_in, cunit, cln, vflag, &
      &               cun_t, lpack)

      !! INPUT :
      !! -------
      !!        vtime  = time array                               [array 1D real8]
      !!        vseries = 1D array containing time-series         [array 1D real4]
      !!        cf_in  = name of the output file                  [character]
      !!        cv_t = name of time                               [character]
      !!        cv_in  = name of the variable                     [character]
      !!        cunit  = unit for treated variable                [character]
      !!        cln = long-name for treated variable              [character]
      !!        vflag = flag value or "0."                        [real]
      !!
      !!        cun_t = unit for time                 |OPTIONAL|  [character]
      !!        lpack = pack/compress data (netcdf4)  |OPTIONAL|  [logical]
      !!
      !!--------------------------------------------------------------------------

      REAL(8), DIMENSION(:),     INTENT(in)   :: vtime
      REAL(4), DIMENSION(:),      INTENT(in)  :: vseries
      CHARACTER(len=*),           INTENT(in)  :: cf_in, cv_t, cv_in, cunit, cln
      REAL(4),                    INTENT(in)  :: vflag
      CHARACTER(len=*), OPTIONAL, INTENT(in)  :: cun_t
      LOGICAL,          OPTIONAL, INTENT(in)  :: lpack
      !!
      INTEGER          :: idf, idv, idtd, idt, nbt, jt
      LOGICAL          :: lp = .FALSE.
      REAL(4)          :: rmin, rmax

      crtn = 'PT_SERIES'

      nbt = size(vseries,1)

      IF ( nbt /= size(vtime,1) ) CALL print_err(crtn, 'Time array and series array dont agree in size!!!')

      IF ( present(lpack) ) THEN
         IF ( lpack ) lp = .TRUE.
      END IF


      IF ( vflag /= 0.) THEN
         rmin =  1.E6 ; rmax = -1.E6
         DO jt = 1, nbt
            IF ((vseries(jt) <= rmin).and.(vseries(jt) /= vflag)) rmin = vseries(jt)
            IF ((vseries(jt) >= rmax).and.(vseries(jt) /= vflag)) rmax = vseries(jt)
         END DO
      ELSE
         rmin = minval(vseries) ; rmax = maxval(vseries)
      END IF

      cu = 'unknown'
      IF ( present(cun_t) ) cu = trim(cun_t)


      vextrema(3,:) = (/minval(vtime),maxval(vtime)/)


      !!           CREATE NETCDF OUTPUT FILE :

      IF ( lp ) THEN
         CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, idf), crtn,cf_in,cv_in)
      ELSE
         CALL sherr( NF90_CREATE(cf_in, NF90_CLOBBER, idf), crtn,cf_in,cv_in)
      END IF

      !! Time
      CALL sherr( NF90_DEF_DIM(idf, trim(cv_t), NF90_UNLIMITED, idtd),      crtn,cf_in,cv_in)
      CALL sherr( NF90_DEF_VAR(idf, trim(cv_t), NF90_DOUBLE,    idtd, idt), crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(idf, idt, 'units',    trim(cun_t)),              crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(idf, idt, 'valid_min', vextrema(3,1)),          crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(idf, idt, 'valid_max', vextrema(3,2)),          crtn,cf_in,cv_in)

      !! Variable
      IF ( lp ) THEN
         CALL sherr( NF90_DEF_VAR(idf, trim(cv_in), NF90_FLOAT, idtd, idv, deflate_level=9), &
            &      crtn,cf_in,cv_in )
      ELSE
         CALL sherr( NF90_DEF_VAR(idf, trim(cv_in), NF90_FLOAT, idtd, idv                 ), &
            &      crtn,cf_in,cv_in )
      END IF
      CALL sherr( NF90_PUT_ATT(idf, idv, 'long_name', trim(cln)),  crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(idf, idv, 'units', trim(cunit) ),   crtn,cf_in,cv_in)
      IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idf, idv,'_FillValue',vflag),  crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(idf, idv,'valid_range', (/rmin,rmax/)),  crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(idf, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_in)


      CALL sherr( NF90_ENDDEF(idf),  crtn,cf_in,cv_in)

      !!       Write time variable :
      CALL sherr( NF90_PUT_VAR(idf, idt, vtime),    crtn,cf_in,cv_in)

      !!      WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(idf, idv, vseries),  crtn,cf_in,cv_in)

      CALL sherr( NF90_CLOSE(idf),  crtn,cf_in,cv_in)

   END SUBROUTINE PT_SERIES







   SUBROUTINE P2D_T(idx_f, idx_v, lt, lct, xlon, xlat, vtime, x2d, cf_in, &
      &           cv_lo, cv_la, cv_t, cv_in, cunit, cln, vflag, cun_t,  &
      &           lpack, cextrainfo)
      !!
      !! INPUT :
      !! -------
      !!        idx_f = ID of the file (takes its value on the first call)
      !!        idx_v = ID of the variable //
      !!        lt    = t dimension of array to plot              [integer]
      !!        lct   = current time step                         [integer]
      !!        xlon  = 2D array of longitude  (nx,ny) or (nx,1)  [double]
      !!        xlat  = 2D array of latitude   (nx,ny) or (ny,1)  [double]
      !!        vtime  = time array                               [array 1D]
      !!        x2d = 2D snap of 2D+T array at time jt to write   [real]
      !!        cf_in  = name of the output file                  [character]
      !!        cv_lo = name of longitude                         [character]
      !!        cv_la = name of latitude                          [character]
      !!        cv_t = name of time                               [character]
      !!        cv_in  = name of the variable                     [character]
      !!        cunit  = unit for treated variable                [character]
      !!        cln = long-name for treated variable              [character]
      !!        vflag = flag value or "0."                        [real]
      !!
      !!        cun_t = unit for time                 |OPTIONAL|  [character]
      !!        lpack = pack/compress data (netcdf4)  |OPTIONAL|  [logical]
      !!        cextrainfo = extra information to go in "Info" of header of netcdf 
      !!
      !!--------------------------------------------------------------------------
      !!
      INTEGER,                    INTENT(inout) :: idx_f, idx_v
      INTEGER,                    INTENT(in)    :: lt, lct
      REAL(8), DIMENSION(:,:),    INTENT(in)    :: xlat, xlon
      REAL(4), DIMENSION(:,:),    INTENT(in)    :: x2d
      REAL(8), DIMENSION(lt),     INTENT(in)    :: vtime
      CHARACTER(len=*),           INTENT(in)    :: cf_in, cv_lo, cv_la, cv_t, cv_in, cunit, cln
      REAL(4),                    INTENT(in)    :: vflag
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: cun_t
      LOGICAL,          OPTIONAL, INTENT(in)    :: lpack
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: cextrainfo
      
      INTEGER          :: lx, ly
      LOGICAL          :: lp = .FALSE.
      REAL(4)          :: rmin, rmax

      crtn = 'P2D_T'


      !! About dimensions of xlon, xlat and x2d:
      CALL ctest_coor(xlon, xlat, x2d, cdt)
      lx = size(x2d,1) ; ly = size(x2d,2)

      IF ( present(lpack) ) THEN
         IF ( lpack ) lp = .TRUE.
      END IF

      IF ( lct == 1 ) THEN

         IF ( vflag /= 0.) THEN
            rmin =  1.E6 ; rmax = -1.E6
            DO jj=1, ly
               DO ji=1, lx
                  IF ((x2d(ji,jj) <= rmin).and.(x2d(ji,jj) /= vflag)) rmin = x2d(ji,jj)
                  IF ((x2d(ji,jj) >= rmax).and.(x2d(ji,jj) /= vflag)) rmax = x2d(ji,jj)
               END DO
            END DO
         ELSE
            rmin = minval(x2d) ; rmax = maxval(x2d)
         END IF

         cu = 'unknown'
         IF ( present(cun_t) ) cu = trim(cun_t)

      END IF ! lct == 1

      IF ( lct == 1 ) THEN

         vextrema(1,:) = (/minval(xlon),maxval(xlon)/); vextrema(2,:) = (/minval(xlat),maxval(xlat)/)
         vextrema(3,:) = (/minval(vtime),maxval(vtime)/)

         !! Opening mesh file for grid quest :
         !! ----------------------------------
         !!
         !!           CREATE NETCDF OUTPUT FILE :

         IF ( lp ) THEN
            CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, idx_f), crtn,cf_in,cv_in)
         ELSE
            CALL sherr( NF90_CREATE(cf_in, NF90_CLOBBER, idx_f), crtn,cf_in,cv_in)
         END IF
         !!
         CALL prepare_nc(idx_f, cdt, lx, ly, cv_lo, cv_la, cv_t, cu, vextrema, &
            &          id_x, id_y, id_t, id_lo, id_la, id_tim, crtn,cf_in,cv_in)
         !!
         !! Variable
         IF ( lp ) THEN
            CALL sherr( NF90_DEF_VAR(idx_f, trim(cv_in), NF90_FLOAT,(/id_x,id_y,id_t/), idx_v, deflate_level=9), &
               &      crtn,cf_in,cv_in )
         ELSE
            CALL sherr( NF90_DEF_VAR(idx_f, trim(cv_in), NF90_FLOAT,(/id_x,id_y,id_t/), idx_v                 ), &
               &      crtn,cf_in,cv_in )
         END IF

         !!  VARIABLE ATTRIBUTES
         !! Long name
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v, 'long_name', trim(cln)),  crtn,cf_in,cv_in)
         !! Units
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v, 'units', trim(cunit) ),   crtn,cf_in,cv_in)
         !!
         IF ( vflag /= 0. ) CALL sherr( NF90_PUT_ATT(idx_f, idx_v,'_FillValue',vflag),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v,'valid_range', (/rmin,rmax/)),  crtn,cf_in,cv_in)

         !! Coordinates
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v,'coordinates', &
            &        TRIM(cv_t)//" "//TRIM(cv_lo)//" "//TRIM(cv_la) ), &
            &  crtn,cf_in,cv_in)

         !! Coordinates
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v,'coordinates', trim(cv_lo)//" "//trim(cv_la) ),   crtn,cf_in,cv_in)
         !! Global attributes
         IF ( PRESENT(cextrainfo) ) &
            CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'Info', TRIM(cextrainfo)),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_in)

         !!           END OF DEFINITION
         CALL sherr( NF90_ENDDEF(idx_f),  crtn,cf_in,cv_in)

         !!       Write longitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_lo, xlon),  crtn,cf_in,cv_in)
         !!
         !!       Write latitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_la, xlat),  crtn,cf_in,cv_in)
         !!
         !!       Write time variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_tim, vtime),  crtn,cf_in,cv_in)
         !!
      END IF

      !!               WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(idx_f, idx_v, x2d,  start=(/1,1,lct/), count=(/lx,ly,1/)),  crtn,cf_in,cv_in)

      IF ( lct == lt ) CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)

   END SUBROUTINE P2D_T





   SUBROUTINE P3D_T(idx_f, idx_v, lt, lct, xlon, xlat, vdpth, vtime, x3d, cf_in, &
      &           cv_lo, cv_la, cv_dpth, cv_t, cv_in, cunit, cln, vflag, &
      &           cun_t, lpack, cun_z, cextrainfo)

      !! INPUT :
      !! -------
      !!        idx_f = ID of the file (takes its value on the first call)
      !!        idx_v = ID of the variable //
      !!        lt    = t dimension of array to plot              [integer]
      !!        lct   = current time step                         [integer]
      !!        xlon  = 2D array of longitude  (nx,ny) or (nx,1)  [double]
      !!        xlat  = 2D array of latitude   (nx,ny) or (ny,1)  [double]
      !!        vdpth = depth array                               [array 1D double]
      !!        vtime  = time array                               [array 1D double]
      !!        x3d = 3D snap of 3D+T array at time jt to write   [real]
      !!        cf_in  = name of the output file                  [character]
      !!        cv_lo = name of longitude                         [character]
      !!        cv_la = name of latitude                          [character]
      !!        cv_dpth = name of depth                           [character]
      !!        cv_t = name of time                               [character]
      !!        cv_in  = name of the variable                     [character]
      !!        cunit  = unit for treated variable                [character]
      !!        cln = long-name for treated variable              [character]
      !!        vflag = flag value or "0."                        [real]
      !!
      !!        cun_t = unit for time                 |OPTIONAL|  [character]
      !!        lpack = pack/compress data (netcdf4)  |OPTIONAL|  [logical]
      !!        cun_z = unit for depth                |OPTIONAL|  [character]
      !!        cextrainfo = extra information to go in "Info" of header of netcdf 
      !!
      !!--------------------------------------------------------------------------

      INTEGER,                    INTENT(inout) :: idx_f, idx_v
      INTEGER,                    INTENT(in)    :: lt, lct
      REAL(4), DIMENSION(:,:,:),  INTENT(in)    :: x3d
      REAL(8), DIMENSION(:,:),    INTENT(in)    :: xlat, xlon
      REAL(8), DIMENSION(:),      INTENT(in)    :: vdpth
      REAL(8), DIMENSION(lt),     INTENT(in)    :: vtime
      CHARACTER(len=*),           INTENT(in)    :: cf_in, cv_lo, cv_la, cv_dpth, cv_t, cv_in, cunit, cln
      REAL(4),                    INTENT(in)    :: vflag
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: cun_t, cun_z
      LOGICAL,          OPTIONAL, INTENT(in)    :: lpack
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: cextrainfo

      INTEGER          :: lx, ly, lz
      LOGICAL          :: lp = .FALSE.
      REAL(4)          :: dr, rmin, rmax

      crtn = 'P3D_T'

      !! About dimensions of xlon, xlat, vdpth and x3d:
      CALL ctest_coor(xlon, xlat, x3d(:,:,1), cdt)
      lx = size(x3d,1) ; ly = size(x3d,2) ; lz = size(vdpth)
      IF ( size(x3d,3) /= lz ) CALL print_err(crtn, 'depth array do not match data')
      !!
      IF ( present(lpack) ) THEN
         IF ( lpack ) lp = .TRUE.
      END IF
      !!
      !!
      IF ( lct == 1 ) THEN
         !!
         IF ( vflag /= 0.) THEN
            rmin =  1.E6 ; rmax = -1.E6
            DO jk=1, lz
               DO jj=1, ly

                  DO ji=1, lx
                     IF ((x3d(ji,jj,jk) <= rmin).and.(x3d(ji,jj,jk) /= vflag)) rmin = x3d(ji,jj,jk)
                     IF ((x3d(ji,jj,jk) >= rmax).and.(x3d(ji,jj,jk) /= vflag)) rmax = x3d(ji,jj,jk)
                  END DO
               END DO
            END DO
         ELSE
            rmin = minval(x3d) ; rmax = maxval(x3d)
         END IF

         dr = (rmax - rmin)/10.0 ; rmin = rmin - dr ; rmax = rmax + dr

         cu = 'unknown'
         IF ( present(cun_t) ) cu = cun_t
         !!
      END IF ! lct == 1

      IF ( lct == 1 ) THEN
         !!
         vextrema(1,:) = (/minval(xlon),maxval(xlon)/); vextrema(2,:) = (/minval(xlat),maxval(xlat)/)
         vextrema(3,:) = (/minval(vtime),maxval(vtime)/)
         !!
         !! Opening mesh file for grid quest :
         !! ----------------------------------
         !!
         !!           CREATE NETCDF OUTPUT FILE :
         IF ( lp ) THEN
            CALL sherr( NF90_CREATE(cf_in, NF90_NETCDF4, idx_f),  crtn,cf_in,cv_in)
         ELSE
            CALL sherr( NF90_CREATE(cf_in, NF90_CLOBBER, idx_f),  crtn,cf_in,cv_in)
         END IF
         !!
         CALL prepare_nc(idx_f, cdt, lx, ly, cv_lo, cv_la, cv_t, cu, vextrema, &
            &          id_x, id_y, id_t, id_lo, id_la, id_tim, crtn,cf_in,cv_in)

         IF ( (trim(cv_dpth) == 'lev').OR.(trim(cv_dpth) == 'depth') ) THEN
            CALL sherr( NF90_DEF_DIM(idx_f, trim(cv_dpth), lz, id_z),  crtn,cf_in,cv_in)
         ELSE
            CALL sherr( NF90_DEF_DIM(idx_f, 'z', lz, id_z),  crtn,cf_in,cv_in)
         END IF

         CALL sherr( NF90_DEF_VAR(idx_f, trim(cv_dpth), NF90_DOUBLE, id_z,id_dpt),  crtn,cf_in,cv_in)
         cu = 'unknown'
         IF ( present(cun_z) ) cu = cun_z
         CALL sherr( NF90_PUT_ATT(idx_f, id_dpt, 'units',     trim(cu)),       crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, id_dpt, 'valid_min', minval(vdpth)),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, id_dpt, 'valid_max', maxval(vdpth)),  crtn,cf_in,cv_in)

         !! Variable
         IF ( lp ) THEN
            CALL sherr( NF90_DEF_VAR(idx_f, TRIM(cv_in), NF90_FLOAT, (/id_x,id_y,id_z,id_t/), idx_v, deflate_level=9), &
               &       crtn,cf_in,cv_in)
         ELSE
            CALL sherr( NF90_DEF_VAR(idx_f, trim(cv_in), NF90_FLOAT, (/id_x,id_y,id_z,id_t/), idx_v),  &
               &       crtn,cf_in,cv_in)
         END IF

         !!  VARIABLE ATTRIBUTES
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v, 'long_name', trim(cln)),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v, 'units', trim(cunit) ),   crtn,cf_in,cv_in)

         IF ( vflag /= 0.) CALL sherr( NF90_PUT_ATT(idx_f, idx_v,'_FillValue',vflag),  crtn,cf_in,cv_in)

         CALL sherr( NF90_PUT_ATT(idx_f, idx_v,'valid_range', (/rmin,rmax/)),  crtn,cf_in,cv_in)

         !! Coordinates
         CALL sherr( NF90_PUT_ATT(idx_f, idx_v,'coordinates', &
            &        TRIM(cv_t)//" "//TRIM(cv_dpth)//" "//TRIM(cv_lo)//" "//TRIM(cv_la) ), &
            &  crtn,cf_in,cv_in)

         !! Global attributes
         IF ( PRESENT(cextrainfo) ) &
            CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'Info', TRIM(cextrainfo)),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_ATT(idx_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_in)
         !!
         !!           END OF DEFINITION
         CALL sherr( NF90_ENDDEF(idx_f),  crtn,cf_in,cv_in)
         !!
         !!
         !!
         !!       Write longitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_lo, xlon),  crtn,cf_in,cv_in)
         !!
         !!       Write latitude variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_la, xlat),  crtn,cf_in,cv_in)
         !!
         !!       Write depth variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_dpt, vdpth),  crtn,cf_in,cv_in)
         !!
         !!       Write time variable :
         CALL sherr( NF90_PUT_VAR(idx_f, id_tim, vtime),  crtn,cf_in,cv_in)
         !!
      END IF

      !!                WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(idx_f, idx_v,  x3d, start=(/1,1,1,lct/), count=(/lx,ly,lz,1/)),  crtn,cf_in,cv_in)

      !! Sync data from buffer to file
      IF ( lct /= lt ) CALL sherr( NF90_SYNC (idx_f),  crtn,cf_in,cv_in)
      IF ( lct == lt ) CALL sherr( NF90_CLOSE(idx_f),  crtn,cf_in,cv_in)

   END SUBROUTINE P3D_T




   SUBROUTINE CHECK_4_MISS(cf_in, cv_in, lmv, rmissval, cmiss)
      !!
      !! o This routine looks for the presence of a missing value attribute
      !!   of variable cv_in into file cf_in
      !!
      !! INPUT :
      !! -------
      !!         * cv_in    = variable                                [character]
      !!         * cf_in    = treated file                              [character]
      !!
      !! OUTPUT :
      !! --------
      !!         * imiss    = 0 -> no missing value, 1 -> missing value found   [integer]
      !!         * rmissval = value of missing value                          [real]
      !!         * [cmiss]  = name of the missing value arg. |OPTIONAL|  [character]
      !!
      !! Author : L. BRODEAU, december 2008
      !!
      !!----------------------------------------------------------------------------
      !!
      CHARACTER(len=*), INTENT(in)  :: cf_in, cv_in
      !!
      LOGICAL,            INTENT(out) :: lmv
      REAL(4),         INTENT(out) :: rmissval
      !!
      CHARACTER(len=*) , OPTIONAL, INTENT(in)  :: cmiss
      !!
      INTEGER :: ierr
      !!
      crtn = 'CHECK_4_MISS'
      !!
      !!
      !! Opening file :
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)
      !!
      !! Chosing variable :
      CALL sherr( NF90_INQ_VARID(id_f, cv_in, id_v),  crtn,cf_in,cv_in)
      !!
      !!
      IF ( present(cmiss) ) THEN
         ierr = NF90_GET_ATT(id_f, id_v, cmiss, rmissval)
      ELSE
         !! Default name for a missing value is "missing_value" :
         ierr = NF90_GET_ATT(id_f, id_v, 'missing_value', rmissval)
      END IF
      !!
      IF ( ierr == -43 ) THEN
         lmv = .FALSE.
         !!
      ELSE
         !!
         IF (ierr ==  NF90_NOERR) THEN
            lmv = .TRUE.
         ELSE
            CALL print_err(crtn, 'problem getting missing_value attribute')
         END IF
      END IF
      !!
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
      !!
   END SUBROUTINE CHECK_4_MISS
   !!
   !!
   SUBROUTINE GET_VAR_INFO(cf_in, cv_in, cunit, clnm)
      !!
      !! o This routine returns the unit and longname of variable if they exist!
      !!
      !! INPUT :
      !! -------
      !!         * cv_in    = variable                                [character]
      !!         * cf_in    = treated file                            [character]
      !!
      !! OUTPUT :
      !! --------
      !!         * cunit = unit of cv_in                              [character]
      !!         * clnm  = name of the missing value arg.            [character]
      !!
      !! Author : L. BRODEAU, 2008
      !!
      !!----------------------------------------------------------------------------
      !!
      CHARACTER(len=*), INTENT(in)  :: cf_in, cv_in
      CHARACTER(len=*) , INTENT(out) :: cunit
      CHARACTER(len=*), INTENT(out) :: clnm
      !!
      INTEGER :: ierr
      CHARACTER(len=400) :: c00
      !!
      crtn = 'GET_VAR_INFO'
      !!
      !!
      !! Opening file :
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)
      !!
      !! Chosing variable :
      CALL sherr( NF90_INQ_VARID(id_f, cv_in, id_v),  crtn,cf_in,cv_in)
      !!
      !!
      c00=''
      ierr = NF90_GET_ATT(id_f, id_v, 'units', c00)
      IF (ierr /= 0) c00 = 'UNKNOWN'
      cunit = trim(c00) ;
      !!
      c00=''
      ierr = NF90_GET_ATT(id_f, id_v, 'long_name', c00)
      IF (ierr /= 0) c00 = 'UNKNOWN'
      clnm = trim(c00)
      !!
      !!
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
      !!
   END SUBROUTINE GET_VAR_INFO



   SUBROUTINE PRTMASK(xmsk, cf_in, cv_in,   xlon, xlat, cv_lo, cv_la)
      !!
      !!-----------------------------------------------------------------------------
      !! This routine prints an integer array of a 2D mask into a netcdf file
      !! ( healthy minds usually agree for earth == 0, sea == 1 )
      !!
      !! INPUT :
      !! -------
      !!        xmsk  = 2D array (lx,ly) contening mask            [real4]
      !!        cf_in  = name of the output file                   [character]
      !!        cv_in  = name of the mask variable                 [character]
      !! OPTIONAL :
      !! ----------
      !!        xlon    = longitude array
      !!        xlat    = latitude array
      !!        cv_lo = longitude name
      !!        cv_la = latitude name
      !!
      !!------------------------------------------------------------------------------
      !!
      !!
      REAL(4),    DIMENSION(:,:), INTENT(in) :: xmsk
      CHARACTER(len=*),           INTENT(in) :: cf_in, cv_in
      !!
      REAL(8), DIMENSION(:,:), INTENT(in), OPTIONAL :: xlon, xlat
      CHARACTER(len=*),        INTENT(in), OPTIONAL :: cv_lo, cv_la
      !!
      INTEGER     :: lx, ly, i01, i02
      !!
      LOGICAL :: lzcoord
      !!
      crtn = 'PRTMASK'
      !!
      lx = size(xmsk,1) ; ly = size(xmsk,2)
      !!
      lzcoord = .FALSE.
      !!
      IF ( present(xlon).AND.present(xlat) ) THEN
         IF ( present(cv_lo).AND.present(cv_la) ) THEN
            lzcoord = .TRUE.
            CALL ctest_coor(xlon, xlat, xmsk, cdt)
         ELSE
            CALL print_err(crtn, 'if you specify xlon and xlat, you must also specify cv_lo and cv_la')
         END IF
         vextrema(1,:) = (/minval(xlon),maxval(xlon)/); vextrema(2,:) = (/minval(xlat),maxval(xlat)/)
      END IF
      !!

      !!
      !!
      !!           CREATE NETCDF OUTPUT FILE :
      !!           ---------------------------
      CALL sherr( NF90_CREATE(CF_IN, NF90_CLOBBER, id_f),  crtn,cf_in,cv_in)
      !!
      !!
      IF ( lzcoord ) THEN
         CALL prepare_nc(id_f, cdt, lx, ly, cv_lo, cv_la, '', '', vextrema, id_x, id_y, i01, id_lo, id_la, i02, &
            &          crtn,cf_in,cv_in)
      ELSE
         CALL prepare_nc(id_f, cdt, lx, ly, '', '', '', '',       vextrema, id_x, id_y, i01, id_lo, id_la, i02, &
            &          crtn,cf_in,cv_in)
      END IF
      !!
      CALL sherr( NF90_DEF_VAR(id_f, trim(cv_in), NF90_FLOAT, (/id_x,id_y/), id_v),       crtn,cf_in,cv_in)
      !!
      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_in,cv_in) ! END OF DEFINITION
      !!
      !!
      !!          WRITE COORDINATES
      IF ( lzcoord ) THEN
         CALL sherr( NF90_PUT_VAR(id_f, id_lo, xlon),  crtn,cf_in,cv_in)
         CALL sherr( NF90_PUT_VAR(id_f, id_la, xlat),  crtn,cf_in,cv_in)
      END IF
      !!
      !!          WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(id_f, id_v, xmsk),  crtn,cf_in,cv_in)
      !!
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
      !!
      !!
   END SUBROUTINE PRTMASK








   SUBROUTINE P2D_MAPPING_AB(cf_out, xlon, xlat, imtrcs, ralfbet, vflag, mproblem)

      CHARACTER(len=*),          INTENT(in) :: cf_out
      REAL(8), DIMENSION(:,:),   INTENT(in) :: xlon, xlat
      INTEGER, DIMENSION(:,:,:), INTENT(in) :: imtrcs
      REAL(8), DIMENSION(:,:,:), INTENT(in) :: ralfbet
      REAL(8),                   INTENT(in) :: vflag
      INTEGER, DIMENSION(:,:),   INTENT(in) :: mproblem

      INTEGER          :: lx, ly, il0, id_n2, id_n3, id_v1, id_v2, id_v3

      crtn = 'P2D_MAPPING_AB'

      lx = size(ralfbet,1) ; ly = size(ralfbet,2)

      il0 = size(ralfbet,3)
      IF ( il0 /= 2 ) THEN
         PRINT *, 'ralfbet in P2D_MAPPING_AB of io_ezcdf.f90 has wrong shape:', lx, ly, il0
         STOP
      END IF

      il0 = size(imtrcs,3)
      IF ( il0 /= 3 ) THEN
         PRINT *, 'imtrcs in P2D_MAPPING_AB of io_ezcdf.f90 has wrong shape:', lx, ly, il0
         STOP
      END IF


      !!           CREATE NETCDF OUTPUT FILE :
      !!           ---------------------------
      !CALL sherr( NF90_CREATE(cf_out, NF90_CLOBBER, id_f),  crtn,cf_out,cdum)
      CALL sherr( NF90_CREATE(cf_out, NF90_NETCDF4, id_f),  crtn,cf_out,cdum)

      CALL sherr( NF90_DEF_DIM(id_f, 'x',  lx, id_x), crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_DIM(id_f, 'y',  ly, id_y), crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_DIM(id_f, 'n2',  2, id_n2), crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_DIM(id_f, 'n3',  3, id_n3), crtn,cf_out,cdum)

      CALL sherr( NF90_DEF_VAR(id_f, 'lon',       NF90_DOUBLE, (/id_x,id_y/),       id_lo, deflate_level=9), &
         &        crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_VAR(id_f, 'lat',       NF90_DOUBLE, (/id_x,id_y/),       id_la, deflate_level=9), &
         &        crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_VAR(id_f, 'metrics',   NF90_INT,    (/id_x,id_y,id_n3/), id_v1, deflate_level=9), &
         &        crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_VAR(id_f, 'alphabeta', NF90_DOUBLE, (/id_x,id_y,id_n2/), id_v2, deflate_level=9), &
         &        crtn,cf_out,cdum)
      CALL sherr( NF90_DEF_VAR(id_f, 'iproblem',  NF90_INT,    (/id_x,id_y/),       id_v3, deflate_level=9), &
         &        crtn,cf_out,cdum)

      IF ( vflag /= 0. ) THEN
         CALL sherr( NF90_PUT_ATT(id_f, id_v1,'_FillValue',INT(vflag)),  crtn,cf_out,cdum)
         CALL sherr( NF90_PUT_ATT(id_f, id_v2,'_FillValue',vflag),       crtn,cf_out,cdum)
      END IF

      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'Info', 'File containing mapping/weight information for bilinear/bicubic interpolation with SOSIE.'), &
         &      crtn,cf_out,cdum)
      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_out,cdum)

      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_out,cdum) ! END OF DEFINITION

      CALL sherr( NF90_PUT_VAR(id_f, id_lo, xlon),     crtn,cf_out,cdum)
      CALL sherr( NF90_PUT_VAR(id_f, id_la, xlat),     crtn,cf_out,cdum)
      CALL sherr( NF90_PUT_VAR(id_f, id_v1,  imtrcs),  crtn,cf_out,cdum)
      CALL sherr( NF90_PUT_VAR(id_f, id_v2, ralfbet),  crtn,cf_out,cdum)
      CALL sherr( NF90_PUT_VAR(id_f, id_v3, mproblem), crtn,cf_out,cdum)

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_out,cdum)

   END SUBROUTINE P2D_MAPPING_AB





   SUBROUTINE  RD_MAPPING_AB(cf_in, imtrcs, ralfbet)

      CHARACTER(len=*),          INTENT(in)  :: cf_in
      INTEGER, DIMENSION(:,:,:), INTENT(out) :: imtrcs
      REAL(8), DIMENSION(:,:,:), INTENT(out) :: ralfbet

      INTEGER :: id_v1, id_v2

      crtn = 'RD_MAPPING_AB'

      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cdum)

      CALL sherr( NF90_INQ_VARID(id_f, 'metrics',   id_v1),  crtn,cf_in,cdum)
      CALL sherr( NF90_INQ_VARID(id_f, 'alphabeta', id_v2),  crtn,cf_in,cdum)


      CALL sherr( NF90_GET_VAR(id_f, id_v1, imtrcs),  crtn,cf_in,cdum)
      CALL sherr( NF90_GET_VAR(id_f, id_v2, ralfbet), crtn,cf_in,cdum)

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cdum)

   END SUBROUTINE RD_MAPPING_AB





   SUBROUTINE PHOVMOLLER(vx, vy, x2d, cf_in, cv_in, cv_x, cv_y, cunit, cln, cuX, cuY)
      !!
      !!-----------------------------------------------------------------------------
      !! This routine prints Hovmoller diagram from a 2D fiels and coordinates given as 2 vectors
      !!
      !! INPUT :
      !! -------
      !!        vx     = X array                                     [real8]
      !!        vy     = Y array                                     [real8]
      !!        x2d    = 2D array (lx,ly) contening mask             [real4]
      !!        cf_in  = name of the output file                   [character]
      !!        cv_in  = name of the mask variable                 [character]
      !!        cf_x   = name of X coordinates                     [character]
      !!        cf_y   = name of Y coordinates                     [character]
      !!        cunit  = unit for treated variable                [character]
      !!        cln = long-name for treated variable              [character]
      !!
      !! OPTIONAL :
      !! ----------
      !!        cuX  = unit X coordinate                          [character]
      !!        cuY  = unit Y coordinate                          [character]
      !!
      !!------------------------------------------------------------------------------

      REAL(8),    DIMENSION(:)  , INTENT(in) :: vx, vy
      REAL(4),    DIMENSION(:,:), INTENT(in) :: x2d
      CHARACTER(len=*),           INTENT(in) :: cf_in, cv_in, cv_x, cv_y, cunit, cln

      CHARACTER(len=*), OPTIONAL, INTENT(in) :: cuX, cuY

      CHARACTER(len=64) :: cuXin, cuYin
      INTEGER :: lx, ly, i01, i02

      crtn = 'PHOVMOLLER'

      lx = size(x2d,1) ; ly = size(x2d,2)

      cuXin = 'unknown'
      cuYin = 'unknown'
      IF ( present(cuX) ) cuXin = trim(cuX)
      IF ( present(cuY) ) cuYin = trim(cuY)



      vextrema(1,:) = (/minval(vx),maxval(vx)/); vextrema(2,:) = (/minval(vy),maxval(vy)/)

      !!           CREATE NETCDF OUTPUT FILE :
      CALL sherr( NF90_CREATE(CF_IN, NF90_CLOBBER, id_f),  crtn,cf_in,cv_in)

      CALL prepare_nc(id_f, '1d', lx, ly, cv_x, cv_y, '', '', vextrema, id_x, id_y, i01, id_lo, id_la, i02, &
         &          crtn,cf_in,cv_in, cu_X=trim(cuXin), cu_Y=trim(cuYin))

      CALL sherr( NF90_DEF_VAR(id_f, trim(cv_in), NF90_FLOAT, (/id_x,id_y/), id_v),       crtn,cf_in,cv_in)

      CALL sherr( NF90_PUT_ATT(id_f, id_v, 'long_name', trim(cln)),  crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_ATT(id_f, id_v, 'units',  trim(cunit) ),  crtn,cf_in,cv_in)
      !! Global attributes
      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_in,cv_in)

      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_in,cv_in) ! END OF DEFINITION


      !!          WRITE COORDINATES
      CALL sherr( NF90_PUT_VAR(id_f, id_lo, vx),  crtn,cf_in,cv_in)
      CALL sherr( NF90_PUT_VAR(id_f, id_la, vy),  crtn,cf_in,cv_in)

      !!          WRITE VARIABLE
      CALL sherr( NF90_PUT_VAR(id_f, id_v, x2d),  crtn,cf_in,cv_in)

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)

   END SUBROUTINE PHOVMOLLER











   SUBROUTINE GET_SF_AO(cf_in, cv_in, rsf, rao)
      !!
      !!-----------------------------------------------------------------------
      !! This routine extracts the 'scale_factor' and 'add_offset' of a given
      !! variable from a netcdf file
      !!
      !! INPUT :
      !! -------
      !!          * cf_in      : name of the input file              (character)
      !!          * cv_in      : name of the variable                (character)
      !!
      !! OUTPUT :
      !! --------
      !!          * rsf       : scale factor                        (real)
      !!          * rao       : add offset                          (real)
      !!
      !!------------------------------------------------------------------------
      !!
      CHARACTER(len=*), INTENT(in) :: cf_in, cv_in
      REAL(4),         INTENT(out) :: rsf, rao
      !!
      !! local :
      INTEGER :: ierr1, ierr2
      !!
      crtn = 'GET_SF_AO'
      !!
      !!
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,cv_in)
      !!
      CALL sherr( NF90_INQ_VARID(id_f, cv_in, id_v),  crtn,cf_in,cv_in)
      !!
      ierr1 = NF90_GET_ATT(id_f, id_v, 'scale_factor', rsf)
      ierr2 = NF90_GET_ATT(id_f, id_v, 'add_offset',   rao)
      !!
      IF ( (ierr1 /= NF90_NOERR).OR.(ierr2 /= NF90_NOERR) ) THEN
         rsf = 1.      ;   rao = 0.
         WRITE(6,*) 'WARNING: variable ', trim(cv_in), ' of file ',trim(cf_in), ' :'
         WRITE(6,*) '       does not have a "scale_factor" and "add_offset" attributes'
         WRITE(6,*) '       => scale_factor =', rsf; WRITE(6,*) '       => add_offset =', rao
         PRINT *, ''
      END IF
      !!
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,cv_in)
      !!
   END SUBROUTINE GET_SF_AO
   !!
   !!
   !!
   SUBROUTINE WHO_IS_MV(cf_in, cv_in, cmv_name, rmissval)
      !!
      !! Who is missing value???
      !!
      !!
      CHARACTER(len=*),  INTENT(in)  :: cf_in, cv_in
      CHARACTER(len=80), INTENT(out) :: cmv_name     !: real name of missing value
      REAL(4),        INTENT(out) :: rmissval  !: value of the missing value
      !!
      LOGICAL    :: lmval
      INTEGER    :: jmval
      !!
      !!
      CALL CHECK_4_MISS(cf_in, cv_in, lmval, rmissval)
      !!
      cmv_name = 'missing_value'
      !!
      !! Maybe the name of missing value is different:
      jmval = 0
      DO WHILE ( .NOT. lmval )
         !!
         jmval = jmval + 1
         cmv_name = trim(c_nm_missing_val(jmval))
         WRITE(6,*) 'Trying "',trim(cmv_name),'"! instead of "missing_value..."'
         CALL CHECK_4_MISS(cf_in, cv_in, lmval, rmissval, cmiss=cmv_name)
         !!
         IF ( ( jmval == nmval ) .AND. ( .NOT. lmval ) ) THEN
            WRITE(6,*) 'Your input file does not contain a missing value!'
            WRITE(6,*) 'Impossible to build a land sea mask array...'
            WRITE(6,*) ' -> you should maybe specify the name of the "missing value"'
            WRITE(6,*) '    found in the netcdf file into module "inter.f90"'
            STOP
         END IF
         !!
      END DO
      !!
      PRINT *, ''; WRITE(6,*) 'Missing value is called "',trim(cmv_name),'" !'
      WRITE(6,*) 'and has the value', rmissval; PRINT *, ''
      !!
      !!
      !!
   END SUBROUTINE WHO_IS_MV
   !!
   !!
   !!
   SUBROUTINE sherr(ierr, croutine, ctf, ctv)
      !!
      !! To handle and display error messages
      !!
      INTEGER,            INTENT(in) :: ierr
      !!
      CHARACTER(len=*) , INTENT(in) :: &
         &            ctf,           &    !: treated file
         &            croutine,      &    !: routine name
         &            ctv                 !: treated varible
      !!
      !!
      IF ( ierr /= NF90_NOERR ) THEN
         PRINT *, ''
         WRITE(6,*) '************************************************'
         WRITE(6,*) 'Error occured in procedure ', trim(croutine),' !'
         PRINT *, ''
         WRITE(6,*) 'Treated file     = ', trim(ctf)
         WRITE(6,*) 'Treated variable = ', trim(ctv)
         PRINT *, ''
         WRITE(6,*) '--> aborting program'
         PRINT *, ''
         WRITE(6,*) 'Netcdf message was :'
         WRITE(6,*) trim(NF90_STRERROR(ierr))
         PRINT *, ''
         WRITE(6,*) '************************************************'
         PRINT *, ''
         STOP
      END IF
      !!
   END SUBROUTINE sherr





   SUBROUTINE ctest_coor(rx, ry, rd, cdm)

      !! Testing if 2D coordinates or 1D, and if match shape of data...

      REAL(8), DIMENSION(:,:), INTENT(in)  :: rx, ry
      REAL(4), DIMENSION(:,:), INTENT(in)  :: rd
      CHARACTER(len=2)       , INTENT(out) :: cdm

      INTEGER :: ix1, ix2, iy1, iy2, id1, id2

      ix1 = size(rx,1) ; ix2 = size(rx,2)
      iy1 = size(ry,1) ; iy2 = size(ry,2)
      id1 = size(rd,1) ; id2 = size(rd,2)

      IF ( (ix2 == 1).AND.(iy2 == 1) ) THEN

         IF ( (ix1 == id1).AND.(iy1 == id2) ) THEN
            cdm = '1d'
         ELSE
            CALL print_err('cdm', 'longitude and latitude array do not match data (1d)')
         END IF

      ELSE

         IF ( (ix1 == id1).AND.(iy1 == id1).AND.(ix2 == id2).AND.(iy2 == id2) ) THEN
            cdm = '2d'
         ELSE
            CALL print_err('cdm', 'longitude and latitude array do not match data (2d)')
         END IF

      END IF

   END SUBROUTINE ctest_coor





   SUBROUTINE prepare_nc(id_file, cdt0, nx, ny, cv_lon, cv_lat, cv_time, cu_t, vxtrm, &
      &                id_ji, id_jj, id_jt, id_lon, id_lat, id_time, cri,cfi,cvi, cu_X, cu_Y)

      INTEGER,                 INTENT(in)  :: id_file, nx, ny
      CHARACTER(len=2),        INTENT(in)  :: cdt0
      CHARACTER(len=*),        INTENT(in)  :: cv_lon, cv_lat, cv_time, cu_t, cri,cfi,cvi
      REAL(8), DIMENSION(3,2), INTENT(in)  :: vxtrm
      INTEGER,                 INTENT(out) :: id_ji, id_jj, id_jt, id_lon, id_lat, id_time

      CHARACTER(len=*),        INTENT(in), OPTIONAL :: cu_X, cu_Y

      CHARACTER(len=64) :: cu_X0, cu_Y0

      cu_X0 = 'degrees_east'
      cu_Y0 = 'degrees_north'
      IF ( present(cu_X) ) cu_X0 = trim(cu_X)
      IF ( present(cu_Y) ) cu_Y0 = trim(cu_Y)




      !!    HORIZONTAL
      IF ( (trim(cv_lon) /= '').AND.(trim(cv_lat) /= '') ) THEN
         !!
         IF ( cdt0 == '2d' ) THEN
            CALL sherr( NF90_DEF_DIM(id_file, 'x', nx, id_ji), cri,cfi,cvi)
            CALL sherr( NF90_DEF_DIM(id_file, 'y', ny, id_jj), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, trim(cv_lon), NF90_DOUBLE, (/id_ji,id_jj/), id_lon), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, trim(cv_lat), NF90_DOUBLE, (/id_ji,id_jj/), id_lat), cri,cfi,cvi)
            !!
         ELSE IF ( cdt0 == '1d' ) THEN
            CALL sherr( NF90_DEF_DIM(id_file, trim(cv_lon), nx, id_ji), cri,cfi,cvi)
            CALL sherr( NF90_DEF_DIM(id_file, trim(cv_lat), ny, id_jj), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, trim(cv_lon), NF90_DOUBLE, id_ji, id_lon), cri,cfi,cvi)
            CALL sherr( NF90_DEF_VAR(id_file, trim(cv_lat), NF90_DOUBLE, id_jj, id_lat), cri,cfi,cvi)
         END IF
         !!
         CALL sherr( NF90_PUT_ATT(id_file, id_lon,  'units', trim(cu_X0)), cri,cfi,cvi)
         CALL sherr( NF90_PUT_ATT(id_file, id_lon, 'valid_min', vxtrm(1,1)), cri,cfi,cvi)
         CALL sherr( NF90_PUT_ATT(id_file, id_lon, 'valid_max', vxtrm(1,2)), cri,cfi,cvi)
         !!
         CALL sherr( NF90_PUT_ATT(id_file, id_lat, 'units', trim(cu_Y0)), cri,cfi,cvi)
         CALL sherr( NF90_PUT_ATT(id_file, id_lat, 'valid_min', vxtrm(2,1)), cri,cfi,cvi)
         CALL sherr( NF90_PUT_ATT(id_file, id_lat, 'valid_max', vxtrm(2,2)), cri,cfi,cvi)
         !!
      ELSE
         CALL sherr( NF90_DEF_DIM(id_file, 'x', nx, id_ji), cri,cfi,cvi)
         CALL sherr( NF90_DEF_DIM(id_file, 'y', ny, id_jj), cri,cfi,cvi)
      END IF
      !!
      !!  TIME
      IF ( trim(cv_time) /= '' ) THEN
         CALL sherr( NF90_DEF_DIM(id_file, trim(cv_time), NF90_UNLIMITED, id_jt), cri,cfi,cvi)
         CALL sherr( NF90_DEF_VAR(id_file, trim(cv_time), NF90_DOUBLE, id_jt, id_time), cri,cfi,cvi)
         CALL sherr( NF90_PUT_ATT(id_file, id_time, 'units',    trim(cu_t)), cri,cfi,cvi)
         CALL sherr( NF90_PUT_ATT(id_file, id_time, 'valid_min',vxtrm(3,1)), cri,cfi,cvi)
         CALL sherr( NF90_PUT_ATT(id_file, id_time, 'valid_max',vxtrm(3,2)), cri,cfi,cvi)
      END IF
      !!
   END SUBROUTINE prepare_nc




   SUBROUTINE print_err(crout, cmess)
      CHARACTER(len=*), INTENT(in) :: crout, cmess
      PRINT *, ''
      WRITE(6,*) 'ERROR in ',trim(crout),' (io_ezcdf.f90): '
      WRITE(6,*) trim(cmess) ; PRINT *, ''
      STOP
   END SUBROUTINE print_err




END MODULE io_ezcdf
