PROGRAM CORR_VECT

   USE mod_conf
   USE mod_init
   USE io_ezcdf

   !!========================================================================
   !! Purpose :  correct vector components 'uraw' and 'vraw' directly
   !!            interpolated on an irregular grid
   !! ---------  Only for ORCA family of grids!!!!
   !!
   !! Author :   Laurent Brodeau, brodeau@gmail.com
   !! --------
   !!
   !!========================================================================

   IMPLICIT NONE

   INTERFACE
      SUBROUTINE angle_dist(xa, xb, xc, xd, xe, xf)
         REAL(8), DIMENSION(:,:), INTENT(in)  :: xa, xc, xb, xd
         REAL(8), DIMENSION(:,:), INTENT(out) :: xe, xf
      END SUBROUTINE angle_dist
   END INTERFACE

   !! Grid :
   CHARACTER(len=80), PARAMETER   :: &
      &    cv_lon_t     = 'glamt',   &   ! input grid longitude name, T-points
      &    cv_lat_t     = 'gphit',   &   ! input grid latitude name,  T-points
      &    cv_lon_u     = 'glamu',   &   ! input grid longitude name, U-points
      &    cv_lat_u     = 'gphiu',   &   ! input grid latitude name,  U-points
      &    cv_lon_v     = 'glamv',   &   ! input grid longitude name, U-points
      &    cv_lat_v     = 'gphiv',   &   ! input grid latitude name,  U-points
      &    cv_depth     = 'deptht'       !  depth at T-points (U-points and V-points too)

   CHARACTER(len=3)    :: cdum
   CHARACTER(len=80)   :: cv_time_0 = 'none', cfext = 'nc'
   CHARACTER(len=800)  :: cr, cf_mm

   CHARACTER(len=80) :: &
      &    cv_out_U    = 'uraw',   &   ! raw U name
      &    cv_out_V    = 'vraw'        ! raw V name

   CHARACTER(len=800)  :: &
      &    cf_out_U,    &  ! file containing u_raw
      &    cf_out_V,    &  ! file containing v_raw
      &    cufilout, cvfilout,    &
      &    cufilin = 'none',  cvfilin = 'none'

   CHARACTER(len=80)  ::  &
      &    cv_rot_U ,  &  ! output name for U corrected
      &    cv_rot_V       ! output name for V corrected

   INTEGER      :: &
      &    jarg, i3d, nbc, &
      &    nlext=3, &
      &    i0, j0, &
      &    ni, nj, nk, nk1, nk2, &
      &    ni1, nj1, Ntr1,       &
      &    ni2, nj2, Ntr2,       &
      &    ni_g, nj_g, nk_g,    &
      &    iargc,         &
      &    idf_u, idv_u, idf_v, idv_v, &
      &    id_f1, id_v1, &
      &    id_f2, id_v2

   INTEGER(2), DIMENSION(:,:,:), ALLOCATABLE :: mask_u, mask_v

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: XCOST, XSINT, U_r, V_r


   REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: U_c, V_c

   REAL(8), DIMENSION(:,:), ALLOCATABLE ::      &
      &    XCOST8, XSINT8, U_r8, V_r8, &
      &    xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v

   REAL(8), DIMENSION(:), ALLOCATABLE ::   vtime, vdepth


   INTEGER :: jt, jk

   LOGICAL :: &
      &    l_inv = .FALSE.,   &
      &    l_3d_inv = .FALSE., &   !: will treat 3d files in inverse mode...
      &    lexist !,  &

   REAL(4), PARAMETER :: zrmv = -9999.

   CHARACTER(LEN=2), DIMENSION(9), PARAMETER :: &
      &            clist_opt = (/ '-I','-h','-m','-p','-x','-y','-f','-i','-t' /)

   WRITE(6,*)''
   WRITE(6,*)'=========================================================='
   WRITE(6,*)'            S  O  S  I  E    version ', trim(csosie_version)
   WRITE(6,*)''
   WRITE(6,*)'            Vector rotation for distorted mapping       '
   WRITE(6,*)'=========================================================='
   WRITE(6,*)''


   !! Getting string arguments :
   !! --------------------------

   !! Some defaults:
   cv_rot_U = 'vectx' ; cv_rot_V = 'vecty'

   jarg = 0

   DO WHILE ( jarg < iargc() )
      !!
      jarg = jarg + 1
      CALL getarg(jarg,cr)
      !!
      !!
      SELECT CASE (trim(cr))
         !!
         !!
      CASE('-h')
         call usage_corr_vect()
         !!
         !!
      CASE('-I')
         l_inv = .TRUE.
         !!
         !!
      CASE('-x')
         IF ( jarg + 1 > iargc() ) THEN
            PRINT *, 'ERROR: Missing zonal component name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: Missing zonal component name!'; call usage_corr_vect()
            ELSE
               cv_rot_U = trim(cr)
            END IF
         END IF
         !!
         !!
      CASE('-y')
         IF ( jarg + 1 > iargc() ) THEN
            PRINT *, 'ERROR: Missing meridional component name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: Missing meridional component name!'; call usage_corr_vect()
            ELSE
               cv_rot_V = trim(cr)
            END IF
         END IF
         !!
         !!
         !CASE('-p')
         !   PRINT *, 'Gonna pack output data in netcdf file'
         !   lpack = .TRUE.
         !!
         !!
      CASE('-f')
         IF ( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
            PRINT *, 'ERROR: Missing namelist name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: ', trim(cr), ' is definitively not the name of the namelist!'
               call usage_corr_vect()
            ELSE
               cf_nml_sosie = trim(cr)
            END IF
         END IF
         !!
      CASE('-m')
         IF ( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
            PRINT *, 'ERROR: Missing mesh_mask file name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: ', trim(cr), ' is definitively not the name of the mesh_mask file!'
               call usage_corr_vect()
            ELSE
               cf_mm = trim(cr)
            END IF
         END IF
         !!
      CASE('-i')
         IF ( jarg + 2 > iargc() ) THEN ! checking that there is at least 2 other arguments following
            PRINT *, 'ERROR: Missing input file names!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: Missing input file name 1 !' ; call usage_corr_vect()
            ELSE
               cufilin = trim(cr)
            END IF
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: Missing input file name 2 !' ; call usage_corr_vect()
            ELSE
               cvfilin = trim(cr)
            END IF
         END IF
         !!
      CASE('-t')
         IF ( jarg + 1 > iargc() ) THEN
            PRINT *, 'ERROR: Missing time name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: Missing time name!'; call usage_corr_vect()
            ELSE
               cv_time_0 = trim(cr)
            END IF
         END IF

         !!
      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
         CALL usage_corr_vect()
         !!
      END SELECT
      !!
      !!
   END DO


   IF ( (.NOT. l_inv).AND.( trim(cufilin) /= 'none' ) ) THEN
      PRINT *, 'ERROR: specify the "-I" switch if you want to perform inverse correction!'
      STOP
   END IF



   PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   PRINT *, ''

   IF ( l_inv ) THEN
      PRINT *, ' * Vector files to unrotate = ', trim(cufilin), ' , ', trim(cvfilin)
      PRINT *, '   => associated variable names = ', trim(cv_rot_U), ' , ', trim(cv_rot_V)
      IF ( trim(cv_time_0) == 'none' ) THEN
         PRINT *, 'ERROR: you must specify the name of time variable with the "-t" switch!'; STOP
      END IF
      PRINT *, '   => time variable name = ', trim(cv_time_0)
   ELSE

      PRINT *, ' * Name for corrected vector components = ', trim(cv_rot_U), ' , ', trim(cv_rot_V)
   END IF

   !PRINT *, ' * Packing output : ', lpack
   PRINT *, ' * mesh_mask file to use = ', trim(cf_mm)
   PRINT *, ' * namelist = ', trim(cf_nml_sosie)
   PRINT *, ''



   INQUIRE(FILE=trim(cf_mm), EXIST=lexist )
   IF ( .NOT. lexist ) THEN
      WRITE(*,'("The mesh_mask file ",a," file was not found!")') trim(cf_mm)
      CALL usage_corr_vect()
   END IF



   IF ( .NOT. l_inv ) THEN

      !! !!     N O R M A L   C O R R E C T I O N

      !! We need the namelist of U or V in that case:
      INQUIRE(FILE=trim(cf_nml_sosie), EXIST=lexist )
      IF ( .NOT. lexist ) THEN
         WRITE(*,'("The namelist file ",a," file was not found!")') trim(cf_nml_sosie)
         CALL usage_corr_vect()
      END IF
      PRINT *, ''
      CALL READ_NMLST(2)


      IF ( lregout ) THEN
         PRINT *, 'Vector correction only makes sense if your target grid is distorded!'
         PRINT *, '  => check "lregout" into the namelist...' ; PRINT *, ''; STOP
      END IF


      IF ( lpcknc4 ) cfext='nc4'

      cf_out_U = trim(cd_out)//'/'//trim(cv_out_U)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.'//trim(cfext)
      cf_out_V = trim(cd_out)//'/'//trim(cv_out_V)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.'//trim(cfext)

      PRINT *, 'The two input pre-interpolated needed files are :'
      PRINT *, trim(cf_out_U) ;   PRINT *, trim(cf_out_V) ;

      cufilout = trim(cd_out)//'/'//trim(cv_rot_U)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.'//trim(cfext)
      cvfilout = trim(cd_out)//'/'//trim(cv_rot_V)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.'//trim(cfext)

      PRINT *, '' ;   PRINT *, 'output files :'
      PRINT *, cufilout ;   PRINT *, cvfilout ; PRINT *, '' ; PRINT *, ''

      PRINT *, 'File containing x and y raw components of vector to be treated :'
      PRINT *, cf_out_U
      PRINT *, cf_out_V ; PRINT *, ''
      PRINT *, 'File containing grid :'
      PRINT *, cf_mm ; PRINT *, ''


      !! Geting array dimension and testing...
      !! -------------------------------------

      CALL DIMS(cf_out_U, cv_out_U, ni1, nj1, nk1, Ntr1)
      CALL DIMS(cf_out_V, cv_out_V, ni2, nj2, nk2, Ntr2)

      CALL DIMS(cf_mm, cv_lon_t, ni_g, nj_g, nk_g, Ntr)

      !! testing ni agreement :
      IF ( (ni1 /= ni2).or.(ni1 /= ni_g).or.(ni2 /= ni_g) ) THEN
         PRINT *, 'Dimension Error! : the 3 files dont agree for x length.' ; STOP
      END IF

      !! testing nj agreement :
      IF ( (nj1 /= nj2).or.(nj1 /= nj_g).or.(nj2 /= nj_g) ) THEN
         PRINT *, 'Dimension Error! : the 3 files dont agree for y length.'; STOP
      END IF

      ni = ni1 ; nj = nj1


      !! testing nk agreement :
      IF ( nk1 /= nk2 ) THEN
         PRINT *, 'Dimension Error! : u and v files dont agree for z length.'
         STOP
      END IF

      nk = nk1


      !! testing nt agreement :
      IF ( Ntr1 /= Ntr2 ) THEN
         PRINT *, 'Dimension Error! : u and v files dont agree for time length.'
         STOP
      END IF

      Ntr = Ntr1


      IF ( nk > 1 ) THEN
         i3d = 1
      ELSE
         nk  = 1
         i3d = 0
      END IF

      WRITE(*,'("Space dimension is : ",i4," x",i4," x",i4)') ni, nj, nk
      WRITE(*,'(i4," time records for u and v")') Ntr
      PRINT *, ''


      !! Allocations :
      !! -------------
      ALLOCATE (XCOST(ni,nj) , XSINT(ni,nj) , U_r(ni,nj) , V_r(ni,nj) ,  &
         &    mask_u(ni,nj,nk) , mask_v(ni,nj,nk) , U_c(ni,nj,nk) , V_c(ni,nj,nk),                              &
         &    xlon_t(ni,nj)  , xlat_t(ni,nj)  , &
         &    xlon_u(ni,nj)  , xlat_u(ni,nj) ,  xlon_v(ni,nj)  , xlat_v(ni,nj) ,  &
         &    vtime(Ntr)   )
      !!
      ALLOCATE (XCOST8(ni,nj) , XSINT8(ni,nj) , U_r8(ni,nj) , V_r8(ni,nj) )

      XCOST = 0. ; XSINT = 0. ; U_r    = 0. ; V_r = 0. ; U_c = 0. ; V_c = 0.
      xlon_t  = 0. ; xlat_t  = 0. ; xlon_u  = 0. ; xlat_u  = 0. ; xlon_v  = 0. ; xlat_v  = 0.
      vtime  = 0.



      IF ( i3d == 1 ) THEN
         ALLOCATE ( vdepth(nk) )
         CALL GETVAR_1D(cf_mm, cv_depth, vdepth)
         CALL GETMASK_3D(cf_lsm_out, 'umask', mask_u(:,:,:))
         CALL GETMASK_3D(cf_lsm_out, 'vmask', mask_v(:,:,:))
      ELSE
         CALL GETMASK_2D(cf_lsm_out, 'umask', mask_u(:,:,1), jlev=1)
         CALL GETMASK_2D(cf_lsm_out, 'vmask', mask_v(:,:,1), jlev=1)
      END IF





      !!  Getting longitude and latitude form grid file :
      !! ------------------------------------------------
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lon_t, 1, 1, 1, xlon_t)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lat_t, 1, 1, 1, xlat_t)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lon_u, 1, 1, 1, xlon_u)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lat_u, 1, 1, 1, xlat_u)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lon_v, 1, 1, 1, xlon_v)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lat_v, 1, 1, 1, xlat_v)


      !!  Getting cos and sin of the grid distorsion angle:
      !! --------------------------------------------------
      CALL angle_dist(xlon_t, xlat_t, xlon_u, xlat_u, XCOST8, XSINT8)


      !!  Getting time from the u_raw file or the namelist :
      !! ---------------------------------------------------
      IF ( lct ) THEN       ! time is being controlled
         DO jt = 1, Ntr
            vtime(jt) = t0 + t_stp*REAL(jt)
         END DO
      ELSE                  ! we use time from input file
         CALL GETVAR_1D(cf_out_U, cv_t_out, vtime)
      END IF




      DO jt = 1, Ntr


         PRINT *, ''; PRINT *, 'Time step =', jt ; PRINT *, ''

         DO jk = 1, nk

            !lolo PRINT *, 'jk =', jk

            !! Getting uncorrected U :
            !! -----------------------
            CALL  GETVAR_2D(idf_u, idv_u, cf_out_U, cv_out_U, Ntr, jk*i3d, jt, U_r, lz=nk)

            !! Getting uncorrected V :
            !! -----------------------
            CALL  GETVAR_2D(idf_v, idv_v, cf_out_V, cv_out_V, Ntr, jk*i3d, jt, V_r, lz=nk)

            U_r8 = U_r ; V_r8 = V_r

            !! Correcting U :
            !! --------------
            U_c(:,:,jk) = REAL(XCOST8*U_r8 + XSINT8*V_r8 , 4)


            IF ( (ni == 1442).and.(nj == 1021) ) THEN
               !!
               !Correction at North-Pole (avoid 2 different vectors)
               !! ecrases :
               U_c(381,1020,jk) = -0.5*U_c(1063,1019,jk) + 0.5*U_c(381,1019,jk)
               U_c(382,1020,jk) = -0.5*U_c(1062,1019,jk) + 0.5*U_c(382,1019,jk)
               U_c(383,1020,jk) = -0.5*U_c(1061,1019,jk) + 0.5*U_c(383,1019,jk)
               U_c(384,1020,jk) = -0.5*U_c(1060,1019,jk) + 0.5*U_c(384,1019,jk)
               !! ecrasant :
               U_c(1063,1020,jk) = 0.5*U_c(1063,1019,jk) - 0.5*U_c(381,1019,jk)
               U_c(1062,1020,jk) = 0.5*U_c(1062,1019,jk) - 0.5*U_c(382,1019,jk)
               U_c(1061,1020,jk) = 0.5*U_c(1061,1019,jk) - 0.5*U_c(383,1019,jk)
               U_c(1060,1020,jk) = 0.5*U_c(1060,1019,jk) - 0.5*U_c(384,1019,jk)
               !!
            END IF



            !! Correcting V :
            !! --------------
            V_c(:,:,jk) = REAL(XCOST8*V_r8 - XSINT8*U_r8 , 4)

            IF ( (ni == 1442).and.(nj == 1021) ) THEN
               ! ORCA025 north-pole correction!!!
               !!
               !Correction at North-Pole (avoid 2 different vectors)
               !! ecrases :
               V_c(381,1020,jk)  = -0.5*V_c(1063,1019,jk) + 0.5*V_c(381,1019,jk)
               V_c(382,1020,jk)  = -0.5*V_c(1062,1019,jk) + 0.5*V_c(382,1019,jk)
               V_c(383,1020,jk)  = -0.5*V_c(1061,1019,jk) + 0.5*V_c(383,1019,jk)
               V_c(384,1020,jk)  = -0.5*V_c(1060,1019,jk) + 0.5*V_c(384,1019,jk)
               !! ecrasant :
               V_c(1063,1020,jk) = 0.5*V_c(1063,1019,jk) - 0.5*V_c(381,1019,jk)
               V_c(1062,1020,jk) = 0.5*V_c(1062,1019,jk) - 0.5*V_c(382,1019,jk)
               V_c(1061,1020,jk) = 0.5*V_c(1061,1019,jk) - 0.5*V_c(383,1019,jk)
               V_c(1060,1020,jk) = 0.5*V_c(1060,1019,jk) - 0.5*V_c(384,1019,jk)
               !!
               !!
            END IF


         END DO ! jk


         IF ( lmout ) THEN
            WHERE ( mask_u == 0 ) U_c = rmaskvalue
            WHERE ( mask_v == 0 ) V_c = rmaskvalue
         ELSE
            rmaskvalue = 0.
         END IF

         IF ( i3d == 1 ) THEN

            CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vdepth, vtime, U_c(:,:,:),  &
               &    cufilout, 'nav_lon_u', 'nav_lat_u', cv_depth, cv_t_out, cv_rot_U, cu_out,       &
               &    cln_out, rmaskvalue, cun_t=cu_t, lpack=lpcknc4)

            CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vdepth, vtime, V_c(:,:,:),  &
               &    cvfilout, 'nav_lon_v', 'nav_lat_v', cv_depth, cv_t_out, cv_rot_V, cu_out,       &
               &    cln_out, rmaskvalue, cun_t=cu_t, lpack=lpcknc4)
         ELSE

            !! Writing file for corrected U and V :
            CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u,         vtime, U_c(:,:,1),     &
               &    cufilout, 'nav_lon_u', 'nav_lat_u', cv_t_out, cv_rot_U, cu_out,       &
               &    cln_out, rmaskvalue, cun_t=cu_t, lpack=lpcknc4)

            CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v,         vtime, V_c(:,:,1), &
               &    cvfilout, 'nav_lon_v', 'nav_lat_v', cv_t_out, cv_rot_V, cu_out,   &
               &    cln_out, rmaskvalue, cun_t=cu_t, lpack=lpcknc4)

         END IF


      END DO ! jt











































   ELSE
      !! !!     I N V E R S E   C O R R E C T I O N

      !! nc4 lolo


      PRINT *, 'Will unrotate vector fields given on an irregular grid!'
      !!
      IF ( lregin ) THEN
         PRINT *, 'Reverse vector correction only makes sense if your source grid is distorded!'
         PRINT *, '  => check "lregin" into the namelist...' ; PRINT *, ''; STOP
      END IF
      !!
      !!
      !!
      cf_out_U = trim(cd_out)//'/'//trim(cv_out_U)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.'//trim(cfext)
      cf_out_V = trim(cd_out)//'/'//trim(cv_out_V)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.'//trim(cfext)
      !!
      PRINT *, 'The two unrotated raw files will be produced :'
      PRINT *, trim(cf_out_U) ;   PRINT *, trim(cf_out_V) ;
      !!
      PRINT *, '' ;   PRINT *, 'Input files :'
      PRINT *, cufilin ;   PRINT *, cvfilin ; PRINT *, '' ; PRINT *, ''
      !!
      PRINT *, 'File containing input grid :'
      PRINT *, cf_mm ; PRINT *, ''


      !! Creating name for unrotated output file:
      nbc = LEN_TRIM(cufilin)
      cdum = cufilin(nbc-2:nbc)
      IF ( cdum /= '.nc' ) THEN
         IF ( cdum == 'nc4' ) THEN
            cfext = 'nc4'
            nlext = 4
            lpcknc4 = .TRUE.
         ELSE
            PRINT *, 'Unknow file extension for ',TRIM(cufilin) ; STOP
         END IF
      END IF

      cf_out_U = cufilin(1:nbc-nlext)
      cf_out_U = trim(cf_out_U)//'_unrotated.'//trim(cfext)

      nbc = LEN_TRIM(cvfilin)
      cf_out_V = cvfilin(1:nbc-nlext)
      cf_out_V = TRIM(cf_out_V)//'_unrotated.'//trim(cfext)

      cv_out_U = trim(cv_rot_U)//'_unrotated'
      cv_out_V = trim(cv_rot_V)//'_unrotated'


      !! Geting array dimension and testing...
      !! -------------------------------------

      CALL DIMS(cufilin, cv_rot_U, ni1, nj1, nk1, Ntr1)
      CALL DIMS(cvfilin, cv_rot_V, ni2, nj2, nk2, Ntr2)

      CALL DIMS(cf_mm,   cv_lon_t, ni_g, nj_g, nk_g, Ntr)

      !! testing ni agreement :
      IF ( (ni1 /= ni2).or.(ni1 /= ni_g).or.(ni2 /= ni_g) ) THEN
         PRINT *, 'Dimension Error! : the 3 files dont agree for x length.'
         STOP
      END IF

      !! testing nj agreement :
      IF ( (nj1 /= nj2).or.(nj1 /= nj_g).or.(nj2 /= nj_g) ) THEN
         PRINT *, 'Dimension Error! : the 3 files dont agree for y length.'
         STOP
      END IF


      nk = 1 ; i3d = 0

      !! testing 3D and nk agreement :

      IF ( (nk1 /= nk2) ) THEN
         PRINT *, 'Dimension Error! : the 2 files dont agree for number of levels.'; STOP
      END IF


      IF ( nk1 /= -1 ) THEN
         l_3d_inv = .TRUE.
         nk = nk1
         i3d = 1
         PRINT *, ''; PRINT *, 'Will perform 3D un-rotating!!!' ; PRINT *, ''
      END IF


      !! testing nt agreement :
      IF ( Ntr1 /= Ntr2 ) THEN
         PRINT *, 'Dimension Error! : u and v files dont agree for time length.'
         STOP
      END IF

      ni = ni1 ; nj = nj1 ; Ntr = Ntr1

      WRITE(*,'("Dimension is : ",i4," x",i4)') ni, nj
      PRINT *, 'Number of levels to treat =>', nk
      WRITE(*,'(i4," time records for u and v")') Ntr
      PRINT *, ''



      !! Allocations :
      !! -------------
      ALLOCATE (XCOST(ni,nj) , XSINT(ni,nj) , U_r(ni,nj) , V_r(ni,nj) ,  &
         &     U_c(ni,nj,nk) , V_c(ni,nj,nk),                              &
         &     xlon_t(ni,nj)  , xlat_t(ni,nj)  , xlon_u(ni,nj) ,  xlat_u(ni,nj), xlon_v(ni,nj) ,  xlat_v(ni,nj) , &
         &     vtime(Ntr) , mask_u(ni,nj,nk) , mask_v(ni,nj,nk)  )
      !!
      ALLOCATE (XCOST8(ni,nj) , XSINT8(ni,nj) , U_r8(ni,nj) , V_r8(ni,nj) )

      XCOST = 0. ; XSINT = 0. ; U_r    = 0. ; V_r = 0. ; U_c = 0. ; V_c = 0.
      xlon_t  = 0. ; xlat_t  = 0. ; xlon_u  = 0. ; xlat_u  = 0. ; xlon_v  = 0. ; xlat_v  = 0.
      vtime  = 0.


      !!  Getting longitude and latitude form grid file :
      !! ------------------------------------------------
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lon_t, 1, 1, 1, xlon_t)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lat_t, 1, 1, 1, xlat_t)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lon_u, 1, 1, 1, xlon_u)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lat_u, 1, 1, 1, xlat_u)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lon_v, 1, 1, 1, xlon_v)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_lat_v, 1, 1, 1, xlat_v)  ; i0=0 ; j0=0

      IF ( l_3d_inv ) THEN
         ALLOCATE ( vdepth(nk) )
         CALL GETVAR_1D(cf_mm, cv_depth, vdepth)
         CALL GETMASK_3D(cf_mm, 'umask', mask_u(:,:,:))
         CALL GETMASK_3D(cf_mm, 'vmask', mask_v(:,:,:))
      ELSE
         CALL GETMASK_2D(cf_mm, 'umask', mask_u(:,:,1))
         CALL GETMASK_2D(cf_mm, 'vmask', mask_v(:,:,1))
      END IF

      !!  Getting cos and sin of the grid distorsion angle:
      !! --------------------------------------------------
      CALL angle_dist(xlon_t, xlat_t, xlon_u, xlat_u, XCOST8, XSINT8)


      !!  Getting time from the u_raw file or the namelist :
      !! ---------------------------------------------------
      IF ( lct ) THEN       ! time is being controlled
         DO jt = 1, Ntr
            vtime(jt) = t0 + t_stp*REAL(jt)
         END DO
      ELSE                  ! we use time from input file
         CALL GETVAR_1D(cufilin, cv_time_0, vtime)
      END IF




      DO jt = 1, Ntr

         PRINT *, ''
         PRINT *, ''
         PRINT *, ''; PRINT *, 'Time step =', jt ; PRINT *, ''

         DO jk = 1, nk



            PRINT *, '  *** Un-rotating level =', jk

            !! Getting U :
            !! -----------
            CALL  GETVAR_2D(idf_u, idv_u, cufilin, cv_rot_U, Ntr, jk*i3d, jt, U_r, lz=nk)

            !! Getting V :
            !! -----------
            CALL  GETVAR_2D(idf_v, idv_v, cvfilin, cv_rot_V, Ntr, jk*i3d, jt, V_r, lz=nk)

            U_r8 = U_r ; V_r8 = V_r


            !! Unrotating U :
            !! --------------
            U_c(:,:,jk) = REAL(XCOST8*U_r8 - XSINT8*V_r8 , 4) ! note the '-' sign --> reverse correction



            !! Unrotating V :
            !! --------------
            V_c(:,:,jk) = REAL(XCOST8*V_r8 + XSINT8*U_r8 , 4) ! note the + sign for reverse correction



         END DO  ! jk


         !lulu
         WHERE ( mask_u == 0 ) U_c(:,:,1:nk) = zrmv
         WHERE ( mask_v == 0 ) V_c(:,:,1:nk) = zrmv

         IF ( l_3d_inv ) THEN

            CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vdepth, vtime, U_c(:,:,:),  &
               &    cf_out_U, 'nav_lon_u', 'nav_lat_u', cv_depth, cv_time_0, cv_out_U, cu_out,       &
               &    cln_out, zrmv, cun_t=cu_t, lpack=lpcknc4, cun_z='m')

            CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vdepth, vtime, V_c(:,:,:),  &
               &    cf_out_V, 'nav_lon_v', 'nav_lat_v', cv_depth, cv_time_0, cv_out_V, cu_out,       &
               &    cln_out, zrmv, cun_t=cu_t, lpack=lpcknc4, cun_z='m')

         ELSE

            CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vtime, U_c(:,:,1),     &
               &    cf_out_U, 'nav_lon_u', 'nav_lat_u', cv_time_0, cv_out_U, cu_out,       &
               &    cln_out, zrmv, cun_t=cu_t, lpack=lpcknc4)

            CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vtime, V_c(:,:,1), &
               &    cf_out_V, 'nav_lon_v', 'nav_lat_v', cv_time_0, cv_out_V, cu_out,   &
               &    cln_out, zrmv, cun_t=cu_t, lpack=lpcknc4)

         END IF


      END DO

      PRINT *, ''
      PRINT *, 'Files created:'; PRINT *, '   ', trim(cf_out_U); PRINT *, '   ', trim(cf_out_V)
      PRINT *, ''

   END IF ! on l_inv


   PRINT *, 'Done!'

END PROGRAM CORR_VECT








SUBROUTINE angle_dist(xlon_t, xlat_t, xlon_u, xlat_u, xcos_t, xsin_t)

   !!========================================================================
   !!
   !!                     *******  ANGLE_DIST  *******
   !!
   !! Given the coordinates at T-points and U-points, this routine computes
   !! sinus and cosinus of the angle of the local grid distorsion at T-points
   !!
   !!
   !! INPUT :     - xlon_t = longitude array at T-points     [REAL]
   !!             - xlat_t = latitude array at T-points      [REAL]
   !!             - xlon_u = longitude array at U-points     [REAL]
   !!             - xlat_t = latitude array at U-points      [REAL]
   !!
   !! OUTPUT :
   !! --------    - xcos_t  = cosinus of the distortion angle [REAL]
   !!             - xsin_t  = sininus of the distortion angle [REAL]
   !!
   !!
   !! Author : Laurent Brodeau
   !!          directly inspired from 'geo2ocean.F90' (NEMOGCM)
   !!
   !!========================================================================

   IMPLICIT NONE

   !! INPUT :
   !! -------
   REAL(8), DIMENSION(:,:), INTENT(in) ::    &
      &       xlon_t, xlon_u,   &  ! latitudes  at point T and U
      &       xlat_t, xlat_u       ! longitudes at point T and U

   !! OUTPUT :
   !! --------
   REAL(8), DIMENSION(:,:), INTENT(out) :: xcos_t, xsin_t

   !! LOCAL :
   !! -------
   INTEGER :: nx, ny, ji

   REAL(8), DIMENSION(:), ALLOCATABLE :: &
      &     zlon, zlat, zxnpt, znnpt, zlan, &
      &     zphh, zxuut, zyuut, zmnput, zynpt, &
      &     ztmp1, ztmp2

   REAL(8), PARAMETER :: &
      &       Pi  = 3.141592653,     &
      &       rad = Pi/180.0

   !! 2D domain shape:
   nx = SIZE(xlon_t,1)
   ny = SIZE(xlon_t,2)

   ALLOCATE ( zlon(ny), zlat(ny), zxnpt(ny), znnpt(ny), zlan(ny), zphh(ny), &
      &       zxuut(ny), zyuut(ny), zmnput(ny), zynpt(ny), ztmp1(ny), ztmp2(ny) )


   DO ji = 1, nx

      !! North pole direction & modulous (at T-point) :
      zlon  = xlon_t(ji,:)
      zlat  = xlat_t(ji,:)

      ztmp1 = TAN(Pi/4. - rad*zlat/2.)
      zxnpt = 0. - 2.*COS( rad*zlon )*ztmp1
      zynpt = 0. - 2.*SIN( rad*zlon )*ztmp1
      znnpt = zxnpt*zxnpt + zynpt*zynpt

      !! "i" direction & modulous (at T-point) :
      !! ---------------------------------------
      !!   ( since we deal with T points we look at U point
      !!     on a V point we would look at the F point )

      zlon = xlon_u(ji,:)
      zlat = xlat_u(ji,:)

      IF ( ji == 1 ) THEN
         zlan = xlon_u(nx-3,:)          !! periodicity of ORCA grid
         zphh = xlat_u(nx-3,:)          !! with overlap of 2 points
      ELSE
         zlan = xlon_u(ji-1,:)
         zphh = xlat_u(ji-1,:)
      END IF

      ztmp1 = TAN(Pi/4. - rad*zlat/2.)
      ztmp2 = TAN(Pi/4. - rad*zphh/2.)
      zxuut  = 2.*COS(rad*zlon)*ztmp1 - 2.*COS(rad*zlan)*ztmp2
      zyuut  = 2.*SIN(rad*zlon)*ztmp1 - 2.*SIN(rad*zlan)*ztmp2
      zmnput = SQRT(znnpt*(zxuut*zxuut + zyuut*zyuut))



      !! Cosinus and sinus using scalar and vectorial products :
      !! -------------------------------------------------------
      !!   (caution, rotation of 90 degres)

      !WHERE ( zmnput < 1.e-14 ) zmnput = 1.e-14
      zmnput = MAX(zmnput, 1.e-14)
      ztmp1 = 1./zmnput
      xsin_t(ji,:) =  ( zxnpt*zxuut + zynpt*zyuut ) * ztmp1
      xcos_t(ji,:) = -( zxnpt*zyuut - zynpt*zxuut ) * ztmp1

   END DO


   DEALLOCATE ( zlon, zlat, zxnpt, znnpt, zlan, zphh, &
      &         zxuut, zyuut, zmnput, zynpt, ztmp1, ztmp2 )


END SUBROUTINE angle_dist






SUBROUTINE usage_corr_vect()
   !!
   PRINT *,''
   PRINT *,'   List of command line options:'
   PRINT *,''
   PRINT *,' -I   => will perform inverse correction, ie un-rotate a ORCA grid'
   PRINT *,''
   PRINT *,''
   PRINT *,'  *** MANDATORY for both normal and inverse mode:'
   PRINT *,''
   PRINT *,' -x  <name U>         => Specify name for x comp. in output file'
   PRINT *,'                         (or input file if inverse mode)'
   PRINT *,''
   PRINT *,' -y  <name V>         => Specify name for y comp. in output file'
   PRINT *,'                         (or input file if inverse mode)'
   PRINT *,''
   PRINT *,' -m  <mesh_mask_file> => Specify which mesh_mask file to use'
   PRINT *,''
   PRINT *,''
   PRINT *,''
   PRINT *,'  ***  MANDATORY for normal mode (no -I switch) :'
   PRINT *,''
   PRINT *,' -f  <namelist_file>  => Specify which namelist file to use'
   PRINT *, '                        No namelist needed when inverse correction'
   PRINT *,''
   PRINT *,''
   PRINT *,'  *** MANDATORY for INVERSE MODE (-I switch) :'
   PRINT *,''
   PRINT *,' -i <x.nc> <y.nc>     =>  unrotate vector fields given in these 2 files'
   PRINT *,'                          to the same grid'
   PRINT *,''
   PRINT *,' -t <time_name>       =>  name of time variable in <x.nc> and <y.nc>'
   PRINT *,''
   PRINT *,''
   PRINT *,'  *** MISC options :'
   PRINT *,''
   !PRINT *,' -p                   => Pack data (short 2-byte) in output netcdf file'
   !PRINT *,'                         (default is float 4-byte)'
   !PRINT *,''
   PRINT *,' -h                   => Show this message'
   PRINT *,''
   STOP
   !!
END SUBROUTINE usage_corr_vect
!!
