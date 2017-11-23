PROGRAM CORR_VECT

   USE mod_conf
   USE mod_init
   USE io_ezcdf
   USE mod_grids,  ONLY: IS_ORCA_NORTH_FOLD
   USE mod_nemotools

   !!========================================================================
   !! Purpose :  correct vector components 'uraw' and 'vraw' directly
   !!            interpolated on an irregular grid
   !! ---------  Only for ORCA family of grids!!!!
   !!
   !! Author :   Laurent Brodeau
   !! --------
   !!
   !!========================================================================

   IMPLICIT NONE

   !! Grid :
   CHARACTER(len=80), PARAMETER   :: &
      &    cv_glamt     = 'glamt',   &   ! input grid longitude name, T-points
      &    cv_gphit     = 'gphit',   &   ! input grid latitude name,  T-points
      &    cv_glamu     = 'glamu',   &   ! input grid longitude name, U-points
      &    cv_gphiu     = 'gphiu',   &   ! input grid latitude name,  U-points
      &    cv_glamv     = 'glamv',   &   ! input grid longitude name, V-points
      &    cv_gphiv     = 'gphiv',   &   ! input grid latitude name,  V-points
      &    cv_glamf     = 'glamf',   &   ! input grid longitude name, F-points
      &    cv_gphif     = 'gphif',   &   ! input grid latitude name,  F-points
      &    cv_depth     = 'deptht'       !  depth at T-points (U-points and V-points too)


   CHARACTER(LEN=400)  :: cextra_x, cextra_y 
   
   CHARACTER(len=3)    :: cdum
   CHARACTER(len=1)    :: cgrid_out='0'
   CHARACTER(len=80)   :: cv_time_0 = 'none', cfext = 'nc'
   CHARACTER(len=800)  :: cr, cf_mm, cnmlst_x, cnmlst_y

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
      &    iorca=0,      &
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

   INTEGER(2), DIMENSION(:,:,:), ALLOCATABLE :: mask_t, mask_u, mask_v

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: Xdum4


   REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: Ut_c, Vt_c, Uu_c, Vv_c

   REAL(8), DIMENSION(:,:), ALLOCATABLE ::      &
      &    XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8, &
      &    U_r8, V_r8,  Xdum8, &
      &    xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f

   REAL(8), DIMENSION(:), ALLOCATABLE ::   vtime, vdepth


   INTEGER :: jt, jk

   LOGICAL :: &
      &    l_inv = .FALSE.,   &
      &    l_anlt = .FALSE.,  & ! analytic input field (U=1, V=0) DEBUG...
      &    l_3d_inv = .FALSE., &   !: will treat 3d files in inverse mode...
      &    lmout_x, lmout_y, &
      &    lexist !,  &
   
   REAL(4), PARAMETER :: zrmv = -9999.

   REAL(8) :: rsgn
   
   CHARACTER(LEN=2), DIMENSION(9), PARAMETER :: &
      &            clist_opt = (/ '-I','-h','-m','-G','-p','-f','-i','-t','-1' /)
      !&            clist_opt = (/ '-I','-h','-m','-G','-p','-x','-y','-f','-i','-t','-1' /)

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
      !CASE('-x')
      !   IF ( jarg + 1 > iargc() ) THEN
      !      PRINT *, 'ERROR: Missing zonal component name!' ; call usage_corr_vect()
      !   ELSE
      !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
      !      IF ( ANY(clist_opt == trim(cr)) ) THEN
      !         PRINT *, 'ERROR: Missing zonal component name!'; call usage_corr_vect()
      !      ELSE
      !         cv_rot_U = trim(cr)
      !      END IF
      !   END IF
         !!
         !!
      !CASE('-y')
      !   IF ( jarg + 1 > iargc() ) THEN
      !      PRINT *, 'ERROR: Missing meridional component name!' ; call usage_corr_vect()
      !   ELSE
      !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
      !      IF ( ANY(clist_opt == trim(cr)) ) THEN
      !         PRINT *, 'ERROR: Missing meridional component name!'; call usage_corr_vect()
      !      ELSE
      !         cv_rot_V = trim(cr)
      !      END IF
      !   END IF
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
      CASE('-G')
         IF ( jarg + 1 > iargc() ) THEN
            PRINT *, 'ERROR: Missing grid type to write to ("T" or "UV"?) !!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == TRIM(cr)) ) THEN
               PRINT *, 'ERROR: Missing grid type to write to ("T" or "UV"?) !!'; CALL usage_corr_vect()
            ELSE
               cgrid_out = TRIM(cr)
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
      CASE('-1')
         l_anlt = .TRUE.
         !!
      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
         CALL usage_corr_vect()
         !!
      END SELECT
      !!
      !!
   END DO


   IF ( (.NOT. l_inv).AND.( TRIM(cufilin) /= 'none' ) ) THEN
      PRINT *, 'ERROR: specify the "-I" switch if you want to perform inverse correction!'
      STOP
   END IF

   IF ( (cgrid_out /= 'T').AND.(cgrid_out /= 'U') ) call usage_corr_vect()

   PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   PRINT *, ''


   IF ( cgrid_out == 'U' ) PRINT *, ' *** Gonna save on grid U and V.'
   IF ( cgrid_out == 'T' ) PRINT *, ' *** Gonna save on grid T.'
   PRINT *, ''


   IF ( l_inv ) THEN
      PRINT *, ' * Vector files to unrotate = ', trim(cufilin), ' , ', trim(cvfilin)
      !PRINT *, '   => associated variable names = ', TRIM(cv_rot_U), ' , ', TRIM(cv_rot_V)
      IF ( trim(cv_time_0) == 'none' ) THEN
         PRINT *, 'ERROR: you must specify the name of time variable with the "-t" switch!'; STOP
      END IF
      PRINT *, '   => time variable name = ', trim(cv_time_0)
   !ELSE
      !PRINT *, ' * Name for corrected vector components = ', TRIM(cv_rot_U), ' , ', TRIM(cv_rot_V)
   END IF

   !PRINT *, ' * Packing output : ', lpack
   PRINT *, ' * mesh_mask file to use = ', TRIM(cf_mm)

   cnmlst_x = TRIM(cf_nml_sosie)//'_x'
   cnmlst_y = TRIM(cf_nml_sosie)//'_y'
   
   PRINT *, ' * namelists we expect => ', TRIM(cnmlst_x)//' and '//TRIM(cnmlst_y)
   PRINT *, ''


   INQUIRE(FILE=trim(cf_mm), EXIST=lexist )
   IF ( .NOT. lexist ) THEN
      WRITE(*,'("The mesh_mask file ",a," file was not found!")') trim(cf_mm)
      CALL usage_corr_vect()
   END IF



   IF ( .NOT. l_inv ) THEN

      !! !!     N O R M A L   C O R R E C T I O N

      !! Namelist of X component:
      INQUIRE(FILE=TRIM(cnmlst_x), EXIST=lexist )
      IF ( .NOT. lexist ) THEN
         WRITE(*,'("The namelist file ",a," file was not found!")') TRIM(cnmlst_x)
         CALL usage_corr_vect()
      END IF
      PRINT *, ''
      cf_nml_sosie = TRIM(cnmlst_x)
      CALL READ_NMLST(2)
      lmout_x  = lmout
      cv_rot_U = cv_out
      cextra_x = cextra

      
      !! Namelist of Y component:
      INQUIRE(FILE=TRIM(cnmlst_y), EXIST=lexist )
      IF ( .NOT. lexist ) THEN
         WRITE(*,'("The namelist file ",a," file was not found!")') TRIM(cnmlst_y)
         CALL usage_corr_vect()
      END IF
      PRINT *, ''
      cf_nml_sosie = TRIM(cnmlst_y)
      CALL READ_NMLST(2)
      lmout_y  = lmout
      cv_rot_V = cv_out
      cextra_y = cextra


      !lolo

      

      IF ( lregout ) THEN
         PRINT *, 'Vector correction only makes sense if your target grid is distorded!'
         PRINT *, '  => check "lregout" into the namelist...' ; PRINT *, ''; STOP
      END IF


      IF ( lpcknc4 ) cfext='nc4'

      cf_out_U = trim(cd_out)//'/'//trim(cv_out_U)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra_x)//'.'//trim(cfext)
      cf_out_V = trim(cd_out)//'/'//trim(cv_out_V)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra_y)//'.'//trim(cfext)

      PRINT *, 'The two input pre-interpolated needed files are :'
      PRINT *, trim(cf_out_U) ;   PRINT *, trim(cf_out_V) ;
      
      cufilout = TRIM(cd_out)//'/'//TRIM(cv_rot_U)//'_'//TRIM(csource)//'-' &
         &   //TRIM(ctarget)//'_'//TRIM(cextra_x)//'.'//TRIM(cfext)
      cvfilout = TRIM(cd_out)//'/'//TRIM(cv_rot_V)//'_'//TRIM(csource)//'-' &
         &   //TRIM(ctarget)//'_'//TRIM(cextra_y)//'.'//TRIM(cfext)
      
      PRINT *, '' ;   PRINT *, 'output files :'
      PRINT *, trim(cufilout) ;   PRINT *, trim(cvfilout) ; PRINT *, '' ; PRINT *, ''

      PRINT *, 'File containing x and y raw components of vector to be treated :'
      PRINT *, trim(cf_out_U)
      PRINT *, trim(cf_out_V) ; PRINT *, ''
      PRINT *, 'File containing grid :'
      PRINT *, trim(cf_mm) ; PRINT *, ''


      !! Geting array dimension and testing...
      !! -------------------------------------

      CALL DIMS(cf_out_U, cv_out_U, ni1, nj1, nk1, Ntr1)
      CALL DIMS(cf_out_V, cv_out_V, ni2, nj2, nk2, Ntr2)

      CALL DIMS(cf_mm, cv_glamt, ni_g, nj_g, nk_g, Ntr)

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
      ALLOCATE ( Xdum4(ni,nj) ,  &
         &    Ut_c(ni,nj,nk) , Vt_c(ni,nj,nk), mask_t(ni,nj,nk), &
         &    xlon_t(ni,nj)  , xlat_t(ni,nj), Xdum8(ni,nj), &
         &    xlon_u(ni,nj)  , xlat_u(ni,nj) , xlon_v(ni,nj)  , xlat_v(ni,nj) , &
         &    xlon_f(ni,nj) , xlat_f(ni,nj), &
         &    vtime(Ntr)   )
      IF ( cgrid_out == 'U' ) ALLOCATE ( Uu_c(ni,nj,nk) , Vv_c(ni,nj,nk) , mask_u(ni,nj,nk) , mask_v(ni,nj,nk)  )

      ALLOCATE (XCOST8(ni,nj) , XSINT8(ni,nj) , XCOSU8(ni,nj) , XSINU8(ni,nj) , XCOSV8(ni,nj) , XSINV8(ni,nj) , &
         &      XCOSF8(ni,nj) , XSINF8(ni,nj) , U_r8(ni,nj) , V_r8(ni,nj) )


      IF ( i3d == 1 ) THEN
         ALLOCATE ( vdepth(nk) )
         CALL GETVAR_1D(cf_mm, cv_depth, vdepth)
         IF ( cgrid_out == 'U' ) THEN
            CALL GETMASK_3D(cf_lsm_out, 'umask', mask_u(:,:,:))
            CALL GETMASK_3D(cf_lsm_out, 'vmask', mask_v(:,:,:))
         ELSE
            CALL GETMASK_3D(cf_lsm_out, 'tmask', mask_t(:,:,:))
         END IF
      ELSE
         IF ( cgrid_out == 'U' ) THEN
            CALL GETMASK_2D(cf_lsm_out, 'umask', mask_u(:,:,1), jlev=1)
            CALL GETMASK_2D(cf_lsm_out, 'vmask', mask_v(:,:,1), jlev=1)
         ELSE
            CALL GETMASK_2D(cf_lsm_out, 'tmask', mask_t(:,:,1), jlev=1)
         END IF
      END IF

      !!  Getting longitude and latitude form grid file :
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamt, 1, 1, 1, xlon_t)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphit, 1, 1, 1, xlat_t)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamu, 1, 1, 1, xlon_u)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphiu, 1, 1, 1, xlat_u)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamv, 1, 1, 1, xlon_v)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphiv, 1, 1, 1, xlat_v)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamf, 1, 1, 1, xlon_f)
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphif, 1, 1, 1, xlat_f)


      PRINT *, ''

      !! Is this a known ORCA grid (just for info now, not used!):
      iorca = IS_ORCA_NORTH_FOLD( xlat_t )
      IF ( iorca == 4 ) PRINT *, ' Grid is an ORCA grid with north-pole T-point folding!'
      IF ( iorca == 6 ) PRINT *, ' Grid is an ORCA grid with north-pole F-point folding!'
      PRINT *, ''


      !!  Getting cosine and sine corresponding to the angle of the local distorsion of the grid:
      CALL ANGLE( iorca, xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f, &
         &        XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8 )

      !CALL PRTMASK(REAL(XCOST8,4), 'cost_angle.nc', 'cos',   xlon_t, xlat_t, cv_glamt, cv_gphit)
      !CALL PRTMASK(REAL(XSINT8,4), 'sint_angle.nc', 'sin',   xlon_t, xlat_t, cv_glamt, cv_gphit)
      !CALL PRTMASK(REAL(XCOSU8,4), 'cosu_angle.nc', 'cos',   xlon_t, xlat_t, cv_glamt, cv_gphit)
      !CALL PRTMASK(REAL(XSINU8,4), 'sinu_angle.nc', 'sin',   xlon_t, xlat_t, cv_glamt, cv_gphit)
      !STOP

      !!  Getting time from the u_raw file or the namelist :
      IF ( lct ) THEN       ! time is being controlled
         DO jt = 1, Ntr
            vtime(jt) = t0 + t_stp*REAL(jt)
         END DO
         nb_att_t = 1
         vatt_info_t(:)%cname = 'null'
         vatt_info_t(1)%cname = 'units'
         vatt_info_t(1)%itype = 2 ! char
         vatt_info_t(1)%val_char = 'unknown'
         vatt_info_t(1)%ilength = LEN('unknown')
      ELSE                  ! we use time from input file
         CALL GETVAR_1D(cf_out_U, cv_t_out, vtime)
         CALL GETVAR_ATTRIBUTES(cf_out_U, cv_t_out, nb_att_t, vatt_info_t)
      END IF



      DO jt = 1, Ntr

         PRINT *, ''; PRINT *, 'Time step =', jt ; PRINT *, ''

         DO jk = 1, nk

            IF ( l_anlt ) THEN
               !lolo: analytical field for debugging purposes...
               U_r8 = 1.
               V_r8 = 0.
            ELSE
               !! Getting uncorrected U on grid T:
               CALL  GETVAR_2D(idf_u, idv_u, cf_out_U, cv_out_U, Ntr, jk*i3d, jt, Xdum4, lz=nk)
               U_r8 = Xdum4

               !! Getting uncorrected V on grid T:
               CALL  GETVAR_2D(idf_v, idv_v, cf_out_V, cv_out_V, Ntr, jk*i3d, jt, Xdum4, lz=nk)
               V_r8 = Xdum4
            END IF

            
            IF ( iorca > 0 ) PRINT *, '   *** goona "lbc_lnk_2d" with iorca =', iorca
            
            rsgn = -1.
            IF ( cgrid_out == 'U' ) THEN
               !! U-V grid:
               !! Correcting U :
               Xdum8 = XCOSU8*U_r8 + XSINU8*V_r8  !lolo: no trcik here???
               IF ( iorca > 0 ) CALL lbc_lnk_2d( iorca, Xdum8, 'U', rsgn )              
               Uu_c(:,:,jk) = REAL(Xdum8 , 4)
               !!
               !! Correcting V :
               Xdum8 = XCOSV8*V_r8 - XSINV8*U_r8  !lolo: no trcik here???
               IF ( iorca > 0 ) CALL lbc_lnk_2d( iorca, Xdum8, 'V', rsgn )               
               Vv_c(:,:,jk) = REAL(Xdum8 , 4)
               !!
            ELSE
               !! T grid:
               !! Correcting U (i-component to east) :
               Xdum8 = XCOST8*U_r8 + XSINT8*V_r8
               IF ( iorca > 0 ) CALL lbc_lnk_2d( iorca, Xdum8, 'T', rsgn )               
               Ut_c(:,:,jk) = REAL(Xdum8 , 4)
               !!
               !! Correcting V (j-component to north):
               Xdum8 = XCOST8*V_r8 - XSINT8*U_r8
               IF ( iorca > 0 ) CALL lbc_lnk_2d( iorca, Xdum8, 'T', rsgn )               
               Vt_c(:,:,jk) = REAL(Xdum8 , 4)
               !!
            END IF

         END DO ! jk

         
         IF ( lmout_x .AND. lmout_y ) THEN
            IF ( cgrid_out == 'U' ) THEN
               WHERE ( mask_u == 0 ) Uu_c = rmaskvalue
               WHERE ( mask_v == 0 ) Vv_c = rmaskvalue
            ELSE
               WHERE ( mask_t == 0 ) Ut_c = rmaskvalue
               WHERE ( mask_t == 0 ) Vt_c = rmaskvalue
            END IF
         ELSE
            rmaskvalue = 0.
         END IF


         IF ( i3d == 1 ) THEN

            !! 3D:
            IF ( cgrid_out == 'U' ) THEN
               CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vdepth, vtime, Uu_c(:,:,:),  &
                  &    cufilout, 'nav_lon_u', 'nav_lat_u', cv_depth, cv_t_out, cv_rot_U,      &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
               CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vdepth, vtime, Vv_c(:,:,:),  &
                  &    cvfilout, 'nav_lon_v', 'nav_lat_v', cv_depth, cv_t_out, cv_rot_V,      &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
            ELSE
               CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_t, xlat_t, vdepth, vtime, Ut_c(:,:,:),  &
                  &    cufilout, 'nav_lon', 'nav_lat', cv_depth, cv_t_out, cv_rot_U,      &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
               CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_t, xlat_t, vdepth, vtime, Vt_c(:,:,:),  &
                  &    cvfilout, 'nav_lon', 'nav_lat', cv_depth, cv_t_out, cv_rot_V,      &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
            END IF

         ELSE

            !! 2D:
            IF ( cgrid_out == 'U' ) THEN
               CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u,         vtime, Uu_c(:,:,1), &
                  &    cufilout, 'nav_lon_u', 'nav_lat_u', cv_t_out, cv_rot_U,       &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
               CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v,         vtime, Vv_c(:,:,1), &
                  &    cvfilout, 'nav_lon_v', 'nav_lat_v', cv_t_out, cv_rot_V,   &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
            ELSE
               CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_t, xlat_t,         vtime, Ut_c(:,:,1), &
                  &    cufilout, 'nav_lon', 'nav_lat', cv_t_out, cv_rot_U,       &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
               CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_t, xlat_t,         vtime, Vt_c(:,:,1), &
                  &    cvfilout, 'nav_lon', 'nav_lat', cv_t_out, cv_rot_V,   &
                  &    rmaskvalue, attr_time=vatt_info_t, lpack=lpcknc4)
            END IF

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
      PRINT *, trim(cufilin) ;   PRINT *, trim(cvfilin) ; PRINT *, '' ; PRINT *, ''
      !!
      PRINT *, 'File containing input grid :'
      PRINT *, trim(cf_mm) ; PRINT *, ''


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

      CALL DIMS(cf_mm,   cv_glamt, ni_g, nj_g, nk_g, Ntr)

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
      ALLOCATE (Xdum4(ni,nj) , &
         &     Ut_c(ni,nj,nk) , Vt_c(ni,nj,nk),                              &
         &     xlon_t(ni,nj) , xlat_t(ni,nj) , xlon_u(ni,nj) , xlat_u(ni,nj), xlon_v(ni,nj) , xlat_v(ni,nj) , &
         &     xlon_f(ni,nj) , xlat_f(ni,nj) , vtime(Ntr) , mask_u(ni,nj,nk) , mask_v(ni,nj,nk)  )
      !!
      ALLOCATE (XCOST8(ni,nj) , XSINT8(ni,nj) , XCOSU8(ni,nj) , XSINU8(ni,nj) , XCOSV8(ni,nj) , XSINV8(ni,nj) , &
         &      XCOSF8(ni,nj) , XSINF8(ni,nj) , U_r8(ni,nj) , V_r8(ni,nj) )


      !!  Getting longitude and latitude form grid file :
      !! ------------------------------------------------
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamt, 1, 1, 1, xlon_t)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphit, 1, 1, 1, xlat_t)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamu, 1, 1, 1, xlon_u)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphiu, 1, 1, 1, xlat_u)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamv, 1, 1, 1, xlon_v)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphiv, 1, 1, 1, xlat_v)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_glamf, 1, 1, 1, xlon_f)  ; i0=0 ; j0=0
      CALL GETVAR_2D_R8(i0, j0, cf_mm, cv_gphif, 1, 1, 1, xlat_f)  ; i0=0 ; j0=0

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
      CALL ANGLE( iorca, xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f, &
         &        XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8 )



      !!  Getting time from the u_raw file or the namelist :
      !! ---------------------------------------------------
      IF ( lct ) THEN       ! time is being controlled
         DO jt = 1, Ntr
            vtime(jt) = t0 + t_stp*REAL(jt)
         END DO
         nb_att_t = 1
         vatt_info_t(:)%cname = 'null'
         vatt_info_t(1)%cname = 'units'
         vatt_info_t(1)%itype = 2 ! char
         vatt_info_t(1)%val_char = 'unknown'
         vatt_info_t(1)%ilength = LEN('unknown')

      ELSE                  ! we use time from input file
         CALL GETVAR_1D(cufilin, cv_time_0, vtime)
         CALL GETVAR_ATTRIBUTES(cufilin, cv_time_0, nb_att_t, vatt_info_t) ; !lolo
      END IF




      DO jt = 1, Ntr

         PRINT *, ''
         PRINT *, ''
         PRINT *, ''; PRINT *, 'Time step =', jt ; PRINT *, ''

         DO jk = 1, nk



            PRINT *, '  *** Un-rotating level =', jk

            !! Getting U :
            !! -----------
            CALL  GETVAR_2D(idf_u, idv_u, cufilin, cv_rot_U, Ntr, jk*i3d, jt, Xdum4, lz=nk)
            U_r8 = Xdum4

            !! Getting V :
            !! -----------
            CALL  GETVAR_2D(idf_v, idv_v, cvfilin, cv_rot_V, Ntr, jk*i3d, jt, Xdum4, lz=nk)
            V_r8 = Xdum4


            !! Unrotating U :
            !! --------------
            Ut_c(:,:,jk) = REAL(XCOST8*U_r8 - XSINT8*V_r8 , 4) ! note the '-' sign --> reverse correction



            !! Unrotating V :
            !! --------------
            Vt_c(:,:,jk) = REAL(XCOST8*V_r8 + XSINT8*U_r8 , 4) ! note the + sign for reverse correction



         END DO  ! jk


         !lulu
         WHERE ( mask_u == 0 ) Ut_c(:,:,1:nk) = zrmv
         WHERE ( mask_v == 0 ) Vt_c(:,:,1:nk) = zrmv

         IF ( l_3d_inv ) THEN

            CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vdepth, vtime, Ut_c(:,:,:),  &
               &    cf_out_U, 'nav_lon_u', 'nav_lat_u', cv_depth, cv_time_0, cv_out_U,       &
               &    zrmv, attr_time=vatt_info_t, lpack=lpcknc4)

            CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vdepth, vtime, Vt_c(:,:,:),  &
               &    cf_out_V, 'nav_lon_v', 'nav_lat_v', cv_depth, cv_time_0, cv_out_V,       &
               &    zrmv, attr_time=vatt_info_t, lpack=lpcknc4)

         ELSE

            CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vtime, Ut_c(:,:,1),     &
               &    cf_out_U, 'nav_lon_u', 'nav_lat_u', cv_time_0, cv_out_U,       &
               &    zrmv, attr_time=vatt_info_t, lpack=lpcknc4)

            CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vtime, Vt_c(:,:,1), &
               &    cf_out_V, 'nav_lon_v', 'nav_lat_v', cv_time_0, cv_out_V,   &
               &    zrmv, attr_time=vatt_info_t, lpack=lpcknc4)

         END IF


      END DO

      PRINT *, ''
      PRINT *, 'Files created:'; PRINT *, '   ', trim(cf_out_U); PRINT *, '   ', trim(cf_out_V)
      PRINT *, ''

   END IF ! on l_inv


   PRINT *, 'Done!'

END PROGRAM CORR_VECT



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
   !PRINT *,' -x  <name U>         => Specify name for x comp. in output file'
   !PRINT *,'                         (or input file if inverse mode)'
   !PRINT *,''
   !PRINT *,' -y  <name V>         => Specify name for y comp. in output file'
   !PRINT *,'                         (or input file if inverse mode)'
   !PRINT *,''
   PRINT *,' -m  <mesh_mask_file> => Specify which mesh_mask file to use'
   PRINT *,''
   PRINT *,''
   PRINT *,'  ***  MANDATORY for normal mode (no -I switch) :'
   PRINT *,''
   PRINT *,' -f  <namelist_file>  => Specify which namelist file to use'
   PRINT *, '                        No namelist needed when inverse correction'
   PRINT *,''
   PRINT *,' -G  <T/U>            => Specify if you want to save rotated vector'
   PRINT *, '                        on T-grid (T) or U- and V-grid (U)'
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
   PRINT *,' -h                   => Show this message'
   PRINT *,''
   STOP
   !!
END SUBROUTINE usage_corr_vect
!!
