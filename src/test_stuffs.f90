PROGRAM TEST_STUFFS

   !USE mod_conf
   !USE mod_init
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


   !CHARACTER(LEN=400)  :: cextra_x, cextra_y

   CHARACTER(len=3)    :: cdum
   CHARACTER(len=1)    :: cgrid_trg='0'
   CHARACTER(len=80)   :: cv_time_0 = 'none', cfext = 'nc'
   CHARACTER(len=800)  :: cr, cf_mm, cnmlst_x, cnmlst_y

   !CHARACTER(len=80) :: &
   !   &    cv_trg_U    = 'uraw',   &   ! raw U name
   !   &    cv_trg_V    = 'vraw'        ! raw V name

   !CHARACTER(len=800)  :: &
   !   &    cf_out_U,    &  ! file containing u_raw
   !   &    cf_out_V,    &  ! file containing v_raw
   !   &    cufilout, cvfilout,    &
   !   &    cufilin = 'none',  cvfilin = 'none'

   !CHARACTER(len=80)  ::  &
   !   &    cv_rot_U ,  &  ! output name for U corrected
   !   &    cv_rot_V       ! output name for V corrected

   INTEGER      :: &
      &    Ntr, &
      &    jarg, i3d, nbc, &
      &    iorca=0,      &
      &    nlext=3, &
      &    i0, j0, &
      &    ji, jj, &
      &    ni, nj, nk, nk1, nk2, &
      &    ni1, nj1, Ntr1,       &
      &    ni2, nj2, Ntr2,       &
      &    ni_g, nj_g, nk_g,    &
      &    iargc,         &
      &    idf_u, idv_u, idf_v, idv_v, &
      &    id_f1, id_v1, &
      &    id_f2, id_v2

   INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: mask_t, mask_u, mask_v

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: Xdum4


   REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: Ut_c, Vt_c, Uu_c, Vv_c

   REAL(8), DIMENSION(:,:), ALLOCATABLE ::      &
      &    XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8, &
      &    U_r8, V_r8,  Xdum8, &
      &    xlamt, xphit, xlamu, xphiu, xlamv, xphiv, xlamf, xphif

   REAL(8), DIMENSION(:), ALLOCATABLE ::   vtime, vdepth


   INTEGER :: jt, jk

   LOGICAL :: &
      &    l_anlt = .FALSE.,  & ! analytic input field (U=1, V=0) DEBUG...
      &    l_3d_inv = .FALSE., &   !: will treat 3d files in inverse mode...
      &    lmout_x, lmout_y, &
      &    lexist !,  &

   REAL(4), PARAMETER :: zrmv = -9999.

   REAL(8) :: rsgn

   CHARACTER(LEN=2), DIMENSION(9), PARAMETER :: &
      &            clist_opt = (/ '-I','-h','-m','-G','-p','-f','-i','-t','-1' /)
   !&            clist_opt = (/ '-I','-h','-m','-G','-p','-x','-y','-f','-i','-t','-1' /)

   !! Getting string arguments :
   !! --------------------------

   !! Some defaults:
   !cv_rot_U = 'vectx' ; cv_rot_V = 'vecty'

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
         call usage_test_stuffs()
         !!
         !!
         !!
         !CASE('-x')
         !   IF ( jarg + 1 > iargc() ) THEN
         !      PRINT *, 'ERROR: Missing zonal component name!' ; call usage_test_stuffs()
         !   ELSE
         !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         !      IF ( ANY(clist_opt == trim(cr)) ) THEN
         !         PRINT *, 'ERROR: Missing zonal component name!'; call usage_test_stuffs()
         !      ELSE
         !         cv_rot_U = trim(cr)
         !      END IF
         !   END IF
         !!
         !!
         !CASE('-y')
         !   IF ( jarg + 1 > iargc() ) THEN
         !      PRINT *, 'ERROR: Missing meridional component name!' ; call usage_test_stuffs()
         !   ELSE
         !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         !      IF ( ANY(clist_opt == trim(cr)) ) THEN
         !         PRINT *, 'ERROR: Missing meridional component name!'; call usage_test_stuffs()
         !      ELSE
         !         cv_rot_V = trim(cr)
         !      END IF
         !   END IF
         !!
         !CASE('-f')
         !   IF ( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
         !      PRINT *, 'ERROR: Missing namelist name!' ; call usage_test_stuffs()
         !   ELSE
         !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         !      IF ( ANY(clist_opt == trim(cr)) ) THEN
         !         PRINT *, 'ERROR: ', trim(cr), ' is definitively not the name of the namelist!'
         !         call usage_test_stuffs()
         !      ELSE
         !         cf_nml_sosie = trim(cr)
         !      END IF
         !   END IF
         !!
      CASE('-m')
         IF ( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
            PRINT *, 'ERROR: Missing mesh_mask file name!' ; call usage_test_stuffs()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *, 'ERROR: ', trim(cr), ' is definitively not the name of the mesh_mask file!'
               call usage_test_stuffs()
            ELSE
               cf_mm = trim(cr)
            END IF
         END IF
         !!
         !CASE('-G')
         !   IF ( jarg + 1 > iargc() ) THEN
         !      PRINT *, 'ERROR: Missing grid type to write to ("T" or "UV"?) !!' ; call usage_test_stuffs()
         !   ELSE
         !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         !      IF ( ANY(clist_opt == TRIM(cr)) ) THEN
         !         PRINT *, 'ERROR: Missing grid type to write to ("T" or "UV"?) !!'; CALL usage_test_stuffs()
         !      ELSE
         !         cgrid_trg = TRIM(cr)
         !      END IF
         !   END IF
         !!
         !CASE('-i')
         !   IF ( jarg + 2 > iargc() ) THEN ! checking that there is at least 2 other arguments following
         !      PRINT *, 'ERROR: Missing input file names!' ; call usage_test_stuffs()
         !   ELSE
         !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         !      IF ( ANY(clist_opt == trim(cr)) ) THEN
         !         PRINT *, 'ERROR: Missing input file name 1 !' ; call usage_test_stuffs()
         !      ELSE
         !         cufilin = trim(cr)
         !      END IF
         !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         !      IF ( ANY(clist_opt == trim(cr)) ) THEN
         !         PRINT *, 'ERROR: Missing input file name 2 !' ; call usage_test_stuffs()
         !      ELSE
         !         cvfilin = trim(cr)
         !      END IF
         !   END IF
         !!
         !CASE('-t')
         !   IF ( jarg + 1 > iargc() ) THEN
         !      PRINT *, 'ERROR: Missing time name!' ; call usage_test_stuffs()
         !   ELSE
         !      jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         !      IF ( ANY(clist_opt == trim(cr)) ) THEN
         !         PRINT *, 'ERROR: Missing time name!'; call usage_test_stuffs()
         !      ELSE
         !         cv_time_0 = trim(cr)
         !      END IF
         !   END IF
         !   !!
         !CASE('-1')
         !   l_anlt = .TRUE.
         !!
      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
         CALL usage_test_stuffs()
         !!
      END SELECT
      !!
      !!
   END DO


   !IF ( (.NOT. l_inv).AND.( TRIM(cufilin) /= 'none' ) ) THEN
   !   PRINT *, 'ERROR: specify the "-I" switch if you want to perform inverse correction!'
   !   STOP
   !END IF

   !IF ( (cgrid_trg /= 'T').AND.(cgrid_trg /= 'U') ) call usage_test_stuffs()

   !PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   !PRINT *, ''


   !IF ( cgrid_trg == 'U' ) PRINT *, ' *** Gonna save on grid U and V.'
   !IF ( cgrid_trg == 'T' ) PRINT *, ' *** Gonna save on grid T.'
   !PRINT *, ''


   !IF ( l_inv ) THEN
   !   PRINT *, ' * Vector files to unrotate = ', trim(cufilin), ' , ', trim(cvfilin)
   !   !PRINT *, '   => associated variable names = ', TRIM(cv_rot_U), ' , ', TRIM(cv_rot_V)
   !   IF ( trim(cv_time_0) == 'none' ) THEN
   !      PRINT *, 'ERROR: you must specify the name of time variable with the "-t" switch!'; STOP
   !   END IF
   !   PRINT *, '   => time variable name = ', trim(cv_time_0)
   !!ELSE
   !   !PRINT *, ' * Name for corrected vector components = ', TRIM(cv_rot_U), ' , ', TRIM(cv_rot_V)
   !END IF

   PRINT *, ' * mesh_mask file to use = ', TRIM(cf_mm)

   !cnmlst_x = TRIM(cf_nml_sosie)//'_x'
   !cnmlst_y = TRIM(cf_nml_sosie)//'_y'

   !PRINT *, ' * namelists we expect => ', TRIM(cnmlst_x)//' and '//TRIM(cnmlst_y)
   !PRINT *, ''


   INQUIRE(FILE=trim(cf_mm), EXIST=lexist )
   IF ( .NOT. lexist ) THEN
      WRITE(*,'("The mesh_mask file ",a," file was not found!")') trim(cf_mm)
      CALL usage_test_stuffs()
   END IF



   CALL DIMS(cf_mm, cv_glamt, ni, nj, nk, Ntr)

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
      &    mask_t(ni,nj,nk), &
      &    xlamt(ni,nj)  , xphit(ni,nj), Xdum8(ni,nj), &
      &    xlamu(ni,nj)  , xphiu(ni,nj) , xlamv(ni,nj)  , xphiv(ni,nj) , &
      &    xlamf(ni,nj) , xphif(ni,nj), &
      &    vtime(Ntr)   )

   ALLOCATE (XCOST8(ni,nj) , XSINT8(ni,nj) , XCOSU8(ni,nj) , XSINU8(ni,nj) , XCOSV8(ni,nj) , XSINV8(ni,nj) , &
      &      XCOSF8(ni,nj) , XSINF8(ni,nj) , U_r8(ni,nj) , V_r8(ni,nj) )


   IF ( i3d == 1 ) THEN
      ALLOCATE ( vdepth(nk) )
      CALL GETVAR_1D(cf_mm, cv_depth, vdepth)
      !IF ( cgrid_trg == 'U' ) THEN
      !   CALL GETMASK_3D(cf_lsm_trg, 'umask', mask_u(:,:,:))
      !   CALL GETMASK_3D(cf_lsm_trg, 'vmask', mask_v(:,:,:))
      !ELSE
      CALL GETMASK_3D(cf_mm, 'tmask', mask_t(:,:,:))
      !END IF
   ELSE
      ! IF ( cgrid_trg == 'U' ) THEN
      !    CALL GETMASK_2D(cf_lsm_trg, 'umask', mask_u(:,:,1), jlev=1)
      !    CALL GETMASK_2D(cf_lsm_trg, 'vmask', mask_v(:,:,1), jlev=1)
      ! ELSE
      CALL GETMASK_2D(cf_mm, 'tmask', mask_t(:,:,1), jlev=1)
      !END IF
   END IF

   !!  Getting longitude and latitude form grid file :
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamt, 1, 1, 1, xlamt)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphit, 1, 1, 1, xphit)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamu, 1, 1, 1, xlamu)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiu, 1, 1, 1, xphiu)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamv, 1, 1, 1, xlamv)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiv, 1, 1, 1, xphiv)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamf, 1, 1, 1, xlamf)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphif, 1, 1, 1, xphif)



   PRINT *, ''

   !! Is this a known ORCA grid (just for info now, not used!):
   iorca = IS_ORCA_NORTH_FOLD( xphit )
   IF ( iorca == 4 ) PRINT *, ' Grid is an ORCA grid with north-pole T-point folding! Type ',iorca
   IF ( iorca == 6 ) PRINT *, ' Grid is an ORCA grid with north-pole F-point folding! Type ',iorca
   PRINT *, ''


   DO jj = 2, nj-1
      DO ji = 2, ni-1
         Xdum8(ji,jj) = COS(4.*xlamt(ji,jj)*rad0)*SIN(4.*xphit(ji,jj)*rad0)
      END DO
   END DO


   CALL lbc_lnk_2d( iorca, Xdum8, 'T', 1.0_8 )

   CALL DUMP_2D_FIELD(REAL(Xdum8,4), 'test.nc', 'test',  xlamt, xphit, cv_glamt, cv_gphit)



   !STOP


   !!  Getting cosine and sine corresponding to the angle of the local distorsion of the grid:
   CALL ANGLE( iorca, xlamt, xphit, xlamu, xphiu, xlamv, xphiv, xlamf, xphif, &
      &        XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8 )


   XSINT8(:,nj) = 0. ; ! cancelling what has be done in ANGLE
   CALL lbc_lnk_2d( iorca, XSINT8, 'T', -1.0_8 )

   XSINU8(:,nj) = 0. ; ! cancelling what has be done in ANGLE
   CALL lbc_lnk_2d( iorca, XSINU8, 'U', -1.0_8 )

   XSINV8(:,nj) = 0. ; ! cancelling what has be done in ANGLE
   CALL lbc_lnk_2d( iorca, XSINV8, 'V', -1.0_8 )


   CALL DUMP_2D_FIELD(REAL(XCOST8,4), 'cost_angle.nc', 'cos',   xlamt, xphit, cv_glamt, cv_gphit)
   CALL DUMP_2D_FIELD(REAL(XSINT8,4), 'sint_angle.nc', 'sin',   xlamt, xphit, cv_glamt, cv_gphit)
   CALL DUMP_2D_FIELD(REAL(XCOSU8,4), 'cosu_angle.nc', 'cos',   xlamt, xphit, cv_glamt, cv_gphit)
   CALL DUMP_2D_FIELD(REAL(XSINU8,4), 'sinu_angle.nc', 'sin',   xlamt, xphit, cv_glamt, cv_gphit)
   CALL DUMP_2D_FIELD(REAL(XCOSV8,4), 'cosv_angle.nc', 'cos',   xlamt, xphit, cv_glamt, cv_gphit)
   CALL DUMP_2D_FIELD(REAL(XSINV8,4), 'sinv_angle.nc', 'sin',   xlamt, xphit, cv_glamt, cv_gphit)
   !STOP

   !CALL P2D_T(id_f1, id_v1, Ntr, jt, xlamu, xphiu,         vtime, Uu_c(:,:,1), &
   !   &    cufilout, 'nav_lon_u', 'nav_lat_u', cv_t_out, cv_rot_U,       &
   !   &    rmiss_val, attr_time=vatt_info_t)


   PRINT *, 'Done!'

END PROGRAM TEST_STUFFS



SUBROUTINE usage_test_stuffs()
   !!
   PRINT *,''
   PRINT *,'   List of command line options:'
   PRINT *,''
   !PRINT *,' -I   => will perform inverse correction, ie un-rotate a ORCA grid'
   !PRINT *,''
   !PRINT *,''
   !PRINT *,'  *** MANDATORY for both normal and inverse mode:'
   !PRINT *,''
   !PRINT *,' -x  <name U>         => Specify name for x comp. in output file'
   !PRINT *,'                         (or input file if inverse mode)'
   !PRINT *,''
   !PRINT *,' -y  <name V>         => Specify name for y comp. in output file'
   !PRINT *,'                         (or input file if inverse mode)'
   !PRINT *,''
   PRINT *,' -m  <mesh_mask_file> => Specify which mesh_mask file to use'
   PRINT *,''
   !PRINT *,''
   !PRINT *,'  ***  MANDATORY for normal mode (no -I switch) :'
   !PRINT *,''
   !PRINT *,' -f  <namelist_file>  => Specify which namelist file to use'
   !PRINT *, '                        No namelist needed when inverse correction'
   !PRINT *,''
   !PRINT *,' -G  <T/U>            => Specify if you want to save rotated vector'
   !PRINT *, '                        on T-grid (T) or U- and V-grid (U)'
   !PRINT *,''
   !PRINT *,''
   !PRINT *,'  *** MANDATORY for INVERSE MODE (-I switch) :'
   !PRINT *,''
   !PRINT *,' -i <x.nc> <y.nc>     =>  unrotate vector fields given in these 2 files'
   !PRINT *,'                          to the same grid'
   !PRINT *,''
   !PRINT *,' -t <time_name>       =>  name of time variable in <x.nc> and <y.nc>'
   !PRINT *,''
   !PRINT *,''
   !PRINT *,'  *** MISC options :'
   !PRINT *,''
   !PRINT *,' -h                   => Show this message'
   !PRINT *,''
   STOP
   !!
END SUBROUTINE usage_test_stuffs
!!
