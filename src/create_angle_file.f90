PROGRAM CREATE_ANGLE_FILE

   !! Purpose of this program is to:
   !! 1/ Read all the `glamX` and `gphiX` into a  NEMO `domain_cfg.nc` or `mesh_mask.nc` file..
   !! 2/ Compute all the rotation angles (actually their sin() and cos())
   !! 3/ Save them into a file..

   !! ==> use `examples/vector_correction/do_mono_angle.sh` to generate a monlithic file from
   !!     all the cos/sin files created by this program...


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

   LOGICAL, PARAMETER :: ldebug = .true.

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

   !CHARACTER(LEN=400)  :: cn_xtr_x, cn_xtr_y, cextra_x, cextra_y

   !CHARACTER(len=3)    :: cdum
   !CHARACTER(len=3)    :: cgrid_trg='T'
   !CHARACTER(len=80)   :: cv_time_0 = 'none'
   CHARACTER(len=800)  :: cr, cf_mm

   !CHARACTER(len=80) :: &
   !   &    cv_u_in='0', cv_v_in='0', &
   !   &    cv_out_U    = 'uraw',     &   ! raw U name
   !   &    cv_out_V    = 'vraw'          ! raw V name

   !CHARACTER(len=800)  :: &
   !   &    cf_raw_U,    &  ! file containing u_raw
   !   &    cf_raw_V,    &  ! file containing v_raw
   !   &    cf_out_U, cf_out_V,    &
   !   &    cufilin = 'none',  cvfilin = 'none'

   !CHARACTER(len=80)  ::  &
   !   &    cv_rot_U ,  &  ! output name for U corrected
   !   &    cv_rot_V       ! output name for V corrected

   TYPE(grid_type) :: gt_orca

   INTEGER      :: &
      &    jarg, &
      &    iorca=0,      &
                                !   &    nlext=3, &
      &    i0, j0, &
      &    ni, nj, nk
   !   &    ni1, nj1, Ntr1,       &
   !   &    ni2, nj2, Ntr2,       &
   !   &    ni_g, nj_g, nk_g,    &
   !   &    iargc,         &
   !   &    idf_u, idv_u, idf_v, idv_v, &
   !   &    id_f1, id_v1, &
   !   &    id_f2, id_v2

   !INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: mask_t, mask_u, mask_v

   !REAL(4), DIMENSION(:,:), ALLOCATABLE :: ztmp4


   !REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: U_c, V_c

   REAL(8), DIMENSION(:,:), ALLOCATABLE ::      &
      &    XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8, &
      &    U_r8, V_r8,  ztmp8, Xdum8, Ydum8, &
      &    xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f

   !REAL(8), DIMENSION(:), ALLOCATABLE ::   vtime, vdepth


   !INTEGER :: jt, jk

   LOGICAL :: &
                                !   &    l_inv    = .FALSE., &
                                !   &    l_anlt   = .FALSE., & !: analytic input field (U=1, V=0) DEBUG...
                                !   &    l_3d_inv = .FALSE., & !: will treat 3d files in inverse mode...
                                !   &    lmout_x, lmout_y,   &
      &    lexist

   CHARACTER(LEN=2), DIMENSION(2), PARAMETER :: &
      &            clist_opt = (/ '-h','-m' /)

   WRITE(6,*)''
   WRITE(6,*)'=========================================================='
   WRITE(6,*)'            S  O  S  I  E    version ', trim(csosie_version)
   WRITE(6,*)''
   WRITE(6,*)'            Rotation angles of distorted NEMO grid        '
   WRITE(6,*)'=========================================================='
   WRITE(6,*)''


   !! Getting string arguments :
   !! --------------------------

   jarg = 0

   DO WHILE ( jarg < iargc() )
      !!
      jarg = jarg + 1
      CALL getarg(jarg,cr)
      !!
      !!
      SELECT CASE (TRIM(cr))
         !!
         !!
      CASE('-h')
         call usage_caf()
         !!
      CASE('-m')
         IF ( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
            PRINT *,'ERROR: Missing `domain_cfg` or `mesh_mask` file name!' ; call usage_caf()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF ( ANY(clist_opt == trim(cr)) ) THEN
               PRINT *,'ERROR: ', trim(cr), ' is definitively not the name of the mesh_mask file!'
               call usage_caf()
            ELSE
               cf_mm = trim(cr)
            END IF
         END IF
         !!
         !!
      CASE DEFAULT
         PRINT *,'Unknown option: ', trim(cr) ; PRINT *,''
         CALL usage_caf()
         !!
      END SELECT
      !!
   END DO

   INQUIRE(FILE=trim(cf_mm), EXIST=lexist )
   IF ( .NOT. lexist ) THEN
      WRITE(*,'("The mesh_mask file ",a," file was not found!")') trim(cf_mm)
      CALL usage_caf()
   END IF

   CALL DIMS(cf_mm, cv_glamt, ni, nj, nk, Ntr)


   WRITE(*,'("Space dimension is : ",i4," x",i4," x",i4)') ni, nj, nk
   WRITE(*,*)

   !! Allocations :
   ALLOCATE (  &
      &    xlon_t(ni,nj)  , xlat_t(ni,nj), &
      &    xlon_u(ni,nj)  , xlat_u(ni,nj) , xlon_v(ni,nj)  , xlat_v(ni,nj) , &
      &    xlon_f(ni,nj) , xlat_f(ni,nj), ztmp8(ni,nj) )

   ALLOCATE ( XCOST8(ni,nj) , XSINT8(ni,nj), &
      &       XCOSU8(ni,nj) , XSINU8(ni,nj), &
      &       XCOSV8(ni,nj) , XSINV8(ni,nj), &
      &       XCOSF8(ni,nj) , XSINF8(ni,nj)  )


   !!  Getting longitude and latitude form grid file :
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamt, 1, 1, 1, xlon_t)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphit, 1, 1, 1, xlat_t)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamu, 1, 1, 1, xlon_u)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiu, 1, 1, 1, xlat_u)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamv, 1, 1, 1, xlon_v)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiv, 1, 1, 1, xlat_v)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_glamf, 1, 1, 1, xlon_f)
   CALL GETVAR_2D(i0, j0, cf_mm, cv_gphif, 1, 1, 1, xlat_f)


   PRINT *,''

   !! Is this a known ORCA grid (just for info now, not used!):
   gt_orca = IS_ORCA_NORTH_FOLD( xlat_t )
   iorca = gt_orca%ifld_nord
   IF ( iorca == 4 ) PRINT *,' Grid is an ORCA grid with north-pole T-point folding!'
   IF ( iorca == 6 ) PRINT *,' Grid is an ORCA grid with north-pole F-point folding!'
   PRINT *,''


   !!  Getting cosine and sine corresponding to the angle of the local distorsion of the grid:
   CALL ANGLE( iorca, xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f, &
      &        XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8 )

   IF ( ldebug ) THEN
      CALL DUMP_FIELD(REAL(XCOST8,4), 'cost_angle.nc', 'cost')
      CALL DUMP_FIELD(REAL(XSINT8,4), 'sint_angle.nc', 'sint')
      CALL DUMP_FIELD(REAL(XCOSU8,4), 'cosu_angle.nc', 'cosu')
      CALL DUMP_FIELD(REAL(XSINU8,4), 'sinu_angle.nc', 'sinu')
      CALL DUMP_FIELD(REAL(XCOSV8,4), 'cosv_angle.nc', 'cosv')
      CALL DUMP_FIELD(REAL(XSINV8,4), 'sinv_angle.nc', 'sinv')
      CALL DUMP_FIELD(REAL(XCOSF8,4), 'cosf_angle.nc', 'cosf')
      CALL DUMP_FIELD(REAL(XSINF8,4), 'sinf_angle.nc', 'sinf')
   END IF


   !! 2D:
   !IF ( cgrid_trg == 'U,V' ) THEN
   !   CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u,         vtime, U_c(:,:,1), &
   !      &    cf_out_U, 'nav_lon_u', 'nav_lat_u', cv_t_out, cv_rot_U,       &
   !      &    rmiss_val, attr_t=vatt_info_t)
   !   CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v,         vtime, V_c(:,:,1), &
   !      &    cf_out_V, 'nav_lon_v', 'nav_lat_v', cv_t_out, cv_rot_V,   &
   !      &    rmiss_val, attr_t=vatt_info_t)
   !ELSE
   !   CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_t, xlat_t,         vtime, U_c(:,:,1), &
   !      &    cf_out_U, 'nav_lon', 'nav_lat', cv_t_out, cv_rot_U,       &
   !      &    rmiss_val, attr_t=vatt_info_t)
   !   CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_t, xlat_t,         vtime, V_c(:,:,1), &
   !      &    cf_out_V, 'nav_lon', 'nav_lat', cv_t_out, cv_rot_V,   &
   !      &    rmiss_val, attr_t=vatt_info_t)
   !END IF


CONTAINS




   SUBROUTINE usage_caf()
      !!
      PRINT *,''
      PRINT *,'   List of command line options:'
      PRINT *,''
      PRINT *,' -m  <mesh_mask_file> => Specify `mesh_mask` or `domain_cfg` file to use'
      PRINT *,''
      PRINT *,' -h                   => Show this message'
      PRINT *,''
      STOP
      !!
   END SUBROUTINE usage_caf






END PROGRAM CREATE_ANGLE_FILE
