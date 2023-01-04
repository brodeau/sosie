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
   CHARACTER(len=6), PARAMETER   :: &
      &    cv_glamt     = 'glamt',   &   ! input grid longitude name, T-points
      &    cv_gphit     = 'gphit',   &   ! input grid latitude name,  T-points
      &    cv_glamu     = 'glamu',   &   ! input grid longitude name, U-points
      &    cv_gphiu     = 'gphiu',   &   ! input grid latitude name,  U-points
      &    cv_glamv     = 'glamv',   &   ! input grid longitude name, V-points
      &    cv_gphiv     = 'gphiv',   &   ! input grid latitude name,  V-points
      &    cv_glamf     = 'glamf',   &   ! input grid longitude name, F-points
      &    cv_gphif     = 'gphif',   &   ! input grid latitude name,  F-points
      &    cv_depth     = 'deptht'       !  depth at T-points (U-points and V-points too)

   CHARACTER(len=800)  :: cr, cf_mm

   TYPE(grid_type) :: gt_orca

   INTEGER      :: &
      &    jarg, &
      &    iorca=0,      &
                                !   &    nlext=3, &
      &    i0, j0, &
      &    ni, nj, nk

   REAL(8), DIMENSION(:,:), ALLOCATABLE ::      &
      &    XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8, &
      &    ztmp8, &
      &    xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f

   LOGICAL :: lexist

   CHARACTER(LEN=2), DIMENSION(2), PARAMETER :: clist_opt = (/ '-h','-m' /)

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
   CALL ANGLE2( iorca, xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f, &
      &         XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8 )

   CALL DUMP_FIELD( XCOST8, 'cost_angle.nc', 'cost')
   CALL DUMP_FIELD( XSINT8, 'sint_angle.nc', 'sint')
   CALL DUMP_FIELD( XCOSU8, 'cosu_angle.nc', 'cosu')
   CALL DUMP_FIELD( XSINU8, 'sinu_angle.nc', 'sinu')
   CALL DUMP_FIELD( XCOSV8, 'cosv_angle.nc', 'cosv')
   CALL DUMP_FIELD( XSINV8, 'sinv_angle.nc', 'sinv')
   CALL DUMP_FIELD( XCOSF8, 'cosf_angle.nc', 'cosf')
   CALL DUMP_FIELD( XSINF8, 'sinf_angle.nc', 'sinf')


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
