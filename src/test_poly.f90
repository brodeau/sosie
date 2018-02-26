PROGRAM TEST_POLY

   !USE mod_conf
   !USE mod_init
   !USE io_ezcdf
   !USE mod_grids,  ONLY: IS_ORCA_NORTH_FOLD
   !USE mod_nemotools
   USE MOD_POLY

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
   !CHARACTER(len=80), PARAMETER   :: &
   !   &    cv_glamt     = 'glamt',   &   ! input grid longitude name, T-points
   !   &    cv_gphit     = 'gphit',   &   ! input grid latitude name,  T-points
   !   &    cv_glamu     = 'glamu',   &   ! input grid longitude name, U-points
   !   &    cv_gphiu     = 'gphiu',   &   ! input grid latitude name,  U-points


   !CHARACTER(LEN=400)  :: cextra_x, cextra_y 
   
   !CHARACTER(len=3)    :: cdum
   !CHARACTER(len=1)    :: cgrid_out='0'
   !CHARACTER(len=80)   :: cv_time_0 = 'none', cfext = 'nc'
   !CHARACTER(len=800)  :: cr, cf_mm, cnmlst_x, cnmlst_y

   !CHARACTER(len=80) :: &
   !   &    cv_out_U    = 'uraw',   &   ! raw U name
   !   &    cv_out_V    = 'vraw'        ! raw V name

   !CHARACTER(len=800)  :: &
   !   &    cf_out_U,    &  ! file containing u_raw
   !   &    cf_out_V,    &  ! file containing v_raw
   !   &    cufilout, cvfilout,    &
   !   &    cufilin = 'none',  cvfilin = 'none'

   !CHARACTER(len=80)  ::  &
   !   &    cv_rot_U ,  &  ! output name for U corrected
   !   &    cv_rot_V       ! output name for V corrected

   !INTEGER      :: &
   !   &    Ntr, &
   !   &    jarg, i3d, nbc, &
  

      !~REAL(4), DIMENSION(:,:), ALLOCATABLE :: Xdum4


      !REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: Ut_c, Vt_c, Uu_c, Vv_c

      !~REAL(8), DIMENSION(:,:), ALLOCATABLE ::      &
   !&    XCOST8, XSINT8, XCOSU8, XSINU8, XCOSV8, XSINV8, XCOSF8, XSINF8, &
   !   &    U_r8, V_r8,  Xdum8, &
   !   &    xlamt, xphit, xlamu, xphiu, xlamv, xphiv, xlamf, xphif

!   REAL(8), DIMENSION(:), ALLOCATABLE ::   vtime, vdepth


 !  INTEGER :: jt, jk

!!

   REAL(4) :: x1, x2, x3, x4, &
      &       y1, y2, y3, y4
   
   REAL(8) :: xp, yp

   LOGICAL :: l_inside
   
   REAL(8), DIMENSION(4) :: vlon, vlat

   

   PRINT *, ' Give the 4 longitudes of the 4 points definingf your mesh:'
   READ(*, *) x1, x2, x3, x4
   PRINT *, ''   
   PRINT *, ' Give the 4 latitudes of the 4 points definingf your mesh:'
   READ(*, *) y1, y2, y3, y4
   PRINT *, ''

   PRINT *, ' Give the "lon lat" of the point to consider:'
   READ(*, *) xp, yp
   PRINT *, ''


   
   !PRINT *, ' x1, x2, x3, x4 =>', x1, x2, x3, x4
   !PRINT *, ' y1, y2, y3, y4 =>', y1, y2, y3, y4
   PRINT *, 'xp, yp          =>', xp, yp


   vlon = (/ x1, x2, x3, x4 /)
   vlat = (/ y1, y2, y3, y4 /)

   
   PRINT *, ' * vlon:', REAL(vlon,4)
   PRINT *, ' * vlat:', REAL(vlat,4)

   l_inside = L_InPoly ( vlon, vlat, xp, yp )
   
   IF ( l_inside ) THEN
      PRINT *, '  *** It is inside !!!'      
   ELSE
      PRINT *, '  *** NOT inside ...'
   END IF
   


END PROGRAM TEST_POLY
