MODULE MOD_MANIP

   !! Misc. manipulations and operations on 2D arrays...

   !! Author: L. Brodeau (brodeau@gmail.com)

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: fill_extra_bands, extra_2_east, extra_2_west, test_xyz, partial_deriv, &
      &      flip_ud_2d, flip_ud_3d, long_reorg_2d, long_reorg_3d, &
      &      distance, find_nearest_point_idiot, find_nearest_point_reg

   LOGICAL, PARAMETER :: ldebug = .FALSE.

   LOGICAL  :: lfirst_dist = .TRUE.



CONTAINS




   SUBROUTINE FILL_EXTRA_BANDS(k_ew, XX, YY, XF, XP4, YP4, FP4)

      !!============================================================================
      !! Extending input arrays with an extraband of two points at north,south,east
      !! and west boundaries.
      !!
      !! The extension is done thanks to Akima's exptrapolation method.
      !!
      !! East-west periodicity of global map is taken into account through 'k_ew' :
      !!
      !!
      !!  k_ew : east-west periodicity on the input file/grid
      !!         k_ew = -1  --> no east-west periodicity (along x)
      !!         k_ew >= 0  --> east-west periodicity with overlap of k_ew points (along x)
      !!
      !!
      !!                       Author : Laurent BRODEAU, 2007
      !!============================================================================

      INTEGER ,                INTENT(in)  :: k_ew
      REAL(8), DIMENSION(:,:), INTENT(in)  :: XX, YY, XF
      REAL(8), DIMENSION(:,:), INTENT(out) :: XP4, YP4, FP4

      !! Local
      INTEGER :: nx, ny, nxp4, nyp4
      INTEGER :: ji, jj

      IF ( (SIZE(XX,1) /= SIZE(YY,1)).OR.(SIZE(XX,2) /= SIZE(YY,2)).OR. &
         & (SIZE(XX,1) /= SIZE(XF,1)).OR.(SIZE(XX,2) /= SIZE(XF,2))) THEN
         PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : size of input coor. and data do not match!!!'; STOP
      END IF

      IF ( (SIZE(XP4,1) /= SIZE(YP4,1)).OR.(SIZE(XP4,2) /= SIZE(YP4,2)).OR. &
         & (SIZE(XP4,1) /= SIZE(FP4,1)).OR.(SIZE(XP4,2) /= SIZE(FP4,2))) THEN
         PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : size of output coor. and data do not match!!!'; STOP
      END IF

      nx = SIZE(XX,1)
      ny = SIZE(XX,2)

      nxp4 = SIZE(XP4,1)
      nyp4 = SIZE(XP4,2)

      IF ( nxp4 /= nx + 4 ) THEN
         PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : target x dim is not ni+4!!!'; STOP
      END IF
      IF ( nyp4 /= ny + 4 ) THEN
         PRINT *, 'ERROR, mod_drown.F90 => FILL_EXTRA_BANDS : target y dim is not nj+4!!!'; STOP
      END IF




      !!   C r e a t i n g   e x t e n d e d   a r r a y s  :
      !!   --------------------------------------------------

      !! Initializing :
      XP4 = 0.
      YP4 = 0.
      FP4 = 0.

      !! Filling center of domain:
      XP4(  3:nxp4-2, 3:nyp4-2) =   XX(:,:)
      YP4(  3:nxp4-2, 3:nyp4-2) =   YY(:,:)
      FP4(3:nxp4-2, 3:nyp4-2) = XF(:,:)



      !! X array :
      !! ---------

      IF (k_ew /= -1) THEN   ! we can use east-west periodicity of input file to
         !!                   ! fill extra bands :
         XP4( 1     , 3:nyp4-2) = XX(nx - 1 - k_ew , :) - 360.
         XP4( 2     , 3:nyp4-2) = XX(nx - k_ew     , :) - 360.
         XP4(nxp4   , 3:nyp4-2) = XX( 2 + k_ew     , :) + 360.
         XP4(nxp4-1 , 3:nyp4-2) = XX( 1 + k_ew     , :) + 360.

      ELSE

         !! Left side :
         XP4(2, 3:nyp4-2) = XX(2,:) - (XX(3,:) - XX(1,:))
         XP4(1, 3:nyp4-2) = XX(1,:) - (XX(3,:) - XX(1,:))

         !! Right side :
         XP4(nxp4-1, 3:nyp4-2) = XX(nx-1,:) + XX(nx,:) - XX(nx-2,:)
         XP4(nxp4  , 3:nyp4-2) = XX(nx,:)   + XX(nx,:) - XX(nx-2,:)

      END IF


      !! Bottom side :
      XP4(:, 2) = XP4(:,4) - (XP4(:,5) - XP4(:,3))
      XP4(:, 1) = XP4(:,3) - (XP4(:,5) - XP4(:,3))

      !! Top side :
      XP4(:,nyp4-1) = XP4(:,nyp4-3) + XP4(:,nyp4-2) - XP4(:,nyp4-4)
      XP4(:,nyp4)   = XP4(:,nyp4-2) + XP4(:,nyp4-2) - XP4(:,nyp4-4)



      !! Y array :
      !! ---------

      !! Top side :
      YP4(3:nxp4-2, nyp4-1) = YY(:, ny-1) + YY(:,ny) - YY(:,ny-2)
      YP4(3:nxp4-2, nyp4)   = YY(:, ny)   + YY(:,ny) - YY(:,ny-2)
      !! Bottom side :
      YP4(3:nxp4-2, 2) = YY(:,2) - (YY(:,3) - YY(:,1))
      YP4(3:nxp4-2, 1) = YY(:,1) - (YY(:,3) - YY(:,1))

      IF (k_ew /= -1) THEN

         YP4( 1     , :) = YP4(nx - 1 - k_ew + 2, :)
         YP4( 2     , :) = YP4(nx - k_ew     + 2, :)
         YP4(nxp4   , :) = YP4( 2 + k_ew     + 2, :)
         YP4(nxp4-1 , :) = YP4( 1 + k_ew     + 2, :)

      ELSE
         !! Left side :
         YP4(2, :) = YP4(4,:) - (YP4(5,:) - YP4(3,:))
         YP4(1, :) = YP4(3,:) - (YP4(5,:) - YP4(3,:))
         !! Right side :
         YP4(nxp4-1,:) = YP4(nxp4-3,:) + YP4(nxp4-2, :) - YP4(nxp4-4, :)
         YP4(nxp4,:)   = YP4(nxp4-2,:) + YP4(nxp4-2,:)  - YP4(nxp4-4, :)

      END IF


      !! Data array :
      !! ------------

      IF (k_ew /= -1) THEN

         FP4( 1     , 3:nyp4-2) = XF(nx - 1 - k_ew , :)
         FP4( 2     , 3:nyp4-2) = XF(nx - k_ew     , :)
         FP4(nxp4   , 3:nyp4-2) = XF( 2 + k_ew     , :)
         FP4(nxp4-1 , 3:nyp4-2) = XF( 1 + k_ew     , :)

      ELSE

         !! Left side :
         DO jj = 3, nyp4-2
            CALL extra_2_east(XP4(nxp4-4,jj),XP4(nxp4-3,jj),XP4(nxp4-2,jj),        &
               &              XP4(nxp4-1,jj),XP4(nxp4,jj),                         &
               &              FP4(nxp4-4,jj),FP4(nxp4-3,jj),FP4(nxp4-2,jj),  &
               &              FP4(nxp4-1,jj),FP4(nxp4,jj) )
         END DO

         !! Right side :
         DO jj = 3, nyp4-2
            CALL extra_2_west(XP4(5,jj),XP4(4,jj),XP4(3,jj),         &
               &              XP4(2,jj),XP4(1,jj),                   &
               &              FP4(5,jj),FP4(4,jj),FP4(3,jj),   &
               &              FP4(2,jj),FP4(1,jj) )
         END DO

      END IF

      !! Top side :
      DO ji = 1, nxp4
         CALL extra_2_east(YP4(ji,nyp4-4),YP4(ji,nyp4-3),YP4(ji,nyp4-2),       &
            &              YP4(ji,nyp4-1),YP4(ji,nyp4),                        &
            &              FP4(ji,nyp4-4),FP4(ji,nyp4-3),FP4(ji,nyp4-2), &
            &              FP4(ji,nyp4-1),FP4(ji,nyp4) )
      END DO

      !! Bottom side :
      DO ji = 1, nxp4
         CALL extra_2_west(YP4(ji,5),YP4(ji,4),YP4(ji,3),       &
            &              YP4(ji,2),YP4(ji,1),                 &
            &              FP4(ji,5),FP4(ji,4),FP4(ji,3), &
            &              FP4(ji,2),FP4(ji,1) )
      END DO

   END SUBROUTINE FILL_EXTRA_BANDS




   SUBROUTINE extra_2_east(x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)
      !!
      !!============================================================================
      !!
      !! Extrapolates 2 extra east (or north) points of a curve with Akima's 1D method
      !!
      !! Input  : x1, x2, x3, x4, x5, y1, y2, y3
      !! Output : y4, y5
      !!
      !!                       Author : Laurent BRODEAU, 2007
      !!============================================================================
      !!
      !!
      REAL(8), INTENT(in)  :: x1, x2, x3, x4, x5, y1, y2, y3
      REAL(8), INTENT(out) :: y4, y5
      !!
      !! Local :
      REAL(8) :: A, B, C, D, ALF, BET
      !!
      !!
      A    = x2 - x1
      B    = x3 - x2
      C    = x4 - x3
      D    = x5 - x4
      !!
      ALF  = y2 - y1
      BET  = y3 - y2
      !!
      IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
         y4 = y3 ; y5 = y3
      ELSE
         y4   = C*(2*BET/B - ALF/A) + y3
         y5   = y4 + y4*D/C + BET*D/B - ALF*D/A - y3*D/C
      END IF
      !!
      !!
   END SUBROUTINE extra_2_east
   !!
   !!
   !!
   !!
   SUBROUTINE extra_2_west(x5, x4, x3, x2, x1, y5, y4, y3, y2, y1)
      !!
      !!============================================================================
      !!
      !! Extrapolates 2 extra west (or south) points of a curve with Akima's 1D method
      !!
      !! Input  : x1, x2, x3, x4, x5, y1, y2, y3
      !! Output : y4, y5
      !!
      !!                       Author : Laurent BRODEAU, 2007
      !!============================================================================
      !!
      !!
      REAL(8), INTENT(in)  :: x1, x2, x3, x4, x5, y5, y4, y3
      REAL(8), INTENT(out) :: y1, y2
      REAL(8) :: A, B, C, D, ALF, BET
      !!
      !! x1 -> x5
      !! x2 -> x4
      !! x3 -> x3
      !! x4 -> x2
      !! x5 -> x1
      !!
      A    = x4 - x5
      B    = x3 - x4
      C    = x2 - x3
      D    = x1 - x2
      !!
      ALF  = y4 - y5
      BET  = y3 - y4
      !!
      IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
         y2 = y3; y1 = y3
      ELSE
         y2   = C*(2*BET/B - ALF/A) + y3
         y1   = y2 + y2*D/C + BET*D/B - ALF*D/A - y3*D/C
      END IF
      !!
      !!
   END SUBROUTINE extra_2_west
   !!
   !!
   FUNCTION TEST_XYZ(rx, ry, rz)
      !!
      !! Testing if 2D coordinates or 1D, and if match shape of data...
      !!
      CHARACTER(len=2) :: TEST_XYZ
      !!
      REAL(8), DIMENSION(:,:), INTENT(in) :: rx, ry
      REAL(4), DIMENSION(:,:), INTENT(in) :: rz
      !!
      INTEGER :: ix1, ix2, iy1, iy2, iz1, iz2
      !!
      ix1 = size(rx,1) ; ix2 = size(rx,2)
      iy1 = size(ry,1) ; iy2 = size(ry,2)
      iz1 = size(rz,1) ; iz2 = size(rz,2)
      !!
      IF ( (ix2 == 1).AND.(iy2 == 1) ) THEN
         !!
         IF ( (ix1 == iz1).AND.(iy1 == iz2) ) THEN
            TEST_XYZ = '1d'
         ELSE
            PRINT *, 'ERROR, mod_drown.F90 = >TEST_XYZ 1 : longitude and latitude array do not match data!'
            PRINT *, ''; STOP
         END IF
         !!
      ELSE
         IF ( (ix1 == iz1).AND.(iy1 == iz1).AND.(ix2 == iz2).AND.(iy2 == iz2) ) THEN
            TEST_XYZ = '2d'
         ELSE
            PRINT *, 'ERROR, mod_drown.F90 = >TEST_XYZ 2 : longitude and latitude array do not match data!'
            PRINT *, ''; STOP
         END IF
      END IF
      !!
      !!
   END FUNCTION TEST_XYZ


   SUBROUTINE PARTIAL_DERIV(k_ew, XX, XY, XF, dFdX, dFdY, d2FdXdY)

      !! Partial derivatives of a field ZF given on a regular gird !!!

      !!  k_ew : east-west periodicity on the input file/grid
      !!         k_ew = -1  --> no east-west periodicity (along x)
      !!         k_ew >= 0  --> east-west periodicity with overlap of k_ew points (along x)

      INTEGER, INTENT(in) :: k_ew

      REAL(8), DIMENSION(:,:), INTENT(in)  :: XX, XY, XF
      REAL(8), DIMENSION(:,:), INTENT(out) :: dFdX, dFdY, d2FdXdY

      !! Local variables :
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: ZX, ZY, ZF
      INTEGER :: nx, ny

      dFdX    = 0.
      dFdY    = 0.
      d2FdXdY = 0.

      nx = SIZE(XF,1) ; ny = SIZE(XF,2)

      !! Extended arrays with a frame of 2 points...
      ALLOCATE ( ZX(nx+4,ny+4), ZY(nx+4,ny+4), ZF(nx+4,ny+4) )

      CALL FILL_EXTRA_BANDS(k_ew, XX, XY, XF, ZX, ZY, ZF)


      !! Second order finite difference:
      !! i+1 => 4:nx+4 / i => 3:nx+2 / i-1 => 2:nx+2
      !! j+1 => 4:ny+4 / j => 3:ny+2 / j-1 => 2:ny+2

      dFdX(:,:) = 0.5*( ZF(4:nx+4,3:ny+2) - ZF(2:nx+2,3:ny+2) ) !/ ( ZX(4:nx+4,3:ny+2) - ZX(2:nx+2,3:ny+2) )
      dFdY(:,:) = 0.5*( ZF(3:nx+2,4:ny+4) - ZF(3:nx+2,2:ny+2) ) !/ ( ZY(3:nx+2,4:ny+4) - ZY(3:nx+2,2:ny+2) )

      !!lolo: Here the denominator is probably wrong:
      d2FdXdY(:,:) = 0.25*( ZF(4:nx+4,4:ny+4) - ZF(4:nx+4,2:ny+2)   -   ZF(2:nx+2,4:ny+4) + ZF(2:nx+2,2:ny+2) ) !&
      !&       / ( ( ZX(4:nx+4,3:ny+2) - ZX(2:nx+2,3:ny+2) ) * ( ZY(3:nx+2,4:ny+4) - ZY(3:nx+2,2:ny+2) ) )

      DEALLOCATE ( ZX, ZY, ZF )

   END SUBROUTINE PARTIAL_DERIV




   SUBROUTINE FLIP_UD_2D(XF)

      REAL(4), DIMENSION(:,:), INTENT(inout) :: XF

      INTEGER :: nx, ny, jj
      REAL(4), DIMENSION(:,:), ALLOCATABLE :: ztmp

      nx = SIZE(XF,1) ; ny = SIZE(XF,2)

      ALLOCATE ( ztmp(nx,ny) )

      ztmp(:,:) = XF(:,:)

      DO jj = 1, ny
         XF(:,jj) =  ztmp(:,ny-jj+1)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE FLIP_UD_2D


   SUBROUTINE LONG_REORG_2D(i_chg_x, XF)

      INTEGER, INTENT(in) :: i_chg_x
      REAL(4), DIMENSION(:,:), INTENT(inout) :: XF

      INTEGER :: nx, ny, ji
      REAL(4), DIMENSION(:,:), ALLOCATABLE :: ztmp

      nx = SIZE(XF,1) ; ny = SIZE(XF,2)

      ALLOCATE ( ztmp(nx,ny) )

      ztmp(:,:) = XF(:,:)

      DO ji = i_chg_x, nx
         XF(ji - i_chg_x + 1 , :) = ztmp(ji,:)
      END DO
      DO ji = 1, i_chg_x - 1
         XF(nx - i_chg_x + 1 + ji , :) = ztmp(ji,:)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE LONG_REORG_2D


   SUBROUTINE FLIP_UD_3D(XF)

      INTEGER(2), DIMENSION(:,:,:), INTENT(inout) :: XF

      INTEGER :: nx, ny, nz, jj
      INTEGER(2), DIMENSION(:,:,:), ALLOCATABLE :: ztmp

      nx = SIZE(XF,1) ; ny = SIZE(XF,2) ; nz = SIZE(XF,3)

      ALLOCATE ( ztmp(nx,ny,nz) )

      ztmp(:,:,:) = XF(:,:,:)

      DO jj = 1, ny
         XF(:,jj,:) =  ztmp(:,ny-jj+1,:)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE FLIP_UD_3D

   SUBROUTINE LONG_REORG_3D(i_chg_x, XF)

      INTEGER, INTENT(in) :: i_chg_x
      INTEGER(2), DIMENSION(:,:,:), INTENT(inout) :: XF

      INTEGER :: nx, ny, nz, ji
      INTEGER(2), DIMENSION(:,:,:), ALLOCATABLE :: ztmp

      nx = SIZE(XF,1) ; ny = SIZE(XF,2) ; nz = SIZE(XF,3)

      ALLOCATE ( ztmp(nx,ny,nz) )

      ztmp(:,:,:) = XF(:,:,:)

      DO ji = i_chg_x, nx
         XF(ji - i_chg_x + 1 , :,:) = ztmp(ji,:,:)
      END DO
      DO ji = 1, i_chg_x - 1
         XF(nx - i_chg_x + 1 + ji , :,:) = ztmp(ji,:,:)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE LONG_REORG_3D


   SUBROUTINE FIND_NEAREST_POINT_IDIOT(rlon, rlat, Xlon, Xlat, iloc, jloc)

      !!----------------------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_POINT_IDIOT  ***
      !!
      !!                Laurent Brodeau, october, 2008
      !!
      !!              Very lame and unefficient but works on all grids!
      !!----------------------------------------------------------------------------

      IMPLICIT NONE

      REAL(8),INTENT(in)                :: rlon, rlat  !: lon and lat of target point to locate
      REAL(8),DIMENSION(:,:),INTENT(in) :: Xlat, Xlon  !: model grid layout
      INTEGER,INTENT (out)              :: iloc, jloc  !: nearest point location

      INTEGER :: nx, ny, ji, jj
      REAL(8) :: zdist, min_dist

      !! This is extended arrays we have
      INTEGER, PARAMETER :: next = 2    ! on each side !LB


      nx = size(Xlat,1)
      ny = size(Xlat,2)

      iloc = -9 ; jloc = -9
      min_dist = 20000.

      DO jj = 1+next, ny-next
         DO ji = 1+next, nx-next

            zdist = distance(rlon, Xlon(ji,jj), rlat, Xlat(ji,jj))

            IF ( zdist < min_dist ) THEN
               min_dist = zdist
               iloc = ji
               jloc = jj
            END IF

         END DO
      END DO


      !! Bug for regular spherical coordinate a la ECMWF: PROBLEM !!!
      IF ( (ABS(Xlon(iloc,jloc) - rlon) > 175.).and.(ABS(Xlon(iloc,jloc) - rlon) < 185.) ) THEN
         CALL locate_point(rlon, rlat, Xlon, Xlat, iloc, jloc)
      END IF


      IF ( ( iloc == -9 ).OR.( iloc == -9 ) ) THEN
         PRINT *, 'ERROR in FIND_NEAREST_POINT_IDIOT of mod_bilin_2d.f90 !'
         PRINT *, 'Point rlon, rlat', rlon, rlat
         PRINT *, 'not found on the source grid!'
         STOP
      END IF

   END SUBROUTINE FIND_NEAREST_POINT_IDIOT




   SUBROUTINE FIND_NEAREST_POINT_REG(rlon, rlat, Xlon, Xlat, iloc, jloc)

      !!----------------------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_POINT_REG  ***
      !!
      !!                Laurent Brodeau, october, 2008
      !!
      !!             Case when source domain (Xlon,Xlat) is actually regular !
      !!
      !!              Very lame and unefficient but works on all grids!
      !!----------------------------------------------------------------------------

      IMPLICIT NONE

      REAL(8),INTENT(in)                :: rlon, rlat  !: lon and lat of target point to locate
      REAL(8),DIMENSION(:,:),INTENT(in) :: Xlat, Xlon  !: model grid layout
      INTEGER,INTENT (out)              :: iloc, jloc  !: nearest point location

      INTEGER :: nx, ny, ji, jj, ip, jp, i_rig, i_lef, j_bot, j_top
      INTEGER, PARAMETER :: iframe = 4
      REAL(8) :: zdist, min_dist
      REAL(8),DIMENSION(:), ALLOCATABLE :: vlat, vlon

      !! This is extended arrays we have
      INTEGER, PARAMETER :: next = 2    ! on each side !LB


      nx = size(Xlat,1)
      ny = size(Xlat,2)

      ALLOCATE ( vlon(nx) , vlat(ny) )

      vlon(:) = Xlon(:,3+next)   ! 3+next => make sure we're inside the domain
      vlat(:) = Xlat(3+next,:)   ! 3+next => make sure we're inside the domain


      !! Finding the rectangular region where to search:
      !! (rectangular region of roughly iframexiframe points...
      
      DO ji = next, nx-next
         IF ( (vlon(ji) <= rlon).AND.(vlon(ji+1) > rlon) ) THEN
            ip = ji
            EXIT
         END IF
      END DO
      
      IF ( vlat(3+next) > vlat(3+next+1) ) THEN
         PRINT *, 'ERROR: FIND_NEAREST_POINT_REG => lat seems to increase...' ; STOP
      END IF
      DO jj = next, ny-next
         IF ( (vlat(jj) <= rlat).AND.(vlat(jj+1) > rlat) ) THEN
            jp = jj
            EXIT
         END IF
      END DO

      i_lef = MAX(ip-(iframe-1) ,  1+next)
      i_rig = MIN(ip+iframe , nx-next)
      
      j_bot = MAX(jp-(iframe-1) ,  1+next)
      j_top = MIN(jp+iframe , ny-next)
      
      !PRINT *, ''
      !PRINT *, 'i_lef, ip, i_rig =>', i_lef, ip, i_rig
      !PRINT *, 'rlon, vlon(i_lef), vlon(i_rig) =>', rlon, vlon(i_lef), vlon(i_rig)
      !PRINT *, 'j_bot, jp, j_top =>', j_bot, jp, j_top
      !PRINT *, 'vlat(j_bot), rlat, vlat(j_top) =>', REAL(vlat(j_bot),4), REAL(rlat,4), REAL(vlat(j_top),4)

      iloc = -9
      jloc = -9
      min_dist = 20000.

      DO jj = j_bot, j_top
         DO ji = i_lef, i_rig

            zdist = distance(rlon, Xlon(ji,jj), rlat, Xlat(ji,jj))
            
            IF ( zdist < min_dist ) THEN
               min_dist = zdist
               iloc = ji
               jloc = jj
            END IF
            
         END DO
      END DO
      
      IF ( ( iloc == -9 ).OR.( iloc == -9 ) ) THEN
         PRINT *, 'ERROR in FIND_NEAREST_POINT_REG of mod_bilin_2d.f90 !'
         PRINT *, 'Point rlon, rlat', rlon, rlat
         PRINT *, 'not found on the source grid!'
         STOP
      END IF
      
      DEALLOCATE ( vlon , vlat )
      
   END SUBROUTINE FIND_NEAREST_POINT_REG


   FUNCTION DISTANCE(plona, plonb, plata, platb)

      !!----------------------------------------------------------
      !!           ***  FUNCTION  DIST  ***
      !!
      !!  ** Purpose : Compute the distance (km) between
      !!          point A (lona, lata) and B(lonb,latb)
      !!
      !!  ** Method : Compute the distance along the orthodromy
      !!
      !! * history : J.M. Molines in CHART, f90, may 2007
      !!----------------------------------------------------------

      IMPLICIT NONE
      ! Argument
      REAL(8), INTENT(in) :: plata, plona, platb, plonb
      REAL(8)             :: distance

      ! Local variables
      REAL(8),SAVE ::  zlatar, zlatbr, zlonar, zlonbr
      REAL(8) ::  zpds
      REAL(8),SAVE :: zux, zuy, zuz
      REAL(8) :: zvx, zvy, zvz

      REAL(8), SAVE :: prevlat=-1000., prevlon=-1000, zr, zpi, zconv


      !! Initialise some values at first call
      IF ( lfirst_dist ) THEN
         lfirst_dist = .FALSE.
         ! constants
         zpi = ACOS(-1._8)
         zconv = zpi/180.  ! for degree to radian conversion
         ! Earth radius
         zr = (6378.137+6356.7523)/2.0 ! km
      ENDIF

      !! Compute these term only if they differ from previous call
      IF ( plata /= prevlat .OR. plona /= prevlon) THEN
         zlatar=plata*zconv
         zlonar=plona*zconv
         zux=COS(zlonar)*COS(zlatar)
         zuy=SIN(zlonar)*COS(zlatar)
         zuz=SIN(zlatar)
         prevlat=plata
         prevlon=plona
      ENDIF

      zlatbr=platb*zconv
      zlonbr=plonb*zconv
      zvx=COS(zlonbr)*COS(zlatbr)
      zvy=SIN(zlonbr)*COS(zlatbr)
      zvz=SIN(zlatbr)

      zpds=zux*zvx+zuy*zvy+zuz*zvz

      IF (zpds >= 1.) THEN
         distance=0.
      ELSE
         distance=zr*ACOS(zpds)
      ENDIF

   END FUNCTION DISTANCE



   SUBROUTINE FIND_NEAREST_POINT_SMART(pplon, pplat, plam, pphi, kpiloc, kpjloc, ldbord)

      !!----------------------------------------------------------------------------
      !!            ***  SUBROUTINE NEAREST  ***
      !!
      !!   ** Purpose:  Computes the positions of the nearest i,j in the grid
      !!                from the given longitudes and latitudes
      !!
      !!   ** Method :  Starts on the middle of the grid, search in a 20x20 box, and move
      !!     the box in the direction where the distance between the box and the
      !!     point is minimum
      !!     Iterates ...
      !!     Stops when the point is outside the grid.
      !!     This algorithm does not work on the Mediteranean grid !
      !!
      !!   * history:
      !!        Anne de Miranda et Pierre-Antoine Darbon Jul. 2000 (CLIPPER)
      !!        Jean-Marc Molines : In NEMO form
      !!----------------------------------------------------------------------------

      IMPLICIT NONE

      !* arguments
      REAL(8), INTENT(in)       ::  pplon, pplat   !: lon and lat of target point
      INTEGER, INTENT (inout)    ::  kpiloc, kpjloc  !: nearest point location
      REAL(8), DIMENSION(:,:),INTENT(in) ::  pphi,plam  !: model grid layout
      LOGICAL                   :: ldbord         !: reach boundary flag

      ! * local variables
      INTEGER :: nx, ny, ji,jj,i0,j0,i1,j1
      INTEGER :: itbl
      REAL(8) ::  zdist,zdistmin,zdistmin0
      LOGICAL ::  lbordcell

      !lulu
      nx = size(pphi,1)
      ny = size(pphi,2)



      ! Initial values
      kpiloc = nx/2 ; kpjloc = ny/2    ! seek from the middle of domain
      itbl = 10                          ! block size for search
      zdistmin=1000000000. ; zdistmin0=1000000000.
      i0=kpiloc ;  j0=kpjloc
      lbordcell=.TRUE.;   ldbord=.FALSE.

      ! loop until found or boundary reach
      DO  WHILE ( lbordcell .AND. .NOT. ldbord)
         i0=kpiloc-itbl ;  i1=kpiloc+itbl
         j0=kpjloc-itbl ;  j1=kpjloc+itbl

         ! search only the inner domain
         IF (i0 <= 1) i0=2
         IF (i1 > nx-1) i1=nx-1
         IF (j0 <= 1) j0=2
         IF( j1 > ny-1) j1=ny-1

         ! within a block itbl+1 x itbl+1:
         DO jj=j0,j1
            DO ji=i0,i1
               ! compute true distance (orthodromy) between target point and grid point
               zdist=distance(pplon,plam(ji,jj),pplat,pphi(ji,jj) )
               zdistmin=MIN(zdistmin,zdist)
               ! update kpiloc, kpjloc if distance decreases
               IF (zdistmin /= zdistmin0 ) THEN
                  kpiloc=ji
                  kpjloc=jj
               ENDIF
               zdistmin0=zdistmin
            END DO
         END DO
         lbordcell=.FALSE.
         ! if kpiloc, kpjloc belong to block boundary proceed to next block, centered on kpiloc, kpjloc
         IF (kpiloc == i0 .OR. kpiloc == i1) lbordcell=.TRUE.
         IF (kpjloc == j0 .OR. kpjloc == j1) lbordcell=.TRUE.
         ! boundary reach ---> not found
         IF (kpiloc == 2  .OR. kpiloc ==nx-1) ldbord=.TRUE.
         IF (kpjloc == 2  .OR. kpjloc ==ny-1) ldbord=.TRUE.
      END DO
      !!
      !!
      !    PRINT *, ''; PRINT *, ''
      !    PRINT *, 'End of NEAREST!'
      !!
      !    PRINT *, 'iloc, jloc', kpiloc, kpjloc
      !
      !
      !!
   END SUBROUTINE  FIND_NEAREST_POINT_SMART



   SUBROUTINE locate_point(rlon_P, rlat_P, Xlon, Xlat, jxfnd, jyfnd)

      IMPLICIT none

      !!===========================================================
      !!
      !! Input :
      !! rlon_P, rlat_P : point to be located
      !! Xlon, Xlat : 2D arrays contening long. and lat.
      !!
      !! Output :
      !! jx, jy     : location
      !!
      !!
      !! L. Brodeau, 2005
      !!===========================================================

      REAL(8),                 INTENT(in)  :: rlon_P, rlat_P
      REAL(8), DIMENSION(:,:), INTENT(in)  :: Xlon, Xlat
      INTEGER,                 INTENT(out) :: jxfnd, jyfnd

      !! Local :
      REAL(8), PARAMETER :: dr = 30.    !! avoid problem at 0. -> 360. for lon
      REAL(8), DIMENSION(4) :: v4x, v4y
      INTEGER :: nx, ny, ji, jj
      REAL(8) :: min_x, max_x, min_y, max_y

      nx = size(Xlon,1)
      ny = size(Xlon,2)

      jxfnd = -1  ; jyfnd = -1

      DO jj=1, ny-1
         DO ji=1, nx-1

            v4x   = (/ Xlon(ji,jj) , Xlon(ji+1,jj) , Xlon(ji+1,jj+1) , Xlon(ji,jj+1) /)
            v4y   = (/ Xlat(ji,jj) , Xlat(ji+1,jj) , Xlat(ji+1,jj+1) , Xlat(ji,jj+1) /)

            min_x = MINVAL(v4x)  ;  max_x = MAXVAL(v4x)
            min_y = MINVAL(v4y)  ;  max_y = MAXVAL(v4y)

            IF ( ((rlon_P >= min_x).and.(rlon_P <= max_x)) .AND. ( (max_x-min_x) < dr)  ) THEN
               IF ( (rlat_P >= min_y).and.(rlat_P <= max_y) ) THEN
                  jxfnd = ji  ; jyfnd = jj
               END IF
            END IF

         END DO
      END DO

   END SUBROUTINE locate_point


END MODULE MOD_MANIP
