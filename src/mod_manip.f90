MODULE MOD_MANIP

   !! Misc. manipulations and operations on 2D arrays...

   !! Author: L. Brodeau


   IMPLICIT NONE

   PRIVATE

   PUBLIC :: fill_extra_bands, fill_extra_north_south, extra_2_east, extra_2_west, test_xyz, partial_deriv, &
      &      flip_ud_1d, flip_ud_1d_double, flip_ud_2d, flip_ud_3d, long_reorg_2d, long_reorg_3d, &
      &      distance, find_nearest_point

   REAL(8), PARAMETER, PUBLIC :: rflg = -9999.

   !LOGICAL, PARAMETER :: ldebug = .TRUE.
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

      USE mod_conf, ONLY: i_orca_in, i_orca_out

      INTEGER ,                INTENT(in)  :: k_ew
      REAL(8), DIMENSION(:,:), INTENT(in)  :: XX, YY, XF
      REAL(8), DIMENSION(:,:), INTENT(out) :: XP4, YP4, FP4

      !! Local
      INTEGER :: nx, ny, nxp4, nyp4
      INTEGER :: ji, jj

      IF ( (SIZE(XX,1) /= SIZE(YY,1)).OR.(SIZE(XX,2) /= SIZE(YY,2)).OR. &
         & (SIZE(XX,1) /= SIZE(XF,1)).OR.(SIZE(XX,2) /= SIZE(XF,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_BANDS : size of input coor. and data do not match!!!'; STOP
      END IF

      IF ( (SIZE(XP4,1) /= SIZE(YP4,1)).OR.(SIZE(XP4,2) /= SIZE(YP4,2)).OR. &
         & (SIZE(XP4,1) /= SIZE(FP4,1)).OR.(SIZE(XP4,2) /= SIZE(FP4,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_BANDS : size of output coor. and data do not match!!!'; STOP
      END IF

      nx = SIZE(XX,1)
      ny = SIZE(XX,2)

      nxp4 = SIZE(XP4,1)
      nyp4 = SIZE(XP4,2)

      IF ( nxp4 /= nx + 4 ) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_BANDS : target x dim is not ni+4!!!'; STOP
      END IF
      IF ( nyp4 /= ny + 4 ) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_BANDS : target y dim is not nj+4!!!'; STOP
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

      IF (k_ew > -1) THEN   ! we can use east-west periodicity of input file to
         !!                   ! fill extra bands :
         XP4( 1     , 3:nyp4-2) = XX(nx - 1 - k_ew , :) - 360.   ! lolo: use or not the 360 stuff???
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
      SELECT CASE( i_orca_in )
         !
      CASE (4)
         XP4(2:nxp4/2             ,nyp4-1) = XP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-5)
         XP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-1) = XP4(2:nxp4/2             ,nyp4-5)
         XP4(2:nxp4/2             ,nyp4)   = XP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-6)
         XP4(nxp4:nxp4-nxp4/2-2:-1,nyp4)   = XP4(2:nxp4/2             ,nyp4-6)
      CASE (6)
         XP4(2:nxp4/2               ,nyp4-1) = XP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-4)
         XP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-1) = XP4(2:nxp4/2               ,nyp4-4)
         XP4(2:nxp4/2               ,nyp4)   = XP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-5)
         XP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4)   = XP4(2:nxp4/2               ,nyp4-5)
      CASE DEFAULT
         XP4(:,nyp4-1) = XP4(:,nyp4-3) + XP4(:,nyp4-2) - XP4(:,nyp4-4)
         XP4(:,nyp4)   = XP4(:,nyp4-2) + XP4(:,nyp4-2) - XP4(:,nyp4-4)

      END SELECT


      !! Y array :
      !! ---------

      !! Top side :
      SELECT CASE( i_orca_in )

      CASE (4)
         YP4(2:nxp4/2             ,nyp4-1) = YP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-5)
         YP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-1) = YP4(2:nxp4/2             ,nyp4-5)
         YP4(2:nxp4/2             ,nyp4)   = YP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-6)
         YP4(nxp4:nxp4-nxp4/2-2:-1,nyp4)   = YP4(2:nxp4/2             ,nyp4-6)
      CASE (6)
         YP4(2:nxp4/2               ,nyp4-1) = YP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-4)
         YP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-1) = YP4(2:nxp4/2               ,nyp4-4)
         YP4(2:nxp4/2               ,nyp4)   = YP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-5)
         YP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4)   = YP4(2:nxp4/2               ,nyp4-5)
      CASE DEFAULT
         YP4(3:nxp4-2, nyp4-1) = YY(:, ny-1) + YY(:,ny) - YY(:,ny-2)
         YP4(3:nxp4-2, nyp4)   = YY(:, ny)   + YY(:,ny) - YY(:,ny-2)
      END SELECT


      !! Bottom side :
      YP4(3:nxp4-2, 2) = YY(:,2) - (YY(:,3) - YY(:,1))
      YP4(3:nxp4-2, 1) = YY(:,1) - (YY(:,3) - YY(:,1))

      IF (k_ew > -1) THEN

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

      IF (k_ew > -1) THEN

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
      SELECT CASE( i_orca_in )

      CASE (4)
         PRINT *, 'ORCA north pole T-point folding type of extrapolation at northern boundary!'
         FP4(2:nxp4/2             ,nyp4-1) = FP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-5)
         FP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-1) = FP4(2:nxp4/2             ,nyp4-5)
         FP4(2:nxp4/2             ,nyp4)   = FP4(nxp4:nxp4-nxp4/2-2:-1,nyp4-6)
         FP4(nxp4:nxp4-nxp4/2-2:-1,nyp4)   = FP4(2:nxp4/2             ,nyp4-6)
      CASE (6)
         PRINT *, 'ORCA north pole F-point folding type of extrapolation at northern boundary!'         
         FP4(2:nxp4/2               ,nyp4-1) = FP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-4)
         FP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-1) = FP4(2:nxp4/2               ,nyp4-4)
         FP4(2:nxp4/2               ,nyp4)   = FP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4-5)
         FP4(nxp4-1:nxp4-nxp4/2+1:-1,nyp4)   = FP4(2:nxp4/2               ,nyp4-5)

      CASE DEFAULT
         DO ji = 1, nxp4
            CALL extra_2_east(YP4(ji,nyp4-4),YP4(ji,nyp4-3),YP4(ji,nyp4-2),       &
               &              YP4(ji,nyp4-1),YP4(ji,nyp4),                        &
               &              FP4(ji,nyp4-4),FP4(ji,nyp4-3),FP4(ji,nyp4-2), &
               &              FP4(ji,nyp4-1),FP4(ji,nyp4) )
         END DO
      END SELECT

      !! Bottom side :
      DO ji = 1, nxp4
         CALL extra_2_west(YP4(ji,5),YP4(ji,4),YP4(ji,3),       &
            &              YP4(ji,2),YP4(ji,1),                 &
            &              FP4(ji,5),FP4(ji,4),FP4(ji,3), &
            &              FP4(ji,2),FP4(ji,1) )
      END DO

   END SUBROUTINE FILL_EXTRA_BANDS





   SUBROUTINE FILL_EXTRA_NORTH_SOUTH(XX, YY, XF, XP4, YP4, FP4)

      !!============================================================================
      !! Extending input arrays with an extraband of two points at northern and
      !! southern boundaries.
      !!
      !! The extension is done thanks to Akima's exptrapolation method or continuity
      !!  / ORCA knowledge...
      !!
      !!============================================================================

      USE mod_conf, ONLY: i_orca_in

      REAL(8), DIMENSION(:,:), INTENT(in)  :: XX, YY, XF
      REAL(8), DIMENSION(:,:), INTENT(out) :: XP4, YP4, FP4

      !! Local
      INTEGER :: nx, ny, nxp4, nyp4
      INTEGER :: ji

      IF ( (SIZE(XX,1) /= SIZE(YY,1)).OR.(SIZE(XX,2) /= SIZE(YY,2)).OR. &
         & (SIZE(XX,1) /= SIZE(XF,1)).OR.(SIZE(XX,2) /= SIZE(XF,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_NORTH_SOUTH : size of input coor. and data do not match!!!'; STOP
      END IF

      IF ( (SIZE(XP4,1) /= SIZE(YP4,1)).OR.(SIZE(XP4,2) /= SIZE(YP4,2)).OR. &
         & (SIZE(XP4,1) /= SIZE(FP4,1)).OR.(SIZE(XP4,2) /= SIZE(FP4,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_NORTH_SOUTH : size of output coor. and data do not match!!!'; STOP
      END IF

      nx = SIZE(XX,1)
      ny = SIZE(XX,2)

      nxp4 = SIZE(XP4,1)
      nyp4 = SIZE(XP4,2)

      IF ( nxp4 /= nx ) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_NORTH_SOUTH : target x dim is not ni!!!'; STOP
      END IF
      IF ( nyp4 /= ny + 4 ) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_NORTH_SOUTH : target y dim is not nj+4!!!'; STOP
      END IF


      !!   C r e a t i n g   e x t e n d e d   a r r a y s  :
      !!   --------------------------------------------------

      !! Initializing :
      XP4 = 0.
      YP4 = 0.
      FP4 = 0.

      !! Filling center of domain:
      XP4(:, 3:nyp4-2) = XX(:,:)
      YP4(:, 3:nyp4-2) = YY(:,:)
      FP4(:, 3:nyp4-2) = XF(:,:)



      !! X array :
      !! ---------

      !! Bottom side :
      XP4(:, 2) = XP4(:,4) - (XP4(:,5) - XP4(:,3))
      XP4(:, 1) = XP4(:,3) - (XP4(:,5) - XP4(:,3))

      !! Top side :
      SELECT CASE( i_orca_in )
         !
      CASE (4)
         XP4(2:nx/2             ,nyp4-1) = XP4(nx:nx-nx/2-2:-1,nyp4-5)
         XP4(nx:nx-nx/2-2:-1,nyp4-1) = XP4(2:nx/2             ,nyp4-5)
         XP4(2:nx/2             ,nyp4)   = XP4(nx:nx-nx/2-2:-1,nyp4-6)
         XP4(nx:nx-nx/2-2:-1,nyp4)   = XP4(2:nx/2             ,nyp4-6)
      CASE (6)
         XP4(2:nx/2               ,nyp4-1) = XP4(nx-1:nx-nx/2+1:-1,nyp4-4)
         XP4(nx-1:nx-nx/2+1:-1,nyp4-1) = XP4(2:nx/2               ,nyp4-4)
         XP4(2:nx/2               ,nyp4)   = XP4(nx-1:nx-nx/2+1:-1,nyp4-5)
         XP4(nx-1:nx-nx/2+1:-1,nyp4)   = XP4(2:nx/2               ,nyp4-5)
      CASE DEFAULT
         XP4(:,nyp4-1) = XP4(:,nyp4-3) + XP4(:,nyp4-2) - XP4(:,nyp4-4)
         XP4(:,nyp4)   = XP4(:,nyp4-2) + XP4(:,nyp4-2) - XP4(:,nyp4-4)

      END SELECT


      !! Y array :
      !! ---------

      !! Top side :
      SELECT CASE( i_orca_in )

      CASE (4)
         YP4(2:nx/2             ,nyp4-1) = YP4(nx:nx-nx/2-2:-1,nyp4-5)
         YP4(nx:nx-nx/2-2:-1,nyp4-1) = YP4(2:nx/2             ,nyp4-5)
         YP4(2:nx/2             ,nyp4)   = YP4(nx:nx-nx/2-2:-1,nyp4-6)
         YP4(nx:nx-nx/2-2:-1,nyp4)   = YP4(2:nx/2             ,nyp4-6)
      CASE (6)
         YP4(2:nx/2               ,nyp4-1) = YP4(nx-1:nx-nx/2+1:-1,nyp4-4)
         YP4(nx-1:nx-nx/2+1:-1,nyp4-1) = YP4(2:nx/2               ,nyp4-4)
         YP4(2:nx/2               ,nyp4)   = YP4(nx-1:nx-nx/2+1:-1,nyp4-5)
         YP4(nx-1:nx-nx/2+1:-1,nyp4)   = YP4(2:nx/2               ,nyp4-5)
      CASE DEFAULT
         YP4(:, nyp4-1) = YY(:, ny-1) + YY(:,ny) - YY(:,ny-2)
         YP4(:, nyp4)   = YY(:, ny)   + YY(:,ny) - YY(:,ny-2)
      END SELECT


      !! Bottom side :
      YP4(:, 2) = YY(:,2) - (YY(:,3) - YY(:,1))
      YP4(:, 1) = YY(:,1) - (YY(:,3) - YY(:,1))



      !! Data array :
      !! ------------

      !! Top side :
      SELECT CASE( i_orca_in )

      CASE (4)
         PRINT *, 'ORCA2 type of extrapolation at northern boundary!'
         FP4(2:nx/2             ,nyp4-1) = FP4(nx:nx-nx/2-2:-1,nyp4-5)
         FP4(nx:nx-nx/2-2:-1,nyp4-1) = FP4(2:nx/2             ,nyp4-5)
         FP4(2:nx/2             ,nyp4)   = FP4(nx:nx-nx/2-2:-1,nyp4-6)
         FP4(nx:nx-nx/2-2:-1,nyp4)   = FP4(2:nx/2             ,nyp4-6)
      CASE (6)
         PRINT *, 'ORCA1 type of extrapolation at northern boundary!'
         FP4(2:nx/2               ,nyp4-1) = FP4(nx-1:nx-nx/2+1:-1,nyp4-4)
         FP4(nx-1:nx-nx/2+1:-1,nyp4-1) = FP4(2:nx/2               ,nyp4-4)
         FP4(2:nx/2               ,nyp4)   = FP4(nx-1:nx-nx/2+1:-1,nyp4-5)
         FP4(nx-1:nx-nx/2+1:-1,nyp4)   = FP4(2:nx/2               ,nyp4-5)
      CASE DEFAULT
         DO ji = 1, nx
            CALL extra_2_east(YP4(ji,nyp4-4),YP4(ji,nyp4-3),YP4(ji,nyp4-2),       &
               &              YP4(ji,nyp4-1),YP4(ji,nyp4),                        &
               &              FP4(ji,nyp4-4),FP4(ji,nyp4-3),FP4(ji,nyp4-2), &
               &              FP4(ji,nyp4-1),FP4(ji,nyp4) )
         END DO
      END SELECT

      !! Bottom side :
      DO ji = 1, nx
         CALL extra_2_west(YP4(ji,5),YP4(ji,4),YP4(ji,3),       &
            &              YP4(ji,2),YP4(ji,1),                 &
            &              FP4(ji,5),FP4(ji,4),FP4(ji,3), &
            &              FP4(ji,2),FP4(ji,1) )
      END DO

   END SUBROUTINE FILL_EXTRA_NORTH_SOUTH



















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
            PRINT *, 'ERROR, mod_manip.f90 = >TEST_XYZ 1: longitude and latitude array do not match data!'
            PRINT *, ''; STOP
         END IF
         !!
      ELSE
         IF ( (ix1 == iz1).AND.(iy1 == iz1).AND.(ix2 == iz2).AND.(iy2 == iz2) ) THEN
            TEST_XYZ = '2d'
         ELSE
            PRINT *, 'ERROR, mod_manip.f90 = >TEST_XYZ 2: longitude and latitude array do not match data!'
            PRINT *, ''; STOP
         END IF
      END IF
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

      DEALLOCATE ( ZX, ZY, ZF )

   END SUBROUTINE PARTIAL_DERIV


   SUBROUTINE FLIP_UD_1D(XF)

      REAL(4), DIMENSION(:), INTENT(inout) :: XF

      INTEGER :: nz, jk
      REAL(4), DIMENSION(:), ALLOCATABLE :: ztmp

      nz = SIZE(XF,1) 

      ALLOCATE ( ztmp(nz) )

      ztmp(:) = XF(:)

      DO jk = 1, nz
         XF(jk) =  ztmp(nz-jk+1)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE FLIP_UD_1D

   SUBROUTINE FLIP_UD_1D_DOUBLE(XF)

      REAL(8), DIMENSION(:), INTENT(inout) :: XF

      INTEGER :: nz, jk
      REAL(8), DIMENSION(:), ALLOCATABLE :: ztmp

      nz = SIZE(XF,1)

      ALLOCATE ( ztmp(nz) )

      ztmp(:) = XF(:)

      DO jk = 1, nz
         XF(jk) =  ztmp(nz-jk+1)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE FLIP_UD_1D_DOUBLE

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






   SUBROUTINE FIND_NEAREST_POINT(Xout, Yout, Xin, Yin, JIpos, JJpos)
      !!
      !!---------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_POINT  ***
      !!
      !!   Source domain: Xin, Yin
      !!   Target domain: Xout, Yout
      !!     => for each point of target domain: i_in,j_in location of nearest point on source domain
      !!        => JIpos & JJpos have same shape as Xout & Yout !
      !!
      !!                Laurent Brodeau, July, 2017
      !!
      !!---------------------------------------------------------------

      !debug: USE io_ezcdf
      USE io_ezcdf
      USE mod_conf,  ONLY: i_orca_out

      IMPLICIT NONE

      INTEGER, PARAMETER :: &
         &                   Nlat_split = 30, &  ! number of latitude bands to split the search work
         &                   nframe_scan = 4  ! domain to scan for nearest point in simple algo => domain of 9x9

      REAL(8), DIMENSION(:,:), INTENT(in)  :: Xout, Yout    !: lon and lat arrays of target domain
      REAL(8), DIMENSION(:,:), INTENT(in)  :: Xin , Yin     !: lon and lat arrays of source domain
      INTEGER, DIMENSION(:,:), INTENT(out) :: JIpos, JJpos  !: nearest point location of point P in Xin,Yin wrt Xout,Yout

      INTEGER :: &
         &    ji, jj, &
         &    nx_in, ny_in, nx_out, ny_out, &
         &    jlat, ji_out, jj_out, ji_in, jj_in, jj_out_old, &
         &    jmin_in, jmax_in, imin_in, imax_in, niter, &
         &    j_strt_out, j_stop_out

      REAL(8) :: rmin_dlat_dj, emax, frac_emax, rlat_low, rlat_hgh, rlon, rlat, rlat_old, zdist, rtmp

      REAL(8) :: y_max_out, y_min_out, y_max_bnd, y_min_bnd, y_max_bnd0, y_min_bnd0, dy, y_max_in, y_min_in
      REAL(8), DIMENSION(:),   ALLOCATABLE :: VLAT_SPLIT_BOUNDS
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: IJ_VLAT_IN

      INTEGER :: jlat_inc    ! 1 if lat increases with j, -1 if decreases with j
      INTEGER, DIMENSION(2) :: jmax_loc, jmin_loc
      INTEGER, DIMENSION(1) :: ip, jp

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xdist, &
         &                                    e1_in, e2_in, &    !: grid layout and metrics
         &                                    ztmp_out
      REAL(8),DIMENSION(:), ALLOCATABLE :: vlat, vlon
      LOGICAL :: l_is_reg_in, l_is_reg_out, lagain


      nx_in  = SIZE(Xin,1)
      ny_in  = SIZE(Xin,2)
      nx_out = SIZE(Xout,1)
      ny_out = SIZE(Xout,2)

      PRINT *, ' Source domain size: ', nx_in, ny_in
      PRINT *, ' Target domain size: ', nx_out, ny_out

      IF ( (SIZE(Yin,1) /= nx_in) .OR. (SIZE(Yin,2) /= ny_in) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Yin dont agree in shape with Xin'
         STOP
      END IF
      IF ( (SIZE(Yout,1) /= nx_out) .OR. (SIZE(Yout,2) /= ny_out) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Yout dont agree in shape with Xout'
         STOP
      END IF
      IF ( (SIZE(JIpos,1) /= nx_out) .OR. (SIZE(JIpos,2) /= ny_out) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): JIpos dont agree in shape with Xout'
         PRINT *, SIZE(JIpos,1), SIZE(JIpos,2), 'vs', nx_out, ny_out
         STOP
      END IF
      IF ( (SIZE(JJpos,1) /= nx_out) .OR. (SIZE(JJpos,2) /= ny_out) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): JJpos dont agree in shape with Xout'
         STOP
      END IF

      !! Checking if source domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_in = L_IS_GRID_REGULAR( Xin , Yin )
      PRINT *, ' *** FIND_NEAREST_POINT => Is source grid regular ??? =>', l_is_reg_in

      !! Checking if target domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_out = L_IS_GRID_REGULAR( Xout , Yout )
      PRINT *, ' *** FIND_NEAREST_POINT => Is target grid regular ??? =>', l_is_reg_out


      ALLOCATE ( Xdist(nx_in,ny_in) )


      IF (ldebug) THEN
         !! Testing the field "distance from point at lon=2./lat=45."
         Xdist(:,:) = DISTANCE_2D(2.0_8, Xin(:,:), 45.0_8, Yin(:,:))
         CALL PRTMASK(REAL(Xdist,4), 'distance_fron_2_45_in.nc', 'dist')
      END IF
      
      !! Going to scan target grid through increasing  (or decreasing) j (latitude)


      
      !! General case:
      j_strt_out = 1
      j_stop_out = ny_out
      jlat_inc   = 1

      
      !! We need to know if the target latitude ONLY keeps on systematically
      !! increasing (or decreasing) as j increases:
      !!
      ALLOCATE ( ztmp_out(nx_out, ny_out) )
      DO jj = 2, ny_out
         ztmp_out(:,jj) = Yout(:,jj) - Yout(:,jj-1)
      END DO

      y_min_in = MINVAL(Yin) ; ! Min and Max latitude of source domain
      y_max_in = MAXVAL(Yin)

      y_min_out = MINVAL(Yout) ; ! Min and Max latitude of target domain
      y_max_out = MAXVAL(Yout)

      IF ( ldebug ) THEN
         PRINT *, ' Min. latitude on target & source domains =>', y_min_out, y_min_in
         PRINT *, ' Max. latitude on target & source domains =>', y_max_out, y_max_in
      END IF


      
      rtmp = SUM(ztmp_out(:,ny_out/2)) ! know if increasing (>0) or decreasing (<0)
      ztmp_out = SIGN(1.0_8 , rtmp)*ztmp_out
      !IF (ldebug) CALL PRTMASK(REAL(ztmp_out,4), 'dlat_dj_out.nc', 'dist')
      rmin_dlat_dj = MINVAL(ztmp_out(:,2:))
      IF (ldebug) PRINT *, ' Minimum dlat_dj_out =>', rmin_dlat_dj



      

      !!  Simplif when [d lat / d j] always has the same sign:
      IF ( (rmin_dlat_dj >= 0.0_8) .OR. l_is_reg_out .OR. (i_orca_out > 0) ) THEN
         !! -> because we need to avoid all the following if target grid is for
         !!    example a polar sterographic projection... (example 5)
         !!
         !! *** Will ignore regions of the TARGET domain that
         !!     are not covered by source domain:
         jmin_loc = MINLOC(Yout, mask=(Yout>=y_min_in))
         jmax_loc = MAXLOC(Yout, mask=(Yout<=y_max_in))
         j_strt_out = jmin_loc(2)  ! smallest j on target source that covers smallest source latitude
         j_stop_out = jmax_loc(2)  ! largest j on target source that covers largest source latitude
         IF ( j_strt_out > j_stop_out ) jlat_inc = -1 ! latitude decreases as j increases (like ECMWF grids...)
         IF (ldebug) THEN
            PRINT *, ' j_strt_out, j_stop_out / nj_out =>', j_strt_out, j_stop_out, '/', ny_out
            PRINT *, ''
         END IF
      END IF ! IF ( (rmin_dlat_dj >= 0.0_8) .OR. l_is_reg_out .OR. (i_orca_out > 0) )


      !! Backround value to spot non-treated regions
      JIpos(:,:) = INT(rflg)
      JJpos(:,:) = INT(rflg)


      !! ---------------------------------------------------------------------------------------
      IF ( l_is_reg_in ) THEN

         PRINT *, '                        => going for simple algorithm !'

         ALLOCATE ( vlon(nx_in) , vlat(ny_in) )

         vlon(:) = Xin(:,3)   ! 3 => make sure we're inside the domain
         vlat(:) = Yin(3,:)   ! 3 => make sure we're inside the domain

         !IF ( vlat(3) > vlat(3+1) ) THEN
         !   PRINT *, 'ERROR: FIND_NEAREST_POINT => lat seems to increase...' ; STOP
         !END IF

         DO jj_out = j_strt_out, j_stop_out, jlat_inc

            IF ( MOD(jj_out,10)==0 ) PRINT *, ' *** Treated j-point of target domain =', jj_out !REAL(rlat,4)

            DO ji_out = 1, nx_out

               rlon = Xout(ji_out,jj_out)
               rlat = Yout(ji_out,jj_out)

               ! Nearest point (in terms of index):
               ip =  MINLOC(ABS(vlon(:)-rlon))
               jp =  MINLOC(ABS(vlat(:)-rlat))

               !! Define the box to scan for shortest distance:
               imin_in = MAX(ip(1)-nframe_scan ,  1)
               imax_in = MIN(ip(1)+nframe_scan , nx_in)
               jmin_in = MAX(jp(1)-nframe_scan ,  1)
               jmax_in = MIN(jp(1)+nframe_scan , ny_in)

               ! Nearest point (in terms of distance):
               Xdist = 1.E12
               Xdist(imin_in:imax_in,jmin_in:jmax_in) = DISTANCE_2D(rlon, Xin(imin_in:imax_in,jmin_in:jmax_in), rlat, Yin(imin_in:imax_in,jmin_in:jmax_in))
               jmin_loc = MINLOC(Xdist(imin_in:imax_in,jmin_in:jmax_in))
               ji_in = jmin_loc(1) + imin_in - 1
               jj_in = jmin_loc(2) + jmin_in - 1
               JIpos(ji_out,jj_out) = ji_in
               JJpos(ji_out,jj_out) = jj_in
               IF ((ji_in==0).OR.(jj_in==0)) THEN
                  PRINT *, ''
                  PRINT *, 'The nearest point was not found!'
                  PRINT *, ' !!! ji_in or jj_in = 0 !!!'
                  PRINT *, ' ** Target point (lon,lat) =>', rlon, rlat
                  PRINT *, ' Nearest point found on source grid:', &
                     &  REAL(Xin(ji_in,jj_in) , 4), &
                     &  REAL(Yin(ji_in,jj_in) , 4)
                  STOP
               END IF

               IF ( ( JIpos(ji_out,jj_out) == INT(rflg) ).OR.( JJpos(ji_out,jj_out) == INT(rflg) ) ) THEN
                  PRINT *, 'ERROR in FIND_NEAREST_POINT of mod_bilin_2d.f90 !'
                  PRINT *, 'Point rlon, rlat', rlon, rlat
                  PRINT *, 'not found on the source grid!'
                  STOP
               END IF

            END DO
         END DO

         DEALLOCATE ( vlon , vlat )


         !! ---------------------------------------------------------------------------------------
      ELSE

         
         !! IRREGULAR CASE !!
         PRINT *, '                        => going for advanced algorithm !'

         ALLOCATE ( e1_in(nx_in,ny_in), e2_in(nx_in,ny_in) )
         !! We need metric of input grid
         e1_in(:,:) = 40000. ;  e2_in(:,:) = 40000.
         DO jj_in=1, ny_in
            DO ji_in=1, nx_in-1
               e1_in(ji_in,jj_in) = distance(Xin(ji_in,jj_in),Xin(ji_in+1,jj_in),Yin(ji_in,jj_in),Yin(ji_in+1,jj_in))*1000. !(m)
            END DO
         END DO
         DO jj_in=1, ny_in-1
            DO ji_in=1, nx_in
               e2_in(ji_in,jj_in) = distance(Xin(ji_in,jj_in),Xin(ji_in,jj_in+1),Yin(ji_in,jj_in),Yin(ji_in,jj_in+1))*1000. !(m)
            END DO
         END DO
         IF (nx_in>1) e1_in(nx_in,:) = e1_in(nx_in-1,:)
         IF (ny_in>1) e2_in(:,ny_in) = e2_in(:,ny_in-1)
         IF (ldebug) THEN
            CALL PRTMASK(REAL(e1_in,4), 'e1_in.nc', 'e1')
            CALL PRTMASK(REAL(e2_in,4), 'e2_in.nc', 'e2')
         END IF

         !! Min and Max latitude to use for binning :
         y_max_bnd = MIN( y_max_out , y_max_in )
         y_max_bnd0 = y_max_bnd
         y_min_bnd = MAX( y_min_out , y_min_in )
         y_min_bnd0 = y_min_bnd
         !PRINT *, ' y_max_bnd #1 => ', y_max_bnd         
         !PRINT *, ' y_min_bnd #1 => ', y_min_bnd
         y_max_bnd = MIN( REAL(INT(y_max_bnd+1),8) ,  90.)
         y_min_bnd = MAX( REAL(INT(y_min_bnd-1),8) , -90.)
         !PRINT *, ' y_max_bnd #2 => ', y_max_bnd         
         !PRINT *, ' y_min_bnd #2 => ', y_min_bnd
         !! Multiple of 0.5:
         y_max_bnd = MIN( NINT(y_max_bnd/0.5)*0.5    ,  90.)
         y_min_bnd = MAX( NINT(y_min_bnd/0.5)*0.5    , -90.)
         !PRINT *, ' y_max_bnd #3 => ', y_max_bnd         
         !PRINT *, ' y_min_bnd #3 => ', y_min_bnd
         
         ALLOCATE ( VLAT_SPLIT_BOUNDS(Nlat_split+1), IJ_VLAT_IN(Nlat_split,2) )
         
         dy = (y_max_bnd - y_min_bnd)/Nlat_split         
         DO jlat=1,Nlat_split+1
            VLAT_SPLIT_BOUNDS(jlat) = y_min_bnd + REAL(jlat-1)*dy
         END DO

         IF ( ldebug ) THEN
            PRINT *, ''
            PRINT *, ' *** Binning between y_min_bnd, y_max_bnd, dy => ', REAL(y_min_bnd,4), REAL(y_max_bnd,4), REAL(dy,4)
            PRINT *, '      * real natural bound =>', REAL(y_min_bnd0,4), REAL(y_max_bnd0,4)
            PRINT *, '     => VLAT_SPLIT_BOUNDS ='
            PRINT *, VLAT_SPLIT_BOUNDS
            PRINT *, ''
         END IF
         
         IF ( (y_min_bnd > y_min_bnd0).OR.(y_max_bnd < y_max_bnd0) ) THEN
            PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Bounds for latitude for VLAT_SPLIT_BOUNDS are bad!'
            PRINT *, ' y_min_bnd, y_max_bnd =', y_min_bnd, y_max_bnd
            PRINT *, ' y_min_bnd0, y_max_bnd0 =', y_min_bnd0, y_max_bnd0
            STOP
         END IF
         
         DO jlat = 1, Nlat_split
            rlat_low = VLAT_SPLIT_BOUNDS(jlat)
            rlat_hgh = VLAT_SPLIT_BOUNDS(jlat+1)
            jmax_loc = MAXLOC(Yin, mask=(Yin<=rlat_hgh))
            jmin_loc = MINLOC(Yin, mask=(Yin>=rlat_low))
            !!
            !! To be sure to include everything, adding 2 extra points below and above:
            IJ_VLAT_IN(jlat,1) = MAX(jmin_loc(2) - 2,   1  )
            IJ_VLAT_IN(jlat,2) = MIN(jmax_loc(2) + 2, ny_in)
            !!
            IF ( ldebug ) THEN
               PRINT *, ' Latitude bin #', jlat
               PRINT *, '     => lat_low, lat_high:', REAL(rlat_low,4), REAL(rlat_hgh,4)
               PRINT *, '     => JJ min and max on input domain =>', IJ_VLAT_IN(jlat,1), IJ_VLAT_IN(jlat,2)
               PRINT *, ''
            END IF
            !!
         END DO

         rlat_old   = rflg
         jj_out_old = -10
         
         DO jj_out = j_strt_out, j_stop_out, jlat_inc
            DO ji_out = 1, nx_out
               
               rlon = Xout(ji_out,jj_out)
               rlat = Yout(ji_out,jj_out)

               !! Display progression in stdout:
               IF ( (ji_out == nx_out/2).AND.(jj_out /= jj_out_old) ) THEN
                  WRITE(*,'("*** Treated latitude of target domain = ",f7.4," (jj_out = ",i5.5,")")') REAL(rlat,4), jj_out
                  jj_out_old = jj_out
               END IF
               
               !! Need to find which jlat of our latitude bins rlat is located in!
               IF ( rlat /= rlat_old ) THEN
                  !IF ( jj_out /= jj_out_old ) WRITE(*,'("*** Treated latitude of target domain = ",f7.4," (jj_out = ",i5.5,")")') REAL(rlat,4), jj_out
                  DO jlat=1,Nlat_split
                     IF (  rlat ==VLAT_SPLIT_BOUNDS(jlat)) EXIT
                     IF ( (rlat > VLAT_SPLIT_BOUNDS(jlat)).AND.(rlat <= VLAT_SPLIT_BOUNDS(jlat+1)) ) EXIT
                  END DO
                  !!
               END IF

               lagain    = .TRUE.
               niter     = 0
               frac_emax = 0.5

               DO WHILE ( lagain )
                  !
                  !! Using band + niter surrounding:
                  !PRINT *, ' jlat, niter =>', jlat, niter 
                  jmin_in = IJ_VLAT_IN(MAX(jlat-niter,1)         , 1)
                  jmax_in = IJ_VLAT_IN(MIN(jlat+niter,Nlat_split), 2)
                  !!
                  IF ( ldebug ) THEN
                     PRINT *, ' *** Treated latitude of target domain =', REAL(rlat,4)
                     PRINT *, '     => bin #', jlat
                     PRINT *, '       => jmin & jmax on source domain =', jmin_in, jmax_in
                  END IF

                  Xdist = 1.E12
                  Xdist(1:nx_in,jmin_in:jmax_in) = DISTANCE_2D(rlon, Xin(1:nx_in,jmin_in:jmax_in), rlat, Yin(1:nx_in,jmin_in:jmax_in))

                  !CALL PRTMASK(REAL(Xdist,4), 'distance_last.nc', 'dist')

                  !! Nearest point is where distance is smallest:
                  jmin_loc = MINLOC(Xdist(1:nx_in,jmin_in:jmax_in))
                  ji_in = jmin_loc(1)
                  jj_in = jmin_loc(2) + jmin_in - 1
                  JIpos(ji_out,jj_out) = ji_in
                  JJpos(ji_out,jj_out) = jj_in

                  IF ((ji_in==0).OR.(jj_in==0)) THEN
                     PRINT *, ''
                     PRINT *, 'The nearest point was not found!'
                     PRINT *, ' !!! ji_in or jj_in = 0 !!!'
                     PRINT *, ' ** Target point (lon,lat) =>', rlon, rlat
                     STOP
                  END IF


                  zdist = Xdist(ji_in,jj_in) ! minimum distance found

                  emax = MAX(e1_in(ji_in,jj_in),e2_in(ji_in,jj_in))/1000.*SQRT(2.)

                  IF (zdist <= frac_emax*emax) THEN

                     lagain = .FALSE.

                  ELSE
                     IF (ldebug) THEN
                        PRINT *, ''
                        PRINT *, ' ** Target point (lon,lat) =>', rlon, rlat
                        PRINT *, ' Nearest point found on source grid:', &
                           &  REAL(Xin(ji_in,jj_in) , 4), &
                           &  REAL(Yin(ji_in,jj_in) , 4)
                        PRINT *, ''
                        PRINT *, 'The nearest point was not found! zdist / frac_emax*emax ', zdist, frac_emax*emax
                        PRINT *, 'Xin(1:6,1) =>', Xin(1:6,1) ; PRINT *, ''
                        PRINT *, 'Xin(nx_in-6:nx_in,1) =>', Xin(nx_in-6:nx_in,1) ; PRINT *, ''
                        PRINT *, ' jmin_in, jmax_in =>', jmin_in, jmax_in
                     END IF

                     IF (niter > Nlat_split/3) THEN
                        !! => increasing max. distance
                        niter = 1
                        !lolo:frac_emax = 1.5*frac_emax
                        frac_emax = 1.25*frac_emax    !lolo
                        !lolo:IF ( frac_emax > 2. ) THEN
                        IF ( frac_emax > 10. ) THEN
                           PRINT *, ' *** WARNING: mod_manip.f90/FIND_NEAREST_POINT: Giving up!!!'
                           PRINT *, '     => did not find nearest point for target coordinates:', &
                              &              REAL(rlon,4), REAL(rlat,4)
                           PRINT *, '     => last tested frac_emax was:', frac_emax
                           PRINT *, ''
                           lagain = .FALSE.
                           !! JIpos(ji_out,jj_out) & JJpos(ji_out,jj_out) should normally contain INT(rflg)
                           !! due to initialization with this value earlier...
                        ELSE
                           IF (ldebug) PRINT *, '   => testing new emax mutiple =>', frac_emax
                        END IF
                     END IF
                     niter = niter + 1
                  END IF

               END DO
               rlat_old = rlat
            END DO
         END DO
         
         DEALLOCATE ( VLAT_SPLIT_BOUNDS, IJ_VLAT_IN, e1_in, e2_in, ztmp_out )
         
      END IF

      
      IF ( ldebug ) THEN
         CALL PRTMASK(REAL(JIpos,4), 'JIpos.nc', 'ji')
         CALL PRTMASK(REAL(JJpos,4), 'JJpos.nc', 'jj')
      END IF

      DEALLOCATE ( Xdist )

   END SUBROUTINE FIND_NEAREST_POINT



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
         zr = (6378.137 + 6356.7523)/2.0 ! km
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




   FUNCTION DISTANCE_2D(plona, Xlonb, plata, Xlatb)

      !!----------------------------------------------------------
      !!           ***  FUNCTION  DIST  ***
      !!
      !!  ** Purpose : Compute the distance (km) between
      !!               point A (lona, lata) and each point of domain
      !!               B (lonb,latb)
      !!
      !!  ** Method : Compute the distance along the orthodromy
      !!
      !! * history : J.M. Molines in CHART, f90, may 2007
      !!----------------------------------------------------------

      IMPLICIT NONE
      ! Argument
      REAL(8),                 INTENT(in) :: plona, plata
      REAL(8), DIMENSION(:,:), INTENT(in) :: Xlonb, Xlatb

      REAL(8), DIMENSION(SIZE(Xlonb,1),SIZE(Xlonb,2)) :: distance_2d

      ! Local variables
      INTEGER :: nx, ny, ji, jj

      REAL(8) :: zlatar, zlatbr, zlonar, zlonbr
      REAL(8) :: zpds
      REAL(8) :: zux, zuy, zuz
      REAL(8) :: zr, zpi, zconv, zvx, zvy, zvz

      nx = SIZE(Xlonb,1)
      ny = SIZE(Xlonb,2)

      !! Initialise some values at first call
      ! constants
      zpi = ACOS(-1._8)
      zconv = zpi/180.  ! for degree to radian conversion
      ! Earth radius
      zr = (6378.137 + 6356.7523)/2.0 ! km

      !! Compute these term only if they differ from previous call
      zlatar=plata*zconv
      zlonar=plona*zconv
      zux=COS(zlonar)*COS(zlatar)
      zuy=SIN(zlonar)*COS(zlatar)
      zuz=SIN(zlatar)

      distance_2d(:,:) = 0.

      DO jj=1,ny
         DO ji=1,nx

            zlatbr=Xlatb(ji,jj)*zconv
            zlonbr=Xlonb(ji,jj)*zconv
            zvx=COS(zlonbr)*COS(zlatbr)
            zvy=SIN(zlonbr)*COS(zlatbr)
            zvz=SIN(zlatbr)

            zpds = zux*zvx + zuy*zvy + zuz*zvz

            IF ( zpds < 1.) distance_2d(ji,jj) = zr*ACOS(zpds)

         END DO
      END DO

   END FUNCTION DISTANCE_2D






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


   FUNCTION L_IS_GRID_REGULAR( Xlon, Xlat )

      !!----------------------------------------------------------
      !! Tell if a grid (1 longitude 2D array and 1 latitude 2D array) is regular
      !!----------------------------------------------------------

      IMPLICIT NONE
      ! Argument
      REAL(8), DIMENSION(:,:), INTENT(in) :: Xlon, Xlat
      LOGICAL                             :: l_is_grid_regular
      INTEGER :: nx, ny, ji, jj

      nx = SIZE(Xlon,1)
      ny = SIZE(Xlon,2)
      IF ( (SIZE(Xlat,1) /= nx) .OR. (SIZE(Xlat,2) /= ny) ) THEN
         PRINT *, ' ERROR (L_IS_GRID_REGULAR of mod_grids.f90): Xlat does not agree in shape with Xlon!'
         STOP
      END IF

      l_is_grid_regular = .TRUE.

      !!  a/ checking on longitude array: (LOLO: use epsilon(Xlon) instead 1.E-12?)
      DO jj = 2, ny
         IF ( SUM( ABS(Xlon(:,jj) - Xlon(:,1)) ) > 1.E-12 ) THEN
            l_is_grid_regular = .FALSE.
            EXIT
         END IF
      END DO
      !!  b/ now on latitude array:
      IF ( l_is_grid_regular ) THEN
         DO ji = 1, nx
            IF ( SUM( ABS(Xlat(ji,:) - Xlat(1,:)) ) > 1.E-12 ) THEN
               l_is_grid_regular = .FALSE.
               EXIT
            END IF
         END DO
      END IF

   END FUNCTION L_IS_GRID_REGULAR


END MODULE MOD_MANIP
