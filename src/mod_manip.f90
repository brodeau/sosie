MODULE MOD_MANIP

   !! Misc. manipulations and operations on 2D arrays...

   !! Author: L. Brodeau


   IMPLICIT NONE

   PRIVATE


   INTERFACE flip_ud
      MODULE PROCEDURE flip_ud_1d_r4, flip_ud_1d_r8, flip_ud_2d_r4, flip_ud_3d_i1
   END INTERFACE flip_ud

   INTERFACE degE_to_degWE
      MODULE PROCEDURE degE_to_degWE_scal, degE_to_degWE_1d, degE_to_degWE_2d
   END INTERFACE degE_to_degWE

   INTERFACE extra_2_east
      MODULE PROCEDURE extra_2_east_r4, extra_2_east_r8
   END INTERFACE extra_2_east

   INTERFACE extra_2_west
      MODULE PROCEDURE extra_2_west_r4, extra_2_west_r8
   END INTERFACE extra_2_west

   INTERFACE long_reorg_3d
      MODULE PROCEDURE long_reorg_3d_i1
   END INTERFACE long_reorg_3d


   PUBLIC :: fill_extra_bands, fill_extra_north_south, extra_2_east, extra_2_west, partial_deriv, &
      &      flip_ud, long_reorg_2d, long_reorg_3d, &
      &      distance, distance_2d, &
      &      find_nearest_point, &
      &      shrink_vector, degE_to_degWE, &
      &      ext_north_to_90_regg

   REAL(8), PARAMETER, PUBLIC :: rflg = -9999.


   !LOGICAL, PARAMETER :: ldebug = .TRUE., l_force_use_of_twisted = .FALSE.
   !LOGICAL, PARAMETER :: ldebug = .FALSE., l_force_use_of_twisted = .TRUE.
   LOGICAL, PARAMETER :: ldebug = .FALSE., l_force_use_of_twisted = .FALSE.

   LOGICAL  :: lfirst_dist = .TRUE.



CONTAINS




   SUBROUTINE FILL_EXTRA_BANDS(k_ew, XX, YY, XF, XP4, YP4, FP4,  is_orca_grid)

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
      INTEGER ,   OPTIONAL,    INTENT(in)  :: is_orca_grid

      INTEGER :: nx, ny, nxp4, nyp4
      INTEGER :: ji, jj, iorca

      iorca = 0
      IF ( PRESENT(is_orca_grid) ) iorca = is_orca_grid


      IF ( (SIZE(XX,1) /= SIZE(YY,1)).OR.(SIZE(XX,2) /= SIZE(YY,2)).OR. &
         & (SIZE(XX,1) /= SIZE(XF,1)).OR.(SIZE(XX,2) /= SIZE(XF,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_BANDS : size of input coor. and data do not match!!!'
         PRINT *, 'SIZE(XX,1), SIZE(YY,1), SIZE(XF,1) =>', SIZE(XX,1), SIZE(YY,1), SIZE(XF,1)
         PRINT *, 'SIZE(XX,2), SIZE(YY,2), SIZE(XF,2) =>', SIZE(XX,2), SIZE(YY,2), SIZE(XF,2)
         STOP
      END IF

      IF ( (SIZE(XP4,1) /= SIZE(YP4,1)).OR.(SIZE(XP4,2) /= SIZE(YP4,2)).OR. &
         & (SIZE(XP4,1) /= SIZE(FP4,1)).OR.(SIZE(XP4,2) /= SIZE(FP4,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => FILL_EXTRA_BANDS : size of output coor. and data do not match!!!'
         PRINT *, 'SIZE(XP4,1), SIZE(YP4,1), SIZE(FP4,1) =>', SIZE(XP4,1), SIZE(YP4,1), SIZE(FP4,1)
         PRINT *, 'SIZE(XP4,2), SIZE(YP4,2), SIZE(FP4,2) =>', SIZE(XP4,2), SIZE(YP4,2), SIZE(FP4,2)
         STOP
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
      XP4(3:nxp4-2, 3:nyp4-2) = XX(:,:)
      YP4(3:nxp4-2, 3:nyp4-2) = YY(:,:)
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

         !! WEST
         XP4(2, 3:nyp4-2) = XX(2,:) - (XX(3,:) - XX(1,:))
         XP4(1, 3:nyp4-2) = XX(1,:) - (XX(3,:) - XX(1,:))

         !! EAST
         XP4(nxp4-1, 3:nyp4-2) = XX(nx-1,:) + XX(nx,:) - XX(nx-2,:)
         XP4(nxp4  , 3:nyp4-2) = XX(nx,:)   + XX(nx,:) - XX(nx-2,:)

      END IF !IF (k_ew > -1)


      
      !! ******************
      !! Southern Extension
      !! ******************
      XP4(:, 2) = XP4(:,4) - (XP4(:,5) - XP4(:,3))
      XP4(:, 1) = XP4(:,3) - (XP4(:,5) - XP4(:,3))

      !! ******************      
      !! Northern Extension
      !! ******************
      
      SELECT CASE( iorca )
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
      SELECT CASE( iorca )

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
      SELECT CASE( iorca )

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





   SUBROUTINE FILL_EXTRA_NORTH_SOUTH(XX, YY, XF, XP4, YP4, FP4,  is_orca_grid)

      !!============================================================================
      !! Extending input arrays with an extraband of two points at northern and
      !! southern boundaries.
      !!
      !! The extension is done thanks to Akima's exptrapolation method or continuity
      !!  / ORCA knowledge...
      !!
      !!============================================================================
      REAL(8), DIMENSION(:,:), INTENT(in)  :: XX, YY, XF
      REAL(8), DIMENSION(:,:), INTENT(out) :: XP4, YP4, FP4
      INTEGER ,   OPTIONAL,    INTENT(in)  :: is_orca_grid

      !! Local
      INTEGER :: nx, ny, nxp4, nyp4
      INTEGER :: ji, iorca

      iorca = 0
      IF ( PRESENT(is_orca_grid) ) iorca = is_orca_grid

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
      SELECT CASE( iorca )
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
      SELECT CASE( iorca )

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
      SELECT CASE( iorca )

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



















   SUBROUTINE extra_2_east_r8(x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)
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
         y4 = y3
         y5 = y3
      ELSE
         y4   = C*(2*BET/B - ALF/A) + y3
         y5   = y4 + y4*D/C + BET*D/B - ALF*D/A - y3*D/C
      END IF
   END SUBROUTINE extra_2_east_r8

   SUBROUTINE extra_2_east_r4(x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)
      REAL(4), INTENT(in)  :: x1, x2, x3, x4, x5, y1, y2, y3
      REAL(4), INTENT(out) :: y4, y5
      REAL(4) :: A, B, C, D, ALF, BET
      A    = x2 - x1
      B    = x3 - x2
      C    = x4 - x3
      D    = x5 - x4
      ALF  = y2 - y1
      BET  = y3 - y2
      IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
         y4 = y3
         y5 = y3
      ELSE
         y4 = C*(2*BET/B - ALF/A) + y3
         y5 = y4 + y4*D/C + BET*D/B - ALF*D/A - y3*D/C
      END IF
   END SUBROUTINE extra_2_east_r4



   SUBROUTINE extra_2_west_r8(x5, x4, x3, x2, x1, y5, y4, y3, y2, y1)
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
         y2 = y3
         y1 = y3
      ELSE
         y2 = C*(2*BET/B - ALF/A) + y3
         y1 = y2 + y2*D/C + BET*D/B - ALF*D/A - y3*D/C
      END IF
   END SUBROUTINE extra_2_west_r8

   SUBROUTINE extra_2_west_r4(x5, x4, x3, x2, x1, y5, y4, y3, y2, y1)
      REAL(4), INTENT(in)  :: x1, x2, x3, x4, x5, y5, y4, y3
      REAL(4), INTENT(out) :: y1, y2
      REAL(4) :: A, B, C, D, ALF, BET
      A    = x4 - x5
      B    = x3 - x4
      C    = x2 - x3
      D    = x1 - x2
      ALF  = y4 - y5
      BET  = y3 - y4
      IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
         y2 = y3
         y1 = y3
      ELSE
         y2 = C*(2*BET/B - ALF/A) + y3
         y1 = y2 + y2*D/C + BET*D/B - ALF*D/A - y3*D/C
      END IF
   END SUBROUTINE extra_2_west_r4




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


   SUBROUTINE FLIP_UD_1D_R4(XF)
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
   END SUBROUTINE FLIP_UD_1D_R4

   SUBROUTINE FLIP_UD_1D_R8(XF)
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
   END SUBROUTINE FLIP_UD_1D_R8


   SUBROUTINE FLIP_UD_2D_R4(XF)

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

   END SUBROUTINE FLIP_UD_2D_R4


   SUBROUTINE FLIP_UD_3D_I1(XF)

      INTEGER(1), DIMENSION(:,:,:), INTENT(inout) :: XF

      INTEGER :: nx, ny, nz, jj
      INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: ztmp

      nx = SIZE(XF,1) ; ny = SIZE(XF,2) ; nz = SIZE(XF,3)

      ALLOCATE ( ztmp(nx,ny,nz) )

      ztmp(:,:,:) = XF(:,:,:)

      DO jj = 1, ny
         XF(:,jj,:) =  ztmp(:,ny-jj+1,:)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE FLIP_UD_3D_I1





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



   SUBROUTINE LONG_REORG_3D_I1(i_chg_x, XF)

      INTEGER, INTENT(in) :: i_chg_x
      INTEGER(1), DIMENSION(:,:,:), INTENT(inout) :: XF

      INTEGER :: nx, ny, nz, ji
      INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: ztmp

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

   END SUBROUTINE LONG_REORG_3D_I1






   SUBROUTINE FIND_NEAREST_POINT(Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp,  mask_domain_trg)
      !!---------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_POINT  ***
      !!
      !!   Source domain: Xsrc, Ysrc
      !!   Target domain: Xtrg, Ytrg
      !!     => for each point of target domain: i_src,j_src location of nearest point on source domain
      !!        => JIp & JJp have same shape as Xtrg & Ytrg !
      !!
      !!                Laurent Brodeau, July, 2017
      !!
      !!
      !! OPTIONAL:
      !!      * mask_domain_trg: ignore (dont't treat) regions of the target domain where mask_domain_trg==0 !
      !!---------------------------------------------------------------
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xtrg, Ytrg    !: lon and lat arrays of target domain
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xsrc , Ysrc     !: lon and lat arrays of source domain
      INTEGER(4), DIMENSION(:,:), INTENT(out) :: JIp, JJp  !: nearest point location of point P in Xsrc,Ysrc wrt Xtrg,Ytrg
      INTEGER(1), OPTIONAL ,DIMENSION(:,:), INTENT(inout) :: mask_domain_trg

      INTEGER :: jj, nx_src, ny_src, nx_trg, ny_trg, j_strt_trg, j_stop_trg, jlat_icr
      REAL(8) :: y_max_src, y_min_src, rmin_dlat_dj, rtmp
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: mask_ignore_trg
      INTEGER,    DIMENSION(:),   ALLOCATABLE :: i1dum
      REAL(8),    DIMENSION(:,:), ALLOCATABLE :: ztmp_trg
      LOGICAL :: l_is_reg_src, l_is_reg_trg

      nx_src  = SIZE(Xsrc,1)
      ny_src  = SIZE(Xsrc,2)
      nx_trg = SIZE(Xtrg,1)
      ny_trg = SIZE(Xtrg,2)

      PRINT *, ' Source domain size: ', nx_src, ny_src
      PRINT *, ' Target domain size: ', nx_trg, ny_trg

      IF ( (SIZE(Ysrc,1) /= nx_src) .OR. (SIZE(Ysrc,2) /= ny_src) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Ysrc dont agree in shape with Xsrc'
         STOP
      END IF
      IF ( (SIZE(Ytrg,1) /= nx_trg) .OR. (SIZE(Ytrg,2) /= ny_trg) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Ytrg dont agree in shape with Xtrg'
         STOP
      END IF
      IF ( (SIZE(JIp,1) /= nx_trg) .OR. (SIZE(JIp,2) /= ny_trg) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): JIp dont agree in shape with Xtrg'
         PRINT *, SIZE(JIp,1), SIZE(JIp,2), 'vs', nx_trg, ny_trg
         STOP
      END IF
      IF ( (SIZE(JJp,1) /= nx_trg) .OR. (SIZE(JJp,2) /= ny_trg) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): JJp dont agree in shape with Xtrg'
         STOP
      END IF

      PRINT *, ''
      !! Checking if source domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_src = L_IS_GRID_REGULAR( Xsrc , Ysrc )
      PRINT *, ' *** FIND_NEAREST_POINT => Is source grid regular ??? =>', l_is_reg_src

      !! Checking if target domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_trg = L_IS_GRID_REGULAR( Xtrg , Ytrg )
      PRINT *, ' *** FIND_NEAREST_POINT => Is target grid regular ??? =>', l_is_reg_trg

      ALLOCATE ( mask_ignore_trg(nx_trg,ny_trg) , i1dum(nx_trg) )
      mask_ignore_trg(:,:) = 1
      IF ( PRESENT( mask_domain_trg ) ) mask_ignore_trg(:,:) = mask_domain_trg(:,:)

      y_min_src  = MINVAL(Ysrc) ; ! Min and Max latitude of source domain
      y_max_src  = MAXVAL(Ysrc)

      !! General case:
      j_strt_trg = 1
      j_stop_trg = ny_trg
      jlat_icr   = 1

      !! We need to know if the target latitude ONLY keeps on systematically
      !! increasing (or decreasing) as j increases:
      rmin_dlat_dj = -100.
      !IF ( REAL(nx_trg)/REAL(ny_trg) > 0.15  ) THEN
      IF ( nx_trg > 2 ) THEN !! We can reasonably say it's a map and not a vector or almost a vectot...
         ALLOCATE ( ztmp_trg(nx_trg, ny_trg) )
         DO jj = 2, ny_trg
            ztmp_trg(:,jj) = Ytrg(:,jj) - Ytrg(:,jj-1)
         END DO
         rtmp = SUM(ztmp_trg(:,ny_trg/2)) ! know if increasing (>0) or decreasing (<0)
         ztmp_trg = SIGN(1.0_8 , rtmp)*ztmp_trg
         !IF (ldebug) CALL DUMP_2D_FIELD(REAL(ztmp_trg,4), 'dlat_dj_trg.nc', 'dist')
         rmin_dlat_dj = MINVAL(ztmp_trg(:,2:))
         IF (ldebug) PRINT *, ' Minimum dlat_dj_trg =>', rmin_dlat_dj
         DEALLOCATE ( ztmp_trg )
      END IF

      !!  Simplif when [d lat / d j] always has the same sign:
      IF ( (rmin_dlat_dj > -1.E-12) .OR. l_is_reg_trg ) THEN    !!!.OR. (i_orca_trg > 0) ) THEN
         !! -> because we need to avoid all the following if target grid is for
         !!    example a polar sterographic projection... (example 5)
         !!
         !! *** Will ignore regions of the TARGET domain that
         !!     are not covered by source domain:
         !ij_min_loc = MINLOC(Ytrg, mask=(Ytrg>=y_min_src))  ! possible bug identified by Max (aka MB)
         !ij_max_loc = MAXLOC(Ytrg, mask=(Ytrg<=y_max_src))
         !j_strt_trg = MINVAL(MINLOC(Ytrg, mask=(Ytrg>=y_min_src), dim=2))  ! smallest j on target source that covers smallest source latitude
         i1dum = MINLOC(Ytrg, mask=(Ytrg>=y_min_src), dim=2)
         j_strt_trg = MINVAL(i1dum, mask=(i1dum>0))
         DEALLOCATE ( i1dum )
         j_stop_trg = MAXVAL(MAXLOC(Ytrg, mask=(Ytrg<=y_max_src), dim=2))  ! largest j on target source that covers largest source latitude
         IF ( j_strt_trg > j_stop_trg ) jlat_icr = -1 ! latitude decreases as j increases (like ECMWF grids...)
         IF (ldebug) THEN
            PRINT *, ' j_strt_trg, j_stop_trg / nj_trg =>', j_strt_trg, j_stop_trg, '/', ny_trg
            PRINT *, ''
         END IF
      END IF ! IF ( (rmin_dlat_dj >= 0.0_8) .OR. l_is_reg_trg .OR. (i_orca_trg > 0) )

      IF ( l_is_reg_src .AND. (.NOT. l_force_use_of_twisted) ) THEN
         PRINT *, '                 => going for simple FIND_NEAREST algorithm !'; PRINT *, ''
         CALL FIND_NEAREST_EASY(    Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp, mask_ignore_trg, &
            &                       j_strt_trg, j_stop_trg, jlat_icr )
      ELSE
         PRINT *, '                  => going for advanced FIND_NEAREST algorithm !'; PRINT *, ''
         CALL FIND_NEAREST_TWISTED( Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp, mask_ignore_trg, &
            &                       j_strt_trg, j_stop_trg, jlat_icr )
      END IF
      PRINT *, ''

      IF ( PRESENT( mask_domain_trg ) ) mask_domain_trg(:,:) = mask_ignore_trg(:,:)
      DEALLOCATE ( mask_ignore_trg )

   END SUBROUTINE FIND_NEAREST_POINT



   SUBROUTINE FIND_NEAREST_EASY( Xtrg, Ytrg, Xsrc, Ysrc, JIpos, JJpos, mask_trg, &
      &                          j_strt_trg, j_stop_trg, jlat_icr )
      !!---------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_EASY  ***
      !!
      !!   Source domain: Xsrc, Ysrc
      !!   Target domain: Xtrg, Ytrg
      !!     => for each point of target domain: i_src,j_src location of nearest point on source domain
      !!        => JIpos & JJpos have same shape as Xtrg & Ytrg !
      !!
      !!                Laurent Brodeau, July, 2017
      !!
      !!
      !!      * mask_trg: ignore (dont't treat) regions of the target domain where mask_trg==0 !
      !!---------------------------------------------------------------
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xtrg, Ytrg    !: lon and lat arrays of target domain
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xsrc , Ysrc     !: lon and lat arrays of source domain
      INTEGER(4), DIMENSION(:,:), INTENT(out) :: JIpos, JJpos  !: nearest point location of point P in Xsrc,Ysrc wrt Xtrg,Ytrg
      INTEGER(1), DIMENSION(:,:), INTENT(in)  :: mask_trg
      INTEGER,                    INTENT(in)  :: j_strt_trg, j_stop_trg, jlat_icr

      !! Important parameters:
      INTEGER, PARAMETER :: nframe_scan = 4  ! domain to scan for nearest point in simple algo => domain of 9x9

      INTEGER :: &
         &    nx_src, ny_src, nx_trg, ny_trg, &
         &    ji_trg, jj_trg, ji_src, jj_src, &
         &    jmin_src, jmax_src, imin_src, imax_src

      REAL(8) :: rlon, rlat

      INTEGER, DIMENSION(2) :: ij_min_loc
      INTEGER, DIMENSION(1) :: ip, jp

      REAL(8), DIMENSION(:),   ALLOCATABLE :: vlat, vlon
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xdist

      nx_src  = SIZE(Xsrc,1)
      ny_src  = SIZE(Xsrc,2)
      nx_trg = SIZE(Xtrg,1)
      ny_trg = SIZE(Xtrg,2)

      JIpos(:,:) = -1
      JJpos(:,:) = -1

      ALLOCATE ( vlon(nx_src) , vlat(ny_src) , Xdist(nx_src,ny_src) )

      vlon(:) = Xsrc(:,3)   ! 3 => make sure we're inside the domain
      vlat(:) = Ysrc(3,:)   ! 3 => make sure we're inside the domain

      DO jj_trg = j_strt_trg, j_stop_trg, jlat_icr
         DO ji_trg = 1, nx_trg

            IF ( mask_trg(ji_trg,jj_trg) == 1 ) THEN

               IF ( (ji_trg==1) .AND. MOD(jj_trg,10)==0 ) PRINT *, ' *** Treated j-point of target domain =', jj_trg !REAL(rlat,4)

               rlon = Xtrg(ji_trg,jj_trg)
               rlat = Ytrg(ji_trg,jj_trg)

               ! Nearest point (in terms of index):
               ip =  MINLOC(ABS(vlon(:)-rlon))
               jp =  MINLOC(ABS(vlat(:)-rlat))

               !! Define the box to scan for shortest distance:
               imin_src = MAX(ip(1)-nframe_scan ,  1)
               imax_src = MIN(ip(1)+nframe_scan , nx_src)
               jmin_src = MAX(jp(1)-nframe_scan ,  1)
               jmax_src = MIN(jp(1)+nframe_scan , ny_src)

               ! Nearest point (in terms of distance):
               Xdist = 1.E12
               Xdist(imin_src:imax_src,jmin_src:jmax_src) = DISTANCE_2D(rlon, Xsrc(imin_src:imax_src,jmin_src:jmax_src), rlat, Ysrc(imin_src:imax_src,jmin_src:jmax_src))
               ij_min_loc = MINLOC(Xdist(imin_src:imax_src,jmin_src:jmax_src))
               ji_src = ij_min_loc(1) + imin_src - 1
               jj_src = ij_min_loc(2) + jmin_src - 1
               JIpos(ji_trg,jj_trg) = ji_src
               JJpos(ji_trg,jj_trg) = jj_src
               IF ((ji_src==0).OR.(jj_src==0)) THEN
                  PRINT *, ''
                  PRINT *, 'The nearest point was not found!'
                  PRINT *, ' !!! ji_src or jj_src = 0 !!!'
                  PRINT *, ' ** Target point (lon,lat) =>', rlon, rlat
                  PRINT *, ' Nearest point found on source grid:', &
                     &  REAL(Xsrc(ji_src,jj_src) , 4), &
                     &  REAL(Ysrc(ji_src,jj_src) , 4)
                  STOP
               END IF

               IF ( ( JIpos(ji_trg,jj_trg) == -1 ).OR.( JJpos(ji_trg,jj_trg) == -1 ) ) THEN
                  PRINT *, 'ERROR in FIND_NEAREST_EASY of mod_bilin_2d.f90 !'
                  PRINT *, 'Point rlon, rlat', rlon, rlat
                  PRINT *, 'not found on the source grid!'
                  STOP
               END IF

            END IF
         END DO
      END DO

      DEALLOCATE ( vlon , vlat , Xdist )

   END SUBROUTINE FIND_NEAREST_EASY



   SUBROUTINE FIND_NEAREST_TWISTED( Xtrg, Ytrg, Xsrc, Ysrc, JIpos, JJpos, mask_trg, &
      &                             j_strt_trg, j_stop_trg, jlat_icr )
      !!---------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_POINT  ***
      !!
      !!   Source domain: Xsrc, Ysrc
      !!   Target domain: Xtrg, Ytrg
      !!     => for each point of target domain: i_src,j_src location of nearest point on source domain
      !!        => JIpos & JJpos have same shape as Xtrg & Ytrg !
      !!
      !!                Laurent Brodeau, July, 2017
      !!
      !!
      !!      * mask_trg: ignore (dont't treat) regions of the target domain where mask_trg==0 !
      !!---------------------------------------------------------------
      REAL(8),    DIMENSION(:,:), INTENT(in)    :: Xtrg, Ytrg    !: lon and lat arrays of target domain
      REAL(8),    DIMENSION(:,:), INTENT(in)    :: Xsrc , Ysrc     !: lon and lat arrays of source domain
      INTEGER(4), DIMENSION(:,:), INTENT(out)   :: JIpos, JJpos  !: nearest point location of point P in Xsrc,Ysrc wrt Xtrg,Ytrg
      INTEGER(1), DIMENSION(:,:), INTENT(inout) :: mask_trg
      INTEGER, INTENT(in)                       :: j_strt_trg, j_stop_trg, jlat_icr
      !!
      !! Important parameters:
      INTEGER, PARAMETER :: Nlat_split = 40   ! number of latitude bands to split the search work (for 180. degree south->north)
      REAL(8), PARAMETER :: frac_emax = 0.51   ! fraction of emax to test if found!
      !!                                        => 0.5 seems to be too small, for example on ORCA1 grid at around 20 deg N... ???
      INTEGER :: &
         &    nx_src, ny_src, nx_trg, ny_trg, &
         &    jlat, ji_trg, jj_trg, ji_src, jj_src, jj_trg_old, &
         &    jmin_src, jmax_src, imin_src, imax_src, niter, &
         &    nsplit, jmax_band, jmin_band

      REAL(8) :: emax, rlat_low, rlat_hgh, rlon, rlat, rlat_old
      REAL(8) :: y_max_trg, y_min_trg, y_max_bnd, y_min_bnd, y_max_bnd0, y_min_bnd0, dy, &
         &       y_max_src, y_min_src
      REAL(8),    DIMENSION(:),   ALLOCATABLE :: VLAT_SPLIT_BOUNDS
      REAL(8),    DIMENSION(:,:), ALLOCATABLE :: Xdist, e1_src, e2_src    !: grid layout and metrics
      INTEGER,    DIMENSION(:,:), ALLOCATABLE :: J_VLAT_SRC
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: mspot_lon, mspot_lat
      INTEGER,    DIMENSION(:),   ALLOCATABLE :: i1dum

      INTEGER, DIMENSION(2) :: ij_min_loc

      LOGICAL :: lagain


      nx_src  = SIZE(Xsrc,1)
      ny_src  = SIZE(Xsrc,2)
      nx_trg = SIZE(Xtrg,1)
      ny_trg = SIZE(Xtrg,2)

      ALLOCATE ( Xdist(nx_src,ny_src) , mspot_lon(nx_src,ny_src) , mspot_lat(nx_src,ny_src) , i1dum(nx_src) )

      y_min_src  = MINVAL(Ysrc) ; ! Min and Max latitude of source domain
      y_max_src  = MAXVAL(Ysrc)
      y_min_trg = MINVAL(Ytrg) ; ! Min and Max latitude of target domain
      y_max_trg = MAXVAL(Ytrg)


      !! Control if source & target domains overlap:
      PRINT *, ' Min. latitude on source & target domains =>', REAL(y_min_src,4), REAL(y_min_trg,4)
      PRINT *, ' Max. latitude on source & target domains =>', REAL(y_max_src,4), REAL(y_max_trg,4)

      IF ( (y_min_trg>=y_max_src).OR.(y_max_trg<=y_min_src) ) THEN
         WRITE(6,*)' ERROR (FIND_NEAREST_TWISTED of mod_manip.f90): Target and source latitudes do not overlap!'
         STOP
      END IF
      !! ---------------------------------------------------------------------------------------
      PRINT *, '                        => going for advanced algorithm !'

      JIpos(:,:) = -1
      JJpos(:,:) = -1

      ALLOCATE ( e1_src(nx_src,ny_src), e2_src(nx_src,ny_src) )

      !! We need metric of input grid
      e1_src(:,:) = 40000. ;  e2_src(:,:) = 40000.
      DO jj_src=1, ny_src
         DO ji_src=1, nx_src-1
            e1_src(ji_src,jj_src) = distance(Xsrc(ji_src,jj_src),Xsrc(ji_src+1,jj_src),Ysrc(ji_src,jj_src),Ysrc(ji_src+1,jj_src))*1000. ! (m)
         END DO
      END DO
      DO jj_src=1, ny_src-1
         DO ji_src=1, nx_src
            e2_src(ji_src,jj_src) = distance(Xsrc(ji_src,jj_src),Xsrc(ji_src,jj_src+1),Ysrc(ji_src,jj_src),Ysrc(ji_src,jj_src+1))*1000. !(m)
         END DO
      END DO
      IF (nx_src>1) e1_src(nx_src,:) = e1_src(nx_src-1,:)
      IF (ny_src>1) e2_src(:,ny_src) = e2_src(:,ny_src-1)
      !IF (ldebug) THEN
      !   CALL DUMP_2D_FIELD(REAL(e1_src,4), 'e1_src.nc', 'e1')
      !   CALL DUMP_2D_FIELD(REAL(e2_src,4), 'e2_src.nc', 'e2')
      !END IF

      !! Min and Max latitude to use for binning :
      PRINT *, ''
      PRINT *, '    *** Latitude binning for more efficient localization *** '
      y_max_bnd = MIN( y_max_trg , y_max_src )  !lolo: why MIN() ??? I don't remember why
      y_max_bnd0 = y_max_bnd
      y_min_bnd = MAX( y_min_trg , y_min_src )  !lolo: same, MAX()???
      y_min_bnd0 = y_min_bnd
      !PRINT *, ' y_max_bnd #1 => ', y_max_bnd ; PRINT *, ' y_min_bnd #1 => ', y_min_bnd
      y_max_bnd = MIN( REAL(INT(y_max_bnd+1),8) ,  90.)
      y_min_bnd = MAX( REAL(INT(y_min_bnd-1),8) , -90.)
      !PRINT *, ' y_max_bnd #2 => ', y_max_bnd ; PRINT *, ' y_min_bnd #2 => ', y_min_bnd
      !! Multiple of 0.5:
      y_max_bnd = MIN( NINT(y_max_bnd/0.5)*0.5    ,  90.)
      y_min_bnd = MAX( NINT(y_min_bnd/0.5)*0.5    , -90.)
      !lolo:
      y_max_bnd = MIN( y_max_bnd , y_max_bnd0)
      y_min_bnd = MAX( y_min_bnd , y_min_bnd0)
      !lolo.

      PRINT *, '   => binning from ', y_min_bnd, ' to ', y_max_bnd
      !!
      nsplit = MAX( INT( REAL(Nlat_split) * (y_max_bnd - y_min_bnd)/180. ) , 1)
      !!
      ALLOCATE ( VLAT_SPLIT_BOUNDS(nsplit+1), J_VLAT_SRC(nsplit,2) )
      dy = (y_max_bnd - y_min_bnd)/REAL(nsplit)
      DO jlat=1,nsplit+1
         VLAT_SPLIT_BOUNDS(jlat) = y_min_bnd + REAL(jlat-1)*dy
      END DO
      !!

      IF ( (y_min_bnd > y_min_bnd0).OR.(y_max_bnd < y_max_bnd0) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Bounds for latitude for VLAT_SPLIT_BOUNDS are bad!'
         PRINT *, ' y_min_bnd, y_max_bnd =', y_min_bnd, y_max_bnd !
         PRINT *, ' y_min_bnd0, y_max_bnd0 =', y_min_bnd0, y_max_bnd0
         STOP
      END IF

      J_VLAT_SRC = 0

      DO jlat = 1, nsplit
         rlat_low = VLAT_SPLIT_BOUNDS(jlat)
         rlat_hgh = VLAT_SPLIT_BOUNDS(jlat+1)
         !! @mbalaro Comment: The two line below can lead to error when working on on small domain ...
         !ij_max_loc = MAXLOC(Ysrc, mask=(Ysrc<=rlat_hgh)) ; jmax_band = ij_max_loc(2)
         !ij_min_loc = MINLOC(Ysrc, mask=(Ysrc>=rlat_low)) ; jmin_band = ij_min_loc(2)
         !! ... it is preferable to look at the min and max value of the ensemble of jj within the range [rlat_low:rlat_hgh]
         !! Largest ever possible j index of the highest latitude in region where Ysrc<=rlat_hgh
         jmax_band = MAXVAL(MAXLOC(Ysrc, mask=(Ysrc<=rlat_hgh), dim=2))  ! MAXLOC(Ysrc, ..., dim=2) returns ni_src values (the max in each column)
         !! Smalles ever possible j index of the smallest latitude in region where Ysrc>=rlat_low
         i1dum = MINLOC(Ysrc, mask=(Ysrc .GT. rlat_low), dim=2)
         jmin_band = MINVAL(i1dum, mask=(i1dum>0))
         DEALLOCATE ( i1dum )
         !!
         !! To be sure to include everything, adding 1 extra points below and above:
         J_VLAT_SRC(jlat,1) = MAX(jmin_band - 1,   1  )
         J_VLAT_SRC(jlat,2) = MIN(jmax_band + 1, ny_src)
         !!
         IF ( ldebug ) THEN
            PRINT *, ' Latitude bin #', jlat
            PRINT *, '     => lat_low, lat_high:', REAL(rlat_low,4), REAL(rlat_hgh,4)
            PRINT *, '     => JJ min and max on input domain =>', J_VLAT_SRC(jlat,1), J_VLAT_SRC(jlat,2)
         END IF
         !!
      END DO

      PRINT *, ''
      PRINT *, '     => VLAT_SPLIT_BOUNDS ='
      PRINT *, VLAT_SPLIT_BOUNDS
      PRINT *, '     => corresponding start values for j_src ='
      PRINT *, J_VLAT_SRC(:,1)
      PRINT *, '     => corresponding stop values for j_src ='
      PRINT *, J_VLAT_SRC(:,2)
      PRINT *, ''

      DO jlat = 1, nsplit
         IF ( J_VLAT_SRC(jlat,2) <= J_VLAT_SRC(jlat,1) ) THEN
            PRINT *, ' ERROR: jj_stop > jj_start ! ', J_VLAT_SRC(jlat,2), J_VLAT_SRC(jlat,1)
            PRINT *, '   => for latitude bin #', jlat
            STOP
         END IF
      END DO

      rlat_old   = rflg
      jj_trg_old = -10

      DO jj_trg = j_strt_trg, j_stop_trg, jlat_icr
         DO ji_trg = 1, nx_trg
            IF ( mask_trg(ji_trg,jj_trg) == 1 ) THEN

               rlon = Xtrg(ji_trg,jj_trg)
               rlat = Ytrg(ji_trg,jj_trg)

               !! Display progression in stdout:
               IF ( (ji_trg == nx_trg/2).AND.(jj_trg /= jj_trg_old) ) THEN
                  WRITE(*,'("*** Treated latitude of target domain = ",f9.4," (jj_trg = ",i5.5,")")') REAL(rlat,4), jj_trg
                  jj_trg_old = jj_trg
               END IF
               IF ( (nx_trg == 1).AND.(MOD(jj_trg,10)==0) ) &
                  & WRITE(*,'("*** Treated point of target domain = ",i7," (ouf of ",i7,")")') jj_trg, ABS(j_stop_trg-j_strt_trg+1) ! in case of trajectory/ephem stuff

               !! Need to find which jlat of our latitude bins rlat is located in!
               IF ( rlat /= rlat_old ) THEN
                  DO jlat=1,nsplit
                     IF (  rlat ==VLAT_SPLIT_BOUNDS(jlat)) EXIT
                     IF ( (rlat > VLAT_SPLIT_BOUNDS(jlat)).AND.(rlat <= VLAT_SPLIT_BOUNDS(jlat+1)) ) EXIT
                  END DO
                  !!
               END IF

               lagain    = .TRUE.
               niter     = -1  ! -1 because first pass is for bluff, we want niter=0 for the first use of latitude binning...
               IF ( rlat > 60. ) niter = 0 ! we skip the bluff part because the grid might be too close to NP boundary cut!

               DO WHILE ( lagain )

                  IF ( niter == -1 ) THEN
                     !! Bluff !
                     !! It's not stupid to assume that the next point to locate is
                     !! pretty near the previously found point (ji_src,jj_src):
                     imin_src = MAX(ji_src - 5 , 1)
                     imax_src = MIN(ji_src + 5 , nx_src)
                     jmin_src = MAX(jj_src - 5 , 1)
                     jmax_src = MIN(jj_src + 5 , ny_src)
                  ELSE
                     imin_src = 1
                     imax_src = nx_src
                     jmin_src = J_VLAT_SRC(MAX(jlat-niter,1)     , 1)  !! MB Comment: Force to 1 if necessary when nearest point not found whereas it exist really
                     jmax_src = J_VLAT_SRC(MIN(jlat+niter,nsplit), 2)  !! MB Comment: Force to ny_src if necessary when nearest point not found whereas it exist really
                     IF ( ldebug ) THEN
                        PRINT *, ' *** Treated latitude of target domain =', REAL(rlat,4), ' iter:', niter, jlat
                        PRINT *, '     => bin #', jlat
                        PRINT *, '       => jmin & jmax on source domain =', jmin_src, jmax_src
                     END IF
                  END IF

                  Xdist = 1.E12
                  Xdist(imin_src:imax_src,jmin_src:jmax_src) = DISTANCE_2D(rlon, Xsrc(imin_src:imax_src,jmin_src:jmax_src), rlat, Ysrc(imin_src:imax_src,jmin_src:jmax_src))
                  !CALL DUMP_2D_FIELD(REAL(Xdist,4), 'distance_last.nc', 'dist')

                  !! Nearest point is where distance is smallest:
                  ij_min_loc = MINLOC(Xdist(imin_src:imax_src,jmin_src:jmax_src))
                  ji_src = ij_min_loc(1) + imin_src - 1
                  jj_src = ij_min_loc(2) + jmin_src - 1

                  IF ((ji_src==0).OR.(jj_src==0)) THEN
                     PRINT *, ''
                     PRINT *, ' W E I R D  !!!'
                     PRINT *, 'The nearest point was not found!'
                     PRINT *, ' !!! ji_src or jj_src = 0 !!!'
                     PRINT *, ' ** Target point (lon,lat) =>', rlon, rlat
                     STOP
                  END IF

                  emax = MAX(e1_src(ji_src,jj_src),e2_src(ji_src,jj_src))/1000.*SQRT(2.)

                  IF ( Xdist(ji_src,jj_src) <= frac_emax*emax) THEN
                     !! Found !
                     lagain = .FALSE.
                     JIpos(ji_trg,jj_trg) = ji_src
                     JJpos(ji_trg,jj_trg) = jj_src
                     IF ( ldebug ) THEN
                        IF ( niter == -1 ) THEN
                           PRINT *, '    --- F O U N D  with bluff !!! ---'
                        ELSE
                           PRINT *, '    --- F O U N D --- niter =', niter
                        END IF
                     END IF
                     !! Found .
                  ELSE
                     !! Not found yet...
                     IF (niter == 0) THEN
                        !! After all the lon,lat couple we are looking for is maybe not part of source domain
                        !! => could do a test and break the iteration if so...
                        mspot_lon = 0 ; mspot_lat = 0
                        WHERE( (Xsrc > MAX(rlon-0.5,  0.)).AND.(Xsrc < MIN(rlon+0.5,360.)) ) mspot_lon = 1
                        WHERE( (Ysrc > MAX(rlat-0.5,-90.)).AND.(Ysrc < MIN(rlat+0.5, 90.)) ) mspot_lat = 1
                        IF ( SUM(mspot_lon*mspot_lat) == 0 ) THEN
                           lagain = .FALSE.
                           IF ( ldebug ) PRINT *, ' *** FIND_NEAREST_POINT: SHORT leave test worked! Aborting search!'
                        END IF
                     END IF

                     IF (niter > nsplit/3) THEN
                        !! We are too far in latitude, giving up...
                        PRINT *, ' *** WARNING: mod_manip.f90/FIND_NEAREST_POINT: Giving up!!!'
                        PRINT *, '     => did not find nearest point for target coordinates:', & !
                           &              REAL(rlon,4), REAL(rlat,4)
                        lagain = .FALSE.
                     END IF

                  END IF    ! IF (Xdist(ji_src,jj_src) <= frac_emax*emax)

                  niter  = niter + 1

               END DO !DO WHILE ( lagain )
               rlat_old = rlat

            END IF ! IF ( mask_trg(ji_trg,jj_trg) == 1 )
         END DO    ! DO ji_trg = 1, nx_trg
      END DO       ! DO jj_trg = j_strt_trg, j_stop_trg, jlat_icr


      WHERE ( JIpos == -1 ) mask_trg = -1
      WHERE ( JJpos == -1 ) mask_trg = -2

      DEALLOCATE ( VLAT_SPLIT_BOUNDS, J_VLAT_SRC, e1_src, e2_src, mspot_lon , mspot_lat , Xdist )

   END SUBROUTINE FIND_NEAREST_TWISTED


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

      ! Argument
      REAL(8),                 INTENT(in) :: plona, plata
      REAL(8), DIMENSION(:,:), INTENT(in) :: Xlonb, Xlatb

      REAL(8), DIMENSION(SIZE(Xlonb,1),SIZE(Xlonb,2)) :: distance_2d

      ! Local variables
      INTEGER :: nx, ny, ji, jj

      REAL(8) :: zlatar, zlatbr, zlonar, zlonbr
      REAL(8) :: zpds
      REAL(8) :: zux, zuy, zuz
      REAL(8) :: zr, zconv, zvx, zvy, zvz

      nx = SIZE(Xlonb,1)
      ny = SIZE(Xlonb,2)

      !! Initialise some values at first call
      ! constants
      zconv = ACOS(-1._8)/180._8  ! for degree to radian conversion
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

            !IF ( zpds < 1.) distance_2d(ji,jj) = zr*ACOS(zpds)

            distance_2d(ji,jj) = zr*ACOS(MIN(zpds,1.))

         END DO
      END DO

   END FUNCTION DISTANCE_2D






   SUBROUTINE locate_point(rlon_P, rlat_P, Xlon, Xlat, jxfnd, jyfnd)


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
      REAL(8), DIMENSION(:,:), INTENT(in) :: Xlon, Xlat
      LOGICAL                             :: l_is_grid_regular
      INTEGER :: nx, ny, ji, jj
      !!-------------------------------------------------------------------------
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





   FUNCTION SHRINK_VECTOR(vect, vmask_ignore, new_size)
      !!----------------------------------------------------------
      !!           ***  FUNCTION  DIST  ***
      !!
      !!  ** Purpose :
      !!
      !!   vect (DATA) and vmask_ignore (points to remove in vect are where
      !!         vmask_ignore == 0) have the same size !
      !!
      !!----------------------------------------------------------
      REAL(8),    DIMENSION(:), INTENT(in) :: vect
      INTEGER(1), DIMENSION(:), INTENT(in) :: vmask_ignore
      INTEGER,                  INTENT(in) :: new_size
      REAL(8), DIMENSION(new_size)         :: SHRINK_VECTOR
      !!
      INTEGER :: jo, jn, nold
      nold = SIZE(vect,1)
      IF ( nold /= SIZE(vmask_ignore,1) ) THEN
         PRINT *, ' ERROR (SHRINK_VECTOR of mod_manip.f90): data vector and mask vector do not agree in length!'
         STOP
      END IF
      IF ( (new_size > nold).OR.(new_size<=0) ) THEN
         PRINT *, ' ERROR (SHRINK_VECTOR of mod_manip.f90): your new_size does not make sense!'
         STOP
      END IF
      jn = 0
      DO jo = 1, nold
         IF ( vmask_ignore(jo) == 1 ) THEN
            jn = jn + 1
            SHRINK_VECTOR(jn) = vect(jo)
         END IF
      END DO
   END FUNCTION SHRINK_VECTOR


   FUNCTION degE_to_degWE_scal( rlong )
      !! From longitude in 0 -- 360 frame to -180 -- +180 frame...
      REAL(8) :: rlong
      REAL(8) :: degE_to_degWE_scal
      degE_to_degWE_scal = SIGN(1._8,180._8-rlong)*MIN(rlong, ABS(rlong-360._8))
   END FUNCTION degE_to_degWE_scal

   FUNCTION degE_to_degWE_1d( vlong )
      !! From longitude in 0 -- 360 frame to -180 -- +180 frame...
      REAL(8), DIMENSION(:) :: vlong
      REAL(8), DIMENSION(SIZE(vlong,1)) :: degE_to_degWE_1d
      degE_to_degWE_1d = SIGN(1._8,180._8-vlong)*MIN(vlong, ABS(vlong-360._8))
   END FUNCTION degE_to_degWE_1d

   FUNCTION degE_to_degWE_2d( xlong )
      !! From longitude in 0 -- 360 frame to -180 -- +180 frame...
      REAL(8), DIMENSION(:,:) :: xlong
      REAL(8), DIMENSION(SIZE(xlong,1),SIZE(xlong,2)) :: degE_to_degWE_2d
      degE_to_degWE_2d = SIGN(1._8,180._8-xlong)*MIN(xlong, ABS(xlong-360._8))
   END FUNCTION degE_to_degWE_2d




   SUBROUTINE EXT_NORTH_TO_90_REGG( XX, YY, XF,  XP, YP, FP )
      !!============================================================================
      !! We accept only regular lat-lon grid (supposed to include the north
      !! pole)!
      !!
      !! XX, YY is a 2D regular lon-lat grid which is supposed to include the
      !! northpole but for which the highest latitude (at j=Nj) is less than
      !! 90. We're going to use "across-North-Pole" continuity to fill the
      !! last upper row where latitude = 90 !
      !!
      !! => XP, YP is the same grid with an extra upper row where latitude = 90
      !!     => Last upper row of FP (on XP,YP) contains interpolated values of
      !!        the field
      !!============================================================================
      REAL(8), DIMENSION(:,:), INTENT(in)  :: XX, YY
      REAL(4), DIMENSION(:,:), INTENT(in)  :: XF
      REAL(8), DIMENSION(:,:), INTENT(out) :: XP, YP
      REAL(4), DIMENSION(:,:), INTENT(out) :: FP

      !! Local
      REAL(8) :: rr, zlon, zlon_m
      INTEGER :: nx, ny, nyp1, nyp2
      INTEGER :: ji, ji_m

      IF ( (SIZE(XX,1) /= SIZE(YY,1)).OR.(SIZE(XX,2) /= SIZE(YY,2)).OR. &
         & (SIZE(XX,1) /= SIZE(XF,1)).OR.(SIZE(XX,2) /= SIZE(XF,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => EXT_NORTH_TO_90_REGG : size of input coor. and data do not match!!!'; STOP
      END IF
      IF ( (SIZE(XP,1) /= SIZE(YP,1)).OR.(SIZE(XP,2) /= SIZE(YP,2)).OR. &
         & (SIZE(XP,1) /= SIZE(FP,1)).OR.(SIZE(XP,2) /= SIZE(FP,2))) THEN
         PRINT *, 'ERROR, mod_manip.f90 => EXT_NORTH_TO_90_REGG : size of output coor. and data do not match!!!'; STOP
      END IF
      nx = SIZE(XX,1)
      ny = SIZE(XX,2)
      nyp2 = SIZE(XP,2)
      IF ( nyp2 /= ny + 2 ) THEN
         PRINT *, 'ERROR, mod_manip.f90 => EXT_NORTH_TO_90_REGG : target y dim is not nj+2!!!'; STOP
      END IF

      nyp1 = ny + 1

      XP = 0.
      YP = 0.
      FP = 0.

      !! Filling center of domain:
      XP(:, 1:ny) = XX(:,:)
      YP(:, 1:ny) = YY(:,:)
      FP(:, 1:ny) = XF(:,:)

      !! Testing if the grid is of the type of what we expect:
      rr = YY(nx/2,ny) ! highest latitude
      IF ( rr == 90.0 ) THEN
         PRINT *, 'ERROR, mod_manip.f90 => EXT_NORTH_TO_90_REGG : mhh well you shouldnt be here I guess, 90 exists!...'
         PRINT *, YY(:,ny)
         STOP
      END IF
      IF ( SUM( (YY(:,ny) - rr)**2 ) > 1.E-12 ) THEN
         PRINT *, 'ERROR, mod_manip.f90 => EXT_NORTH_TO_90_REGG : mhh well you shouldnt be here I guess, grid doesnt seem to be regular!...'
         STOP
      END IF

      !! Longitude points for the extra upper row are just the same!
      XP(:,nyp1) = XX(:,ny)
      XP(:,nyp2) = XX(:,ny)

      !! For latitude it's easy:
      YP(:,nyp1) = 90.0
      YP(:,nyp2) = 90.0 + (YY(nx/2,ny) - YY(nx/2,ny-1)) ! 90. + dlon ! lolo bad???

      DO ji=1, nx
         zlon = XX(ji,ny) ! ji => zlon
         !! at what ji do we arrive when crossing the northpole => ji_m!
         zlon_m = MOD(zlon+180.,360.)
         !PRINT *, ' zlon, zlon_m =>', zlon, zlon_m
         ji_m = MINLOC(ABS(XX(:,ny)-zlon_m), dim=1)
         !PRINT *, '  ji_m =', ji_m
         !PRINT *, 'XX(ji_m,ny) =', XX(ji_m,ny)
         !! Well so the northpole is righ in between so:
         FP(ji,nyp1) = 0.5*(XF(ji,ny) + XF(ji_m,ny)) ! lolo fix!
         FP(ji,nyp2) =  XF(ji_m,ny-1) ! lolo bad???
      END DO

      !PRINT *, 'LOLO EXT_NORTH_TO_90_REGG: YP =', YP(nx/2,:)
      !PRINT *, ''
      !PRINT *, 'LOLO EXT_NORTH_TO_90_REGG: FP =', FP(nx/2,:)
      !STOP 'boo'

   END SUBROUTINE EXT_NORTH_TO_90_REGG





END MODULE MOD_MANIP
