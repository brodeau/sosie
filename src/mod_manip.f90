MODULE MOD_MANIP

   !! Misc. manipulations and operations on 2D arrays...

   !! Author: L. Brodeau

   USE io_ezcdf, ONLY: DUMP_FIELD  ! debug
   USE mod_conf, ONLY: STOP_THIS, iverbose, rd2rad, rradE, rmissval, l_reg_src, l_reg_trg, idb, jdb

   IMPLICIT NONE

   PRIVATE


   INTERFACE flip_ud
      MODULE PROCEDURE flip_ud_1d_r4, flip_ud_1d_r8, flip_ud_2d_r4, flip_ud_3d_i1
   END INTERFACE flip_ud

   INTERFACE to_degE
      MODULE PROCEDURE to_degE_scal, to_degE_1d, to_degE_2d
   END INTERFACE to_degE

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

   PUBLIC :: EXTEND_ARRAY_2D_COOR, EXTEND_ARRAY_2D_DATA, fill_extra_north_south, &
      &      extra_2_east, extra_2_west, partial_deriv, &
      &      flip_ud, long_reorg_2d, long_reorg_3d, &
      &      distance, distance_2d, &
      &      find_nearest_point, &
      &      shrink_vector, to_degE, degE_to_degWE, &
      &      APPLY_EW_PRDCT

   !      &      ext_beyond_90_reg, l_ext2np_reg

   LOGICAL, PARAMETER :: l_force_use_of_twisted = .FALSE.


CONTAINS

   SUBROUTINE EXTEND_ARRAY_2D_COOR(k_ew, pX, pY, pXx, pYx,  is_orca_grid)
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
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pX, pY
      REAL(8), DIMENSION(:,:), INTENT(out) :: pXx, pYx
      INTEGER ,      OPTIONAL, INTENT(in)  :: is_orca_grid
      !!
      INTEGER :: nx, ny, nxx, nyx
      INTEGER :: iorca=0
      !!============================================================================
      IF( PRESENT(is_orca_grid) ) iorca = is_orca_grid

      nx = SIZE(pX ,1)
      ny = SIZE(pX ,2)
      IF( (nx  /= SIZE(pY ,1)).OR.(ny  /= SIZE(pY ,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_2D_COOR] => size of input longitude do not match!!!')

      nxx = SIZE(pXx,1)
      nyx = SIZE(pXx,2)

      IF( nxx /= nx + 4 )                                CALL STOP_THIS('[EXTEND_ARRAY_2D_COOR] => target x dim is not ni+4!!!')
      IF( nyx /= ny + 4 )                                CALL STOP_THIS('[EXTEND_ARRAY_2D_COOR] => target y dim is not nj+4!!!')
      IF( (nxx /= SIZE(pYx,1)).OR.(nyx /= SIZE(pYx,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_2D_COOR] => size of input latitude does not match longitude !!!')

      !! Initializing :
      pXx(:,:) = 0.
      pYx(:,:) = 0.

      !! Filling center of domain:
      pXx(3:nxx-2, 3:nyx-2) = pX(:,:)
      pYx(3:nxx-2, 3:nyx-2) = pY(:,:)


      !! Northern Extension
      SELECT CASE( iorca )
      CASE (4)
         pXx(2:nxx/2           ,nyx-1) = pXx(nxx:nxx-nxx/2-2:-1,nyx-5)
         pXx(nxx:nxx-nxx/2-2:-1,nyx-1) = pXx(2:nxx/2             ,nyx-5)
         pXx(2:nxx/2             ,nyx) = pXx(nxx:nxx-nxx/2-2:-1,nyx-6)
         pXx(nxx:nxx-nxx/2-2:-1,nyx)   = pXx(2:nxx/2             ,nyx-6)
         !!
         pYx(2:nxx/2           ,nyx-1) = pYx(nxx:nxx-nxx/2-2:-1,nyx-5)
         pYx(nxx:nxx-nxx/2-2:-1,nyx-1) = pYx(2:nxx/2             ,nyx-5)
         pYx(2:nxx/2             ,nyx) = pYx(nxx:nxx-nxx/2-2:-1,nyx-6)
         pYx(nxx:nxx-nxx/2-2:-1,nyx)   = pYx(2:nxx/2             ,nyx-6)
         !!
      CASE (6)
         pXx(2:nxx/2             ,nyx-1) = pXx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pXx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pXx(2:nxx/2               ,nyx-4)
         pXx(2:nxx/2               ,nyx) = pXx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pXx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pXx(2:nxx/2               ,nyx-5)
         !!
         pYx(2:nxx/2             ,nyx-1) = pYx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pYx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pYx(2:nxx/2               ,nyx-4)
         pYx(2:nxx/2               ,nyx) = pYx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pYx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pYx(2:nxx/2               ,nyx-5)
         !!
      CASE DEFAULT
         !!
         pXx(:,nyx-1) = pXx(:,nyx-3) + pXx(:,nyx-2) - pXx(:,nyx-4)
         pXx(:,nyx)   = pXx(:,nyx-2) + pXx(:,nyx-2) - pXx(:,nyx-4)
         !!
         pYx(3:nxx-2, nyx-1) = pY(:, ny-1) + pY(:,ny) - pY(:,ny-2)
         pYx(3:nxx-2, nyx)   = pY(:, ny)   + pY(:,ny) - pY(:,ny-2)
         !!
      END SELECT

      !! Southern Extension
      pXx(:, 2) = pXx(:,4) - (pXx(:,5) - pXx(:,3))
      pXx(:, 1) = pXx(:,3) - (pXx(:,5) - pXx(:,3))
      !!
      pYx(3:nxx-2, 2) = pY(:,2) - (pY(:,3) - pY(:,1))
      pYx(3:nxx-2, 1) = pY(:,1) - (pY(:,3) - pY(:,1))

      !! Western/Easten Extensions
      IF(k_ew > -1) THEN   ! we can use east-west periodicity of input file to
         CALL APPLY_EW_PRDCT( k_ew, 2, pXx(:,:), l_is_longitude=.TRUE. )
         CALL APPLY_EW_PRDCT( k_ew, 2, pYx(:,:) )
      ELSE
         !!
         !!LOLO: in the fowlowing the 4 corners are not treated?
         pXx(2, 3:nyx-2) = pX(2,:) - (pX(3,:) - pX(1,:))
         pXx(1, 3:nyx-2) = pX(1,:) - (pX(3,:) - pX(1,:))
         pYx(2, :) = pYx(4,:) - (pYx(5,:) - pYx(3,:))
         pYx(1, :) = pYx(3,:) - (pYx(5,:) - pYx(3,:))
         !!
         pXx(nxx-1, 3:nyx-2) = pX(nx-1,:) + pX(nx,:) - pX(nx-2,:)
         pXx(nxx  , 3:nyx-2) = pX(nx,:)   + pX(nx,:) - pX(nx-2,:)
         pYx(nxx-1,:) = pYx(nxx-3,:) + pYx(nxx-2, :) - pYx(nxx-4, :)
         pYx(nxx,:)   = pYx(nxx-2,:) + pYx(nxx-2,:)  - pYx(nxx-4, :)
         !!
      END IF !IF(k_ew > -1)

   END SUBROUTINE EXTEND_ARRAY_2D_COOR


   SUBROUTINE EXTEND_ARRAY_2D_DATA( k_ew, pXx, pYx, pF, pFx,  is_orca_grid, l_smart_NP )
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
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pXx, pYx ! extended coordinates !
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pF
      REAL(8), DIMENSION(:,:), INTENT(out) :: pFx
      INTEGER,       OPTIONAL, INTENT(in)  :: is_orca_grid
      LOGICAL,       OPTIONAL, INTENT(in)  :: l_smart_NP
      !!
      INTEGER :: nx, ny, nxx, nyx
      INTEGER :: ji, jj, iext_np
      INTEGER :: iorca=0
      INTEGER :: lsNoP=.false.
      REAL(8), DIMENSION(:,:), ALLOCATABLE  :: ztmp
      !!============================================================================

      IF( PRESENT(is_orca_grid) ) iorca = is_orca_grid
      IF( PRESENT( l_smart_NP)  ) lsNoP = l_smart_NP

      nx  = SIZE(pF ,1) ; ny  = SIZE(pF ,2)
      nxx = SIZE(pFx,1) ; nyx = SIZE(pFx,2)

      IF( nxx /= nx + 4 )                            CALL STOP_THIS('[EXTEND_ARRAY_2D_DATA] => target x dim is not ni+4!!!')
      IF( nyx /= ny + 4 )                            CALL STOP_THIS('[EXTEND_ARRAY_2D_DATA] => target y dim is not nj+4!!!')
      IF( (nxx /= SIZE(pXx,1)).OR.(nyx /= SIZE(pXx,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_2D_DATA] => size of input longitude do not match!!!')
      IF( (nxx /= SIZE(pYx,1)).OR.(nyx /= SIZE(pYx,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_2D_DATA] => size of input latitude  do not match!!!')

      !! Initializing :
      pFx(:,:) = 0.

      !! Filling center of domain:
      pFx(3:nxx-2, 3:nyx-2)             = pF(:,:)

      IF(k_ew > -1) THEN
         CALL APPLY_EW_PRDCT( k_ew, 2, pFx(:,3:nyx-2) )
      ELSE
         !! East:
         DO jj = 3, nyx-2
            CALL extra_2_east(pXx(nxx-4,jj),pXx(nxx-3,jj),pXx(nxx-2,jj),        &
               &              pXx(nxx-1,jj),pXx(nxx,jj),                         &
               &              pFx(nxx-4,jj),pFx(nxx-3,jj),pFx(nxx-2,jj),  &
               &              pFx(nxx-1,jj),pFx(nxx,jj) )
         END DO
         !! West:
         DO jj = 3, nyx-2
            CALL extra_2_west(pXx(5,jj),pXx(4,jj),pXx(3,jj),         &
               &              pXx(2,jj),pXx(1,jj),                   &
               &              pFx(5,jj),pFx(4,jj),pFx(3,jj),   &
               &              pFx(2,jj),pFx(1,jj) )
         END DO
      END IF
      !!
      !! Top:
      SELECT CASE( iorca )
      CASE (4)
         IF(iverbose>0) PRINT *, 'ORCA north pole T-point folding type of extrapolation at northern boundary!'
         pFx(2:nxx/2             ,nyx-1) = pFx(nxx:nxx-nxx/2-2:-1,nyx-5)
         pFx(nxx:nxx-nxx/2-2:-1,nyx-1) = pFx(2:nxx/2             ,nyx-5)
         pFx(2:nxx/2             ,nyx)   = pFx(nxx:nxx-nxx/2-2:-1,nyx-6)
         pFx(nxx:nxx-nxx/2-2:-1,nyx)   = pFx(2:nxx/2             ,nyx-6)
         !!
      CASE (6)
         IF(iverbose>0) PRINT *, 'ORCA north pole F-point folding type of extrapolation at northern boundary!'
         pFx(2:nxx/2               ,nyx-1) = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pFx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pFx(2:nxx/2               ,nyx-4)
         pFx(2:nxx/2               ,nyx)   = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pFx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pFx(2:nxx/2               ,nyx-5)
         !!
      CASE DEFAULT

         !! Before defaulting to extrapolations, there is the case when a regular source grid ends at the North-Pole:
         iext_np = I_EXT2NP_REG(pXx(3:nxx-2,3:nyx-2),pYx(3:nxx-2,3:nyx-2))
         IF( (lsNoP).AND.(iext_np>0) ) THEN
            IF(iverbose>0) PRINT *, 'Going for smart `EXT_BEYOND_90_REG` NorthPole extrapolation!'
            ALLOCATE( ztmp(nxx,ny+2) )
            ztmp(:,:) = pFx(:,:nyx-2)
            CALL EXT_BEYOND_90_REG( iext_np, pXx(:,:nyx-2), pYx(:,:nyx-2), ztmp(:,:), pFx(:,:) )
            DEALLOCATE( ztmp )
            !!
         ELSE
            !! Defaulting to extrapolation:
            IF(iverbose>0) PRINT *, 'Going for default dumb NorthPole extrapolation!'
            DO ji = 1, nxx
               CALL extra_2_east(pYx(ji,nyx-4),pYx(ji,nyx-3),pYx(ji,nyx-2),       &
                  &              pYx(ji,nyx-1),pYx(ji,nyx),                        &
                  &              pFx(ji,nyx-4),pFx(ji,nyx-3),pFx(ji,nyx-2), &
                  &              pFx(ji,nyx-1),pFx(ji,nyx) )
            END DO
         END IF

      END SELECT
      !!
      !! Bottom:
      DO ji = 1, nxx
         CALL extra_2_west(pYx(ji,5),pYx(ji,4),pYx(ji,3),       &
            &              pYx(ji,2),pYx(ji,1),                 &
            &              pFx(ji,5),pFx(ji,4),pFx(ji,3), &
            &              pFx(ji,2),pFx(ji,1) )
      END DO
      !!
   END SUBROUTINE EXTEND_ARRAY_2D_DATA




   SUBROUTINE FILL_EXTRA_NORTH_SOUTH(pX, YY, pF, pXx, pYx, pFx,  is_orca_grid)

      !!============================================================================
      !! Extending input arrays with an extraband of two points at northern and
      !! southern boundaries.
      !!
      !! The extension is done thanks to Akima's exptrapolation method or continuity
      !!  / ORCA knowledge...
      !!
      !!============================================================================
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pX, YY, pF
      REAL(8), DIMENSION(:,:), INTENT(out) :: pXx, pYx, pFx
      INTEGER ,   OPTIONAL,    INTENT(in)  :: is_orca_grid

      !! Local
      INTEGER :: nx, ny, nxx, nyx
      INTEGER :: ji, iorca

      iorca = 0
      IF( PRESENT(is_orca_grid) ) iorca = is_orca_grid

      nx = SIZE(pX,1)
      ny = SIZE(pX,2)
      nxx = SIZE(pXx,1)
      nyx = SIZE(pXx,2)

      IF( (nx /= SIZE(YY,1)).OR.(ny /= SIZE(YY,2)).OR. &
         & (nx /= SIZE(pF,1)).OR.(ny /= SIZE(pF,2))) THEN
         CALL STOP_THIS('[FILL_EXTRA_NORTH_SOUTH] => FILL_EXTRA_NORTH_SOUTH : size of input coor. and data do not match!!!')
      END IF

      IF( (nxx /= SIZE(pYx,1)).OR.(nyx /= SIZE(pYx,2)).OR. &
         & (nxx /= SIZE(pFx,1)).OR.(nyx /= SIZE(pFx,2))) THEN
         CALL STOP_THIS('[FILL_EXTRA_NORTH_SOUTH] => FILL_EXTRA_NORTH_SOUTH : size of output coor. and data do not match!!!')
      END IF


      IF( nxx /= nx   ) CALL STOP_THIS('[FILL_EXTRA_NORTH_SOUTH] => target x dim is not ni!!!')
      IF( nyx /= ny+4 ) CALL STOP_THIS('[FILL_EXTRA_NORTH_SOUTH] => target y dim is not nj+4!!!')


      !!   C r e a t i n g   e x t e n d e d   a r r a y s  :
      !!   --------------------------------------------------

      !! Initializing :
      pXx = 0.
      pYx = 0.
      pFx = 0.

      !! Filling center of domain:
      pXx(:, 3:nyx-2) = pX(:,:)
      pYx(:, 3:nyx-2) = YY(:,:)
      pFx(:, 3:nyx-2) = pF(:,:)



      !! X array :
      !! ---------

      !! Bottom side :
      pXx(:, 2) = pXx(:,4) - (pXx(:,5) - pXx(:,3))
      pXx(:, 1) = pXx(:,3) - (pXx(:,5) - pXx(:,3))

      !! Top side :
      SELECT CASE( iorca )
         !
      CASE (4)
         pXx(2:nx/2             ,nyx-1) = pXx(nx:nx-nx/2-2:-1,nyx-5)
         pXx(nx:nx-nx/2-2:-1,nyx-1) = pXx(2:nx/2             ,nyx-5)
         pXx(2:nx/2             ,nyx)   = pXx(nx:nx-nx/2-2:-1,nyx-6)
         pXx(nx:nx-nx/2-2:-1,nyx)   = pXx(2:nx/2             ,nyx-6)
      CASE (6)
         pXx(2:nx/2               ,nyx-1) = pXx(nx-1:nx-nx/2+1:-1,nyx-4)
         pXx(nx-1:nx-nx/2+1:-1,nyx-1) = pXx(2:nx/2               ,nyx-4)
         pXx(2:nx/2               ,nyx)   = pXx(nx-1:nx-nx/2+1:-1,nyx-5)
         pXx(nx-1:nx-nx/2+1:-1,nyx)   = pXx(2:nx/2               ,nyx-5)
      CASE DEFAULT
         pXx(:,nyx-1) = pXx(:,nyx-3) + pXx(:,nyx-2) - pXx(:,nyx-4)
         pXx(:,nyx)   = pXx(:,nyx-2) + pXx(:,nyx-2) - pXx(:,nyx-4)

      END SELECT


      !! Y array :
      !! ---------

      !! Top side :
      SELECT CASE( iorca )

      CASE (4)
         pYx(2:nx/2             ,nyx-1) = pYx(nx:nx-nx/2-2:-1,nyx-5)
         pYx(nx:nx-nx/2-2:-1,nyx-1) = pYx(2:nx/2             ,nyx-5)
         pYx(2:nx/2             ,nyx)   = pYx(nx:nx-nx/2-2:-1,nyx-6)
         pYx(nx:nx-nx/2-2:-1,nyx)   = pYx(2:nx/2             ,nyx-6)
      CASE (6)
         pYx(2:nx/2               ,nyx-1) = pYx(nx-1:nx-nx/2+1:-1,nyx-4)
         pYx(nx-1:nx-nx/2+1:-1,nyx-1) = pYx(2:nx/2               ,nyx-4)
         pYx(2:nx/2               ,nyx)   = pYx(nx-1:nx-nx/2+1:-1,nyx-5)
         pYx(nx-1:nx-nx/2+1:-1,nyx)   = pYx(2:nx/2               ,nyx-5)
      CASE DEFAULT
         pYx(:, nyx-1) = YY(:, ny-1) + YY(:,ny) - YY(:,ny-2)
         pYx(:, nyx)   = YY(:, ny)   + YY(:,ny) - YY(:,ny-2)
      END SELECT


      !! Bottom side :
      pYx(:, 2) = YY(:,2) - (YY(:,3) - YY(:,1))
      pYx(:, 1) = YY(:,1) - (YY(:,3) - YY(:,1))



      !! Data array :
      !! ------------

      !! Top side :
      SELECT CASE( iorca )

      CASE (4)
         PRINT *, 'ORCA2 type of extrapolation at northern boundary!'
         pFx(2:nx/2             ,nyx-1) = pFx(nx:nx-nx/2-2:-1,nyx-5)
         pFx(nx:nx-nx/2-2:-1,nyx-1) = pFx(2:nx/2             ,nyx-5)
         pFx(2:nx/2             ,nyx)   = pFx(nx:nx-nx/2-2:-1,nyx-6)
         pFx(nx:nx-nx/2-2:-1,nyx)   = pFx(2:nx/2             ,nyx-6)
      CASE (6)
         PRINT *, 'ORCA1 type of extrapolation at northern boundary!'
         pFx(2:nx/2               ,nyx-1) = pFx(nx-1:nx-nx/2+1:-1,nyx-4)
         pFx(nx-1:nx-nx/2+1:-1,nyx-1) = pFx(2:nx/2               ,nyx-4)
         pFx(2:nx/2               ,nyx)   = pFx(nx-1:nx-nx/2+1:-1,nyx-5)
         pFx(nx-1:nx-nx/2+1:-1,nyx)   = pFx(2:nx/2               ,nyx-5)
      CASE DEFAULT
         DO ji = 1, nx
            CALL extra_2_east(pYx(ji,nyx-4),pYx(ji,nyx-3),pYx(ji,nyx-2),       &
               &              pYx(ji,nyx-1),pYx(ji,nyx),                        &
               &              pFx(ji,nyx-4),pFx(ji,nyx-3),pFx(ji,nyx-2), &
               &              pFx(ji,nyx-1),pFx(ji,nyx) )
         END DO
      END SELECT

      !! Bottom side :
      DO ji = 1, nx
         CALL extra_2_west(pYx(ji,5),pYx(ji,4),pYx(ji,3),       &
            &              pYx(ji,2),pYx(ji,1),                 &
            &              pFx(ji,5),pFx(ji,4),pFx(ji,3), &
            &              pFx(ji,2),pFx(ji,1) )
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
      IF( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
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
      IF( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
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
      IF( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
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
      IF( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
         y2 = y3
         y1 = y3
      ELSE
         y2 = C*(2*BET/B - ALF/A) + y3
         y1 = y2 + y2*D/C + BET*D/B - ALF*D/A - y3*D/C
      END IF
   END SUBROUTINE extra_2_west_r4



   SUBROUTINE PARTIAL_DERIV(k_ew, pX, pY, pF, dFdX, dFdY, d2FdXdY)

      !! Partial derivatives of a field ZF given on a regular gird !!!

      !!  k_ew : east-west periodicity on the input file/grid
      !!         k_ew = -1  --> no east-west periodicity (along x)
      !!         k_ew >= 0  --> east-west periodicity with overlap of k_ew points (along x)

      INTEGER, INTENT(in) :: k_ew

      REAL(8), DIMENSION(:,:), INTENT(in)  :: pX, pY, pF
      REAL(8), DIMENSION(:,:), INTENT(out) :: dFdX, dFdY, d2FdXdY

      !! Local variables :
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: ZX, ZY, ZF
      INTEGER :: nx, ny

      dFdX    = 0.
      dFdY    = 0.
      d2FdXdY = 0.

      nx = SIZE(pF,1) ; ny = SIZE(pF,2)

      !! Extended arrays with a frame of 2 points...
      ALLOCATE ( ZX(nx+4,ny+4), ZY(nx+4,ny+4), ZF(nx+4,ny+4) )

      CALL EXTEND_ARRAY_2D_COOR(k_ew, pX, pY, ZX, ZY)
      CALL EXTEND_ARRAY_2D_DATA(k_ew, ZX, ZY, pF, ZF)

      !! Second order finite difference:
      !! i+1 => 4:nx+4 / i => 3:nx+2 / i-1 => 2:nx+2
      !! j+1 => 4:ny+4 / j => 3:ny+2 / j-1 => 2:ny+2


      dFdX(:,:) = 0.5*( ZF(4:nx+4,3:ny+2) - ZF(2:nx+2,3:ny+2) ) !/ ( ZX(4:nx+4,3:ny+2) - ZX(2:nx+2,3:ny+2) )
      dFdY(:,:) = 0.5*( ZF(3:nx+2,4:ny+4) - ZF(3:nx+2,2:ny+2) ) !/ ( ZY(3:nx+2,4:ny+4) - ZY(3:nx+2,2:ny+2) )

      !!lolo: Here the denominator is probably wrong:
      d2FdXdY(:,:) = 0.25*( ZF(4:nx+4,4:ny+4) - ZF(4:nx+4,2:ny+2)   -   ZF(2:nx+2,4:ny+4) + ZF(2:nx+2,2:ny+2) ) !&

      DEALLOCATE ( ZX, ZY, ZF )

   END SUBROUTINE PARTIAL_DERIV


   SUBROUTINE FLIP_UD_1D_R4(pF)
      REAL(4), DIMENSION(:), INTENT(inout) :: pF
      INTEGER :: nz, jk
      REAL(4), DIMENSION(:), ALLOCATABLE :: ztmp
      nz = SIZE(pF,1)
      ALLOCATE ( ztmp(nz) )
      ztmp(:) = pF(:)
      DO jk = 1, nz
         pF(jk) =  ztmp(nz-jk+1)
      END DO
      DEALLOCATE ( ztmp )
   END SUBROUTINE FLIP_UD_1D_R4

   SUBROUTINE FLIP_UD_1D_R8(pF)
      REAL(8), DIMENSION(:), INTENT(inout) :: pF
      INTEGER :: nz, jk
      REAL(8), DIMENSION(:), ALLOCATABLE :: ztmp
      nz = SIZE(pF,1)
      ALLOCATE ( ztmp(nz) )
      ztmp(:) = pF(:)
      DO jk = 1, nz
         pF(jk) =  ztmp(nz-jk+1)
      END DO
      DEALLOCATE ( ztmp )
   END SUBROUTINE FLIP_UD_1D_R8


   SUBROUTINE FLIP_UD_2D_R4(pF)

      REAL(4), DIMENSION(:,:), INTENT(inout) :: pF

      INTEGER :: nx, ny, jj
      REAL(4), DIMENSION(:,:), ALLOCATABLE :: ztmp

      nx = SIZE(pF,1) ; ny = SIZE(pF,2)

      ALLOCATE ( ztmp(nx,ny) )

      ztmp(:,:) = pF(:,:)

      DO jj = 1, ny
         pF(:,jj) =  ztmp(:,ny-jj+1)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE FLIP_UD_2D_R4


   SUBROUTINE FLIP_UD_3D_I1(pF)

      INTEGER(1), DIMENSION(:,:,:), INTENT(inout) :: pF

      INTEGER :: nx, ny, nz, jj
      INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: ztmp

      nx = SIZE(pF,1) ; ny = SIZE(pF,2) ; nz = SIZE(pF,3)

      ALLOCATE ( ztmp(nx,ny,nz) )

      ztmp(:,:,:) = pF(:,:,:)

      DO jj = 1, ny
         pF(:,jj,:) =  ztmp(:,ny-jj+1,:)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE FLIP_UD_3D_I1





   SUBROUTINE LONG_REORG_2D(i_chg_x, pF)

      INTEGER, INTENT(in) :: i_chg_x
      REAL(4), DIMENSION(:,:), INTENT(inout) :: pF

      INTEGER :: nx, ny, ji
      REAL(4), DIMENSION(:,:), ALLOCATABLE :: ztmp

      nx = SIZE(pF,1) ; ny = SIZE(pF,2)

      ALLOCATE ( ztmp(nx,ny) )

      ztmp(:,:) = pF(:,:)

      DO ji = i_chg_x, nx
         pF(ji - i_chg_x + 1 , :) = ztmp(ji,:)
      END DO
      DO ji = 1, i_chg_x - 1
         pF(nx - i_chg_x + 1 + ji , :) = ztmp(ji,:)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE LONG_REORG_2D



   SUBROUTINE LONG_REORG_3D_I1(i_chg_x, pF)

      INTEGER, INTENT(in) :: i_chg_x
      INTEGER(1), DIMENSION(:,:,:), INTENT(inout) :: pF

      INTEGER :: nx, ny, nz, ji
      INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: ztmp

      nx = SIZE(pF,1) ; ny = SIZE(pF,2) ; nz = SIZE(pF,3)

      ALLOCATE ( ztmp(nx,ny,nz) )

      ztmp(:,:,:) = pF(:,:,:)

      DO ji = i_chg_x, nx
         pF(ji - i_chg_x + 1 , :,:) = ztmp(ji,:,:)
      END DO
      DO ji = 1, i_chg_x - 1
         pF(nx - i_chg_x + 1 + ji , :,:) = ztmp(ji,:,:)
      END DO

      DEALLOCATE ( ztmp )

   END SUBROUTINE LONG_REORG_3D_I1






   SUBROUTINE FIND_NEAREST_POINT(Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp,  pmsk_dom_trg)
      !!---------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_POINT  ***
      !!
      !!   Source domain: Xsrc, Ysrc
      !!   Target domain: Xtrg, Ytrg
      !!     => for each point of target domain: i_s,j_s location of nearest point on source domain
      !!        => JIp & JJp have same shape as Xtrg & Ytrg !
      !!
      !!                Laurent Brodeau, July, 2017
      !!
      !!
      !! OPTIONAL:
      !!      * pmsk_dom_trg: ignore (dont't treat) regions of the target domain where pmsk_dom_trg==0 !
      !!---------------------------------------------------------------
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xtrg, Ytrg    !: lon and lat arrays of target domain
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xsrc , Ysrc     !: lon and lat arrays of source domain
      INTEGER(4), DIMENSION(:,:), INTENT(out) :: JIp, JJp  !: nearest point location of point P in Xsrc,Ysrc wrt Xtrg,Ytrg
      !INTEGER,    OPTIONAL,       INTENT(in)  :: ithread
      INTEGER(1), OPTIONAL ,DIMENSION(:,:), INTENT(inout) :: pmsk_dom_trg

      INTEGER :: jj, nx_s, ny_s, nx_t, ny_t, j_strt_t, j_stop_t, jlat_icr
      REAL(8) :: y_max_s, y_min_s, rmin_dlat_dj, rtmp
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: mask_ignore_t
      INTEGER,    DIMENSION(:),   ALLOCATABLE :: i1dum
      REAL(8),    DIMENSION(:,:), ALLOCATABLE :: ztmp_t
      LOGICAL :: l_is_reg_s, l_is_reg_t
      CHARACTER(len=128) :: cmsg

      !ithrd = 1 ! no OpenMP !
      !IF( PRESENT(ithread) ) ithrd = ithread

      nx_s  = SIZE(Xsrc,1)
      ny_s  = SIZE(Xsrc,2)
      nx_t = SIZE(Xtrg,1)
      ny_t = SIZE(Xtrg,2)

      PRINT *, ' Source domain size: ', nx_s, ny_s
      PRINT *, ' Target domain size: ', nx_t, ny_t

      cmsg = '`FIND_NEAREST_POINT` of mod_manip.f90 =>'
      IF((SIZE(Ysrc,1)/= nx_s).OR.(SIZE(Ysrc,2)/= ny_s)) CALL STOP_THIS(cmsg//' Ysrc dont agree in shape with Xsrc')
      IF((SIZE(Ytrg,1)/= nx_t).OR.(SIZE(Ytrg,2)/= ny_t)) CALL STOP_THIS(cmsg//' Ytrg dont agree in shape with Xtrg')
      IF((SIZE(JIp,1) /= nx_t).OR.(SIZE(JIp,2) /= ny_t)) CALL STOP_THIS(cmsg//' JIp dont agree in shape with Xtrg')
      IF((SIZE(JJp,1) /= nx_t).OR.(SIZE(JJp,2) /= ny_t)) CALL STOP_THIS(cmsg//' JJp dont agree in shape with Xtrg')

      PRINT *, ''
      !! Checking if source domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_s = L_IS_GRID_REGULAR( Xsrc , Ysrc )
      PRINT *, ' *** FIND_NEAREST_POINT => Is source grid regular ??? =>', l_is_reg_s
      IF( l_reg_src .AND. (.NOT. l_is_reg_s) )  CALL STOP_THIS( ' ==> BUT namelist says `l_reg_src=T`!')

      !! Checking if target domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_t = L_IS_GRID_REGULAR( Xtrg , Ytrg )
      PRINT *, ' *** FIND_NEAREST_POINT => Is target grid regular ??? =>', l_is_reg_t
      IF( l_reg_trg .AND. (.NOT. l_is_reg_t) )  CALL STOP_THIS( ' ==> BUT namelist says `l_reg_trg=T`!')

      ALLOCATE ( mask_ignore_t(nx_t,ny_t) , i1dum(nx_t) )
      mask_ignore_t(:,:) = 1
      IF( PRESENT( pmsk_dom_trg ) ) mask_ignore_t(:,:) = pmsk_dom_trg(:,:)

      y_min_s  = MINVAL(Ysrc) ; ! Min and Max latitude of source domain
      y_max_s  = MAXVAL(Ysrc)

      !! General case:
      j_strt_t = 1
      j_stop_t = ny_t
      jlat_icr   = 1

      !! We need to know if the target latitude ONLY keeps on systematically
      !! increasing (or decreasing) as j increases:
      rmin_dlat_dj = -100.
      !IF( REAL(nx_t)/REAL(ny_t) > 0.15  ) THEN
      IF( nx_t > 2 ) THEN !! We can reasonably say it's a map and not a vector or almost a vectot...
         ALLOCATE ( ztmp_t(nx_t, ny_t) )
         DO jj = 2, ny_t
            ztmp_t(:,jj) = Ytrg(:,jj) - Ytrg(:,jj-1)
         END DO
         rtmp = SUM(ztmp_t(:,ny_t/2)) ! know if increasing (>0) or decreasing (<0)
         ztmp_t = SIGN(1.0_8 , rtmp)*ztmp_t
         !IF(iverbose==2) CALL DUMP_FIELD(REAL(ztmp_t,4), 'dlat_dj_t.nc', 'dist')
         rmin_dlat_dj = MINVAL(ztmp_t(:,2:))
         IF(iverbose==2) PRINT *, ' Minimum dlat_dj_t =>', rmin_dlat_dj
         DEALLOCATE ( ztmp_t )
      END IF

      !!  Simplif when [d lat / d j] always has the same sign:
      IF( (rmin_dlat_dj > -1.E-12) .OR. l_is_reg_t ) THEN    !!!.OR. (i_orca_t > 0) ) THEN
         !! -> because we need to avoid all the following if target grid is for
         !!    example a polar sterographic projection... (example 5)
         !!
         !! *** Will ignore regions of the TARGET domain that
         !!     are not covered by source domain:
         !ij_min_loc = MINLOC(Ytrg, mask=(Ytrg>=y_min_s))  ! possible bug identified by Max (aka MB)
         !ij_max_loc = MAXLOC(Ytrg, mask=(Ytrg<=y_max_s))
         !j_strt_t = MINVAL(MINLOC(Ytrg, mask=(Ytrg>=y_min_s), dim=2))  ! smallest j on target source that covers smallest source latitude
         i1dum = MINLOC(Ytrg, mask=(Ytrg>=y_min_s), dim=2)
         j_strt_t = MINVAL(i1dum, mask=(i1dum>0))
         DEALLOCATE ( i1dum )
         j_stop_t = MAXVAL(MAXLOC(Ytrg, mask=(Ytrg<=y_max_s), dim=2))  ! largest j on target source that covers largest source latitude
         IF( j_strt_t > j_stop_t ) jlat_icr = -1 ! latitude decreases as j increases (like ECMWF grids...)
         IF(iverbose==2) THEN
            PRINT *, ' j_strt_t, j_stop_t / nj_t =>', j_strt_t, j_stop_t, '/', ny_t
            PRINT *, ''
         END IF
      END IF ! IF( (rmin_dlat_dj >= 0.0_8) .OR. l_is_reg_t .OR. (i_orca_t > 0) )

      IF( l_is_reg_s .AND. (.NOT. l_force_use_of_twisted) ) THEN
         PRINT *, '                 => going for simple FIND_NEAREST algorithm !'; PRINT *, ''
         CALL FIND_NEAREST_EASY(    Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp, mask_ignore_t, &
            &                       j_strt_t, j_stop_t, jlat_icr )
      ELSE
         PRINT *, '                  => going for advanced FIND_NEAREST algorithm !'; PRINT *, ''
         CALL FIND_NEAREST_TWISTED( Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp, mask_ignore_t, &
            &                       j_strt_t, j_stop_t, jlat_icr )
      END IF
      PRINT *, ''

      IF( PRESENT( pmsk_dom_trg ) ) pmsk_dom_trg(:,:) = mask_ignore_t(:,:)
      DEALLOCATE ( mask_ignore_t )

   END SUBROUTINE FIND_NEAREST_POINT







   SUBROUTINE FIND_NEAREST_EASY( Xtrg, Ytrg, Xsrc, Ysrc, JIpos, JJpos, mask_t, &
      &                          j_strt_t, j_stop_t, jlat_icr )
      !!---------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_EASY  ***
      !!
      !!   Source domain: Xsrc, Ysrc
      !!   Target domain: Xtrg, Ytrg
      !!     => for each point of target domain: i_s,j_s location of nearest point on source domain
      !!        => JIpos & JJpos have same shape as Xtrg & Ytrg !
      !!
      !!                Laurent Brodeau, July, 2017
      !!
      !!
      !!      * mask_t: ignore (dont't treat) regions of the target domain where mask_t==0 !
      !!---------------------------------------------------------------
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xtrg, Ytrg    !: lon and lat arrays of target domain
      REAL(8),    DIMENSION(:,:), INTENT(in)  :: Xsrc , Ysrc     !: lon and lat arrays of source domain
      INTEGER(4), DIMENSION(:,:), INTENT(out) :: JIpos, JJpos  !: nearest point location of point P in Xsrc,Ysrc wrt Xtrg,Ytrg
      INTEGER(1), DIMENSION(:,:), INTENT(in)  :: mask_t
      INTEGER,                    INTENT(in)  :: j_strt_t, j_stop_t, jlat_icr

      !! Important parameters:
      INTEGER, PARAMETER :: nframe_scan = 4  ! domain to scan for nearest point in simple algo => domain of 9x9

      INTEGER :: &
         &    nx_s, ny_s, nx_t, ny_t, &
         &    ji_t, jj_t, ji_s, jj_s, &
         &    j1s, j2s, i1s, i2s

      REAL(8) :: rlon, rlat

      INTEGER, DIMENSION(2) :: ij_min_loc
      INTEGER, DIMENSION(1) :: ip, jp

      REAL(8), DIMENSION(:),   ALLOCATABLE :: vlat, vlon
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: Xdist

      nx_s  = SIZE(Xsrc,1)
      ny_s  = SIZE(Xsrc,2)
      nx_t = SIZE(Xtrg,1)
      ny_t = SIZE(Xtrg,2)

      JIpos(:,:) = -1
      JJpos(:,:) = -1

      ALLOCATE ( vlon(nx_s) , vlat(ny_s) , Xdist(nx_s,ny_s) )

      vlon(:) = Xsrc(:,3)   ! 3 => make sure we're inside the domain
      vlat(:) = Ysrc(3,:)   ! 3 => make sure we're inside the domain

      DO jj_t = j_strt_t, j_stop_t, jlat_icr
         DO ji_t = 1, nx_t

            IF( mask_t(ji_t,jj_t) == 1 ) THEN

               IF( (ji_t==1) .AND. MOD(jj_t,10)==0 ) &
                  &  WRITE(6,'(" *** Treated j-point of target domain = ",i4.4," (thread # ",i2.2,")")') jj_t, 0

               rlon = Xtrg(ji_t,jj_t)
               rlat = Ytrg(ji_t,jj_t)

               ! Nearest point (in terms of index):
               ip =  MINLOC(ABS(vlon(:)-rlon))
               jp =  MINLOC(ABS(vlat(:)-rlat))

               !! Define the box to scan for shortest distance:
               i1s = MAX(ip(1)-nframe_scan ,  1)
               i2s = MIN(ip(1)+nframe_scan , nx_s)
               j1s = MAX(jp(1)-nframe_scan ,  1)
               j2s = MIN(jp(1)+nframe_scan , ny_s)

               ! Nearest point (in terms of distance):
               Xdist = 1.E12
               Xdist(i1s:i2s,j1s:j2s) = DISTANCE_2D(rlon, Xsrc(i1s:i2s,j1s:j2s), rlat, Ysrc(i1s:i2s,j1s:j2s))
               ij_min_loc = MINLOC(Xdist(i1s:i2s,j1s:j2s))
               ji_s = ij_min_loc(1) + i1s - 1
               jj_s = ij_min_loc(2) + j1s - 1
               JIpos(ji_t,jj_t) = ji_s
               JJpos(ji_t,jj_t) = jj_s
               IF((ji_s==0).OR.(jj_s==0)) THEN
                  PRINT *, ''
                  PRINT *, 'The nearest point was not found!'
                  PRINT *, ' !!! ji_s or jj_s = 0 !!!'
                  PRINT *, ' ** Target point (lon,lat) =>', rlon, rlat
                  PRINT *, ' Nearest point found on source grid:', &
                     &  REAL(Xsrc(ji_s,jj_s) , 4), &
                     &  REAL(Ysrc(ji_s,jj_s) , 4)
                  STOP
               END IF

               IF( ( JIpos(ji_t,jj_t) == -1 ).OR.( JJpos(ji_t,jj_t) == -1 ) ) THEN
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



   SUBROUTINE FIND_NEAREST_TWISTED( Xtrg, Ytrg, Xsrc, Ysrc, JIpos, JJpos, mask_t, &
      &                             j_strt_t, j_stop_t, jlat_icr )
      !!---------------------------------------------------------------
      !!            ***  SUBROUTINE FIND_NEAREST_POINT  ***
      !!
      !!   Source domain: Xsrc, Ysrc
      !!   Target domain: Xtrg, Ytrg
      !!     => for each point of target domain: i_s,j_s location of nearest point on source domain
      !!        => JIpos & JJpos have same shape as Xtrg & Ytrg !
      !!
      !!                Laurent Brodeau, July, 2017
      !!
      !!
      !!      * mask_t: ignore (dont't treat) regions of the target domain where mask_t==0 !
      !!---------------------------------------------------------------
      REAL(8),    DIMENSION(:,:), INTENT(in)    :: Xtrg, Ytrg    !: lon and lat arrays of target domain
      REAL(8),    DIMENSION(:,:), INTENT(in)    :: Xsrc , Ysrc     !: lon and lat arrays of source domain
      INTEGER(4), DIMENSION(:,:), INTENT(out)   :: JIpos, JJpos  !: nearest point location of point P in Xsrc,Ysrc wrt Xtrg,Ytrg
      INTEGER(1), DIMENSION(:,:), INTENT(inout) :: mask_t
      INTEGER, INTENT(in)                       :: j_strt_t, j_stop_t, jlat_icr
      !!
      !! Important parameters:
      INTEGER, PARAMETER :: Nlat_split = 40   ! number of latitude bands to split the search work (for 180. degree south->north)
      REAL(8), PARAMETER :: frac_emax = 0.51  ! fraction of emax to test if found!
      !!                                        => 0.5 seems to be too small, for example on ORCA1 grid at around 20 deg N... ???
      INTEGER :: &
         &    nx_s, ny_s, nx_t, ny_t, &
         &    jlat, ji_t, jj_t, ji_s, jj_s, jj_t_old, &
         &    j1s, j2s, i1s, i2s, niter, &
         &    nsplit, jmax_band, jmin_band

      REAL(8) :: emax, rlat_low, rlat_hgh, rlon, rlat, rlat_old
      REAL(8) :: y_max_t, y_min_t, y_max_bnd, y_min_bnd, y_max_bnd0, y_min_bnd0, dy, &
         &       y_max_s, y_min_s
      REAL(8),    DIMENSION(:),   ALLOCATABLE :: VLAT_SPLIT_BOUNDS
      REAL(8),    DIMENSION(:,:), ALLOCATABLE :: Xdist, e1_s, e2_s    !: grid layout and metrics
      INTEGER,    DIMENSION(:,:), ALLOCATABLE :: J_VLAT_S
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: mspot_lon, mspot_lat
      INTEGER,    DIMENSION(:),   ALLOCATABLE :: i1dum

      INTEGER, DIMENSION(2) :: ij_min_loc

      LOGICAL :: lagain, lPdbg


      nx_s  = SIZE(Xsrc,1)
      ny_s  = SIZE(Xsrc,2)
      nx_t = SIZE(Xtrg,1)
      ny_t = SIZE(Xtrg,2)

      ALLOCATE ( Xdist(nx_s,ny_s) , mspot_lon(nx_s,ny_s) , mspot_lat(nx_s,ny_s) , i1dum(nx_s) )

      y_min_s  = MINVAL(Ysrc) ; ! Min and Max latitude of source domain
      y_max_s  = MAXVAL(Ysrc)
      y_min_t = MINVAL(Ytrg) ; ! Min and Max latitude of target domain
      y_max_t = MAXVAL(Ytrg)


      !! Control if source & target domains overlap:
      PRINT *, ' Min. latitude on source & target domains =>', REAL(y_min_s,4), REAL(y_min_t,4)
      PRINT *, ' Max. latitude on source & target domains =>', REAL(y_max_s,4), REAL(y_max_t,4)

      IF( (y_min_t>=y_max_s).OR.(y_max_t<=y_min_s) ) &
         & CALL STOP_THIS( '`FIND_NEAREST_TWISTED@mod_manip.f90`: Target and source latitudes do not overlap!' )
      !! ---------------------------------------------------------------------------------------
      PRINT *, '       => going for advanced algorithm !' ;! (thread #', ithrd, ')'

      JIpos(:,:) = -1
      JJpos(:,:) = -1

      ALLOCATE ( e1_s(nx_s,ny_s), e2_s(nx_s,ny_s) )

      !! We need metric of input grid
      e1_s(:,:) = 40000. ;  e2_s(:,:) = 40000.
      DO jj_s=1, ny_s
         DO ji_s=1, nx_s-1
            e1_s(ji_s,jj_s) = distance(Xsrc(ji_s,jj_s),Xsrc(ji_s+1,jj_s),Ysrc(ji_s,jj_s),Ysrc(ji_s+1,jj_s))*1000. ! (m)
         END DO
      END DO
      DO jj_s=1, ny_s-1
         DO ji_s=1, nx_s
            e2_s(ji_s,jj_s) = distance(Xsrc(ji_s,jj_s),Xsrc(ji_s,jj_s+1),Ysrc(ji_s,jj_s),Ysrc(ji_s,jj_s+1))*1000. !(m)
         END DO
      END DO
      IF(nx_s>1) e1_s(nx_s,:) = e1_s(nx_s-1,:)
      IF(ny_s>1) e2_s(:,ny_s) = e2_s(:,ny_s-1)
      !CALL DUMP_FIELD(REAL(e1_s,4), 'e1_s.nc', 'e1')
      !CALL DUMP_FIELD(REAL(e2_s,4), 'e2_s.nc', 'e2')

      !! Min and Max latitude to use for binning :
      PRINT *, ''
      PRINT *, '    *** Latitude binning for more efficient localization *** '
      y_max_bnd = MIN( y_max_t , y_max_s )  !lolo: why MIN() ??? I don't remember why
      y_max_bnd0 = y_max_bnd
      y_min_bnd = MAX( y_min_t , y_min_s )  !lolo: same, MAX()???
      y_min_bnd0 = y_min_bnd
      !PRINT *, 'LOLO y_max_bnd #1 => ', y_max_bnd ; PRINT *, ' y_min_bnd #1 => ', y_min_bnd
      y_max_bnd = MIN( REAL(INT(y_max_bnd+2.),8) ,  90.)
      y_min_bnd = MAX( REAL(INT(y_min_bnd-2.),8) , -90.)
      !PRINT *, 'LOLO y_max_bnd #2 => ', y_max_bnd ; PRINT *, ' y_min_bnd #2 => ', y_min_bnd
      !! Multiple of 0.5:
      y_max_bnd = MIN( NINT(y_max_bnd/0.5)*0.5    ,  90.)
      y_min_bnd = MAX( NINT(y_min_bnd/0.5)*0.5    , -90.)
      !PRINT *, 'LOLO y_max_bnd #3 => ', y_max_bnd ; PRINT *, ' y_min_bnd #3 => ', y_min_bnd

      !lolo: WHY????
      y_max_bnd = MAX( y_max_bnd , y_max_bnd0)
      y_min_bnd = MIN( y_min_bnd , y_min_bnd0)
      !lolo.

      PRINT *, '   => binning from ', y_min_bnd, ' to ', y_max_bnd
      !!
      nsplit = MAX( INT( REAL(Nlat_split) * (y_max_bnd - y_min_bnd)/180. ) , 1)
      !!
      ALLOCATE ( VLAT_SPLIT_BOUNDS(nsplit+1), J_VLAT_S(nsplit,2) )
      dy = (y_max_bnd - y_min_bnd)/REAL(nsplit)
      DO jlat=1,nsplit+1
         VLAT_SPLIT_BOUNDS(jlat) = y_min_bnd + REAL(jlat-1)*dy
      END DO
      !!

      IF( (y_min_bnd > y_min_bnd0).OR.(y_max_bnd < y_max_bnd0) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Bounds for latitude for VLAT_SPLIT_BOUNDS are bad!'
         PRINT *, ' y_min_bnd, y_max_bnd =', y_min_bnd, y_max_bnd !
         PRINT *, ' y_min_bnd0, y_max_bnd0 =', y_min_bnd0, y_max_bnd0
         STOP
      END IF

      J_VLAT_S = 0

      DO jlat = 1, nsplit
         rlat_low = VLAT_SPLIT_BOUNDS(jlat)
         rlat_hgh = VLAT_SPLIT_BOUNDS(jlat+1)
         !! @mbalaro Comment: The two line below can lead to error when working on on small domain ...
         !ij_max_loc = MAXLOC(Ysrc, mask=(Ysrc<=rlat_hgh)) ; jmax_band = ij_max_loc(2)
         !ij_min_loc = MINLOC(Ysrc, mask=(Ysrc>=rlat_low)) ; jmin_band = ij_min_loc(2)
         !! ... it is preferable to look at the min and max value of the ensemble of jj within the range [rlat_low:rlat_hgh]
         !! Largest ever possible j index of the highest latitude in region where Ysrc<=rlat_hgh
         jmax_band = MAXVAL(MAXLOC(Ysrc, mask=(Ysrc<=rlat_hgh), dim=2))  ! MAXLOC(Ysrc, ..., dim=2) returns ni_s values (the max in each column)
         !! Smalles ever possible j index of the smallest latitude in region where Ysrc>=rlat_low
         i1dum = MINLOC(Ysrc, mask=(Ysrc .GT. rlat_low), dim=2)
         jmin_band = MINVAL(i1dum, mask=(i1dum>0))
         DEALLOCATE ( i1dum )
         !!
         !! To be sure to include everything, adding 1 extra points below and above:
         J_VLAT_S(jlat,1) = MAX(jmin_band - 1,   1  )
         J_VLAT_S(jlat,2) = MIN(jmax_band + 1, ny_s)
         !!
         IF(iverbose==2) THEN
            PRINT *, ' Latitude bin #', jlat
            PRINT *, '     => lat_low, lat_high:', REAL(rlat_low,4), REAL(rlat_hgh,4)
            PRINT *, '     => JJ min and max on input domain =>', J_VLAT_S(jlat,1), J_VLAT_S(jlat,2)
         END IF
         !!
      END DO

      PRINT *, ''
      PRINT *, '     => VLAT_SPLIT_BOUNDS ='
      PRINT *, VLAT_SPLIT_BOUNDS
      PRINT *, '     => corresponding start values for j_s ='
      PRINT *, J_VLAT_S(:,1)
      PRINT *, '     => corresponding stop values for j_s ='
      PRINT *, J_VLAT_S(:,2)
      PRINT *, ''

      DO jlat = 1, nsplit
         IF( J_VLAT_S(jlat,2) <= J_VLAT_S(jlat,1) ) THEN
            PRINT *, ' ERROR: jj_stop > jj_start ! ', J_VLAT_S(jlat,2), J_VLAT_S(jlat,1)
            PRINT *, '   => for latitude bin #', jlat
            STOP
         END IF
      END DO

      rlat_old = rmissval
      jj_t_old = -10

      DO jj_t = j_strt_t, j_stop_t, jlat_icr
         DO ji_t = 1, nx_t
            IF( mask_t(ji_t,jj_t) == 1 ) THEN

               lPdbg = ( (iverbose>0).AND.(ji_t==idb).AND.(jj_t==jdb) )
               !lilo
               IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED):'
               IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): WILL DEBUG for target point ji, jj =', idb, jdb

               rlon = Xtrg(ji_t,jj_t)
               rlat = Ytrg(ji_t,jj_t)

               IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): Target lon,lat =', REAL(rlon,4), REAL(rlat,4)

               !! Display progression in stdout:
               IF( (ji_t == nx_t/2).AND.(jj_t /= jj_t_old) ) THEN
                  WRITE(*,'("*** Treated latitude of target domain = ",f9.4," (jj_t = ",i5.5,")")') REAL(rlat,4), jj_t
                  jj_t_old = jj_t
               END IF
               IF( (nx_t == 1).AND.(MOD(jj_t,10)==0) ) &
                  & WRITE(*,'("*** Treated point of target domain = ",i7," (ouf of ",i7,")")') jj_t, ABS(j_stop_t-j_strt_t+1) ! in case of trajectory/ephem stuff

               !! Need to find which jlat of our latitude bins rlat is located in!
               IF( rlat /= rlat_old ) THEN
                  DO jlat=1,nsplit
                     IF(  rlat ==VLAT_SPLIT_BOUNDS(jlat)) EXIT
                     IF( (rlat > VLAT_SPLIT_BOUNDS(jlat)).AND.(rlat <= VLAT_SPLIT_BOUNDS(jlat+1)) ) EXIT
                  END DO
                  !!
               END IF

               lagain    = .TRUE.
               niter     = -1  ! -1 because first pass is for bluff, we want niter=0 for the first use of latitude binning...
               IF( rlat > 60. ) niter = 0 ! we skip the bluff part because the grid might be too close to NP boundary cut!

               DO WHILE ( lagain )

                  IF( niter == -1 ) THEN
                     !! Bluff !
                     !! It's not stupid to assume that the next point to locate is
                     !! pretty near the previously found point (ji_s,jj_s):
                     i1s = MAX(ji_s - 5 , 1)
                     i2s = MIN(ji_s + 5 , nx_s)
                     j1s = MAX(jj_s - 5 , 1)
                     j2s = MIN(jj_s + 5 , ny_s)
                  ELSE
                     i1s = 1
                     i2s = nx_s
                     j1s = J_VLAT_S(MAX(jlat-niter,1)     , 1)  !! MB Comment: Force to 1 if necessary when nearest point not found whereas it exist really
                     j2s = J_VLAT_S(MIN(jlat+niter,nsplit), 2)  !! MB Comment: Force to ny_s if necessary when nearest point not found whereas it exist really
                     IF(iverbose==2) THEN
                        PRINT *, ' *** Treated latitude of target domain =', REAL(rlat,4), ' iter:', niter, jlat
                        PRINT *, '     => bin #', jlat
                        PRINT *, '       => jmin & jmax on source domain =', j1s, j2s
                     END IF
                  END IF

                  Xdist = 1.E12
                  Xdist(i1s:i2s,j1s:j2s) = DISTANCE_2D(rlon, Xsrc(i1s:i2s,j1s:j2s), rlat, Ysrc(i1s:i2s,j1s:j2s))
                  !CALL DUMP_FIELD(REAL(Xdist,4), 'distance_last.nc', 'dist')

                  !! Nearest point is where distance is smallest:
                  ij_min_loc = MINLOC(Xdist(i1s:i2s,j1s:j2s))
                  IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): `ij_min_loc` = ', ij_min_loc
                  ji_s = ij_min_loc(1) + i1s - 1
                  jj_s = ij_min_loc(2) + j1s - 1

                  IF((ji_s==0).OR.(jj_s==0)) THEN
                     PRINT *, ''
                     PRINT *, ' W E I R D  !!!'
                     PRINT *, 'The nearest point was not found!'
                     PRINT *, ' !!! ji_s or jj_s = 0 !!!'
                     PRINT *, ' ** Target point (lon,lat) =>', rlon, rlat
                     STOP
                  END IF

                  emax = MAX( e1_s(ji_s,jj_s),e2_s(ji_s,jj_s) )/1000.*SQRT(2.)
                  IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): `emax` = ', emax, ' km'
                  
                  IF( (Xdist(ji_s,jj_s) <= frac_emax*emax).AND.(ji_s<nx_s) ) THEN
                     !!LOLO: the `.AND.(ji_s<nx_s)` has been the miracle fix :)                     
                     IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): FOUND THE EASY WAY!!!'
                     IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): Distance =', REAL(Xdist(ji_s,jj_s),4), ' km'
                     !! Found !
                     lagain = .FALSE.
                     JIpos(ji_t,jj_t) = ji_s
                     JJpos(ji_t,jj_t) = jj_s
                     !!
                     IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): ji_s, jj_s =', ji_s, jj_s
                     IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): nx_s, ny_s =', nx_s, ny_s
                     !IF(lPdbg) STOP'LOLOX'
                     !!
                  ELSE
                     IF(lPdbg) WRITE(6,*)'#DBG(FIND_NEAREST_TWISTED): GOING FOR REFINED SEARCHING!'
                     !! Not found yet...
                     IF(niter == 0) THEN
                        !! After all the lon,lat couple we are looking for is maybe not part of source domain
                        !! => could do a test and break the iteration if so...
                        mspot_lon = 0 ; mspot_lat = 0
                        WHERE( (Xsrc > MAX(rlon-0.5,  0.)).AND.(Xsrc < MIN(rlon+0.5,360.)) ) mspot_lon = 1
                        WHERE( (Ysrc > MAX(rlat-0.5,-90.)).AND.(Ysrc < MIN(rlat+0.5, 90.)) ) mspot_lat = 1
                        IF( SUM(mspot_lon*mspot_lat) == 0 ) THEN
                           lagain = .FALSE.
                           IF(iverbose==2) PRINT *, ' *** FIND_NEAREST_POINT: SHORT leave test worked! Aborting search!'
                        END IF
                     END IF

                     IF(niter > nsplit/3) THEN
                        !! We are too far in latitude, giving up...
                        IF(iverbose>0) THEN
                           PRINT *, ' *** WARNING: mod_manip.f90/FIND_NEAREST_POINT: Giving up!!!'
                           PRINT *, '     => did not find nearest point for target coordinates:', & !
                              &              REAL(rlon,4), REAL(rlat,4)
                        END IF
                        lagain = .FALSE.
                     END IF

                  END IF    ! IF(Xdist(ji_s,jj_s) <= frac_emax*emax)

                  niter  = niter + 1

               END DO !DO WHILE ( lagain )
               rlat_old = rlat

            END IF ! IF( mask_t(ji_t,jj_t) == 1 )
         END DO    ! DO ji_t = 1, nx_t
      END DO       ! DO jj_t = j_strt_t, j_stop_t, jlat_icr


      WHERE ( JIpos == -1 ) mask_t = -1
      WHERE ( JJpos == -1 ) mask_t = -2

      DEALLOCATE ( VLAT_SPLIT_BOUNDS, J_VLAT_S, e1_s, e2_s, mspot_lon , mspot_lat , Xdist )

   END SUBROUTINE FIND_NEAREST_TWISTED


   FUNCTION DISTANCE(plona, plonb, plata, platb)

      !!----------------------------------------------------------
      !!           ***  FUNCTION  DIST  ***
      !!
      !!  ** Purpose : Compute the distance (km) between
      !!          point A (lona, lata) and B (lonb,latb)
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

      REAL(8), SAVE :: prevlat=-1000., prevlon=-1000

      !! Compute these term only if they differ from previous call
      IF( plata /= prevlat .OR. plona /= prevlon) THEN
         zlatar=plata*rd2rad
         zlonar=plona*rd2rad
         zux=COS(zlonar)*COS(zlatar)
         zuy=SIN(zlonar)*COS(zlatar)
         zuz=SIN(zlatar)
         prevlat=plata
         prevlon=plona
      ENDIF

      zlatbr=platb*rd2rad
      zlonbr=plonb*rd2rad
      zvx=COS(zlonbr)*COS(zlatbr)
      zvy=SIN(zlonbr)*COS(zlatbr)
      zvz=SIN(zlatbr)

      zpds=zux*zvx+zuy*zvy+zuz*zvz

      IF(zpds >= 1.) THEN
         distance = 0.
      ELSE
         distance = rradE*ACOS(zpds)
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
      REAL(8) :: zvx, zvy, zvz

      nx = SIZE(Xlonb,1)
      ny = SIZE(Xlonb,2)

      !! Compute these term only if they differ from previous call
      zlatar=plata*rd2rad
      zlonar=plona*rd2rad
      zux=COS(zlonar)*COS(zlatar)
      zuy=SIN(zlonar)*COS(zlatar)
      zuz=SIN(zlatar)

      distance_2d(:,:) = 0.

      DO jj=1,ny
         DO ji=1,nx

            zlatbr=Xlatb(ji,jj)*rd2rad
            zlonbr=Xlonb(ji,jj)*rd2rad
            zvx=COS(zlonbr)*COS(zlatbr)
            zvy=SIN(zlonbr)*COS(zlatbr)
            zvz=SIN(zlatbr)

            zpds = zux*zvx + zuy*zvy + zuz*zvz

            !IF( zpds < 1.) distance_2d(ji,jj) = rradE*ACOS(zpds)

            distance_2d(ji,jj) = rradE*ACOS(MIN(zpds,1.))

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

            IF( ((rlon_P >= min_x).and.(rlon_P <= max_x)) .AND. ( (max_x-min_x) < dr)  ) THEN
               IF( (rlat_P >= min_y).and.(rlat_P <= max_y) ) THEN
                  jxfnd = ji  ; jyfnd = jj
               END IF
            END IF

         END DO
      END DO

   END SUBROUTINE locate_point


   FUNCTION L_IS_GRID_REGULAR( pX, pY )
      !!-------------------------------------------------------------------------
      !! Tell if a grid (1 longitude 2D array and 1 latitude 2D array) is regular
      !!-------------------------------------------------------------------------
      REAL(8), DIMENSION(:,:), INTENT(in) :: pX, pY
      LOGICAL                             :: L_IS_GRID_REGULAR
      INTEGER :: nx, ny, ji, jj
      !!-------------------------------------------------------------------------
      nx = SIZE(pX,1)
      ny = SIZE(pX,2)
      IF((SIZE(pY,1)/=nx).OR.(SIZE(pY,2)/=ny)) CALL STOP_THIS('`L_IS_GRID_REGULAR@mod_manip` => shape(pY) /= shape(pX)!')
      L_IS_GRID_REGULAR = .TRUE.
      !!  a/ checking on longitude array: (LOLO: use epsilon(pX) instead 1.E-12?)
      DO jj = 2, ny
         IF( SUM( ABS(pX(:,jj) - pX(:,1)) ) > 1.E-9 ) THEN
            L_IS_GRID_REGULAR = .FALSE.
            EXIT
         END IF
      END DO
      !!  b/ now on latitude array:
      IF( L_IS_GRID_REGULAR ) THEN
         DO ji = 2, nx
            IF( SUM( ABS(pY(ji,:) - pY(1,:)) ) > 1.E-9 ) THEN
               L_IS_GRID_REGULAR = .FALSE.
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
      IF( nold /= SIZE(vmask_ignore,1) )     CALL STOP_THIS('`SHRINK_VECTOR` of mod_manip.f90 => data vector and mask vector do not agree in length!')
      IF((new_size > nold).OR.(new_size<=0)) CALL STOP_THIS('`SHRINK_VECTOR` of mod_manip.f90 => your new_size does not make sense!')
      jn = 0
      DO jo = 1, nold
         IF( vmask_ignore(jo) == 1 ) THEN
            jn = jn + 1
            SHRINK_VECTOR(jn) = vect(jo)
         END IF
      END DO
   END FUNCTION SHRINK_VECTOR



   FUNCTION to_degE_scal( rlong )
      !! From any longitude to something between 0 and 360 !
      REAL(8), INTENT(in) :: rlong
      REAL(8)             :: to_degE_scal
      to_degE_scal = MOD( rlong + 360._8 , 360._8 )
   END FUNCTION to_degE_scal
   !!
   FUNCTION to_degE_1d( vlong )
      !! From any longitude to something between 0 and 360 !
      REAL(8), DIMENSION(:), INTENT(in) :: vlong
      REAL(8), DIMENSION(SIZE(vlong,1)) :: to_degE_1d
      to_degE_1d(:) = MOD( vlong(:) + 360._8 , 360._8 )
   END FUNCTION to_degE_1d
   !!
   FUNCTION to_degE_2d( xlong )
      !! From any longitude to something between 0 and 360 !
      REAL(8), DIMENSION(:,:), INTENT(in) :: xlong
      REAL(8), DIMENSION(SIZE(xlong,1),SIZE(xlong,2)) :: to_degE_2d
      to_degE_2d(:,:) = MOD( xlong(:,:) + 360._8 , 360._8 )
   END FUNCTION to_degE_2d


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
      degE_to_degWE_1d(:) = SIGN(1._8,180._8-vlong(:))*MIN(vlong(:), ABS(vlong(:)-360._8))
   END FUNCTION degE_to_degWE_1d

   FUNCTION degE_to_degWE_2d( xlong )
      !! From longitude in 0 -- 360 frame to -180 -- +180 frame...
      REAL(8), DIMENSION(:,:) :: xlong !
      REAL(8), DIMENSION(SIZE(xlong,1),SIZE(xlong,2)) :: degE_to_degWE_2d
      degE_to_degWE_2d(:,:) = SIGN(1._8,180._8-xlong(:,:))*MIN(xlong(:,:), ABS(xlong(:,:)-360._8))
   END FUNCTION degE_to_degWE_2d


   SUBROUTINE EXT_BEYOND_90_REG( kr90, pX, pY, pF, pFe )
      !!============================================================================
      !! We accept only regular lat-lon grid (supposed to include the north
      !! pole)!
      !!
      !! pX, pY is a 2D regular lon-lat grid which is supposed to include the
      !! northpole but for which the highest latitude (at j=Nj) is less than or EQUAL to
      !! 90. We're going to use "across-North-Pole" continuity to fill these 2
      !! extra upper rows!
      !!============================================================================
      INTEGER,                 INTENT(in)  :: kr90 ! 1 -> longitude does not reach 90N, 2 -> it does!
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pX, pY
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pF
      REAL(8), DIMENSION(:,:), INTENT(out) :: pFe
      !!============================================================================
      REAL(8) :: zr, zlon, zlon_m
      INTEGER :: nx, ny, nyp1, nyp2, ji, ji_m
      !!============================================================================
      nx = SIZE(pF,1)
      ny = SIZE(pF,2)

      IF( (kr90<1).OR.(kr90>2) ) CALL STOP_THIS('[EXT_BEYOND_90_REG] => wrong value for kr90!')

      nyp2 = SIZE(pFe,2)
      IF( nyp2 /= ny + 2 ) CALL STOP_THIS('[EXT_BEYOND_90_REG] => target y dim is not nj+2!!!')
      IF( (nx /= SIZE(pY,1)).OR.(ny /= SIZE(pY,2)) ) THEN
         CALL STOP_THIS('[EXT_BEYOND_90_REG] => size of input coor. and data do not match!!!')
      END IF

      nyp1 = ny + 1

      pFe(:,:)    = 0.
      pFe(:,1:ny) = pF(:,:) ! filling center of the domain

      DO ji=1, nx
         zlon = pX(ji,ny) ! ji => zlon
         !! at what ji do we arrive when crossing the northpole => ji_m!
         zlon_m = MOD(zlon+180.,360.)
         ji_m = MINLOC(ABS(pX(:,ny)-zlon_m), dim=1)  ! corresponding i index
         !!                                doesn't reach 90 | reaches 90
         pFe(ji,nyp1) =  pF(ji_m,ny-kr90+1)  !     ny       |   ny-1
         pFe(ji,nyp2) =  pF(ji_m,ny-kr90  )  !    ny-1      |   ny-2
      END DO
   END SUBROUTINE EXT_BEYOND_90_REG


   FUNCTION I_EXT2NP_REG( pX, pY )
      !!
      !! Tells if a grid is regular and stops at the North Pole by returning:
      !!  * 0 => no
      !!  * 1 => yes, doesn't reach longitude = 90N
      !!  * 2 => yes, and reaches longitude = 90N
      !!
      REAL(8), DIMENSION(:,:), INTENT(in) :: pX, pY
      INTEGER                             :: I_EXT2NP_REG
      !!
      INTEGER :: nx, ny
      REAL(8) :: zlat_max, zlat_mm1

      nx = SIZE(pY,1)
      ny = SIZE(pY,2)

      I_EXT2NP_REG = 0

      !! 1st, make sure it's a a map not a trajectory or something else:
      IF( (nx>20).AND.(ny>20) ) THEN
         !! 2nd, it must be regular:
         IF( L_IS_GRID_REGULAR( pX, pY ) ) THEN
            !! 3rd, it must "stop" at the North Pole:
            zlat_max = pY(nx/2,ny)
            zlat_mm1 = pY(nx/2,ny-1)
            IF ( (zlat_max <  90._8).AND.( 2.*zlat_max - zlat_mm1 > 90._8 ) )  I_EXT2NP_REG = 1
            IF (  zlat_max == 90._8 )                                          I_EXT2NP_REG = 2
         END IF
      END IF
   END FUNCTION I_EXT2NP_REG



   SUBROUTINE APPLY_EW_PRDCT( kewp, nbp, pX,  l_is_longitude )
      !!=====================================================================================
      !! Fills the `nbp` first and `nbp` last columns of an array taking
      !! advantage of East-West periodicity with an overlap of `kewp` points
      !!=====================================================================================
      INTEGER,                 INTENT(in)    :: kewp, nbp
      REAL(8), DIMENSION(:,:), INTENT(inout) :: pX
      LOGICAL, OPTIONAL,       INTENT(in)    :: l_is_longitude ! is this a longitude array?
      INTEGER :: nx, ny, kk
      REAL(8) :: zcorr
      !!=====================================================================================
      IF( nbp < 1 ) CALL STOP_THIS( '`APPLY_EW_PRDCT@mod_manip` => `nbp` cannot be < 1 !')
      !!
      zcorr = 0._8
      IF( PRESENT(l_is_longitude) ) THEN
         IF( l_is_longitude ) zcorr = 360._8
      END IF
      !!
      IF( kewp >= 0 ) THEN
         nx = SIZE(pX,1)
         ny = SIZE(pX,2)
         IF( nbp > nx/2 ) CALL STOP_THIS( '`APPLY_EW_PRDCT@mod_manip` => `nbp` too big!')
         !!
         DO kk = 0, nbp-1
            !pX( 1+kk,:) = pX(nx-(kewp-1-kk),:) - zcorr
            !pX(nx-kk,:) = pX( 1+(kewp-1-kk),:) + zcorr
            pX( 1+kk,:) = pX(nx-(kewp+nbp+1-kk),:) - zcorr
            pX(nx-kk,:) = pX( 1+(kewp+nbp+1-kk),:) + zcorr
         END DO
      END IF
   END SUBROUTINE APPLY_EW_PRDCT


END MODULE MOD_MANIP
