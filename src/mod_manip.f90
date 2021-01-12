MODULE MOD_MANIP

   !! Misc. manipulations and operations on 2D arrays...

   !! Author: L. Brodeau

   USE io_ezcdf, ONLY: DUMP_FIELD  ! debug
   USE mod_conf, ONLY: STOP_THIS, iverbose, rd2rad, rradE, rmissval

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

   PUBLIC :: extend_2d_arrays, EXTEND_ARRAY_COOR_4, EXTEND_ARRAY_DATA_4, fill_extra_north_south, extra_2_east, extra_2_west, partial_deriv, &
      &      flip_ud, long_reorg_2d, long_reorg_3d, &
      &      distance, distance_2d, &
      &      find_nearest_point, &
      &      shrink_vector, to_degE, degE_to_degWE, &
      &      ext_north_to_90_regg

   LOGICAL, PARAMETER :: l_force_use_of_twisted = .FALSE.


CONTAINS



   SUBROUTINE EXTEND_2D_ARRAYS(k_ew, pX, pY, pXx, pYx,   pF, pFx, is_orca_grid)

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
      !!
      REAL(8), DIMENSION(:,:), OPTIONAL, INTENT(in)  :: pF
      REAL(8), DIMENSION(:,:), OPTIONAL, INTENT(out) :: pFx
      INTEGER ,                OPTIONAL, INTENT(in)  :: is_orca_grid
      !!
      LOGICAL :: l_extend_F=.FALSE.
      INTEGER :: nx, ny, nxx, nyx
      INTEGER :: ji, jj, iorca
      !!============================================================================

      l_extend_F=.FALSE.
      IF( PRESENT(pF).AND.PRESENT(pFx) ) l_extend_F=.TRUE.
      IF( (PRESENT(pF).AND.(.NOT. PRESENT(pFx))) .OR. (PRESENT(pFx).AND.(.NOT. PRESENT(pF))) ) &
         &  CALL STOP_THIS('[EXTEND_2D_ARRAYS] => you must provide both "pF" and "pFx" as opt. arguments!')

      iorca = 0
      IF( PRESENT(is_orca_grid) ) iorca = is_orca_grid

      nx = SIZE(pX,1)
      ny = SIZE(pX,2)
      nxx = SIZE(pXx,1)
      nyx = SIZE(pXx,2)

      IF( nxx /= nx + 4 ) CALL STOP_THIS('[EXTEND_2D_ARRAYS] => target x dim is not ni+4!!!')
      IF( nyx /= ny + 4 ) CALL STOP_THIS('[EXTEND_2D_ARRAYS] => target y dim is not nj+4!!!')

      IF( (nx /= SIZE(pY,1)).OR.(ny /= SIZE(pY,2)) )     CALL STOP_THIS('[EXTEND_2D_ARRAYS] => size of input coor. do not match!!!')
      IF( (nxx /= SIZE(pYx,1)).OR.(nyx /= SIZE(pYx,2)) ) CALL STOP_THIS('[EXTEND_2D_ARRAYS] => size of output coor. do not match!!!')

      IF( l_extend_F ) THEN
         IF( (nx /= SIZE(pF,1)).OR.(ny /= SIZE(pF,2)) )     CALL STOP_THIS('[EXTEND_2D_ARRAYS] => size of input  coor. and data do not match!!!')
         IF( (nxx /= SIZE(pFx,1)).OR.(nyx /= SIZE(pFx,2)) ) CALL STOP_THIS('[EXTEND_2D_ARRAYS] => size of output coor. and data do not match!!!')
      END IF



      !!   C r e a t i n g   e x t e n d e d   a r r a y s  :
      !!   --------------------------------------------------

      !! Initializing :
      pXx = 0.
      pYx = 0.
      IF( l_extend_F ) pFx = 0.

      !! Filling center of domain:
      pXx(3:nxx-2, 3:nyx-2) = pX(:,:)
      pYx(3:nxx-2, 3:nyx-2) = pY(:,:)
      IF( l_extend_F ) pFx(3:nxx-2, 3:nyx-2) = pF(:,:)



      !! X array :
      !! ---------

      IF(k_ew > -1) THEN   ! we can use east-west periodicity of input file to
         !!                   ! fill extra bands :
         pXx( 1     , 3:nyx-2) = pX(nx - 1 - k_ew , :) - 360.   ! lolo: use or not the 360 stuff???
         pXx( 2     , 3:nyx-2) = pX(nx - k_ew     , :) - 360.
         pXx(nxx   , 3:nyx-2) = pX( 2 + k_ew     , :) + 360.
         pXx(nxx-1 , 3:nyx-2) = pX( 1 + k_ew     , :) + 360.

      ELSE

         !! WEST
         pXx(2, 3:nyx-2) = pX(2,:) - (pX(3,:) - pX(1,:))
         pXx(1, 3:nyx-2) = pX(1,:) - (pX(3,:) - pX(1,:))

         !! EAST
         pXx(nxx-1, 3:nyx-2) = pX(nx-1,:) + pX(nx,:) - pX(nx-2,:)
         pXx(nxx  , 3:nyx-2) = pX(nx,:)   + pX(nx,:) - pX(nx-2,:)

      END IF !IF(k_ew > -1)



      !! ******************
      !! Southern Extension
      !! ******************
      pXx(:, 2) = pXx(:,4) - (pXx(:,5) - pXx(:,3))
      pXx(:, 1) = pXx(:,3) - (pXx(:,5) - pXx(:,3))

      !! ******************
      !! Northern Extension
      !! ******************

      SELECT CASE( iorca )
         !
      CASE (4)
         pXx(2:nxx/2             ,nyx-1) = pXx(nxx:nxx-nxx/2-2:-1,nyx-5)
         pXx(nxx:nxx-nxx/2-2:-1,nyx-1) = pXx(2:nxx/2             ,nyx-5)
         pXx(2:nxx/2             ,nyx)   = pXx(nxx:nxx-nxx/2-2:-1,nyx-6)
         pXx(nxx:nxx-nxx/2-2:-1,nyx)   = pXx(2:nxx/2             ,nyx-6)
      CASE (6)
         pXx(2:nxx/2               ,nyx-1) = pXx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pXx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pXx(2:nxx/2               ,nyx-4)
         pXx(2:nxx/2               ,nyx)   = pXx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pXx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pXx(2:nxx/2               ,nyx-5)
      CASE DEFAULT
         pXx(:,nyx-1) = pXx(:,nyx-3) + pXx(:,nyx-2) - pXx(:,nyx-4)
         pXx(:,nyx)   = pXx(:,nyx-2) + pXx(:,nyx-2) - pXx(:,nyx-4)

      END SELECT


      !! Y array :
      !! ---------

      !! Top side :
      SELECT CASE( iorca )

      CASE (4)
         pYx(2:nxx/2             ,nyx-1) = pYx(nxx:nxx-nxx/2-2:-1,nyx-5)
         pYx(nxx:nxx-nxx/2-2:-1,nyx-1) = pYx(2:nxx/2             ,nyx-5)
         pYx(2:nxx/2             ,nyx)   = pYx(nxx:nxx-nxx/2-2:-1,nyx-6)
         pYx(nxx:nxx-nxx/2-2:-1,nyx)   = pYx(2:nxx/2             ,nyx-6)
      CASE (6)
         pYx(2:nxx/2               ,nyx-1) = pYx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pYx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pYx(2:nxx/2               ,nyx-4)
         pYx(2:nxx/2               ,nyx)   = pYx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pYx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pYx(2:nxx/2               ,nyx-5)
      CASE DEFAULT
         pYx(3:nxx-2, nyx-1) = pY(:, ny-1) + pY(:,ny) - pY(:,ny-2)
         pYx(3:nxx-2, nyx)   = pY(:, ny)   + pY(:,ny) - pY(:,ny-2)
      END SELECT


      !! Bottom side :
      pYx(3:nxx-2, 2) = pY(:,2) - (pY(:,3) - pY(:,1))
      pYx(3:nxx-2, 1) = pY(:,1) - (pY(:,3) - pY(:,1))

      IF(k_ew > -1) THEN

         pYx( 1     , :) = pYx(nx - 1 - k_ew + 2, :)
         pYx( 2     , :) = pYx(nx - k_ew     + 2, :)
         pYx(nxx   , :) = pYx( 2 + k_ew     + 2, :)
         pYx(nxx-1 , :) = pYx( 1 + k_ew     + 2, :)

      ELSE
         !! Left side :
         pYx(2, :) = pYx(4,:) - (pYx(5,:) - pYx(3,:))
         pYx(1, :) = pYx(3,:) - (pYx(5,:) - pYx(3,:))
         !! Right side :
         pYx(nxx-1,:) = pYx(nxx-3,:) + pYx(nxx-2, :) - pYx(nxx-4, :)
         pYx(nxx,:)   = pYx(nxx-2,:) + pYx(nxx-2,:)  - pYx(nxx-4, :)

      END IF


      !! Data array :
      !! ------------
      IF( l_extend_F ) THEN

         IF(k_ew > -1) THEN

            pFx( 1     , 3:nyx-2) = pF(nx - 1 - k_ew , :)
            pFx( 2     , 3:nyx-2) = pF(nx - k_ew     , :)
            pFx(nxx   , 3:nyx-2) = pF( 2 + k_ew     , :)
            pFx(nxx-1 , 3:nyx-2) = pF( 1 + k_ew     , :)

         ELSE

            !! Left side :
            DO jj = 3, nyx-2
               CALL extra_2_east(pXx(nxx-4,jj),pXx(nxx-3,jj),pXx(nxx-2,jj),        &
                  &              pXx(nxx-1,jj),pXx(nxx,jj),                         &
                  &              pFx(nxx-4,jj),pFx(nxx-3,jj),pFx(nxx-2,jj),  &
                  &              pFx(nxx-1,jj),pFx(nxx,jj) )
            END DO

            !! Right side :
            DO jj = 3, nyx-2
               CALL extra_2_west(pXx(5,jj),pXx(4,jj),pXx(3,jj),         &
                  &              pXx(2,jj),pXx(1,jj),                   &
                  &              pFx(5,jj),pFx(4,jj),pFx(3,jj),   &
                  &              pFx(2,jj),pFx(1,jj) )
            END DO

         END IF

         !! Top :
         SELECT CASE( iorca )
         CASE (4)
            IF(iverbose>0) PRINT *, 'ORCA north pole T-point folding type of extrapolation at northern boundary!'
            pFx(2:nxx/2             ,nyx-1) = pFx(nxx:nxx-nxx/2-2:-1,nyx-5)
            pFx(nxx:nxx-nxx/2-2:-1,nyx-1) = pFx(2:nxx/2             ,nyx-5)
            pFx(2:nxx/2             ,nyx)   = pFx(nxx:nxx-nxx/2-2:-1,nyx-6)
            pFx(nxx:nxx-nxx/2-2:-1,nyx)   = pFx(2:nxx/2             ,nyx-6)
         CASE (6)
            IF(iverbose>0) PRINT *, 'ORCA north pole F-point folding type of extrapolation at northern boundary!'
            pFx(2:nxx/2               ,nyx-1) = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
            pFx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pFx(2:nxx/2               ,nyx-4)
            pFx(2:nxx/2               ,nyx)   = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
            pFx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pFx(2:nxx/2               ,nyx-5)

         CASE DEFAULT
            DO ji = 1, nxx
               CALL extra_2_east(pYx(ji,nyx-4),pYx(ji,nyx-3),pYx(ji,nyx-2),       &
                  &              pYx(ji,nyx-1),pYx(ji,nyx),                        &
                  &              pFx(ji,nyx-4),pFx(ji,nyx-3),pFx(ji,nyx-2), &
                  &              pFx(ji,nyx-1),pFx(ji,nyx) )
            END DO
         END SELECT

         !! Bottom side :
         DO ji = 1, nxx
            CALL extra_2_west(pYx(ji,5),pYx(ji,4),pYx(ji,3),       &
               &              pYx(ji,2),pYx(ji,1),                 &
               &              pFx(ji,5),pFx(ji,4),pFx(ji,3), &
               &              pFx(ji,2),pFx(ji,1) )
         END DO

      END IF !IF( l_extend_F )

   END SUBROUTINE EXTEND_2D_ARRAYS







   SUBROUTINE EXTEND_ARRAY_4(k_ew, pX, pY, pF, pFx, ctype,  is_orca_grid)

      !! => needs pX, pY to interpolate the data if it's a data field other than lon or lat...

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
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pF
      REAL(8), DIMENSION(:,:), INTENT(out) :: pFx
      CHARACTER(len=1),        INTENT(in)  :: ctype
      INTEGER ,      OPTIONAL, INTENT(in)  :: is_orca_grid
      !!
      LOGICAL :: l_data=.FALSE., l_lon=.FALSE., l_lat=.FALSE.
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: zXx, zYx
      INTEGER :: nx, ny, nxx, nyx
      INTEGER :: ji, jj, iorca
      !!============================================================================
      iorca = 0
      IF( PRESENT(is_orca_grid) ) iorca = is_orca_grid

      SELECT CASE( ctype )
      CASE ('d')
         l_data = .TRUE.
         l_lon  = .TRUE.
         l_lat  = .TRUE.
      CASE ('x')
         l_lon  = .TRUE.
      CASE ('y')
         l_lat  = .TRUE.
      CASE DEFAULT
         CALL STOP_THIS('[EXTEND_ARRAY_4] => only "d", "x", or "y" allowed for "ctype" !')
      END SELECT

      nx  = SIZE(pF ,1) ; ny  = SIZE(pF ,2)
      nxx = SIZE(pFx,1) ; nyx = SIZE(pFx,2)

      IF( nxx /= nx + 4 )                            CALL STOP_THIS('[EXTEND_ARRAY_4] => target x dim is not ni+4!!!')
      IF( nyx /= ny + 4 )                            CALL STOP_THIS('[EXTEND_ARRAY_4] => target y dim is not nj+4!!!')
      IF( (nx /= SIZE(pX,1)).OR.(ny /= SIZE(pX,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_4] => size of input longitude do not match!!!')
      IF( (nx /= SIZE(pY,1)).OR.(ny /= SIZE(pY,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_4] => size of input latitude  do not match!!!')

      !! Initializing :
      IF( l_lon ) ALLOCATE( zXx(nxx,nyx) )
      IF( l_lat ) ALLOCATE( zYx(nxx,nyx) )
      pFx(:,:) = 0.

      !! Filling center of domain:
      IF( l_lon ) zXx(3:nxx-2, 3:nyx-2) = pX(:,:)
      IF( l_lat ) zYx(3:nxx-2, 3:nyx-2) = pY(:,:)
      pFx(3:nxx-2, 3:nyx-2)             = pF(:,:)


      IF( l_lon ) THEN
         !! X array :
         !! ---------
         IF(k_ew > -1) THEN   ! we can use east-west periodicity of input file to
            !!                   ! fill extra bands :
            zXx( 1     , 3:nyx-2) = pX(nx - 1 - k_ew , :) - 360.   ! lolo: use or not the 360 stuff???
            zXx( 2     , 3:nyx-2) = pX(nx - k_ew     , :) - 360.
            zXx(nxx   , 3:nyx-2) = pX( 2 + k_ew     , :) + 360.
            zXx(nxx-1 , 3:nyx-2) = pX( 1 + k_ew     , :) + 360.
         ELSE
            !! WEST:
            zXx(2, 3:nyx-2) = pX(2,:) - (pX(3,:) - pX(1,:))
            zXx(1, 3:nyx-2) = pX(1,:) - (pX(3,:) - pX(1,:))
            !! EAST:
            zXx(nxx-1, 3:nyx-2) = pX(nx-1,:) + pX(nx,:) - pX(nx-2,:)
            zXx(nxx  , 3:nyx-2) = pX(nx,:)   + pX(nx,:) - pX(nx-2,:)
         END IF !IF(k_ew > -1)
         !!
         !! Southern Extension
         zXx(:, 2) = zXx(:,4) - (zXx(:,5) - zXx(:,3))
         zXx(:, 1) = zXx(:,3) - (zXx(:,5) - zXx(:,3))
         !!
         !! Northern Extension
         SELECT CASE( iorca )
         CASE (4)
            zXx(2:nxx/2             ,nyx-1) = zXx(nxx:nxx-nxx/2-2:-1,nyx-5)
            zXx(nxx:nxx-nxx/2-2:-1,nyx-1) = zXx(2:nxx/2             ,nyx-5)
            zXx(2:nxx/2             ,nyx)   = zXx(nxx:nxx-nxx/2-2:-1,nyx-6)
            zXx(nxx:nxx-nxx/2-2:-1,nyx)   = zXx(2:nxx/2             ,nyx-6)
         CASE (6)
            zXx(2:nxx/2               ,nyx-1) = zXx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
            zXx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = zXx(2:nxx/2               ,nyx-4)
            zXx(2:nxx/2               ,nyx)   = zXx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
            zXx(nxx-1:nxx-nxx/2+1:-1,nyx)   = zXx(2:nxx/2               ,nyx-5)
         CASE DEFAULT
            zXx(:,nyx-1) = zXx(:,nyx-3) + zXx(:,nyx-2) - zXx(:,nyx-4)
            zXx(:,nyx)   = zXx(:,nyx-2) + zXx(:,nyx-2) - zXx(:,nyx-4)
         END SELECT
      END IF !IF( l_lon )


      IF( l_lat ) THEN
         !! Y array :
         !! ---------
         !! Top:
         SELECT CASE( iorca )
         CASE (4)
            zYx(2:nxx/2             ,nyx-1) = zYx(nxx:nxx-nxx/2-2:-1,nyx-5)
            zYx(nxx:nxx-nxx/2-2:-1,nyx-1) = zYx(2:nxx/2             ,nyx-5)
            zYx(2:nxx/2             ,nyx)   = zYx(nxx:nxx-nxx/2-2:-1,nyx-6)
            zYx(nxx:nxx-nxx/2-2:-1,nyx)   = zYx(2:nxx/2             ,nyx-6)
         CASE (6)
            zYx(2:nxx/2               ,nyx-1) = zYx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
            zYx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = zYx(2:nxx/2               ,nyx-4)
            zYx(2:nxx/2               ,nyx)   = zYx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
            zYx(nxx-1:nxx-nxx/2+1:-1,nyx)   = zYx(2:nxx/2               ,nyx-5)
         CASE DEFAULT
            zYx(3:nxx-2, nyx-1) = pY(:, ny-1) + pY(:,ny) - pY(:,ny-2)
            zYx(3:nxx-2, nyx)   = pY(:, ny)   + pY(:,ny) - pY(:,ny-2)
         END SELECT
         !!
         !! Bottom:
         zYx(3:nxx-2, 2) = pY(:,2) - (pY(:,3) - pY(:,1))
         zYx(3:nxx-2, 1) = pY(:,1) - (pY(:,3) - pY(:,1))
         !!
         !! East-West:
         IF(k_ew > -1) THEN
            zYx( 1     , :) = zYx(nx - 1 - k_ew + 2, :)
            zYx( 2     , :) = zYx(nx - k_ew     + 2, :)
            zYx(nxx   , :) = zYx( 2 + k_ew     + 2, :)
            zYx(nxx-1 , :) = zYx( 1 + k_ew     + 2, :)
         ELSE
            !! West:
            zYx(2, :) = zYx(4,:) - (zYx(5,:) - zYx(3,:))
            zYx(1, :) = zYx(3,:) - (zYx(5,:) - zYx(3,:))
            !! East:
            zYx(nxx-1,:) = zYx(nxx-3,:) + zYx(nxx-2, :) - zYx(nxx-4, :)
            zYx(nxx,:)   = zYx(nxx-2,:) + zYx(nxx-2,:)  - zYx(nxx-4, :)
         END IF
         !!
      END IF !IF( l_lat )


      !! Data array :
      !! ------------
      IF( l_data ) THEN

         IF(k_ew > -1) THEN
            pFx( 1     , 3:nyx-2) = pF(nx - 1 - k_ew , :)
            pFx( 2     , 3:nyx-2) = pF(nx - k_ew     , :)
            pFx(nxx   , 3:nyx-2) = pF( 2 + k_ew     , :)
            pFx(nxx-1 , 3:nyx-2) = pF( 1 + k_ew     , :)
         ELSE
            !! East:
            DO jj = 3, nyx-2
               CALL extra_2_east(zXx(nxx-4,jj),zXx(nxx-3,jj),zXx(nxx-2,jj),        &
                  &              zXx(nxx-1,jj),zXx(nxx,jj),                         &
                  &              pFx(nxx-4,jj),pFx(nxx-3,jj),pFx(nxx-2,jj),  &
                  &              pFx(nxx-1,jj),pFx(nxx,jj) )
            END DO
            !! West:
            DO jj = 3, nyx-2
               CALL extra_2_west(zXx(5,jj),zXx(4,jj),zXx(3,jj),         &
                  &              zXx(2,jj),zXx(1,jj),                   &
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
         CASE (6)
            IF(iverbose>0) PRINT *, 'ORCA north pole F-point folding type of extrapolation at northern boundary!'
            pFx(2:nxx/2               ,nyx-1) = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
            pFx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pFx(2:nxx/2               ,nyx-4)
            pFx(2:nxx/2               ,nyx)   = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
            pFx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pFx(2:nxx/2               ,nyx-5)
         CASE DEFAULT
            DO ji = 1, nxx
               CALL extra_2_east(zYx(ji,nyx-4),zYx(ji,nyx-3),zYx(ji,nyx-2),       &
                  &              zYx(ji,nyx-1),zYx(ji,nyx),                        &
                  &              pFx(ji,nyx-4),pFx(ji,nyx-3),pFx(ji,nyx-2), &
                  &              pFx(ji,nyx-1),pFx(ji,nyx) )
            END DO
         END SELECT
         !!
         !! Bottom:
         DO ji = 1, nxx
            CALL extra_2_west(zYx(ji,5),zYx(ji,4),zYx(ji,3),       &
               &              zYx(ji,2),zYx(ji,1),                 &
               &              pFx(ji,5),pFx(ji,4),pFx(ji,3), &
               &              pFx(ji,2),pFx(ji,1) )
         END DO
         !!
      ELSE

         IF( l_lon ) pFx(:,:) = zXx(:,:)
         IF( l_lat ) pFx(:,:) = zYx(:,:)

      END IF !IF( l_data )


      IF( l_lon ) DEALLOCATE( zXx )
      IF( l_lat ) DEALLOCATE( zYx )

   END SUBROUTINE EXTEND_ARRAY_4




   
   SUBROUTINE EXTEND_ARRAY_COOR_4(k_ew, pX, pY, pXx, pYx,  is_orca_grid)
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
      INTEGER :: iorca
      !!============================================================================
      iorca = 0
      IF( PRESENT(is_orca_grid) ) iorca = is_orca_grid

      nx  = SIZE(pX ,1) ; ny  = SIZE(pX ,2)
      IF( (nx  /= SIZE(pY ,1)).OR.(ny  /= SIZE(pY ,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_COOR_4] => size of input longitude do not match!!!')

      nxx = SIZE(pXx,1) ; nyx = SIZE(pXx,2)
      IF( nxx /= nx + 4 )                                CALL STOP_THIS('[EXTEND_ARRAY_COOR_4] => target x dim is not ni+4!!!')
      IF( nyx /= ny + 4 )                                CALL STOP_THIS('[EXTEND_ARRAY_COOR_4] => target y dim is not nj+4!!!')
      IF( (nxx /= SIZE(pYx,1)).OR.(nyx /= SIZE(pYx,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_COOR_4] => size of input latitude does not match longitude !!!')

      !! Initializing :
      pXx(:,:) = 0.
      pYx(:,:) = 0.

      !! Filling center of domain:
      pXx(3:nxx-2, 3:nyx-2) = pX(:,:)
      pYx(3:nxx-2, 3:nyx-2) = pY(:,:)

      !! X array :
      !! ---------
      IF(k_ew > -1) THEN   ! we can use east-west periodicity of input file to
         !!                   ! fill extra bands :
         pXx( 1     , 3:nyx-2) = pX(nx - 1 - k_ew , :) - 360.   ! lolo: use or not the 360 stuff???
         pXx( 2     , 3:nyx-2) = pX(nx - k_ew     , :) - 360.
         pXx(nxx   , 3:nyx-2) = pX( 2 + k_ew     , :) + 360.
         pXx(nxx-1 , 3:nyx-2) = pX( 1 + k_ew     , :) + 360.
      ELSE
         !! WEST:
         pXx(2, 3:nyx-2) = pX(2,:) - (pX(3,:) - pX(1,:))
         pXx(1, 3:nyx-2) = pX(1,:) - (pX(3,:) - pX(1,:))
         !! EAST:
         pXx(nxx-1, 3:nyx-2) = pX(nx-1,:) + pX(nx,:) - pX(nx-2,:)
         pXx(nxx  , 3:nyx-2) = pX(nx,:)   + pX(nx,:) - pX(nx-2,:)
      END IF !IF(k_ew > -1)
      !!
      !! Southern Extension
      pXx(:, 2) = pXx(:,4) - (pXx(:,5) - pXx(:,3))
      pXx(:, 1) = pXx(:,3) - (pXx(:,5) - pXx(:,3))
      !!
      !! Northern Extension
      SELECT CASE( iorca )
      CASE (4)
         pXx(2:nxx/2             ,nyx-1) = pXx(nxx:nxx-nxx/2-2:-1,nyx-5)
         pXx(nxx:nxx-nxx/2-2:-1,nyx-1) = pXx(2:nxx/2             ,nyx-5)
         pXx(2:nxx/2             ,nyx)   = pXx(nxx:nxx-nxx/2-2:-1,nyx-6)
         pXx(nxx:nxx-nxx/2-2:-1,nyx)   = pXx(2:nxx/2             ,nyx-6)
      CASE (6)
         pXx(2:nxx/2               ,nyx-1) = pXx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pXx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pXx(2:nxx/2               ,nyx-4)
         pXx(2:nxx/2               ,nyx)   = pXx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pXx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pXx(2:nxx/2               ,nyx-5)
      CASE DEFAULT
         pXx(:,nyx-1) = pXx(:,nyx-3) + pXx(:,nyx-2) - pXx(:,nyx-4)
         pXx(:,nyx)   = pXx(:,nyx-2) + pXx(:,nyx-2) - pXx(:,nyx-4)
      END SELECT


      !! Y array :
      !! ---------
      !! Top:
      SELECT CASE( iorca )
      CASE (4)
         pYx(2:nxx/2             ,nyx-1) = pYx(nxx:nxx-nxx/2-2:-1,nyx-5)
         pYx(nxx:nxx-nxx/2-2:-1,nyx-1) = pYx(2:nxx/2             ,nyx-5)
         pYx(2:nxx/2             ,nyx)   = pYx(nxx:nxx-nxx/2-2:-1,nyx-6)
         pYx(nxx:nxx-nxx/2-2:-1,nyx)   = pYx(2:nxx/2             ,nyx-6)
      CASE (6)
         pYx(2:nxx/2               ,nyx-1) = pYx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pYx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pYx(2:nxx/2               ,nyx-4)
         pYx(2:nxx/2               ,nyx)   = pYx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pYx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pYx(2:nxx/2               ,nyx-5)
      CASE DEFAULT
         pYx(3:nxx-2, nyx-1) = pY(:, ny-1) + pY(:,ny) - pY(:,ny-2)
         pYx(3:nxx-2, nyx)   = pY(:, ny)   + pY(:,ny) - pY(:,ny-2)
      END SELECT
      !!
      !! Bottom:
      pYx(3:nxx-2, 2) = pY(:,2) - (pY(:,3) - pY(:,1))
      pYx(3:nxx-2, 1) = pY(:,1) - (pY(:,3) - pY(:,1))
      !!
      !! East-West:
      IF(k_ew > -1) THEN
         pYx( 1     , :) = pYx(nx - 1 - k_ew + 2, :)
         pYx( 2     , :) = pYx(nx - k_ew     + 2, :)
         pYx(nxx   , :) = pYx( 2 + k_ew     + 2, :)
         pYx(nxx-1 , :) = pYx( 1 + k_ew     + 2, :)
      ELSE
         !! West:
         pYx(2, :) = pYx(4,:) - (pYx(5,:) - pYx(3,:))
         pYx(1, :) = pYx(3,:) - (pYx(5,:) - pYx(3,:))
         !! East:
         pYx(nxx-1,:) = pYx(nxx-3,:) + pYx(nxx-2, :) - pYx(nxx-4, :)
         pYx(nxx,:)   = pYx(nxx-2,:) + pYx(nxx-2,:)  - pYx(nxx-4, :)
      END IF
      !!
   END SUBROUTINE EXTEND_ARRAY_COOR_4
   

   SUBROUTINE EXTEND_ARRAY_DATA_4(k_ew, pXx, pYx, pF, pFx,  is_orca_grid)
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
      INTEGER ,      OPTIONAL, INTENT(in)  :: is_orca_grid
      !!
      INTEGER :: nx, ny, nxx, nyx
      INTEGER :: ji, jj, iorca
      !!============================================================================
      iorca = 0
      IF( PRESENT(is_orca_grid) ) iorca = is_orca_grid

      nx  = SIZE(pF ,1) ; ny  = SIZE(pF ,2)
      nxx = SIZE(pFx,1) ; nyx = SIZE(pFx,2)

      IF( nxx /= nx + 4 )                            CALL STOP_THIS('[EXTEND_ARRAY_DATA_4] => target x dim is not ni+4!!!')
      IF( nyx /= ny + 4 )                            CALL STOP_THIS('[EXTEND_ARRAY_DATA_4] => target y dim is not nj+4!!!')
      IF( (nxx /= SIZE(pXx,1)).OR.(nyx /= SIZE(pXx,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_DATA_4] => size of input longitude do not match!!!')
      IF( (nxx /= SIZE(pYx,1)).OR.(nyx /= SIZE(pYx,2)) ) CALL STOP_THIS('[EXTEND_ARRAY_DATA_4] => size of input latitude  do not match!!!')

      !! Initializing :
      pFx(:,:) = 0.

      !! Filling center of domain:
      pFx(3:nxx-2, 3:nyx-2)             = pF(:,:)

      !! Data array :
      !! ------------
      IF(k_ew > -1) THEN
         pFx( 1     , 3:nyx-2) = pF(nx - 1 - k_ew , :)
         pFx( 2     , 3:nyx-2) = pF(nx - k_ew     , :)
         pFx(nxx   , 3:nyx-2) = pF( 2 + k_ew     , :)
         pFx(nxx-1 , 3:nyx-2) = pF( 1 + k_ew     , :)
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
      CASE (6)
         IF(iverbose>0) PRINT *, 'ORCA north pole F-point folding type of extrapolation at northern boundary!'
         pFx(2:nxx/2               ,nyx-1) = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-4)
         pFx(nxx-1:nxx-nxx/2+1:-1,nyx-1) = pFx(2:nxx/2               ,nyx-4)
         pFx(2:nxx/2               ,nyx)   = pFx(nxx-1:nxx-nxx/2+1:-1,nyx-5)
         pFx(nxx-1:nxx-nxx/2+1:-1,nyx)   = pFx(2:nxx/2               ,nyx-5)
      CASE DEFAULT
         DO ji = 1, nxx
            CALL extra_2_east(pYx(ji,nyx-4),pYx(ji,nyx-3),pYx(ji,nyx-2),       &
               &              pYx(ji,nyx-1),pYx(ji,nyx),                        &
               &              pFx(ji,nyx-4),pFx(ji,nyx-3),pFx(ji,nyx-2), &
               &              pFx(ji,nyx-1),pFx(ji,nyx) )
         END DO
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
   END SUBROUTINE EXTEND_ARRAY_DATA_4












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

      CALL EXTEND_2D_ARRAYS(k_ew, pX, pY, ZX, ZY,  pF=pF, pFx=ZF)


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






   SUBROUTINE FIND_NEAREST_POINT(Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp,  ithread, pmsk_dom_trg)
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
      INTEGER,    OPTIONAL,       INTENT(in)  :: ithread
      INTEGER(1), OPTIONAL ,DIMENSION(:,:), INTENT(inout) :: pmsk_dom_trg

      INTEGER :: ithrd, jj, nx_s, ny_s, nx_t, ny_t, j_strt_t, j_stop_t, jlat_icr
      REAL(8) :: y_max_s, y_min_s, rmin_dlat_dj, rtmp
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: mask_ignore_t
      INTEGER,    DIMENSION(:),   ALLOCATABLE :: i1dum
      REAL(8),    DIMENSION(:,:), ALLOCATABLE :: ztmp_t
      LOGICAL :: l_is_reg_s, l_is_reg_t

      ithrd = 0 ! no OpenMP !
      IF( PRESENT(ithread) ) ithrd = ithread

      nx_s  = SIZE(Xsrc,1)
      ny_s  = SIZE(Xsrc,2)
      nx_t = SIZE(Xtrg,1)
      ny_t = SIZE(Xtrg,2)

      PRINT *, ' Source domain size: ', nx_s, ny_s
      PRINT *, ' Target domain size: ', nx_t, ny_t

      IF( (SIZE(Ysrc,1) /= nx_s) .OR. (SIZE(Ysrc,2) /= ny_s) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Ysrc dont agree in shape with Xsrc'
         STOP
      END IF
      IF( (SIZE(Ytrg,1) /= nx_t) .OR. (SIZE(Ytrg,2) /= ny_t) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): Ytrg dont agree in shape with Xtrg'
         STOP
      END IF
      IF( (SIZE(JIp,1) /= nx_t) .OR. (SIZE(JIp,2) /= ny_t) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): JIp dont agree in shape with Xtrg'
         PRINT *, SIZE(JIp,1), SIZE(JIp,2), 'vs', nx_t, ny_t
         STOP
      END IF
      IF( (SIZE(JJp,1) /= nx_t) .OR. (SIZE(JJp,2) /= ny_t) ) THEN
         PRINT *, ' ERROR (FIND_NEAREST_POINT of mod_manip.f90): JJp dont agree in shape with Xtrg'
         STOP
      END IF

      PRINT *, ''
      !! Checking if source domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_s = L_IS_GRID_REGULAR( Xsrc , Ysrc )
      PRINT *, ' *** FIND_NEAREST_POINT => Is source grid regular ??? =>', l_is_reg_s

      !! Checking if target domain is regular or not.  => will allow later to
      !!  decide the level of complexity of the algorithm that find nearest
      !!  points...
      l_is_reg_t = L_IS_GRID_REGULAR( Xtrg , Ytrg )
      PRINT *, ' *** FIND_NEAREST_POINT => Is target grid regular ??? =>', l_is_reg_t

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
            &                       j_strt_t, j_stop_t, jlat_icr, ithrd )
      ELSE
         PRINT *, '                  => going for advanced FIND_NEAREST algorithm !'; PRINT *, ''
         CALL FIND_NEAREST_TWISTED( Xtrg, Ytrg, Xsrc, Ysrc, JIp, JJp, mask_ignore_t, &
            &                       j_strt_t, j_stop_t, jlat_icr, ithrd )
      END IF
      PRINT *, ''

      IF( PRESENT( pmsk_dom_trg ) ) pmsk_dom_trg(:,:) = mask_ignore_t(:,:)
      DEALLOCATE ( mask_ignore_t )

   END SUBROUTINE FIND_NEAREST_POINT







   SUBROUTINE FIND_NEAREST_EASY( Xtrg, Ytrg, Xsrc, Ysrc, JIpos, JJpos, mask_t, &
      &                          j_strt_t, j_stop_t, jlat_icr, ithrd )
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
      INTEGER,                    INTENT(in)  :: ithrd, j_strt_t, j_stop_t, jlat_icr

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
                  &  WRITE(6,'(" *** Treated j-point of target domain = ",i4.4," (thread # ",i2.2,")")') jj_t, ithrd

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
      &                             j_strt_t, j_stop_t, jlat_icr, ithrd )
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
      INTEGER, INTENT(in)                       :: j_strt_t, j_stop_t, jlat_icr, ithrd
      !!
      !! Important parameters:
      INTEGER, PARAMETER :: Nlat_split = 40   ! number of latitude bands to split the search work (for 180. degree south->north)
      REAL(8), PARAMETER :: frac_emax = 0.51   ! fraction of emax to test if found!
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

      LOGICAL :: lagain


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

      IF( (y_min_t>=y_max_s).OR.(y_max_t<=y_min_s) ) THEN
         WRITE(6,*)' ERROR (FIND_NEAREST_TWISTED of mod_manip.f90): Target and source latitudes do not overlap!'
         STOP
      END IF
      !! ---------------------------------------------------------------------------------------
      PRINT *, '       => going for advanced algorithm ! (thread #', ithrd, ')'

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
      !IF(iverbose==2) THEN
      !   CALL DUMP_FIELD(REAL(e1_s,4), 'e1_s.nc', 'e1')
      !   CALL DUMP_FIELD(REAL(e2_s,4), 'e2_s.nc', 'e2')
      !END IF

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

               rlon = Xtrg(ji_t,jj_t)
               rlat = Ytrg(ji_t,jj_t)

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
                  !PRINT *, 'LOLO STOP! mod_manip.f90 !!!'; STOP

                  !! Nearest point is where distance is smallest:
                  ij_min_loc = MINLOC(Xdist(i1s:i2s,j1s:j2s))
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

                  emax = MAX(e1_s(ji_s,jj_s),e2_s(ji_s,jj_s))/1000.*SQRT(2.)

                  IF( Xdist(ji_s,jj_s) <= frac_emax*emax) THEN
                     !! Found !
                     lagain = .FALSE.
                     JIpos(ji_t,jj_t) = ji_s
                     JJpos(ji_t,jj_t) = jj_s
                     IF(iverbose==2) THEN
                        IF( niter == -1 ) THEN
                           PRINT *, '    --- F O U N D  with bluff !!! ---'
                        ELSE
                           PRINT *, '    --- F O U N D --- niter =', niter
                        END IF
                     END IF
                     !! Found .
                  ELSE
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
      IF( (SIZE(Xlat,1) /= nx) .OR. (SIZE(Xlat,2) /= ny) ) THEN
         PRINT *, ' ERROR (L_IS_GRID_REGULAR of mod_grids.f90): Xlat does not agree in shape with Xlon!'
         STOP
      END IF
      l_is_grid_regular = .TRUE.
      !!  a/ checking on longitude array: (LOLO: use epsilon(Xlon) instead 1.E-12?)
      DO jj = 2, ny
         IF( SUM( ABS(Xlon(:,jj) - Xlon(:,1)) ) > 1.E-12 ) THEN
            l_is_grid_regular = .FALSE.
            EXIT
         END IF
      END DO
      !!  b/ now on latitude array:
      IF( l_is_grid_regular ) THEN
         DO ji = 1, nx
            IF( SUM( ABS(Xlat(ji,:) - Xlat(1,:)) ) > 1.E-12 ) THEN
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
      IF( nold /= SIZE(vmask_ignore,1) ) THEN
         PRINT *, ' ERROR (SHRINK_VECTOR of mod_manip.f90): data vector and mask vector do not agree in length!'
         STOP
      END IF
      IF( (new_size > nold).OR.(new_size<=0) ) THEN
         PRINT *, ' ERROR (SHRINK_VECTOR of mod_manip.f90): your new_size does not make sense!'
         STOP
      END IF
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
      REAL(8), DIMENSION(:,:) :: xlong
      REAL(8), DIMENSION(SIZE(xlong,1),SIZE(xlong,2)) :: degE_to_degWE_2d
      degE_to_degWE_2d(:,:) = SIGN(1._8,180._8-xlong(:,:))*MIN(xlong(:,:), ABS(xlong(:,:)-360._8))
   END FUNCTION degE_to_degWE_2d




   SUBROUTINE EXT_NORTH_TO_90_REGG( pX, YY, pF,  XP, YP, FP )
      !!============================================================================
      !! We accept only regular lat-lon grid (supposed to include the north
      !! pole)!
      !!
      !! pX, YY is a 2D regular lon-lat grid which is supposed to include the
      !! northpole but for which the highest latitude (at j=Nj) is less than
      !! 90. We're going to use "across-North-Pole" continuity to fill the
      !! last upper row where latitude = 90 !
      !!
      !! => XP, YP is the same grid with an extra upper row where latitude = 90
      !!     => Last upper row of FP (on XP,YP) contains interpolated values of
      !!        the field
      !!============================================================================
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pX, YY
      REAL(4), DIMENSION(:,:), INTENT(in)  :: pF
      REAL(8), DIMENSION(:,:), INTENT(out) :: XP, YP
      REAL(4), DIMENSION(:,:), INTENT(out) :: FP

      !! Local
      REAL(8) :: rr, zlon, zlon_m
      INTEGER :: nx, ny, nyp1, nyp2
      INTEGER :: ji, ji_m

      nx = SIZE(pX,1)
      ny = SIZE(pX,2)

      nyp2 = SIZE(XP,2)
      IF( nyp2 /= ny + 2 ) CALL STOP_THIS('[EXTEND_2D_ARRAYS] => EXT_NORTH_TO_90_REGG : target y dim is not nj+2!!!')

      IF( (nx /= SIZE(YY,1)).OR.(ny /= SIZE(YY,2)).OR. &
         & (nx /= SIZE(pF,1)).OR.(ny /= SIZE(pF,2))) THEN
         CALL STOP_THIS('[EXTEND_2D_ARRAYS] => EXT_NORTH_TO_90_REGG : size of input coor. and data do not match!!!')
      END IF
      IF( (SIZE(XP,1) /= SIZE(YP,1)).OR.(SIZE(XP,2) /= SIZE(YP,2)).OR. &
         & (SIZE(XP,1) /= SIZE(FP,1)).OR.(SIZE(XP,2) /= SIZE(FP,2))) THEN
         CALL STOP_THIS('[EXTEND_2D_ARRAYS] => EXT_NORTH_TO_90_REGG : size of output coor. and data do not match!!!')
      END IF

      nyp1 = ny + 1

      XP = 0.
      YP = 0.
      FP = 0.

      !! Filling center of domain:
      XP(:, 1:ny) = pX(:,:)
      YP(:, 1:ny) = YY(:,:)
      FP(:, 1:ny) = pF(:,:)

      !! Testing if the grid is of the type of what we expect:
      rr = YY(nx/2,ny) ! highest latitude
      IF( rr == 90.0 ) CALL STOP_THIS('[EXTEND_2D_ARRAYS] => EXT_NORTH_TO_90_REGG : mhh well you shouldnt be here I guess, 90 exists!...')
      IF( SUM( (YY(:,ny) - rr)**2 ) > 1.E-12 ) CALL STOP_THIS('[EXTEND_2D_ARRAYS] => EXT_NORTH_TO_90_REGG : mhh well you shouldnt be here I guess, grid doesnt seem to be regular!...')

      !! Longitude points for the extra upper row are just the same!
      XP(:,nyp1) = pX(:,ny)
      XP(:,nyp2) = pX(:,ny)

      !! For latitude it's easy:
      YP(:,nyp1) = 90.0
      YP(:,nyp2) = 90.0 + (YY(nx/2,ny) - YY(nx/2,ny-1)) ! 90. + dlon ! lolo bad???

      DO ji=1, nx
         zlon = pX(ji,ny) ! ji => zlon
         !! at what ji do we arrive when crossing the northpole => ji_m!
         zlon_m = MOD(zlon+180.,360.)
         !PRINT *, ' zlon, zlon_m =>', zlon, zlon_m
         ji_m = MINLOC(ABS(pX(:,ny)-zlon_m), dim=1)
         !PRINT *, '  ji_m =', ji_m
         !PRINT *, 'pX(ji_m,ny) =', pX(ji_m,ny)
         !! Well so the northpole is righ in between so:
         FP(ji,nyp1) = 0.5*(pF(ji,ny) + pF(ji_m,ny)) ! lolo fix!
         FP(ji,nyp2) =  pF(ji_m,ny-1) ! lolo bad???
      END DO

      !PRINT *, 'LOLO EXT_NORTH_TO_90_REGG: YP =', YP(nx/2,:)
      !PRINT *, ''
      !PRINT *, 'LOLO EXT_NORTH_TO_90_REGG: FP =', FP(nx/2,:)
      !STOP 'boo'

   END SUBROUTINE EXT_NORTH_TO_90_REGG





END MODULE MOD_MANIP
