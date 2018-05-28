MODULE MOD_AKIMA_1D
   !!
   !! Author:
   !! =======
   !!           Coded from scratch by Laurent BRODEAU, 2008
   !!
   !!
   !!
   USE mod_manip, ONLY: extra_2_west, extra_2_east
   !!
   IMPLICIT NONE

   
   INTERFACE AKIMA_1D
      MODULE PROCEDURE AKIMA_1D_1D    !, AKIMA_1D_3D
   END INTERFACE AKIMA_1D

   PRIVATE


   PUBLIC :: AKIMA_1D


CONTAINS



   SUBROUTINE AKIMA_1D_1D(vx1, vy1, vx2, vy2)
      !!
      REAL(4), DIMENSION(:), INTENT(in)  :: vx1, vy1, vx2
      REAL(4), DIMENSION(:), INTENT(out) :: vy2
      !!___________________________________________________
      !!
      !! Local :
      INTEGER :: n1, n2
      !!
      REAL(8), DIMENSION(:), ALLOCATABLE  :: vx1p4, vy1p4, vslp
      !!
      REAL(8) :: a0, a1, a2, a3, rr, dx
      INTEGER  :: n4, j1, j2, p1, p2, j1p4
      LOGICAL  :: lfnd
      !!
      !!
      IF ( size(vx1) /= size(vy1) ) THEN
         PRINT *, 'ERROR, mod_akima_1d.f90 => AKIMA_1D_1D: vx1 and vy1 do not have the same size!'; STOP
      END IF
      IF ( size(vx2) /= size(vy2) ) THEN
         PRINT *, 'ERROR, mod_akima_1d.f90 => AKIMA_1D_1D: vx2 and vy2 do not have the same size!'; STOP
      END IF
      !!
      n1 = size(vx1)
      n2 = size(vx2)
      !!
      n4 = n1 + 4
      !!
      ALLOCATE ( vx1p4(n4), vy1p4(n4), vslp(n4) )
      !!
      vx1p4(3:n1+2) = REAL(vx1(:), 8)
      vy1p4(3:n1+2) = REAL(vy1(:), 8)
      !!
      !! Extending input X array :
      !! =========================
      !! - if vx is regularly spaced, it's not a big deal, otherwise we use
      !!   what's been proposed by Akima (1970) :
      !!
      !! Bottom (or West) :
      vx1p4(2) =  REAL(vx1p4(4) - (vx1p4(5) - vx1p4(3)), 8)
      vx1p4(1) =  REAL(vx1p4(3) - (vx1p4(5) - vx1p4(3)), 8)
      !!
      !! Top (or East) :
      vx1p4(n4-1) =  REAL(vx1p4(n4-3) + vx1p4(n4-2) - vx1p4(n4-4), 8)
      vx1p4(n4)   =  REAL(vx1p4(n4-2) + vx1p4(n4-2) - vx1p4(n4-4), 8)
      !!
      !!
      !! Now extrapolating input Y values on these 4 extra points :
      !! ==========================================================
      !!
      !! Bottom (or West) :
      CALL extra_2_west(vx1p4(5), vx1p4(4), vx1p4(3), vx1p4(2), vx1p4(1), &
         &            vy1p4(5), vy1p4(4), vy1p4(3), vy1p4(2), vy1p4(1) )
      !!
      !! Top (or East) :
      CALL extra_2_east(vx1p4(n4-4), vx1p4(n4-3), vx1p4(n4-2), vx1p4(n4-1), vx1p4(n4), &
         &            vy1p4(n4-4), vy1p4(n4-3), vy1p4(n4-2), vy1p4(n4-1), vy1p4(n4) )
      !!
      !!
      !! Computing slopes :
      !! ==================
      !!
      CALL SLOPES_1D(vx1p4, vy1p4, vslp)
      !!
      !!
      !!
      !! Ok! Now, each point of the target grid must be interpolated :
      !! =============================================================
      !!
      !! We INTERPOLATE and don't extrapolate so checking the bounds !
      !! input X array is supposed to be orgnised so :
      !    xstart = vx1(1)  !  xstart = minval(x1) ;   xstop  = maxval(x1)
      !    xstop  = vx1(n1)
      !!
      !!
      !! Treating each target point :
      !! ============================
      !!
      !! !LB : we assume so far that x2 is totally in-organised !!!
      !!        -> should be fixed !!!
      !!        x1 is totally ORGANISED (increasing!)
      !!
      !!
      DO j2 = 1, n2
         !!
         j1   = 1
         lfnd = .FALSE.
         !!
         !! Persistence: if point of output grid is shallower
         !!              than first point of input grid
         !!              (should occur only when j2 = 1)
         IF ( vx2(j2) < vx1(1) ) THEN
            vy2(j2) = vy1(1)
            GOTO 10
         END IF
         !!
         !!
         DO WHILE ( (j1 < n1).and.(.not. lfnd) )
            !!
            !!
            IF ( (vx2(j2) > vx1(j1)).and.(vx2(j2) < vx1(j1+1)) ) THEN
               j1p4  = j1 + 2
               p1    = j1p4
               p2    = j1p4 + 1
               lfnd = .TRUE.
               !!
            ELSE
               !!
               IF ( vx2(j2) == vx1(j1) ) THEN
                  vy2(j2) = vy1(j1)
                  GOTO 10
               ELSE
                  IF ( vx2(j2) == vx1(j1+1) ) THEN
                     vy2(j2) = vy1(j1+1)
                     GOTO 10
                  END IF
               END IF
               !!
            END IF
            !!
            j1 = j1 + 1
            !!
         END DO
         !!  (PM) Take bottom value below the last data
         IF (( j1 == n1) .and. (.not. lfnd)) THEN
            vy2(j2) = vy1(n1)
            GOTO 10
         END IF
         !!
         IF ( .not. lfnd ) THEN
            vy2(j2) = -6666.0
            GOTO 10
         END IF
         !!
         !!
         !!
         !! Estimating vy2(j2) with a third order polynome :
         !! -----------------------------------------------
         !! Coefficients of the polynome
         !!
         !! MIND ! : p1 and p2 are given in 'vp4'
         a0 = vy1p4(p1)
         a1 = vslp(p1)
         a2 = (3*(vy1p4(p2) - vy1p4(p1))/(vx1p4(p2) - vx1p4(p1)) - 2*vslp(p1) - vslp(p2)) &
            & /(vx1p4(p2) - vx1p4(p1))
         a3 = (vslp(p1) + vslp(p2) - 2*(vy1p4(p2)-vy1p4(p1))/(vx1p4(p2)-vx1p4(p1)))  &
            &         /(vx1p4(p2) - vx1p4(p1))**2
         !!
         dx = REAL(vx2(j2), 8) - vx1p4(p1)
         !!
         rr = a0 + a1*dx + a2*dx**2 + a3*dx**3
         vy2(j2) = REAL(rr, 4)
         !!
         !!
10       CONTINUE
         !!
      END DO
      !!
      !!
   END SUBROUTINE AKIMA_1D_1D







!  SUBROUTINE AKIMA_1D_3D(vx1, vy1, vx2, vy2)
!     !!
!     REAL(4), DIMENSION(:,:,:), INTENT(in)  :: vx1, vy1, vx2
!     REAL(4), DIMENSION(:,:,:), INTENT(out) :: vy2
!     !!___________________________________________________
!     !!
!     !! Local :
!     INTEGER :: nx, ny, nz1, nz2
!     !!
!     REAL(4), DIMENSION(:,:,:), ALLOCATABLE  :: vx1p4, vy1p4, vslp
!     !!
!     REAL(4), DIMENSION(:,:), ALLOCATABLE  :: a0, a1, a2, a3, dx
!     INTEGER  :: n4, j1, j2, p1, p2, j1p4
!     LOGICAL, DIMENSION(:,:), ALLOCATABLE  :: lfnd
!     !!
!     !!
!      !IF ( SIZE(vx1) /= SIZE(vy1) ) THEN
!      !   PRINT *, 'ERROR, mod_akima_1d.f90 => AKIMA_1D_3D: vx1 and vy1 do not have the same size!'; STOP
!      !END IF
!      !IF ( size(vx2) /= size(vy2) ) THEN
!      !   PRINT *, 'ERROR, mod_akima_1d.f90 => AKIMA_1D_3D: vx2 and vy2 do not have the same size!'; STOP
!      !END IF
!      !!
!      nx =  SIZE(vx1,1)
!      ny =  SIZE(vx1,2)
!      nz1 = SIZE(vx1,3)
!      nz2 = SIZE(vx2,3)
!      !!
!      n4 = nz1 + 4
!      !!
!      ALLOCATE ( vx1p4(nx,ny,n4), vy1p4(nx,ny,n4), vslp(nx,ny,n4), lfnd(nx,ny) )
!      ALLOCATE ( a0(nx,ny), a1(nx,ny), a2(nx,ny), a3(nx,ny), dx(nx,ny) )
!      !!
!      vx1p4(:,:,3:nz1+2) = vx1(:,:,:)
!      vy1p4(:,:,3:nz1+2) = vy1(:,:,:)
!      !!
!      !! Extending input X array :
!      !! =========================
!      !! - if vx is regularly spaced, it's not a big deal, otherwise we use
!      !!   what's been proposed by Akima (1970) :
!      !!
!      !! Bottom (or West) :
!      vx1p4(:,:,2) =  vx1p4(:,:,4) - (vx1p4(:,:,5) - vx1p4(:,:,3))
!      vx1p4(:,:,1) =  vx1p4(:,:,3) - (vx1p4(:,:,5) - vx1p4(:,:,3))
!      !!
!      !! Top (or East) :
!      vx1p4(:,:,n4-1) =  vx1p4(:,:,n4-3) + vx1p4(:,:,n4-2) - vx1p4(:,:,n4-4)
!      vx1p4(:,:,n4)   =  vx1p4(:,:,n4-2) + vx1p4(:,:,n4-2) - vx1p4(:,:,n4-4)
!      !!
!      !!
!      !! Now extrapolating input Y values on these 4 extra points :
!      !! ==========================================================
!      !!
!      !! Bottom :
!      CALL extra_2_bottom(vx1p4(:,:,5), vx1p4(:,:,4), vx1p4(:,:,3), vx1p4(:,:,2), vx1p4(:,:,1), &
!         &            vy1p4(:,:,5), vy1p4(:,:,4), vy1p4(:,:,3), vy1p4(:,:,2), vy1p4(:,:,1) )
!      !!
!      !! Top :
!      CALL extra_2_top(vx1p4(:,:,n4-4), vx1p4(:,:,n4-3), vx1p4(:,:,n4-2), vx1p4(:,:,n4-1), vx1p4(:,:,n4), &
!         &            vy1p4(:,:,n4-4), vy1p4(:,:,n4-3), vy1p4(:,:,n4-2), vy1p4(:,:,n4-1), vy1p4(:,:,n4) )
!      !!
!      !!
!      !! Computing slopes :
!      !! ==================
!      !!
!      CALL SLOPES_3D(vx1p4, vy1p4, vslp)
!      !!
!      !!
!      !!
!      !! Ok! Now, each point of the target grid must be interpolated :
!      !! =============================================================
!      DO j2 = 1, nz2
!         !!
!         j1   = 1
!         lfnd(:,:) = .FALSE.
!         !!
!         !IF ( vx2(j2) < vx1(1) ) THEN
!         !   vy2(j2) = vy1(1)
!         !   GOTO 10
!         !END IF
!         !!
!         !!
!         DO WHILE ( j1 < nz1 )
!            WHERE ( .NOT. lfnd )
!               !!
!               !!
!               WHERE ( (vx2(:,:,j2) > vx1(:,:,j1)).AND.(vx2(:,:,j2) < vx1(:,:,j1+1)) )
!                  j1p4  = j1 + 2
!                  p1    = j1p4
!                  p2    = j1p4 + 1
!                  lfnd = .TRUE.
!                  !!
!               ELSEWHERE
!                  !!
!                  WHERE ( vx2(:,:,j2) == vx1(:,:,j1) )
!                     vy2(:,:,j2) = vy1(:,:,j1)
!                     !GOTO 10
!                  ELSEWHERE
!                     WHERE ( vx2(:,:,j2) == vx1(:,:,j1+1) )
!                        vy2(:,:,j2) = vy1(:,:,j1+1)
!                        !GOTO 10
!                     END WHERE
!                  END WHERE
!                  !!
!               END WHERE
!               !!
!               j1 = j1 + 1
!               !!
!            END WHERE
!         END DO
!
!         IF ( j1 == nz1) THEN
!            WHERE (.NOT. lfnd)
!               vy2(:,:,j2) = vy1(:,:,nz1)
!               !GOTO 10
!            END WHERE
!         END IF
!
!         !IF ( .NOT. lfnd ) THEN
!         !   vy2(:,:,j2) = -6666.0
!         !   GOTO 10
!         !END IF
!
!
!
!         !!
!         !! Estimating vy2(j2) with a third order polynome :
!         !! -----------------------------------------------
!         !! Coefficients of the polynome
!         !!
!         !! MIND ! : p1 and p2 are given in 'vp4'
!         a0(:,:) = vy1p4(:,:,p1)
!         a1(:,:) = vslp(:,:,p1)
!         a2(:,:) = (3*(vy1p4(:,:,p2) - vy1p4(:,:,p1))/(vx1p4(:,:,p2) - vx1p4(:,:,p1)) - 2*vslp(:,:,p1) - vslp(:,:,p2)) &
!            & /(vx1p4(:,:,p2) - vx1p4(:,:,p1))
!         a3(:,:) = (vslp(:,:,p1) + vslp(:,:,p2) - 2*(vy1p4(:,:,p2)-vy1p4(:,:,p1))/(vx1p4(:,:,p2)-vx1p4(:,:,p1)))  &
!            &         /(vx1p4(:,:,p2) - vx1p4(:,:,p1))**2
!         !!
!         dx(:,:) = vx2(j2) - vx1p4(:,:,p1)
!         !!
!         vy2(j2) = a0(:,:) + a1(:,:)*dx(:,:) + a2(:,:)*dx(:,:)**2 + a3(:,:)*dx(:,:)**3
!         !!
!         !!
!10       CONTINUE
!        !!
!     END DO
!     !!
!     DEALLOCATE ( vx1p4, vy1p4, vslp, lfnd )
!     DEALLOCATE ( a0, a1, a2, a3, dx )
!     !!
!  END SUBROUTINE AKIMA_1D_3D



   SUBROUTINE slopes_1d(x, y, slope)
      !!
      !!   C o m p u t a t i o n   o f   s l o p e s  :
      !!   --------------------------------------------
      !! * For each given point, we compute its local slope, from the four
      !!   surrounding points and itself; following Akima 1970, p.590, eq.(1)
      !! * The slopes are stored into an array 'slp(:)' of length n
      !!
      !! INPUT :
      !! -------
      !!          - x       = 1D array (n) contening input data coordinates
      !!          - y       = 1D array (n) contains values y(x)
      !!
      !! OUTPUT :
      !! --------
      !!          - slope  = 1D array (n) contening the slopes
      !!
      !!
      !! Author : Laurent Brodeau, dec. 2004
      !!
      !!-------------------------------------------------------------------
      REAL(8), DIMENSION(:), INTENT(in)  :: x, y
      REAL(8), DIMENSION(:), INTENT(out) :: slope
      !!_______________________________________________________
      !!
      !! Local :
      INTEGER :: nz, k
      REAL(8) :: m1, m2, m3, m4
      !!
      nz = SIZE(x)
      !!
      !! Treating middle of array ( at least to points away from the bordures ) :
      !!-------------------------------------------------------------------------
      DO k=3, nz-2
         m1 =  ( y(k-1) - y(k-2) ) / ( x(k-1) - x(k-2) )
         m2 =  ( y(k)   - y(k-1) ) / ( x(k)   - x(k-1) )
         m3 =  ( y(k+1) - y(k)   ) / ( x(k+1) - x(k)   )
         m4 =  ( y(k+2) - y(k+1) ) / ( x(k+2) - x(k+1) )
         !!
         IF ( (m1 == m2).and.(m3 == m4) ) THEN
            slope(k) = 0.5*(m2 + m3)
         ELSE
            slope(k) =   ( ABS(m4-m3)*m2 + ABS(m2-m1)*m3 ) &
               &     / ( ABS(m4-m3)    + ABS(m2-m1) )
         END IF
         !!
      END DO
      !!
      !!
      !! Treating 'second' and 'before last' points :
      !! --------------------------------------------
      !! Point k=2
      m2 =  ( y(2)   - y(1)   ) / ( x(2)   - x(1)   )
      m3 =  ( y(3)   - y(2)   ) / ( x(3)   - x(2)   )
      slope(2)   = 0.5*( m2 + m3 )
      !! Point k=n-1
      m2 =  ( y(nz-1) - y(nz-2) ) / ( x(nz-1) - x(nz-2) )
      m3 =  ( y(nz)   - y(nz-1) ) / ( x(nz)   - x(nz-1) )
      slope(nz-1) = 0.5*( m2 + m3 )
      !!
      !! Treating 'first' and 'last' points :
      !! --------------------------------------------
      !! Point k = 1
      m3 =  ( y(2)   - y(1)   ) / ( x(2)   - x(1)   )
      slope(1)   = m3
      !! Point k = nz
      m2 =  ( y(nz) - y(nz-1) ) / ( x(nz) - x(nz-1) )
      slope(nz)   = m2
      !!
      !!
   END SUBROUTINE slopes_1d



   SUBROUTINE slopes_3d(x, y, slope)
      REAL(4), DIMENSION(:,:,:), INTENT(in)  :: x, y
      REAL(4), DIMENSION(:,:,:), INTENT(out) :: slope
      !! Local :
      INTEGER :: k, nx, ny, nz
      REAL(4), DIMENSION(:,:), ALLOCATABLE :: m1, m2, m3, m4
      !!
      nx = SIZE(x,1)
      ny = SIZE(x,2)
      nz = SIZE(x,3)
      !!
      ALLOCATE ( m1(nx,ny), m2(nx,ny), m3(nx,ny), m4(nx,ny) )
      !!
      DO k=3, nz-2
         m1(:,:) =  ( y(:,:,k-1) - y(:,:,k-2) ) / ( x(:,:,k-1) - x(:,:,k-2) )
         m2(:,:) =  ( y(:,:,k)   - y(:,:,k-1) ) / ( x(:,:,k)   - x(:,:,k-1) )
         m3(:,:) =  ( y(:,:,k+1) - y(:,:,k)   ) / ( x(:,:,k+1) - x(:,:,k)   )
         m4(:,:) =  ( y(:,:,k+2) - y(:,:,k+1) ) / ( x(:,:,k+2) - x(:,:,k+1) )
         !!
         WHERE ( (m1 == m2).AND.(m3 == m4) )
            slope(:,:,k) = 0.5*(m2(:,:) + m3(:,:))
         ELSEWHERE
            slope(:,:,k) =  ( ABS(m4(:,:)-m3(:,:))*m2(:,:) + ABS(m2(:,:)-m1(:,:))*m3(:,:) ) &
               &     / ( ABS(m4(:,:)-m3(:,:))    + ABS(m2(:,:)-m1(:,:)) )
         END WHERE
         !!
      END DO
      !!
      !!
      !! Treating 'second' and 'before last' points :
      !! --------------------------------------------
      !! Point k=2
      m2(:,:) =  ( y(:,:,2)   - y(:,:,1)   ) / ( x(:,:,2)   - x(:,:,1)   )
      m3(:,:) =  ( y(:,:,3)   - y(:,:,2)   ) / ( x(:,:,3)   - x(:,:,2)   )
      slope(:,:,2)   = 0.5*( m2(:,:) + m3(:,:) )
      !! Point k=n-1
      m2(:,:) =  ( y(:,:,nz-1) - y(:,:,nz-2) ) / ( x(:,:,nz-1) - x(:,:,nz-2) )
      m3(:,:) =  ( y(:,:,nz)   - y(:,:,nz-1) ) / ( x(:,:,nz)   - x(:,:,nz-1) )
      slope(:,:,nz-1) = 0.5*( m2(:,:) + m3(:,:) )
      !!
      !! Treating 'first' and 'last' points :
      !! --------------------------------------------
      !! Point k = 1
      m3(:,:) = ( y(:,:,2)   - y(:,:,1)   ) / ( x(:,:,2)   - x(:,:,1)   )
      slope(:,:,1) = m3(:,:)
      !! Point k = nz
      m2(:,:) =  ( y(:,:,nz) - y(:,:,nz-1) ) / ( x(:,:,nz) - x(:,:,nz-1) )
      slope(:,:,nz) = m2(:,:)
      !!
      DEALLOCATE ( m1, m2, m3, m4 )
      !!
   END SUBROUTINE slopes_3d





   SUBROUTINE extra_2_top_3d(x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)
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
      REAL(8), DIMENSION(:,:),                   INTENT(in)   :: x1, x2, x3, x4, x5, y1, y2, y3
      REAL(8), DIMENSION(SIZE(x1,1),SIZE(x1,2)), INTENT(out)  :: y4, y5
      !!
      !! Local :
      REAL(8),  DIMENSION(:,:), ALLOCATABLE :: A, B, C, D, ALF, BET
      !!
      INTEGER :: nx, ny
      !!
      nx = SIZE(x1,1)
      ny = SIZE(x1,2)
      !!
      ALLOCATE ( A(nx,ny), B(nx,ny), C(nx,ny), D(nx,ny), ALF(nx,ny), BET(nx,ny) )
      !!
      A    = x2 - x1
      B    = x3 - x2
      C    = x4 - x3
      D    = x5 - x4
      !!
      ALF  = y2 - y1
      BET  = y3 - y2
      !!
      WHERE ( (A == 0.).OR.(B == 0.).OR.(C == 0.) )
         y4 = y3
         y5 = y3
      ELSEWHERE
         y4   = C*(2*BET/B - ALF/A) + y3
         y5   = y4 + y4*D/C + BET*D/B - ALF*D/A - y3*D/C
      END WHERE
      !!
      DEALLOCATE ( A, B, C, D, ALF, BET )
      !!
   END SUBROUTINE extra_2_top_3d


   SUBROUTINE extra_2_bottom_3d(x5, x4, x3, x2, x1, y5, y4, y3, y2, y1)
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
      REAL(8), DIMENSION(:,:),                    INTENT(in) :: x1, x2, x3, x4, x5, y5, y4, y3
      REAL(8), DIMENSION(SIZE(x5,1),SIZE(x5,2)), INTENT(out) :: y1, y2
      REAL(8),  DIMENSION(:,:), ALLOCATABLE                  :: A, B, C, D, ALF, BET
      !!
      INTEGER :: nx, ny
      !!
      nx = SIZE(x5,1)
      ny = SIZE(x5,2)
      !!
      ALLOCATE ( A(nx,ny), B(nx,ny), C(nx,ny), D(nx,ny), ALF(nx,ny), BET(nx,ny) )
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
      WHERE ( (A == 0.).OR.(B == 0.).OR.(C == 0.) )
         y2 = y3
         y1 = y3
      ELSEWHERE
         y2   = C*(2*BET/B - ALF/A) + y3
         y1   = y2 + y2*D/C + BET*D/B - ALF*D/A - y3*D/C
      END WHERE
      !!
      DEALLOCATE ( A, B, C, D, ALF, BET )
      !!
   END SUBROUTINE extra_2_bottom_3d


END MODULE MOD_AKIMA_1D
!!
