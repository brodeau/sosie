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
      MODULE PROCEDURE AKIMA_1D_1D, AKIMA_1D_3D
   END INTERFACE AKIMA_1D

   PRIVATE


   PUBLIC :: AKIMA_1D


CONTAINS



   SUBROUTINE AKIMA_1D_1D( vz1, vf1,  vz2, vf2)
      !!
      REAL(4), DIMENSION(:), INTENT(in)  :: vz1, vf1, vz2
      REAL(4), DIMENSION(:), INTENT(out) :: vf2
      !!___________________________________________________
      INTEGER :: nk1, nk2, nk1e, jk1, jk2, jp1, jp2, jk1e
      REAL(4)  :: a0, a1, a2, a3, dz, dz2
      LOGICAL  :: lfnd, l_do_interp
      REAL(4), DIMENSION(:), ALLOCATABLE  :: vz1e, vf1e, vslp
      !!___________________________________________________
      !!
      !!
      IF ( size(vz1) /= size(vf1) ) THEN
         PRINT *, 'ERROR, mod_akima_1d.f90 => AKIMA_1D_1D: vz1 and vf1 do not have the same size!'; STOP
      END IF
      IF ( size(vz2) /= size(vf2) ) THEN
         PRINT *, 'ERROR, mod_akima_1d.f90 => AKIMA_1D_1D: vz2 and vf2 do not have the same size!'; STOP
      END IF
      !!
      nk1 = size(vz1)
      nk2 = size(vz2)
      !!
      nk1e = nk1 + 4
      !!
      ALLOCATE ( vz1e(nk1e), vf1e(nk1e), vslp(nk1e) )
      !!
      vz1e(3:nk1+2) = vz1(:)
      vf1e(3:nk1+2) = vf1(:)
      !!
      !! Extending input X array :
      !! =========================
      !! - if vz is regularly spaced, it's not a big deal, otherwise we use
      !!   what's been proposed by Akima (1970) :
      !!
      !! Bottom (or West) :
      vz1e(2) =  vz1e(4) - (vz1e(5) - vz1e(3))
      vz1e(1) =  vz1e(3) - (vz1e(5) - vz1e(3))
      !!
      !! Top (or East) :
      vz1e(nk1e-1) =  vz1e(nk1e-3) + vz1e(nk1e-2) - vz1e(nk1e-4)
      vz1e(nk1e)   =  vz1e(nk1e-2) + vz1e(nk1e-2) - vz1e(nk1e-4)
      !!
      !!
      !! Now extrapolating input Y values on these 4 extra points :
      !! ==========================================================
      !!
      !! Bottom (or West) :
      CALL extra_2_west(vz1e(5), vz1e(4), vz1e(3), vz1e(2), vz1e(1), &
         &              vf1e(5), vf1e(4), vf1e(3), vf1e(2), vf1e(1) )
      !!
      !! Top (or East) :
      CALL extra_2_east(vz1e(nk1e-4), vz1e(nk1e-3), vz1e(nk1e-2), vz1e(nk1e-1), vz1e(nk1e), &
         &              vf1e(nk1e-4), vf1e(nk1e-3), vf1e(nk1e-2), vf1e(nk1e-1), vf1e(nk1e) )
      !!
      !!
      !! Computing slopes :
      !! ==================
      !!
      CALL SLOPES_1D(vz1e, vf1e, vslp)
      !!
      !!
      !!
      !! Ok! Now, each point of the target grid must be interpolated :
      !! =============================================================
      !!
      !! We INTERPOLATE and don't extrapolate so checking the bounds !
      !! input X array is supposed to be orgnised so :
      !    xstart = vz1(1)  !  xstart = minval(x1) ;   xstop  = maxval(x1)
      !    xstop  = vz1(nk1)
      !!
      !!
      !! Treating each target point :
      !! ============================
      !!
      !! !LB : we assume so far that x2 is totally in-organised !!!
      !!        -> should be fixed !!!
      !!        x1 is totally ORGANISED (increasing!)
      !!
      !jk1_o = 1
      !!
      DO jk2 = 1, nk2
         !!
         !jk1  = MAX(jk1_o,1)
         l_do_interp = .TRUE.
         jk1 = 1
         lfnd = .FALSE.
         !!
         !! Persistence: if point of output grid is shallower
         !!              than first point of input grid
         !!              (should occur only when jk2 = 1)
         IF ( vz2(jk2) < vz1(1) ) THEN
            vf2(jk2) = vf1(1)
            l_do_interp = .FALSE.
         END IF
         !!
         !!
         DO WHILE ( (.NOT. lfnd).AND.(l_do_interp) )
            !!
            IF ( jk1 == nk1 ) EXIT
            !!
            ! JMM fix in case of tgr depth == src depth : w/o fix put missing_value on this level
            !            IF ( (vz2(jk2) > vz1(jk1)).and.(vz2(jk2) < vz1(jk1+1)) ) THEN
            IF ( (vz2(jk2) > vz1(jk1)).and.(vz2(jk2) <= vz1(jk1+1)) ) THEN
               jk1e  = jk1 + 2
               jp1    = jk1e
               jp2    = jk1e + 1
               lfnd = .TRUE.
               !!
            ELSE
               !!
               IF ( vz2(jk2) == vz1(jk1) ) THEN
                  vf2(jk2) = vf1(jk1)
                  l_do_interp = .FALSE.
                  EXIT
               ELSE
                  IF ( vz2(jk2) == vz1(jk1+1) ) THEN
                     vf2(jk2) = vf1(jk1+1)
                     l_do_interp = .FALSE.
                     EXIT
                  END IF
               END IF
               !!
            END IF
            !!
            jk1 = jk1 + 1
            !!
         END DO ! DO WHILE ( .NOT. lfnd )

         !!  (PM) Take bottom value below the last data
         IF (( jk1 == nk1) .and. (.not. lfnd)) THEN
            vf2(jk2) = vf1(nk1)
            l_do_interp = .FALSE.
         END IF
         !!
         IF ( .not. lfnd ) THEN
            vf2(jk2) = -6666.0
            l_do_interp = .FALSE.
         END IF
         !!
         !!
         !!
         IF ( l_do_interp ) THEN
            !! Estimating vf2(jk2) with a third order polynome :
            !! -----------------------------------------------
            !! Coefficients of the polynome
            !!
            !! MIND ! : jp1 and jp2 are given in 've' (n+4 extended arrays)
            a0 = vf1e(jp1)
            a1 = vslp(jp1)
            dz = 1./(vz1e(jp2) - vz1e(jp1))
            dz2 = dz*dz
            a2 = ( 3*(vf1e(jp2) - vf1e(jp1))*dz - 2*vslp(jp1) - vslp(jp2) ) * dz
            a3 = ( vslp(jp1) + vslp(jp2) - 2*(vf1e(jp2)-vf1e(jp1))/(vz1e(jp2)-vz1e(jp1)) ) * dz2
            !!
            dz = vz2(jk2) - vz1e(jp1)
            dz2 = dz*dz
            !!
            vf2(jk2) = a0 + a1*dz + a2*dz2 + a3*dz2*dz
            !!
         END IF
         !!
      END DO  !DO jk2 = 1, nk2
      !!
      DEALLOCATE ( vz1e, vf1e, vslp )
      !!
   END SUBROUTINE AKIMA_1D_1D










   SUBROUTINE AKIMA_1D_3D( vz1, xf1, vz2, xf2, rfill_val)
      !!
      REAL(4), DIMENSION(:),     INTENT(in)  :: vz1, vz2
      REAL(4), DIMENSION(:,:,:), INTENT(in)  :: xf1
      REAL(4), DIMENSION(:,:,:), INTENT(out) :: xf2
      REAL(4), INTENT(in)  :: rfill_val
      !!___________________________________________________
      INTEGER :: nk1, nk2, nk1e, nx, ny, jk1, jk2, jp1, jp2, jk1e
      REAL(4) :: dz, dz2
      LOGICAL :: lfnd, l_do_interp
      REAL(4), DIMENSION(:,:),   ALLOCATABLE :: a0, a1, a2, a3
      REAL(4), DIMENSION(:),     ALLOCATABLE :: vz1e
      REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: xf1e, xslp
      !!___________________________________________________

      !! Vertical dimensions and n+4 extension:
      nk1 = SIZE(vz1)
      nk2 = SIZE(vz2)
      nk1e = nk1 + 4

      !! Horizontal dimensions:
      nx =  SIZE(xf1,1)
      ny =  SIZE(xf1,2)

      IF ( (SIZE(xf2,1)/=nx).OR.(SIZE(xf2,2)/=ny) ) THEN
         PRINT *, 'ERROR, mod_akima_1d.f90 => AKIMA_1D_3D: xf1 and xf2 do not have the same horizontal shape!!!'
         STOP
      END IF

      ALLOCATE ( vz1e(nk1e), xf1e(nx,ny,nk1e), xslp(nx,ny,nk1e) )
      ALLOCATE ( a0(nx,ny), a1(nx,ny), a2(nx,ny), a3(nx,ny) )

      vz1e(    3:nk1+2) =     vz1(:)
      xf1e(:,:,3:nk1+2) = xf1(:,:,:)

      !! Bottom (or West) :
      vz1e(2) =  vz1e(4) - (vz1e(5) - vz1e(3))
      vz1e(1) =  vz1e(3) - (vz1e(5) - vz1e(3))

      !! Top (or East) :
      vz1e(nk1e-1) =  vz1e(nk1e-3) + vz1e(nk1e-2) - vz1e(nk1e-4)
      vz1e(nk1e)   =  vz1e(nk1e-2) + vz1e(nk1e-2) - vz1e(nk1e-4)
      !!
      !!
      !! Now extrapolating input Y values on these 4 extra points :
      !! ==========================================================
      !!
      !! Bottom (or West) :
      CALL extra_2_bottom_3d(vz1e(5), vz1e(4), vz1e(3), vz1e(2), vz1e(1), &
         &                   xf1e(:,:,5), xf1e(:,:,4), xf1e(:,:,3), xf1e(:,:,2), xf1e(:,:,1) )
      !!
      !! Top (or East) :
      CALL extra_2_top_3d(vz1e(nk1e-4), vz1e(nk1e-3), vz1e(nk1e-2), vz1e(nk1e-1), vz1e(nk1e), &
         &                xf1e(:,:,nk1e-4), xf1e(:,:,nk1e-3), xf1e(:,:,nk1e-2), xf1e(:,:,nk1e-1), xf1e(:,:,nk1e) )

      CALL slopes_3d(vz1e, xf1e, xslp)

      DO jk2 = 1, nk2

         l_do_interp = .TRUE.
         jk1 = 1
         lfnd = .FALSE.

         IF ( vz2(jk2) < vz1(1) ) THEN       !! Persistence: if point of output grid is shallower
            xf2(:,:,jk2) = xf1(:,:,1)                !!              than first point of input grid
            l_do_interp = .FALSE.            !!              (should occur only when jk2 = 1)
         END IF

         DO WHILE ( (.NOT. lfnd).AND.(l_do_interp) )
            IF ( jk1 == nk1 ) EXIT
            ! JMM fix in case of tgr depth == src depth : w/o fix put missing_value on this level
            !            IF ( (vz2(jk2) > vz1(jk1)).AND.(vz2(jk2) < vz1(jk1+1)) ) THEN
            IF ( (vz2(jk2) > vz1(jk1)).AND.(vz2(jk2) <= vz1(jk1+1)) ) THEN
               jk1e  = jk1 + 2
               jp1    = jk1e
               jp2    = jk1e + 1
               lfnd = .TRUE.
            ELSE
               IF ( vz2(jk2) == vz1(jk1) ) THEN
                  xf2(:,:,jk2) = xf1(:,:,jk1)
                  l_do_interp = .FALSE.
                  EXIT
               ELSE
                  IF ( vz2(jk2) == vz1(jk1+1) ) THEN
                     xf2(:,:,jk2) = xf1(:,:,jk1+1)
                     l_do_interp = .FALSE.
                     EXIT
                  END IF
               END IF
            END IF
            jk1 = jk1 + 1
         END DO ! DO WHILE ( .NOT. lfnd )

         IF ( (jk1 == nk1).AND.(.NOT. lfnd) ) THEN  !!  (PM) Take bottom value below the last data
            xf2(:,:,jk2) = xf1(:,:,nk1)
            l_do_interp = .FALSE.
         END IF

         IF ( .NOT. lfnd ) THEN
            xf2(:,:,jk2) = rfill_val
            l_do_interp = .FALSE.
         END IF

         IF ( l_do_interp ) THEN
            !! Estimating xf2(:,:,jk2) with a third order polynome :
            !! -----------------------------------------------
            !! Coefficients of the polynome
            !!
            !! MIND ! : jp1 and jp2 are given in 've' (n+4 extended arrays)
            a0 = xf1e(:,:,jp1)
            a1 = xslp(:,:,jp1)
            dz = 1./(vz1e(jp2) - vz1e(jp1))
            dz2 = dz*dz
            a2 = ( 3*(xf1e(:,:,jp2) - xf1e(:,:,jp1))*dz - 2*xslp(:,:,jp1) - xslp(:,:,jp2) ) * dz
            a3 = ( xslp(:,:,jp1) + xslp(:,:,jp2) - 2*(xf1e(:,:,jp2)-xf1e(:,:,jp1))/(vz1e(jp2)-vz1e(jp1)) ) * dz2
            !!
            dz = vz2(jk2) - vz1e(jp1)
            dz2 = dz*dz
            !!
            xf2(:,:,jk2) = a0 + a1*dz + a2*dz2 + a3*dz2*dz
            !!
         END IF
         !!
      END DO  !DO jk2 = 1, nk2
      !!
      DEALLOCATE ( vz1e, xf1e, xslp )
      DEALLOCATE ( a0, a1, a2, a3   )
      !!

      ! SURFACE: persistence when first layer is shallower in target domain than
      ! in source domain:
      jk2 = 1
      DO WHILE ( (vz2(jk2) < vz1(1)).AND.(jk2<nk2) )
         xf2(:,:,jk2) = xf1(:,:,1)
         jk2 = jk2 + 1
      END DO

      ! BOTTOM: persistence when last layer is deeper in target domain than in
      ! source domain:
      jk2 = nk2
      DO WHILE ( (vz2(jk2) > vz1(nk1)).AND.(jk2>0) )
         xf2(:,:,jk2) = xf1(:,:,nk1)
         jk2 = jk2 - 1
      END DO

   END SUBROUTINE AKIMA_1D_3D




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
      REAL(4), DIMENSION(:), INTENT(in)  :: x, y
      REAL(4), DIMENSION(:), INTENT(out) :: slope
      !!_______________________________________________________
      !!
      !! Local :
      INTEGER :: nz, k
      REAL(4) :: m1, m2, m3, m4
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



   SUBROUTINE slopes_3d(vz, Xf, slope)
      REAL(4), DIMENSION(:),     INTENT(in)  :: vz
      REAL(4), DIMENSION(:,:,:), INTENT(in)  :: Xf
      REAL(4), DIMENSION(:,:,:), INTENT(out) :: slope
      !! Local :
      INTEGER :: k, nx, ny, nz
      REAL(4), DIMENSION(:,:), ALLOCATABLE :: m1, m2, m3, m4
      !!
      nx = SIZE(Xf,1)
      ny = SIZE(Xf,2)
      nz = SIZE(Xf,3)

      IF ( SIZE(vz)/= nz) THEN
         PRINT *, 'ERROR, mod_akima_1d.f90 => slopes_3d: Xf and vz do not have the same number of levels!!!'
         STOP
      END IF

      ALLOCATE ( m1(nx,ny), m2(nx,ny), m3(nx,ny), m4(nx,ny) )
      !!
      DO k=3, nz-2
         m1(:,:) =  ( Xf(:,:,k-1) - Xf(:,:,k-2) ) / ( vz(k-1) - vz(k-2) )
         m2(:,:) =  ( Xf(:,:,k)   - Xf(:,:,k-1) ) / ( vz(k)   - vz(k-1) )
         m3(:,:) =  ( Xf(:,:,k+1) - Xf(:,:,k)   ) / ( vz(k+1) - vz(k)   )
         m4(:,:) =  ( Xf(:,:,k+2) - Xf(:,:,k+1) ) / ( vz(k+2) - vz(k+1) )
         !!
         WHERE ( (m1 == m2).AND.(m3 == m4) )
            slope(:,:,k) = 0.5*(m2(:,:) + m3(:,:))
         ELSEWHERE
            slope(:,:,k) =  ( ABS(m4(:,:)-m3(:,:))*m2(:,:) + ABS(m2(:,:)-m1(:,:))*m3(:,:) ) &
               &          / ( ABS(m4(:,:)-m3(:,:))         + ABS(m2(:,:)-m1(:,:)) )
         END WHERE
         !!
      END DO
      !!
      !!
      !! Treating 'second' and 'before last' points :
      !! --------------------------------------------
      !! Point k=2
      m2(:,:) =  ( Xf(:,:,2)   - Xf(:,:,1)   ) / ( vz(2)   - vz(1)   )
      m3(:,:) =  ( Xf(:,:,3)   - Xf(:,:,2)   ) / ( vz(3)   - vz(2)   )
      slope(:,:,2)   = 0.5*( m2(:,:) + m3(:,:) )
      !! Point k=n-1
      m2(:,:) =  ( Xf(:,:,nz-1) - Xf(:,:,nz-2) ) / ( vz(nz-1) - vz(nz-2) )
      m3(:,:) =  ( Xf(:,:,nz)   - Xf(:,:,nz-1) ) / ( vz(nz)   - vz(nz-1) )
      slope(:,:,nz-1) = 0.5*( m2(:,:) + m3(:,:) )
      !!
      !! Treating 'first' and 'last' points :
      !! --------------------------------------------
      !! Point k = 1
      m3(:,:) = ( Xf(:,:,2)   - Xf(:,:,1)   )  / ( vz(2)   - vz(1)   )
      slope(:,:,1) = m3(:,:)
      !! Point k = nz
      m2(:,:) =  ( Xf(:,:,nz) - Xf(:,:,nz-1) ) / ( vz(nz) - vz(nz-1) )
      slope(:,:,nz) = m2(:,:)
      !!
      DEALLOCATE ( m1, m2, m3, m4 )
      !!
   END SUBROUTINE slopes_3d





   SUBROUTINE extra_2_top_3d(x1, x2, x3, x4, x5, Xf1, Xf2, Xf3, Xf4, Xf5)
      !!============================================================================
      !!
      !! Extrapolates 2 extra east (or north) points of a curve with Akima's 1D method
      !!
      !! Input  : x1, x2, x3, x4, x5, Xf1, Xf2, Xf3
      !! Output : Xf4, Xf5
      !!
      !!                       Author : Laurent BRODEAU, 2007
      !!============================================================================
      REAL(4),                                     INTENT(in)   :: x1, x2, x3, x4, x5
      REAL(4), DIMENSION(:,:),                     INTENT(in)   :: Xf1, Xf2, Xf3
      REAL(4), DIMENSION(SIZE(Xf1,1),SIZE(Xf1,2)), INTENT(out)  :: Xf4, Xf5
      !! Local :
      INTEGER :: nx, ny
      REAL(4) :: A, B, C, D
      REAL(4),  DIMENSION(:,:), ALLOCATABLE :: ALF, BET

      nx = SIZE(Xf1,1)
      ny = SIZE(Xf1,2)
      ALLOCATE ( ALF(nx,ny), BET(nx,ny) )

      A    = x2 - x1
      B    = x3 - x2
      C    = x4 - x3
      D    = x5 - x4

      ALF  = Xf2 - Xf1
      BET  = Xf3 - Xf2

      IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
         Xf4 = Xf3
         Xf5 = Xf3
      ELSE
         Xf4   = C*(2*BET/B - ALF/A) + Xf3
         Xf5   = Xf4 + Xf4*D/C + BET*D/B - ALF*D/A - Xf3*D/C
      END IF

      DEALLOCATE ( ALF, BET )
   END SUBROUTINE extra_2_top_3d


   SUBROUTINE extra_2_bottom_3d(x5, x4, x3, x2, x1, Xf5, Xf4, Xf3, Xf2, Xf1)
      !!============================================================================
      REAL(4),                                    INTENT(in) :: x1, x2, x3, x4, x5
      REAL(4), DIMENSION(:,:),                    INTENT(in) :: Xf5, Xf4, Xf3
      REAL(4), DIMENSION(SIZE(Xf5,1),SIZE(Xf5,2)), INTENT(out) :: Xf1, Xf2
      !! Local:
      INTEGER :: nx, ny
      REAL(4) :: A, B, C, D
      REAL(4), DIMENSION(:,:), ALLOCATABLE :: ALF, BET

      nx = SIZE(Xf5,1)
      ny = SIZE(Xf5,2)
      ALLOCATE ( ALF(nx,ny), BET(nx,ny) )

      A    = x4 - x5
      B    = x3 - x4
      C    = x2 - x3
      D    = x1 - x2

      ALF  = Xf4 - Xf5
      BET  = Xf3 - Xf4

      IF ( (A == 0.).OR.(B == 0.).OR.(C == 0.) ) THEN
         Xf2 = Xf3
         Xf1 = Xf3
      ELSE
         Xf2   = C*(2*BET/B - ALF/A) + Xf3
         Xf1   = Xf2 + Xf2*D/C + BET*D/B - ALF*D/A - Xf3*D/C
      END IF

      DEALLOCATE ( ALF, BET )
   END SUBROUTINE extra_2_bottom_3d









END MODULE MOD_AKIMA_1D
!!
