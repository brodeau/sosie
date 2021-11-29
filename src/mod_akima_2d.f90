MODULE MOD_AKIMA_2D

   !!-----------------------------------------------------------------
   !!         A Method Of Bivariate Interpolation And Smooth
   !!            Surface Fitting Based On Local Procedures
   !!
   !!                        Hiroshi AKIMA, 1974
   !!
   !!  author={Akima, H.},
   !!  title={A Method of Bivariate Interpolation and Smooth Surface Fitting
   !!         Based on Local Procedures},
   !!  journal={Commun. ACM},
   !!  year={1974},
   !!  pages={18-20},
   !!  volume={17},
   !!  number={1}
   !!
   !!
   !!   AUTHORS:
   !!            Coded from scratch by Laurent BRODEAU, 2007
   !!
   !!            Stylish and more efficient method to solve the 16x16 linear system:
   !!            Jean-Michel Brankart, October 2010
   !!
   !!            Last update: Laurent Brodeau, June 2014
   !!
   !!            Contact: https://github.com/brodeau/sosie
   !!
   !!-----------------------------------------------------------------

   USE mod_conf
   USE mod_manip, ONLY: FILL_EXTRA_BANDS
   USE io_ezcdf,  ONLY: TEST_XYZ

   IMPLICIT NONE


   LOGICAL, PUBLIC, SAVE :: &
      &    l_always_first_call  = .FALSE.

   INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: ixy_pos !: table storing source/target grids mapping

   PRIVATE

   PUBLIC :: AKIMA_2D

   REAL(8), PARAMETER :: repsilon = 1.E-9

   INTEGER, PARAMETER  :: nsys = 16 !: Dimmension of the linear sytem to solve

CONTAINS







   SUBROUTINE AKIMA_2D(k_ew_per, X10, Y10, Z1, X20, Y20, Z2,    icall)

      !!================================================================
      !!
      !! INPUT :     k_ew_per : east-west periodicity
      !!                        k_ew_per = -1  --> no periodicity
      !!                        k_ew_per >= 0  --> periodicity with overlap of k_ew_per points
      !!             X10   : 2D source longitude array (ni*nj) or (ni*1) (must be regular!)
      !!             Y10   : 2D source latitude  array (ni*nj) or (nj*1) (must be regular!)
      !!             Z1    : source field on source grid
      !!
      !!             X20   : 2D target longitude array (ni*nj) or (ni*1) (can be irregular)
      !!             Y20   : 2D target latitude  array (ni*nj) or (nj*1) (can be irregular)
      !!
      !! OUTPUT :
      !!             Z2    : input field on target grid
      !!
      !! input (optional)  : icall : if icall=1, will always force 'l_first_call_interp_routine' to .TRUE.
      !!
      !!================================================================


      !! Input/Output arguments
      INTEGER,                 INTENT(in)  :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X10, Y10
      REAL(4), DIMENSION(:,:), INTENT(in)  :: Z1
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X20, Y20
      REAL(4), DIMENSION(:,:), INTENT(out) :: Z2
      INTEGER,       OPTIONAL, INTENT(in)  :: icall


      !! Local variables
      INTEGER :: nx1, ny1, nx2, ny2, ji, jj, kewp

      INTEGER, PARAMETER :: n_extd = 4    ! source grid extension

      INTEGER :: &
         &     ji1, jj1, ji2, jj2,   &
         &     ni1, nj1

      REAL(8), DIMENSION(nsys) ::  vpl

      REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::  poly

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: &
         &    X1, Y1, X2, Y2,   &
         &    zFs , zXs , zYs

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: slpx, slpy, slpxy

      REAL(8) :: &
         &  px2, py2,  &
         &  min_lon1, max_lon1, min_lat1, max_lat1,   &
         &  min_lon2, max_lon2, min_lat2, max_lat2

      REAL(8), DIMENSION(4) :: xy_range_src

      CHARACTER(len=2) :: ctype

      IF ( present(icall) ) THEN
         IF ( icall == 1 ) THEN
            l_first_call_interp_routine = .TRUE.
            l_always_first_call  = .TRUE.
         END IF
      END IF



      !! Create 2D (ni*nj) arrays out of 1d (ni*1 and nj*1) arrays if needed:
      !! => TEST_XYZ tests if a 2D array is a true 2D array (NxM) or a fake (Nx1)
      !!    and returns '1d' if it is a fake 2D array

      ctype = TEST_XYZ(X10, Y10, Z1)
      nx1 = SIZE(Z1,1)
      ny1 = SIZE(Z1,2)
      ALLOCATE ( X1(nx1,ny1) , Y1(nx1,ny1) )
      IF ( ctype == '1d' ) THEN
         FORALL (jj = 1:ny1) X1(:,jj) = X10(:,1)
         FORALL (ji = 1:nx1) Y1(ji,:) = Y10(:,1)
      ELSE
         X1 = X10 ; Y1 = Y10
      END IF

      ctype = '00'
      ctype = TEST_XYZ(X20, Y20, Z2)
      nx2 = SIZE(Z2,1)
      ny2 = SIZE(Z2,2)
      ALLOCATE ( X2(nx2,ny2) , Y2(nx2,ny2) )
      IF ( ctype == '1d' ) THEN
         FORALL (jj=1:ny2) X2(:,jj) = X20(:,1)
         FORALL (ji=1:nx2) Y2(ji,:) = Y20(:,1)
      ELSE
         X2 = X20 ; Y2 = Y20
      END IF


      !!                       S T A R T

      !! Extending the source 2D domain with a frame of 2 points:
      !!    We extend initial 2D array with a frame, adding n_extd points in each
      !!    dimension This is really needed specially for preserving good east-west
      !!    perdiodicity...

      ni1 = nx1 + n_extd  ;   nj1 = ny1 + n_extd

      ALLOCATE ( zFs(ni1,nj1), zXs(ni1,nj1), zYs(ni1,nj1), &
         &       slpx(ni1,nj1),   slpy(ni1,nj1),  slpxy(ni1,nj1), &
         &       poly(ni1-1,nj1-1,nsys)    )

      CALL FILL_EXTRA_BANDS(k_ew_per, X1, Y1, REAL(Z1,8), zXs, zYs, zFs,  is_orca_grid=i_orca_src)

      DEALLOCATE (X1, Y1)



      !! Computation of partial derivatives:
      !! * since we use extended arrays (2-point wide frame) we need to adapt E-W periodicity
      kewp = -1
      IF( k_ew_per > -1 ) kewp = k_ew_per + 4
      CALL SLOPES_AKIMA(kewp, zXs, zYs, zFs, slpx, slpy, slpxy)

      !! Polynome:
      CALL build_pol(zXs, zYs, zFs, slpx, slpy, slpxy, poly)

      DEALLOCATE ( slpx, slpy, slpxy )

      !! Checking if the target grid does not overlap source grid :
      min_lon1 = minval(zXs) ;  max_lon1 = maxval(zXs)
      min_lat1 = minval(zYs) ;  max_lat1 = maxval(zYs)
      xy_range_src(:) = (/ min_lon1,max_lon1 , min_lat1,max_lat1 /)

      min_lon2 = MINVAL(X2)     ;  max_lon2 = MAXVAL(X2)
      min_lat2 = minval(Y2)     ;  max_lat2 = maxval(Y2)

      !! Doing the mapping once for all and saving into ixy_pos:
      IF ( l_first_call_interp_routine ) THEN
         ALLOCATE ( ixy_pos(nx2, ny2, 2) )
         ixy_pos(:,:,:) = 0
         CALL find_nearest_akima( zXs, zYs, xy_range_src, X2, Y2, ixy_pos )
      END IF


      !! Loop on target domain:
      DO jj2 = 1, ny2
         DO ji2 = 1, nx2

            !! The coordinates of current target point are (px2,py2) :
            px2 = X2(ji2,jj2) ;  py2 = Y2(ji2,jj2)

            !! Checking if this belongs to source domain :
            IF ( ((px2>=min_lon1).AND.(px2<=max_lon1)).AND.((py2>=min_lat1).AND.(py2<=max_lat1)) ) THEN

               !! We know the right location from time = 1 :
               ji1 = ixy_pos(ji2,jj2,1)
               jj1 = ixy_pos(ji2,jj2,2)

               !! It's time to interpolate:
               px2 = px2 - zXs(ji1,jj1)
               py2 = py2 - zYs(ji1,jj1)
               vpl = poly(ji1,jj1,:)

               Z2(ji2,jj2) = REAL( pol_val(px2, py2, vpl) , 4)  ! back to real(4)

            ELSE
               Z2(ji2,jj2) = 0.  ! point is not on source domain!
            END IF

         END DO
      END DO

      !! Deallocation :
      DEALLOCATE ( zFs , zXs , zYs, poly, X2, Y2 )

      l_first_call_interp_routine = .FALSE.

      IF ( l_always_first_call ) THEN
         DEALLOCATE ( ixy_pos )
         l_first_call_interp_routine = .TRUE.
      END IF

   END SUBROUTINE AKIMA_2D




   !! ########################
   !! LOCAL PRIVATE ROUTINES :
   !! ########################


   SUBROUTINE build_pol(ZX, ZY, ZF, sx, sy, sxy, XPLNM)

      !!==================================================================
      !!
      !! Compute the 16 coefficients of the polynomes used to interpolate
      !!
      !!==================================================================

      REAL(8), DIMENSION(:,:),   INTENT(in)  :: ZX, ZY, ZF, sx, sy, sxy
      REAL(8), DIMENSION(:,:,:), INTENT(out) :: XPLNM

      !! Local variables :
      !! -----------------
      REAL(8), DIMENSION(:), ALLOCATABLE :: VX
      INTEGER                  :: nx, ny, ji, jj

      REAL(8) :: &
         &   x, x2, x3, y, y2, y3, xy, &
         &   b1, b2, b3, b4, b5, b6, b7, b8, &
         &   b9, b10, b11, b12, b13, b14, b15, b16, &
         
         &   c1, c2, c3, c4, c5, c6, c7, c8, &
         &   c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, &
         &   d1, d2, d3, d4, d5, d6, d7, d8, d9, &
         &   f1, f2, f3, f4, f5, f6

      nx = SIZE(ZF,1)
      ny = SIZE(ZF,2)

      ALLOCATE( VX(nsys) )

      DO jj=1, ny-1
         DO ji=1, nx-1

            VX(:) = 0.

            !! Local dx and dy :
            x = ZX(ji+1,jj) - ZX(ji,jj)
            y = ZY(ji,jj+1) - ZY(ji,jj)

            x2 = x*x ; x3 = x2*x
            y2 = y*y ; y3 = y2*y
            xy = x*y

            !! Vector B, value at each points, d/dx, d/dy and d2/dxdy :
            b1  =  ZF(ji,jj) ; b2  =  ZF(ji+1,jj) ;  b3  =  ZF(ji+1,jj+1) ; b4  =  ZF(ji,jj+1)
            b5  =  sx(ji,jj) ; b6  =  sx(ji+1,jj) ;  b7  =  sx(ji+1,jj+1) ; b8  =  sx(ji,jj+1)
            b9  =  sy(ji,jj) ; b10 =  sy(ji+1,jj) ;  b11 =  sy(ji+1,jj+1) ; b12 =  sy(ji,jj+1)
            b13 = sxy(ji,jj) ; b14 = sxy(ji+1,jj) ;  b15 = sxy(ji+1,jj+1) ; b16 = sxy(ji,jj+1)



            ! Linear System 16x16 to be solved ( A.X = B )
            ! ============================================
            !
            !     (/ 0.    0.    0.   0.   0.    0.    0.   0.  0.   0.   0.  0. 0.  0. 0. 1. /)
            !     (/ 0.    0.    0.   x^3  0.    0.    0.   x^2 0.   0.   0.  x  0.  0. 0. 1. /)
            !     (/ x^3*y^3 x^3*y^2 x^3*y x^3  x^2*y^3 x^2*y^2 x^2*y x^2 x*y^3 x*y^2 x*y x  y^3  y^2 y 1. /)
            !     (/ 0.    0.    0.   0.  0.    0.    0.   0. 0.   0.   0.  0. y^3  y^2 y  1. /)
            !     (/ 0.      0.      0.     0.   0.     0.     0.    0.  0. 0. 0. 1. 0. 0. 0. 0. /)
            !     (/ 0.      0.      0.     3*x^2 0.     0.     0.    2*x 0. 0. 0. 1. 0. 0. 0. 0. /)
            !     (/ 3*x^2*y^3 3*x^2*y^2 3*x^2*y 3*x^2 2*x*y^3 2*x*y^2 2*x*y 2*x y^3 y^2 y  1. 0. 0. 0. 0. /)
            ! A = (/ 0.      0.      0.     0.   0.     0.     0.    0.  y^3 y^2 y  1. 0. 0. 0. 0. /)
            !     (/ 0.      0.     0.  0. 0.      0.     0. 0. 0.     0.    0. 0. 0.   0.  1. 0. /)
            !     (/ 0.      0.     x^3  0. 0.      0.     x^2 0. 0.     0.    x  0. 0.   0.  1. 0. /)
            !     (/ 3*x^3*y^2 2*x^3*y x^3  0. 3*x^2*y^2 2*x^2*y x^2 0. 3*x*y^2 2*x*y x  0. 3*y^2 2*y 1. 0. /)
            !     (/ 0.      0.     0.   0. 0.     0.    0. 0. 0.    0.  0. 0. 3*y^2 2*y 1. 0. /)
            !     (/ 0.      0.     0.   0. 0.     0.    0.  0. 0.   0.  1. 0. 0. 0. 0. 0. /)
            !     (/ 0.      0.   3*x^2  0. 0.     0.    2*x 0. 0.   0.  1. 0. 0. 0. 0. 0. /)
            !     (/ 9*x^2*y^2 6*x^2*y 3*x^2 0. 6*x*y^2 4*x*y 2*x 0. 3*y^2 2*y 1. 0. 0. 0. 0. 0. /)
            !     (/ 0.      0.     0.   0. 0.     0.    0.  0. 3*y^2 2*y 1. 0. 0. 0. 0. 0. /)
            !
            ! X = (/ a33, a32, a31, a30, a23, a22, a21, a20, a13, a12, a11, a10, a03, a02, a01, a00 /)
            !
            ! B = (/  b1,  b2,  b3,  b4,  b5,  b6,  b7,  b8,  b9, b10, b11, b12, b13, b14, b15, b16 /)


            !! I keep Jean-Michel's comments in french cause I just love them,
            !! they are the soul of this code! They are followed my my english
            !! translation. Note that Jean-Michel's english is excellent but he
            !! just thought this would stay between us.
            !! /laurent

            !! 1) D'abord, calculer les 4 inconnues (du système à 16 éqs) qui
            !! dépendent directement des second membres bj (avant de faire la
            !! modification des bj):

            !! 1) First, calculate the 4 unknowns (of the 16-eq. system) that
            !! directly depend on second members bj, and do this before modifying
            !! bj :

            VX(11) = b13 ; VX(12) = b5 ; VX(15) = b9 ; VX(16) = b1



            !! 2) Ensuite, mettre à échelle les seconds membres bj
            !! (b1 à b4 restent inchangés):

            !! 2) Then, scale bj second members (b1 to b4 remain unchanged):

            b5  = x*b5   ; b6  = x*b6   ; b7  = x*b7   ; b8  = x*b8
            b9  = y*b9   ; b10 = y*b10  ; b11 = y*b11  ; b12 = y*b12
            b13 = xy*b13 ; b14 = xy*b14 ; b15 = xy*b15 ; b16 = xy*b16


            !! 3) Puis, résoudre le système avec x=1 et y=1 (ou bien utiliser ta
            !! solution symbolique avec x=1 et y=1), en remarquant que, dans ce
            !! cas particulier, certaines combinaisons apparaissent souvent, et
            !! qu'il vaut mieux les calculer d'abord:

            !! 3) Then, solve the system with x=1 and y=1, taking advantage of the
            !! fact, that in that particular case, some combinations are often
            !! appearing and that it is therefore better to calculate them before:

            !! a) Premierement:
            c1  = b1-b2     ; c2  = b3-b4     ; c3  = b5+b6
            c4  = b7+b8     ; c5  = b9-b10    ; c6  = b11-b12
            c7  = b13+b14   ; c8  = b15+b16   ; c9  = 2*b5+b6
            c10 = b7+2*b8   ; c11 = 2*b13+b14 ; c12 = b15+2*b16
            c13 = b5-b8     ; c14 = b1-b4     ; c15 = b13+b16
            c16 =  2*b13 + b16 ; c17 = b9+b12    ; c18 = 2*b9+b12

            !! b) Deuxiemement:
            d1 = c1+c2   ; d2 = c3-c4   ; d3 = c5-c6
            d4 = c7+c8   ; d5 = c9-c10  ; d6 = 2*c5-c6
            d7 = 2*c7+c8 ; d8 = c11+c12 ; d9 = 2*c11+c12

            !! c) Troisiemement:
            f1 = 2*d1+d2 ; f2 = 2*d3+d4 ; f3 = 2*d6+d7
            f4 = 3*d1+d5 ; f5 = 3*d3+d8 ; f6 = 3*d6+d9

            !! d) De sorte que la solution s'écrit simplement:
            !! d) So that the solution simply writes:
            VX(1)  = 2*f1+f2      ; VX(2)  = -(3*f1+f3) ; VX(3)  = 2*c5+c7
            VX(4)  = 2*c1+c3      ; VX(5)  = -(2*f4+f5) ; VX(6)  = 3*f4+f6
            VX(7)  = -(3*c5+c11)  ; VX(8)  = -(3*c1+c9) ; VX(9)  = 2*c13+c15
            VX(10) = -(3*c13+c16) ; VX(13) = 2*c14+c17  ; VX(14) = -(3*c14+c18)

            !! On remarque même que les seules mulitplications qui apparaissent
            !! sont des doublements ou des triplements, qu'on pourrait donc
            !! remplacer par des additions, par ex.: cc14=c14+c14, puis ccc14=cc14+c14,
            !! et donc résoudre le systeme simplifié sans aucune multiplication
            !! ni division, ce qui est tout de même assez émouvant :-)

            !! One can even notice that the only multiplications that show up are
            !! multiplications by 2 or by 3. It would therefore be possible to
            !! replace them by additions, for example: cc14=c14+c14, then
            !! ccc14=cc14+c14, and therefore solve the simplified system without
            !! any multiplication and division, which I found rather touching :-)


            !! 4) Et finalement, il ne reste plus qu'à remettre à échelle la solution:

            !! 4) Finally, all that remains to be done is to scale back the solution:

            VX(1)=VX(1)/(x3*y3) ; VX(2)=VX(2)/(x3*y2) ; VX(3)=VX(3)/(x3*y) ; VX(4)=VX(4)/x3
            VX(5)=VX(5)/(x2*y3) ; VX(6)=VX(6)/(x2*y2) ; VX(7)=VX(7)/(x2*y) ; VX(8)=VX(8)/x2
            VX(9)=VX(9)/(x*y3)  ; VX(10)=VX(10)/(x*y2); VX(13)=VX(13)/y3   ; VX(14)=VX(14)/y2


            !! Bien sûr, je n'ai pas vérifié et j'ai pu faire des erreurs de calcul
            !! ou de recopiage, mais le principe me semble bon, et on aboutit
            !! à un nombre d'opérations:
            !! Etape 2: 9x
            !! Etape 3: 45+, 26x (ou même simplement 69+)
            !! Etape 4: 8x, 12/
            !! Soit un total de: 45+, 43x, 12/ pour résoudre le systeme,
            !! au lieu de: 84+, 140x, 16/ pour la solution symbolique de SOSIE.
            !!
            !! Je suis conscient que ça ne sert sans doute pas à grand chose
            !! car le calcul n'est fait que pour chaque point de la grille de départ
            !! (si je comprends bien) et ce coût n'était probablement déjà plus dominant.
            !! Il était déjà largement dépassé (j'imagine) par le coût d'évaluation
            !! des polynômes aux noeuds de la grille d'arrivée, du moins quand
            !! sa résolution est beaucoup plus fine que celle de la grille de départ.


            !! Number of operations:
            !! Point 2: 9x
            !! Point 3: 45+, 26x (or even possibly 69+)
            !! Point 4: 8x, 12/
            !! A total of 45+,  43x and 12/ to solve the system
            !! Instead of 84+, 140x and 16/ for the original symbolic solution

            !! Storing the 16 coefficients of the polynome for point [ji,jj] :
            XPLNM(ji,jj,:) = VX(:)

         END DO
      END DO

      DEALLOCATE( VX )

   END SUBROUTINE build_pol



   FUNCTION pol_val(x, y, V)

      !!======================================================
      !!
      !! Give value of polynome pol_val(x,y), polynome coefficients
      !! are stored in to vector V
      !!
      !! Optimizing by using Horner's scheme (http://en.wikipedia.org/wiki/Horner_scheme)
      !!
      !!  "It has been shown that the Horner scheme is optimal,
      !!  in the sense that any algorithm to evaluate an arbitrary polynomial
      !!  must use at least as many operations.
      !!  That the number of additions required is minimal was shown
      !!  by Alexander Ostrowski in 1954; that the number of multiplications
      !!  is minimal by Victor Pan."
      !!
      !! => we perform the minimum number of multiplications possible: 15 !
      !!
      !! Big thanks to Jean-Michel Brankart (Jean-Michel.Brankart@hmg.inpg.fr) for
      !! adapting this scheme to Sosie.
      !!
      !!=====================================================================

      REAL(8)                              :: pol_val
      REAL(8), INTENT(in)                  :: x, y
      REAL(8), DIMENSION(nsys), INTENT(in) :: V

      REAL(8) :: p1, p2, p3, p4

      p1      = ( ( V(1)  * y + V(2)  ) * y + V(3)  ) * y + V(4)
      p2      = ( ( V(5)  * y + V(6)  ) * y + V(7)  ) * y + V(8)
      p3      = ( ( V(9)  * y + V(10) ) * y + V(11) ) * y + V(12)
      p4      = ( ( V(13) * y + V(14) ) * y + V(15) ) * y + V(16)
      pol_val = ( ( p1    * x + p2    ) * x + p3    ) * x + p4

   END FUNCTION pol_val


   SUBROUTINE SLOPES_AKIMA(k_ew, XX, XY, XF, dFdX, dFdY, d2FdXdY)

      !! Slopes ~ partial derivatives of a field ZF according to Akima method
      !! given on a regular gird !!

      !!  k_ew : east-west periodicity on the source file/grid
      !!         k_ew = -1  --> no east-west periodicity (along x)
      !!         k_ew >= 0  --> east-west periodicity with overlap of k_ew points (along x)
      !! lulu

      INTEGER, INTENT(in) :: k_ew

      REAL(8), DIMENSION(:,:), INTENT(in)  :: XX, XY, XF
      REAL(8), DIMENSION(:,:), INTENT(out) :: dFdX, dFdY, d2FdXdY

      !! Local variables :
      REAL(8), DIMENSION(:,:), ALLOCATABLE :: ZX, ZY, ZF
      INTEGER :: nx, ny, ji, jj, ip1, ip2, jp1, jp2
      REAL(8) :: m1, m2, m3, m4, rx
      REAL(8) :: Wx2, Wx3, Wy2, Wy3
      REAL(8) :: d22, e22, d23, e23, d42, e32, d43, e33

      dFdX    = 0.
      dFdY    = 0.
      d2FdXdY = 0.

      nx = SIZE(XF,1)
      ny = SIZE(XF,2)

      !! Extended arrays with a frame of 2 points...
      ALLOCATE ( ZX(nx+4,ny+4), ZY(nx+4,ny+4), ZF(nx+4,ny+4) )
      ! We don't want any smart extrapolation in the north because it has already be done...
      ! (and `k_ew` should take into account that we are dealing with arrays extended with 2-point frame)
      CALL FILL_EXTRA_BANDS(k_ew, XX, XY, XF, ZX, ZY, ZF,  is_orca_grid=0)
      !debug:
      !CALL DUMP_FIELD(REAL(ZX,4), 'slopes_akima_extended_LON.nc', 'lon') !lolo: bug corners!!!
      !CALL DUMP_FIELD(REAL(ZY,4), 'slopes_akima_extended_LAT.nc', 'lat')
      !CALL DUMP_FIELD(REAL(ZF,4), 'slopes_akima_extended_FLD.nc', 'field')
      !STOP'SLOPES_AKIMA'
      !debug.


      !! Treating middle of array ( at least 2 points away from the bordures ) :
      !!------------------------------------------------------------------------

      DO jj=1, ny
         DO ji=1, nx

            ! Initialisation MB Problem when they are close to zero but not zero
            m1=0. ; m2=0. ; m3=0. ; m4=0. ; Wx2=0. ; Wx3=0. ; Wy2=0. ; Wy3=0.

            ip1 = ji+1
            ip2 = ji+2
            jp1 = jj+1
            jp2 = jj+2

            !!   SLOPE / X :
            !!   ***********

            rx = ZX(ip1,jp2) - ZX(ji,  jp2)
            m1 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m1 = (ZF(ip1,jp2) - ZF(ji,  jp2)) / m1

            rx = ZX(ip2,jp2) - ZX(ip1,jp2)
            m2 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m2 = (ZF(ip2,jp2) - ZF(ip1,jp2)) / m2

            rx = ZX(ji+3,jp2) - ZX(ip2,jp2)
            m3 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m3 = (ZF(ji+3,jp2) - ZF(ip2,jp2)) / m3

            rx = ZX(ji+4,jp2) - ZX(ji+3,jp2)
            m4 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m4 = (ZF(ji+4,jp2) - ZF(ji+3,jp2)) / m4

            IF ( (m1 == m2).and.(m3 == m4) ) THEN
               dFdX(ji,jj) = 0.5*(m2 + m3)
            ELSE
               Wx2 = ABS(m4 - m3)
               Wx3 = ABS(m2 - m1)
               dFdX(ji,jj) = ( Wx2*m2 + Wx3*m3 ) / ( Wx2 + Wx3 )
            END IF

            !!   SLOPE / Y :
            !!   ***********
            rx = ZY(ip2,jp1) - ZY(ip2,jj  )
            m1 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m1 = (ZF(ip2,jp1) - ZF(ip2,jj  )) / m1

            rx = ZY(ip2,jp2) - ZY(ip2,jp1)
            m2 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m2 = (ZF(ip2,jp2) - ZF(ip2,jp1)) / m2

            rx = ZY(ip2,jj+3) - ZY(ip2,jp2)
            m3 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m3 = (ZF(ip2,jj+3) - ZF(ip2,jp2)) / m3

            rx = ZY(ip2,jj+4) - ZY(ip2,jj+3)
            m4 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            m4 = (ZF(ip2,jj+4) - ZF(ip2,jj+3)) / m4


            IF ( (m1 == m2).and.(m3 == m4) ) THEN
               dFdY(ji,jj) = 0.5*(m2 + m3)
            ELSE
               Wy2 =  ABS(m4 - m3)
               Wy3 =  ABS(m2 - m1)
               dFdY(ji,jj) =   ( Wy2*m2 + Wy3*m3 ) / ( Wy2 + Wy3 )
            END IF


            !!   CROSS DERIVATIVE /XY :
            !!   **********************
            !! d22 = d(i-1,j-1) = [ z(i-1,j)-z(i-1,j-1) ] / [ y(j) - y(j-1) ]
            rx  = ZY(ip1,jp2) - ZY(ip1,jp1)
            d22 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            d22 = (ZF(ip1,jp2) - ZF(ip1,jp1)) / d22

            !! d23 = d(i-1 , j) = [ z(i-1,j+1)-z(i-1,j) ] / [ y(j+1) - y(j) ]
            rx  = ZY(ip1,jj+3) - ZY(ip1,jp2)
            d23 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            d23 = (ZF(ip1,jj+3) - ZF(ip1,jp2)) / d23

            !! d42 = d(i+1 , j-1) = [ z(i+1 , j) - z(i+1 , j-1) ] / [ y(j) - y(j-1) ]
            rx  = ZY(ji+3,jp2) - ZY(ji+3,jp1)
            d42 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            d42 = (ZF(ji+3,jp2) - ZF(ji+3,jp1)) / d42

            !! d43 = d(i+1 , j) = [ z(i+1 , j+1)-z(i+1 , j) ] / [ y(j+1) - y(j) ]
            rx = ZY(ji+3,jj+3) - ZY(ji+3,jp2)
            d43 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            d43 = (ZF(ji+3,jj+3) - ZF(ji+3,jp2)) / d43


            !! e22  = [ m2 - d22 ] / [ x(i) - x(i-1) ]
            rx  = ZX(ip2,jp1) - ZX(ip1,jp1)
            e22 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            e22 = ( m2 - d22 ) / e22

            !! e23  = [ m3 - d23 ] / [ x(i) - x(i-1) ]
            rx  = ZX(ip2,jp2) - ZX(ip1,jp2)
            e23 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            e23 = ( m3 - d23 ) / e23

            !! e32  = [ d42 - m2 ] / [ x(i+1) - x(i) ]
            rx  = ZX(ji+3,jp2) - ZX(ip2,jp2)
            e32 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            e32 = ( d42 - m2 ) / e32

            !! e33  = [ d43 - m3 ] / [ x(i+1) - x(i) ]
            rx  = ZX(ji+3,jp2) - ZX(ip2,jp2)
            e33 = SIGN(1.d0 , rx) * MAX(ABS(rx),repsilon)
            e33 = ( d43 - m3 ) / e33


            IF ( ((Wx2 == 0).and.(Wx3 == 0)).or.((Wy2 == 0).and.(Wy3 == 0)) ) THEN
               IF ( (Wx2 == 0).and.(Wx3 == 0) ) THEN
                  Wx2 = 1.  ; Wx3 = 1.
               END IF

               IF ( (Wy2 == 0).and.(Wy3 == 0) ) THEN
                  Wy2 = 1.  ; Wy3 = 1.
               END IF

            END IF

            d2FdXdY(ji,jj) = ( Wx2*(Wy2*e22 + Wy3*e23) + Wx3*(Wy2*e32 + Wy3*e33) )  &
               &           / ( (Wx2 + Wx3) * (Wy2 + Wy3) )

         END DO
      END DO

      DEALLOCATE ( ZX, ZY, ZF )

      !debug:
      !CALL DUMP_FIELD(REAL(dFdX,4), 'slopes_akima_extended_dFdX.nc', 'dFdX')
      !CALL DUMP_FIELD(REAL(dFdY,4), 'slopes_akima_extended_dFdY.nc', 'dFdX')
      !CALL DUMP_FIELD(REAL(d2FdXdY,4), 'slopes_akima_extended_d2FdXdY.nc', 'd2FdXdY')
      !STOP'SLOPES_AKIMA'
      !debug.

   END SUBROUTINE SLOPES_AKIMA





   SUBROUTINE find_nearest_akima( plon_src, plat_src, pxyr_src, plon_trg, plat_trg, ixyp_trg )
      !!---------------------------------------------------------------------------------------
      !! Laboriuously scanning the entire source grid to find location of treated point
      !!---------------------------------------------------------------------------------------
      REAL(8), DIMENSION(:,:)  , INTENT(in)  :: plon_src, plat_src
      REAL(8), DIMENSION(4)    , INTENT(in)  :: pxyr_src       ! (/ lon_min,lon_max, lat_min,lat_max /)
      REAL(8), DIMENSION(:,:)  , INTENT(in)  :: plon_trg, plat_trg
      INTEGER, DIMENSION(:,:,:), INTENT(out) :: ixyp_trg

      INTEGER :: jis, jjs, jit, jjt, nxs, nys, nxt, nyt
      REAL(8) :: pxt, pyt
      LOGICAL :: l_x_found, l_y_found

      nxs = SIZE(plon_src,1)
      nys = SIZE(plon_src,2)
      nxt = SIZE(ixyp_trg,1)
      nyt = SIZE(ixyp_trg,2)

      !! Loop on target domain:
      DO jjt = 1, nyt
         DO jit = 1, nxt

            jis = 1
            jjs = 1

            l_x_found = .FALSE.
            l_y_found = .FALSE.

            pxt = plon_trg(jit,jjt)
            pyt = plat_trg(jit,jjt)

            !! Only searching if inside source domain:
            IF ( ((pxt>=pxyr_src(1)).AND.(pxt<=pxyr_src(2))).AND.((pyt>=pxyr_src(3)).AND.(pyt<=pxyr_src(4))) ) THEN

               DO WHILE ( .NOT. (l_x_found .AND. l_y_found) )

                  l_x_found = .FALSE.
                  DO WHILE ( .NOT. l_x_found )
                     IF (jis < nxs) THEN
                        IF ((plon_src(jis,jjs) <= pxt).and.(plon_src(jis+1,jjs) > pxt)) THEN
                           l_x_found = .TRUE.
                        ELSE
                           jis = jis+1
                        END IF
                     ELSE   ! jis = nxs
                        jis = jis-1  ! we are at the top need to use former pol.
                        l_x_found = .TRUE.
                     END IF
                  END DO

                  l_y_found = .FALSE.
                  DO WHILE ( .NOT. l_y_found )
                     IF ( jjs < nys ) THEN
                        IF ((plat_src(jis,jjs) <= pyt).and.(plat_src(jis,jjs+1) > pyt)) THEN
                           l_y_found = .TRUE.
                        ELSE
                           jjs = jjs + 1
                           l_x_found = .FALSE.
                           l_y_found = .TRUE. ! just so that we exit the loop on l_y_found
                        END IF
                     ELSE   ! jjs == nys
                        jjs = nys-1        ! we are using polynome at (ji,nys-1)
                        l_y_found = .TRUE. ! for extreme right boundary
                     END IF
                  END DO

               END DO !DO WHILE ( .NOT. (l_x_found .AND. l_y_found) )

               ixyp_trg(jit,jjt,:) = (/ jis, jjs /)

            END IF !IF ( ((pxt>=pxyr_src(1)).AND.(pxt<=pxyr_src(2))).AND.((pyt>=pxyr_src(3)).AND.(pyt<=pxyr_src(4))) )

         END DO !DO jit = 1, nxt
      END DO !DO jjt = 1, nyt

   END SUBROUTINE find_nearest_akima

END MODULE MOD_AKIMA_2D
