MODULE mod_nemotools

   !! A SET OF ROUTINES DIRECTLY "TAKEN" FROM NEMO GCM...




   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! Ocean        : lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)  Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment
   !!            3.5  ! 2012     (S.Mocavero, I. Epicoco) optimization of BDY comm. via lbc_bdy_lnk and lbc_obc_lnk
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie, G. Reffray)  add a C1D case
   !!            3.6  ! 2015-06  (O. Tintó and M. Castrillo)  add lbc_lnk_multi
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option                              shared memory computing
   !!----------------------------------------------------------------------
   !!   lbc_lnk_2d    : set the lateral boundary condition on a 2D variable on ocean mesh
   !!----------------------------------------------------------------------
   !USE oce             ! ocean dynamics and tracers
   !USE dom_oce         ! ocean space and time domain
   !USE in_out_manager  ! I/O manager
   !USE lbcnfd          ! north fold

   IMPLICIT NONE
   PRIVATE


   INTERFACE lbc_lnk
      MODULE PROCEDURE lbc_lnk_2d_r8, lbc_lnk_2d_r4
   END INTERFACE lbc_lnk




   PUBLIC   lbc_lnk, angle



   REAL(8), PARAMETER, PUBLIC :: &
      &       rPi0 = 3.141592653,     &
      &       rad0 = rPi0/180.0

   !TYPE arrayptr
   !   REAL , DIMENSION (:,:),  POINTER :: pt2d
   !END TYPE arrayptr
   !PUBLIC   arrayptr



   INTEGER, PARAMETER :: wp=8


   !! For lbc_nfd_2d:
   INTEGER, PUBLIC, PARAMETER            ::   jpmaxngh = 3               !:
   INTEGER, PUBLIC                       ::   nsndto, nfsloop, nfeloop   !:
   INTEGER, PUBLIC, DIMENSION (jpmaxngh) ::   isendto                    !: processes to which communicate


CONTAINS

   SUBROUTINE lbc_lnk_2d_r8( nperio, pt2d, cd_type, psgn, pval )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d_r8  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 2D array (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   )           ::   nperio
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(wp), DIMENSION(:,:)    , INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign
      REAL(wp)                    , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !!
      REAL(wp) ::   zland
      !!----------------------------------------------------------------------

      INTEGER :: jpi, jpj, jpim1

      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value (zero by default)
      ELSE                         ;   zland = 0._wp
      ENDIF

      jpi = SIZE(pt2d,1)
      jpj = SIZE(pt2d,2)
      jpim1 = jpi-1

      !
      !                                     ! East-West boundaries
      !                                     ! ====================
      SELECT CASE ( nperio )
         !
      CASE ( 1 , 4 , 6 )                       !** cyclic east-west
         pt2d( 1 ,:) = pt2d(jpim1,:)               ! all points
         pt2d(jpi,:) = pt2d(  2  ,:)
         !
      CASE DEFAULT                             !** East closed  --  West closed
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'U' , 'V' , 'W' )            ! T-, U-, V-, W-points
            pt2d( 1 ,:) = zland
            pt2d(jpi,:) = zland
         CASE ( 'F' )                              ! F-point
            pt2d(jpi,:) = zland
         END SELECT
         !
      END SELECT
      !                                     ! North-South boundaries
      !                                     ! ======================
      SELECT CASE ( nperio )
         !
      CASE ( 2 )                               !**  South symmetric  --  North closed
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'U' , 'W' )                   ! T-, U-, W-points
            pt2d(:, 1 ) = pt2d(:,3)
            pt2d(:,jpj) = zland
         CASE ( 'V' , 'F' )                         ! V-, F-points
            pt2d(:, 1 ) = psgn * pt2d(:,2)
            pt2d(:,jpj) = zland
         END SELECT
         !
      CASE ( 3 , 4 , 5 , 6 )                   !**  North fold  T or F-point pivot  --  South closed
         SELECT CASE ( cd_type )                    ! South : closed
         CASE ( 'T' , 'U' , 'V' , 'W' , 'I' )             ! all points except F-point
            pt2d(:, 1 ) = zland
         END SELECT
         !                                          ! North fold
         CALL lbc_nfd_2d( nperio, pt2d(:,:), cd_type, psgn )
         !
      CASE DEFAULT                             !**  North closed  --  South closed
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
            pt2d(:, 1 ) = zland
            pt2d(:,jpj) = zland
         CASE ( 'F' )                               ! F-point
            pt2d(:,jpj) = zland
         END SELECT
         !
      END SELECT
      !
      !
   END SUBROUTINE lbc_lnk_2d_r8



   SUBROUTINE lbc_lnk_2d_r4( nperio, pt2d, cd_type, psgn, pval )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE lbc_lnk_2d_r8  ***
      !!
      !! ** Purpose :   set lateral boundary conditions on a 2D array (non mpp case)
      !!
      !! ** Method  :   psign = -1 :    change the sign across the north fold
      !!                      =  1 : no change of the sign across the north fold
      !!                      =  0 : no change of the sign across the north fold and
      !!                             strict positivity preserved: use inner row/column
      !!                             for closed boundaries.
      !!----------------------------------------------------------------------
      INTEGER                     , INTENT(in   )           ::   nperio
      CHARACTER(len=1)            , INTENT(in   )           ::   cd_type   ! nature of pt3d grid-points
      REAL(4), DIMENSION(:,:)     , INTENT(inout)           ::   pt2d      ! 2D array on which the lbc is applied
      REAL(wp)                    , INTENT(in   )           ::   psgn      ! control of the sign
      REAL(wp)                    , INTENT(in   ), OPTIONAL ::   pval      ! background value (for closed boundaries)
      !!
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: pt2dr8
      !!----------------------------------------------------------------------

      INTEGER :: jpi, jpj

      jpi = SIZE(pt2d,1)
      jpj = SIZE(pt2d,2)

      ALLOCATE ( pt2dr8(jpi,jpj) )
      pt2dr8(:,:) = REAL(pt2d , wp)
      IF( PRESENT( pval ) ) THEN
         CALL lbc_lnk_2d_r8( nperio, pt2dr8, cd_type, psgn, pval=pval )
      ELSE
         CALL lbc_lnk_2d_r8( nperio, pt2dr8, cd_type, psgn, pval=pval )
      END IF

      pt2d(:,:) = REAL(pt2dr8 , 4)

      DEALLOCATE ( pt2dr8 )

   END SUBROUTINE lbc_lnk_2d_r4







   SUBROUTINE lbc_nfd_2d( npolj, pt2d, cd_type, psgn)
      !!----------------------------------------------------------------------
      !!                  ***  routine lbc_nfd_2d  ***
      !!
      !! ** Purpose :   2D lateral boundary condition : North fold treatment
      !!       without processor exchanges.
      !!
      !! ** Method  :
      !!
      !! ** Action  :   pt2d with updated values along the north fold
      !!----------------------------------------------------------------------
      INTEGER                 , INTENT(in   ) ::   npolj     !: north fold mark (0, 3 or 4)
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type   ! define the nature of ptab array grid-points
      !                                                      ! = T , U , V , F , W points
      REAL(wp)                , INTENT(in   ) ::   psgn      ! control of the sign change
      !                                                      !   = -1. , the sign is changed if north fold boundary
      !                                                      !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d      ! 2D array on which the boundary condition is applied
      !INTEGER , OPTIONAL      , INTENT(in   ) ::   pr2dj     !OPTIONAL! number of additional halos
      !
      INTEGER  ::   ji, jl, ipr2dj
      INTEGER  ::   ijt, iju, ijpj, ijpjm1
      !!----------------------------------------------------------------------
      INTEGER :: jpi, jpj, jpiglo, jpjglo
      !!----------------------------------------------------------------------

      jpi = SIZE(pt2d,1) ; jpiglo = jpi
      jpj = SIZE(pt2d,2) ; jpjglo = jpj


      !SELECT CASE ( jpni )
      !CASE ( 1 )     ;   ijpj = nlcj      ! 1 proc only  along the i-direction
      !CASE DEFAULT   ;   ijpj = 4         ! several proc along the i-direction
      !END SELECT
      !
      !IF( PRESENT(pr2dj) ) THEN           ! use of additional halos
      !   ipr2dj = pr2dj
      !   IF( jpni > 1 )   ijpj = ijpj + ipr2dj
      !ELSE
      !   ipr2dj = 0
      !ENDIF
      !

      !! 1 proc (jpni=1):
      ijpj   = jpj
      ipr2dj = 0


      ijpjm1 = ijpj-1

      !! ORCA2 => npolj => 4
      !! ORCA1 => npolj => 6

      SELECT CASE ( npolj )
         !
      CASE ( 3, 4 )                       ! *  North fold  T-point pivot
         !
         SELECT CASE ( cd_type )
            !
         CASE ( 'T' , 'W' )                               ! T- , W-points
            DO jl = 0, ipr2dj
               DO ji = 2, jpiglo
                  ijt=jpiglo-ji+2
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-2-jl)
               END DO
            END DO
            pt2d(1,ijpj)   = psgn * pt2d(3,ijpj-2)
            DO ji = jpiglo/2+1, jpiglo
               ijt=jpiglo-ji+2
               pt2d(ji,ijpj-1) = psgn * pt2d(ijt,ijpj-1)
            END DO
         CASE ( 'U' )                                     ! U-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-2-jl)
               END DO
            END DO
            pt2d(   1  ,ijpj  ) = psgn * pt2d(    2   ,ijpj-2)
            pt2d(jpiglo,ijpj  ) = psgn * pt2d(jpiglo-1,ijpj-2)
            pt2d(1     ,ijpj-1) = psgn * pt2d(jpiglo  ,ijpj-1)
            DO ji = jpiglo/2, jpiglo-1
               iju = jpiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'V' )                                     ! V-point
            DO jl = -1, ipr2dj
               DO ji = 2, jpiglo
                  ijt = jpiglo-ji+2
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-3-jl)
               END DO
            END DO
            pt2d( 1 ,ijpj)   = psgn * pt2d( 3 ,ijpj-3)
         CASE ( 'F' )                                     ! F-point
            DO jl = -1, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-3-jl)
               END DO
            END DO
            pt2d(   1  ,ijpj)   = psgn * pt2d(    2   ,ijpj-3)
            pt2d(jpiglo,ijpj)   = psgn * pt2d(jpiglo-1,ijpj-3)
            pt2d(jpiglo,ijpj-1) = psgn * pt2d(jpiglo-1,ijpj-2)
            pt2d(   1  ,ijpj-1) = psgn * pt2d(    2   ,ijpj-2)
         CASE ( 'I' )                                     ! ice U-V point (I-point)
            DO jl = 0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'J' )                                     ! first ice U-V point
            DO jl =0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'K' )                                     ! second ice U-V point
            DO jl =0, ipr2dj
               pt2d(2,ijpj+jl) = psgn * pt2d(3,ijpj-1+jl)
               DO ji = 3, jpiglo
                  iju = jpiglo - ji + 3
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
         END SELECT
         !
      CASE ( 5, 6 )                        ! *  North fold  F-point pivot
         !
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'W' )                               ! T-, W-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'U' )                                     ! U-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-1-jl)
               END DO
            END DO
            pt2d(jpiglo,ijpj) = psgn * pt2d(1,ijpj-1)
         CASE ( 'V' )                                     ! V-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo
                  ijt = jpiglo-ji+1
                  pt2d(ji,ijpj+jl) = psgn * pt2d(ijt,ijpj-2-jl)
               END DO
            END DO
            DO ji = jpiglo/2+1, jpiglo
               ijt = jpiglo-ji+1
               pt2d(ji,ijpjm1) = psgn * pt2d(ijt,ijpjm1)
            END DO
         CASE ( 'F' )                               ! F-point
            DO jl = 0, ipr2dj
               DO ji = 1, jpiglo-1
                  iju = jpiglo-ji
                  pt2d(ji,ijpj+jl) = psgn * pt2d(iju,ijpj-2-jl)
               END DO
            END DO
            pt2d(jpiglo,ijpj) = psgn * pt2d(1,ijpj-2)
            DO ji = jpiglo/2+1, jpiglo-1
               iju = jpiglo-ji
               pt2d(ji,ijpjm1) = psgn * pt2d(iju,ijpjm1)
            END DO
         CASE ( 'I' )                                  ! ice U-V point (I-point)
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = 0.e0
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= 0.5 * ( pt2d(ji,ijpj-1-jl) + psgn * pt2d(ijt,ijpj-1-jl) )
               END DO
            END DO
         CASE ( 'J' )                                  ! first ice U-V point
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = 0.e0
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= pt2d(ji,ijpj-1-jl)
               END DO
            END DO
         CASE ( 'K' )                                  ! second ice U-V point
            pt2d( 2 ,ijpj:ijpj+ipr2dj) = 0.e0
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  ijt = jpiglo - ji + 2
                  pt2d(ji,ijpj+jl)= pt2d(ijt,ijpj-1-jl)
               END DO
            END DO
         END SELECT
         !
      CASE DEFAULT                           ! *  closed : the code probably never go through
         !
         SELECT CASE ( cd_type)
         CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'F' )                                   ! F-point
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'I' )                                   ! ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'J' )                                   ! first ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         CASE ( 'K' )                                   ! second ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,ijpj:ijpj+ipr2dj) = 0.e0
         END SELECT
         !
      END SELECT
      !
   END SUBROUTINE lbc_nfd_2d



   SUBROUTINE angle( nperio, glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif,  &
      &                      gcost, gsint, gcosu, gsinu, gcosv, gsinv, gcosf, gsinf     )

      !!========================================================================
      !!
      !!                     *******  ANGLE  *******
      !!
      !! Given the coordinates at T-points and U-points, this routine computes
      !! sinus and cosinus of the angle of the local grid distorsion at T-points
      !!
      !!
      !! INPUT :     - glamt = longitude array at T-points     [REAL]
      !!             - gphit = latitude array at T-points      [REAL]
      !!             - glamu = longitude array at U-points     [REAL]
      !!             - gphiu = latitude array at U-points      [REAL]
      !!             - glamv = longitude array at V-points     [REAL]
      !!             - gphiv = latitude array at V-points      [REAL]
      !!             - glamf = longitude array at F-points     [REAL]
      !!             - gphif = latitude array at F-points      [REAL]
      !!
      !! OUTPUT :
      !! --------    - gcost  = cosinus of the distortion angle at T-points  [REAL]
      !!             - gsint  = sininus of the distortion angle at T-points  [REAL]
      !!             - gcosu  = cosinus of the distortion angle at U-points  [REAL]
      !!             - gsinu  = sininus of the distortion angle at U-points  [REAL]
      !!             - gcosv  = cosinus of the distortion angle at V-points  [REAL]
      !!             - gsinv  = sininus of the distortion angle at V-points  [REAL]
      !!             - gcosf  = cosinus of the distortion angle at F-points  [REAL]
      !!             - gsinf  = sininus of the distortion angle at F-points  [REAL]
      !!
      !!
      !! Author : Laurent Brodeau
      !!           => adapted 'geo2ocean.F90' (NEMOGCM)
      !!
      !!========================================================================

      IMPLICIT NONE

      !! INPUT :
      !! -------
      INTEGER                     , INTENT(in   )   ::   nperio
      REAL(8), DIMENSION(:,:), INTENT(in) ::    &
         &       glamt, gphit, &
         &       glamu, gphiu, &
         &       glamv, gphiv, &
         &       glamf, gphif

      !! OUTPUT :
      !! --------
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcost, gsint, gcosu, gsinu, gcosv, gsinv, gcosf, gsinf

      !! LOCAL :
      !! -------
      INTEGER :: nx, ny, ji, jj

      REAL(8) :: zlam, zphi

      !   REAL(8), DIMENzyuutION(:), ALLOCATABLE :: &
      REAL(8) :: &
         &     zxnpt, znnpt, zxnpu, znnpu, zxnpv, znnpv, zlan, &
         &     zxvvt, zyvvt, znvvt, znuuf, znnpf, znffv, znffu, &
         &     zynpu, zynpv, zxffu, zyffu, zyffv, zxuuf, zyuuf, zxffv, &
         &     zphh, zynpt, zxnpf, zynpf


      !! 2D domain shape:
      nx = SIZE(glamt,1)
      ny = SIZE(glamt,2)

      !loloDO jj = 1, ny
      DO jj = 2, ny-1
         !loloDO ji = 1, nx
         DO ji = 2, nx

            zlam = glamt(ji,jj)     ! north pole direction & modulous (at t-point)
            zphi = gphit(ji,jj)
            zxnpt = 0. - 2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            zynpt = 0. - 2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            znnpt = zxnpt*zxnpt + zynpt*zynpt
            !
            zlam = glamu(ji,jj)     ! north pole direction & modulous (at u-point)
            zphi = gphiu(ji,jj)
            zxnpu = 0. - 2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            zynpu = 0. - 2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            znnpu = zxnpu*zxnpu + zynpu*zynpu

            zlam = glamv(ji,jj)     ! north pole direction & modulous (at v-point)
            zphi = gphiv(ji,jj)
            zxnpv = 0. - 2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            zynpv = 0. - 2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            znnpv = zxnpv*zxnpv + zynpv*zynpv

            zlam = glamf(ji,jj)     ! north pole direction & modulous (at f-point)
            zphi = gphif(ji,jj)
            zxnpf = 0. - 2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            zynpf = 0. - 2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )
            znnpf = zxnpf*zxnpf + zynpf*zynpf
            !
            zlam = glamv(ji,jj  )   ! j-direction: v-point segment direction (around t-point)
            zphi = gphiv(ji,jj  )
            zlan = glamv(ji,jj-1)
            zphh = gphiv(ji,jj-1)
            zxvvt =  2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * COS( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            zyvvt =  2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * SIN( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            znvvt = SQRT( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
            znvvt = MAX( znvvt, 1.e-14 )

            zlam = glamf(ji,jj  )   ! j-direction: f-point segment direction (around u-point)
            zphi = gphif(ji,jj  )
            zlan = glamf(ji,jj-1)
            zphh = gphif(ji,jj-1)
            zxffu =  2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * COS( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            zyffu =  2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * SIN( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            znffu = SQRT( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            znffu = MAX( znffu, 1.e-14 )

            zlam = glamf(ji  ,jj)   ! i-direction: f-point segment direction (around v-point)
            zphi = gphif(ji  ,jj)
            zlan = glamf(ji-1,jj)
            zphh = gphif(ji-1,jj)
            zxffv =  2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * COS( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            zyffv =  2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * SIN( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            znffv = SQRT( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            znffv = MAX( znffv, 1.e-14 )

            zlam = glamu(ji,jj+1)   ! j-direction: u-point segment direction (around f-point)
            zphi = gphiu(ji,jj+1)
            zlan = glamu(ji,jj  )
            zphh = gphiu(ji,jj  )
            zxuuf =  2. * COS( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * COS( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            zyuuf =  2. * SIN( rad0*zlam ) * TAN( rPi0/4. - rad0*zphi/2. )   &
               &  -  2. * SIN( rad0*zlan ) * TAN( rPi0/4. - rad0*zphh/2. )
            znuuf = SQRT( znnpf * ( zxuuf*zxuuf + zyuuf*zyuuf )  )
            znuuf = MAX( znuuf, 1.e-14 )

            ! cosinus and sinus using dot and cross products
            !
            gsint(ji,jj) = ( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt
            gcost(ji,jj) = ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt
            !
            gsinu(ji,jj) = ( zxnpu*zyffu - zynpu*zxffu ) / znffu
            gcosu(ji,jj) = ( zxnpu*zxffu + zynpu*zyffu ) / znffu
            !
            gsinf(ji,jj) = ( zxnpf*zyuuf - zynpf*zxuuf ) / znuuf
            gcosf(ji,jj) = ( zxnpf*zxuuf + zynpf*zyuuf ) / znuuf
            !
            gsinv(ji,jj) = ( zxnpv*zxffv + zynpv*zyffv ) / znffv
            gcosv(ji,jj) =-( zxnpv*zyffv - zynpv*zxffv ) / znffv     ! (caution, rotation of 90 degres)



         END DO
      END DO


      DO jj = 2, ny-1
         DO ji = 2, nx
            IF( MOD( ABS( glamv(ji,jj) - glamv(ji,jj-1) ), 360. ) < 1.e-8 ) THEN
               gsint(ji,jj) = 0.
               gcost(ji,jj) = 1.
            ENDIF
            IF( MOD( ABS( glamf(ji,jj) - glamf(ji,jj-1) ), 360. ) < 1.e-8 ) THEN
               gsinu(ji,jj) = 0.
               gcosu(ji,jj) = 1.
            ENDIF
            IF(      ABS( gphif(ji,jj) - gphif(ji-1,jj) )         < 1.e-8 ) THEN
               gsinv(ji,jj) = 0.
               gcosv(ji,jj) = 1.
            ENDIF
            IF( MOD( ABS( glamu(ji,jj) - glamu(ji,jj+1) ), 360. ) < 1.e-8 ) THEN
               gsinf(ji,jj) = 0.
               gcosf(ji,jj) = 1.
            ENDIF
         END DO
      END DO

      !! If NEMO grid (nperio > 0):
      IF ( nperio > 0 ) THEN
         !
         PRINT *, '  *** ANGLE (mod_nemotools.f90) => using LBC_LNK!'
         !           ! lateral boundary cond.: T-, U-, V-, F-pts, sgn
         CALL lbc_lnk( nperio, gcost, 'T', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsint, 'T', -1.0_8 )
         CALL lbc_lnk( nperio, gcosu, 'U', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsinu, 'U', -1.0_8 )
         CALL lbc_lnk( nperio, gcosv, 'V', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsinv, 'V', -1.0_8 )
         CALL lbc_lnk( nperio, gcosf, 'F', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsinf, 'F', -1.0_8 )
      END IF

   END SUBROUTINE angle


   !!======================================================================
END MODULE mod_nemotools
