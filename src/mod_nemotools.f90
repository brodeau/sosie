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
   USE mod_manip, ONLY : extra_1_east, extra_1_west
   
   IMPLICIT NONE
   PRIVATE


   INTERFACE lbc_lnk
      MODULE PROCEDURE lbc_lnk_2d_r8, lbc_lnk_2d_r4
   END INTERFACE lbc_lnk




   PUBLIC   lbc_lnk, angle2



   REAL(8), PARAMETER, PUBLIC :: &
      &       rpi = 3.141592653, &
      &       rad = rpi/180.0,   &
      &       rtiny8 = 1.E-8

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



   SUBROUTINE angle( nperio, plamt, pphit, plamu, pphiu, plamv, pphiv, plamf, pphif, &
      &                      gcost, gsint, gcosu, gsinu, gcosv, gsinv, gcosf, gsinf   )
      !!
      !!              ***************************************
      !!              *  Taken from "OCE/SBC/geo2ocean.F90" *
      !!              ***************************************
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE angle  ***
      !!
      !! ** Purpose :   Compute angles between model grid lines and the North direction
      !!
      !! ** Method  :   sinus and cosinus of the angle between the north-south axe
      !!              and the j-direction at t, u, v and f-points
      !!                dot and cross products are used to obtain cos and sin, resp.
      !!
      !! ** Action  : - gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! INPUT :
      !! -------
      INTEGER                , INTENT(in) :: nperio
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamt, pphit
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamu, pphiu
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamv, pphiv
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamf, pphif

      !! OUTPUT :
      !! --------
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcost, gsint
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcosu, gsinu
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcosv, gsinv
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcosf, gsinf

      !! LOCAL :
      !! -------
      INTEGER :: nx, ny, ji, jj
      REAL(8) :: zlam, zphi
      !!
      REAL(8) :: &
         &     zxnpt, znnpt, zxnpu, znnpu, zxnpv, znnpv, zlan, &
         &     zxvvt, zyvvt, znvvt, znuuf, znnpf, znffv, znffu, &
         &     zynpu, zynpv, zxffu, zyffu, zyffv, zxuuf, zyuuf, zxffv, &
         &     zphh, zynpt, zxnpf, zynpf

      !! 2D domain shape:
      nx = SIZE(plamt,1)
      ny = SIZE(plamt,2)

      ! ============================= !
      ! Compute the cosinus and sinus !
      ! ============================= !
      ! (computation done on the north stereographic polar plane)
      !
      DO jj = 2, ny-1
         DO ji = 2, nx
            !
            zlam = plamt(ji,jj)     ! north pole direction & modulous (at t-point)
            zphi = pphit(ji,jj)
            zxnpt = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpt = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpt = zxnpt*zxnpt + zynpt*zynpt
            !
            zlam = plamu(ji,jj)     ! north pole direction & modulous (at u-point)
            zphi = pphiu(ji,jj)
            zxnpu = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpu = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpu = zxnpu*zxnpu + zynpu*zynpu
            !
            zlam = plamv(ji,jj)     ! north pole direction & modulous (at v-point)
            zphi = pphiv(ji,jj)
            zxnpv = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpv = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpv = zxnpv*zxnpv + zynpv*zynpv
            !
            zlam = plamf(ji,jj)     ! north pole direction & modulous (at f-point)
            zphi = pphif(ji,jj)
            zxnpf = 0. - 2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            zynpf = 0. - 2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )
            znnpf = zxnpf*zxnpf + zynpf*zynpf
            !
            zlam = plamv(ji,jj  )   ! j-direction: v-point segment direction (around t-point)
            zphi = pphiv(ji,jj  )
            zlan = plamv(ji,jj-1)
            zphh = pphiv(ji,jj-1)
            zxvvt =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyvvt =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znvvt = SQRT( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
            znvvt = MAX( znvvt, 1.e-14 )
            !
            zlam = plamf(ji,jj  )   ! j-direction: f-point segment direction (around u-point)
            zphi = pphif(ji,jj  )
            zlan = plamf(ji,jj-1)
            zphh = pphif(ji,jj-1)
            zxffu =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyffu =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znffu = SQRT( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            znffu = MAX( znffu, 1.e-14 )
            !
            zlam = plamf(ji  ,jj)   ! i-direction: f-point segment direction (around v-point)
            zphi = pphif(ji  ,jj)
            zlan = plamf(ji-1,jj)
            zphh = pphif(ji-1,jj)
            zxffv =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyffv =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znffv = SQRT( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            znffv = MAX( znffv, 1.e-14 )
            !
            zlam = plamu(ji,jj+1)   ! j-direction: u-point segment direction (around f-point)
            zphi = pphiu(ji,jj+1)
            zlan = plamu(ji,jj  )
            zphh = pphiu(ji,jj  )
            zxuuf =  2. * COS( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            zyuuf =  2. * SIN( rad*zlam ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlan ) * TAN( rpi/4. - rad*zphh/2. )
            znuuf = SQRT( znnpf * ( zxuuf*zxuuf + zyuuf*zyuuf )  )
            znuuf = MAX( znuuf, 1.e-14 )
            !
            !                       ! cosinus and sinus using dot and cross products
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
            !
         END DO
      END DO

      ! =============== !
      ! Geographic mesh !
      ! =============== !

      DO jj = 2, ny-1
         DO ji = 2, nx
            IF( MOD( ABS( plamv(ji,jj) - plamv(ji,jj-1) ), 360. ) < rtiny8 ) THEN
               gsint(ji,jj) = 0.
               gcost(ji,jj) = 1.
            ENDIF
            IF( MOD( ABS( plamf(ji,jj) - plamf(ji,jj-1) ), 360. ) < rtiny8 ) THEN
               gsinu(ji,jj) = 0.
               gcosu(ji,jj) = 1.
            ENDIF
            IF(      ABS( pphif(ji,jj) - pphif(ji-1,jj) )         < rtiny8 ) THEN
               gsinv(ji,jj) = 0.
               gcosv(ji,jj) = 1.
            ENDIF
            IF( MOD( ABS( plamu(ji,jj) - plamu(ji,jj+1) ), 360. ) < rtiny8 ) THEN
               gsinf(ji,jj) = 0.
               gcosf(ji,jj) = 1.
            ENDIF
         END DO
      END DO

      ! =========================== !
      ! Lateral boundary conditions !
      ! =========================== !

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


   

   SUBROUTINE angle2( nperio, plamt, pphit, plamu, pphiu, plamv, pphiv, plamf, pphif, &
      &                       gcost, gsint, gcosu, gsinu, gcosv, gsinv, gcosf, gsinf   )
      !!
      !!              ***************************************
      !!              *  Taken from "OCE/SBC/geo2ocean.F90" *
      !!              ***************************************
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE angle2  ***
      !!
      !! ** Purpose :   Compute angles between model grid lines and the North direction
      !!
      !! ** Method  :   sinus and cosinus of the angle2 between the north-south axe
      !!              and the j-direction at t, u, v and f-points
      !!                dot and cross products are used to obtain cos and sin, resp.
      !!
      !! ** Action  : - gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! INPUT :
      !! -------
      INTEGER                , INTENT(in) :: nperio
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamt, pphit
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamu, pphiu
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamv, pphiv
      REAL(8), DIMENSION(:,:), INTENT(in) :: plamf, pphif

      !! OUTPUT :
      !! --------
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcost, gsint
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcosu, gsinu
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcosv, gsinv
      REAL(8), DIMENSION(:,:), INTENT(out) :: gcosf, gsinf

      !! LOCAL :
      !! -------
      LOGICAL :: l_cut_180
      INTEGER :: nx, ny, ji, jj
      REAL(8) :: zphi, zlamt, zlamu, zlamv, zlamf, &
         &       zlamvjm1, zlamfjm1, zlamfim1, zlamujp1

      !!
      REAL(8) :: &
         &     zxnpt, znnpt, zxnpu, znnpu, zxnpv, znnpv, zlan, &
         &     zxvvt, zyvvt, znvvt, znuuf, znnpf, znffv, znffu, &
         &     zynpu, zynpv, zxffu, zyffu, zyffv, zxuuf, zyuuf, zxffv, &
         &     zphh, zynpt, zxnpf, zynpf

      !! 2D domain shape:
      nx = SIZE(plamt,1)
      ny = SIZE(plamt,2)
      
      l_cut_180 = (      (MAXVAL(plamt)<= 180.).AND.(MAXVAL(plamt)> 175.) &
         &          .AND.(MINVAL(plamt)>=-180.).AND.(MINVAL(plamt)<-175.) )


      PRINT *, 'LOLO AAA: l_cut_180 =', l_cut_180
      
      ! ============================= !
      ! Compute the cosinus and sinus !
      ! ============================= !
      ! (computation done on the north stereographic polar plane)
      !
      DO jj = 2, ny-1
         DO ji = 2, nx
            !
            ! Longitudes, must mind if near the 180 -> -180 cut line, then must work in the 0:360 frame:            
            zlamt = plamt(ji,jj)     ! north pole direction & modulous (at t-point)
            zlamu = plamu(ji,jj)     ! north pole direction & modulous (at u-point)
            zlamv = plamv(ji,jj)     ! north pole direction & modulous (at v-point)
            zlamf = plamf(ji,jj)     ! north pole direction & modulous (at f-point)
            !!
            zlamvjm1 = plamv(ji,jj-1)
            zlamfjm1 = plamf(ji,jj-1)
            zlamfim1 = plamf(ji-1,jj)
            zlamujp1 = plamu(ji,jj+1)

            IF( l_cut_180 ) THEN
               !! Are we close to this cut line?
               IF( (ABS(zlamt)>170.).OR.(ABS(zlamu)>170.).OR.(ABS(zlamv)>170.).OR.(ABS(zlamf)>170.) ) THEN
                  !! => move to the 0:360 frame:
                  zlamt = MOD( 360. + plamt(ji,jj) , 360. )
                  zlamu = MOD( 360. + plamu(ji,jj) , 360. )
                  zlamv = MOD( 360. + plamv(ji,jj) , 360. )
                  zlamf = MOD( 360. + plamf(ji,jj) , 360. )
                  !!
                  zlamvjm1 = MOD( 360. + plamv(ji,jj-1) , 360. )
                  zlamfjm1 = MOD( 360. + plamf(ji,jj-1) , 360. )
                  zlamfim1 = MOD( 360. + plamf(ji-1,jj) , 360. )
                  zlamujp1 = MOD( 360. + plamu(ji,jj+1) , 360. )
               END IF
            END IF
            
            zphi  = pphit(ji,jj)
            zxnpt = 0. - 2. * COS( rad*zlamt ) * TAN( rpi/4. - rad*zphi/2. )
            zynpt = 0. - 2. * SIN( rad*zlamt ) * TAN( rpi/4. - rad*zphi/2. )
            znnpt = zxnpt*zxnpt + zynpt*zynpt

            zphi  = pphiu(ji,jj)
            zxnpu = 0. - 2. * COS( rad*zlamu ) * TAN( rpi/4. - rad*zphi/2. )
            zynpu = 0. - 2. * SIN( rad*zlamu ) * TAN( rpi/4. - rad*zphi/2. )
            znnpu = zxnpu*zxnpu + zynpu*zynpu

            zphi  = pphiv(ji,jj)
            zxnpv = 0. - 2. * COS( rad*zlamv ) * TAN( rpi/4. - rad*zphi/2. )
            zynpv = 0. - 2. * SIN( rad*zlamv ) * TAN( rpi/4. - rad*zphi/2. )
            znnpv = zxnpv*zxnpv + zynpv*zynpv

            zphi  = pphif(ji,jj)
            zxnpf = 0. - 2. * COS( rad*zlamf ) * TAN( rpi/4. - rad*zphi/2. )
            zynpf = 0. - 2. * SIN( rad*zlamf ) * TAN( rpi/4. - rad*zphi/2. )
            znnpf = zxnpf*zxnpf + zynpf*zynpf

            zphi  = pphiv(ji,jj  )
            zphh  = pphiv(ji,jj-1)
            zxvvt =  2. * COS( rad*zlamv ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlamvjm1 ) * TAN( rpi/4. - rad*zphh/2. )
            zyvvt =  2. * SIN( rad*zlamv ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlamvjm1 ) * TAN( rpi/4. - rad*zphh/2. )
            znvvt = SQRT( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
            znvvt = MAX( znvvt, 1.e-14 )

            zphi  = pphif(ji,jj  )
            zphh  = pphif(ji,jj-1)
            zxffu =  2. * COS( rad*zlamf ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlamfjm1 ) * TAN( rpi/4. - rad*zphh/2. )
            zyffu =  2. * SIN( rad*zlamf ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlamfjm1 ) * TAN( rpi/4. - rad*zphh/2. )
            znffu = SQRT( znnpu * ( zxffu*zxffu + zyffu*zyffu )  )
            znffu = MAX( znffu, 1.e-14 )
            
            zphi  = pphif(ji  ,jj)
            zphh  = pphif(ji-1,jj)
            zxffv =  2. * COS( rad*zlamf ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlamfim1 ) * TAN( rpi/4. - rad*zphh/2. )
            zyffv =  2. * SIN( rad*zlamf ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlamfim1 ) * TAN( rpi/4. - rad*zphh/2. )
            znffv = SQRT( znnpv * ( zxffv*zxffv + zyffv*zyffv )  )
            znffv = MAX( znffv, 1.e-14 )

            zphi    = pphiu(ji,jj+1)
            zphh    = pphiu(ji,jj  )
            zxuuf =  2. * COS( rad*zlamujp1 ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * COS( rad*zlamu ) * TAN( rpi/4. - rad*zphh/2. )
            zyuuf =  2. * SIN( rad*zlamujp1 ) * TAN( rpi/4. - rad*zphi/2. )   &
               &  -  2. * SIN( rad*zlamu ) * TAN( rpi/4. - rad*zphh/2. )
            znuuf = SQRT( znnpf * ( zxuuf*zxuuf + zyuuf*zyuuf )  )
            znuuf = MAX( znuuf, 1.e-14 )
            !
            !                       ! cosinus and sinus using dot and cross products
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

      ! =========================== !
      ! Lateral boundary conditions !
      ! =========================== !

      !! If NEMO grid (nperio > 0):
      IF ( nperio > 0 ) THEN
         !
         PRINT *, '  *** ANGLE2 (mod_nemotools.f90) => using LBC_LNK!'
         !           ! lateral boundary cond.: T-, U-, V-, F-pts, sgn
         CALL lbc_lnk( nperio, gcost, 'T', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsint, 'T', -1.0_8 )
         CALL lbc_lnk( nperio, gcosu, 'U', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsinu, 'U', -1.0_8 )
         CALL lbc_lnk( nperio, gcosv, 'V', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsinv, 'V', -1.0_8 )
         CALL lbc_lnk( nperio, gcosf, 'F', -1.0_8 )   ;   CALL lbc_lnk( nperio, gsinf, 'F', -1.0_8 )

      ELSE

         !DO jj = 2, ny-1

         ! Northernmost row:
         DO ji = 1, nx
            CALL extra_1_east(pphit(ji,ny-3), pphit(ji,ny-2), pphit(ji,ny-1), pphit(ji,ny), gcost(ji,ny-3), gcost(ji,ny-2), gcost(ji,ny-1),  gcost(ji,ny))
            CALL extra_1_east(pphit(ji,ny-3), pphit(ji,ny-2), pphit(ji,ny-1), pphit(ji,ny), gsint(ji,ny-3), gsint(ji,ny-2), gsint(ji,ny-1),  gsint(ji,ny))
            CALL extra_1_east(pphiu(ji,ny-3), pphiu(ji,ny-2), pphiu(ji,ny-1), pphiu(ji,ny), gcosu(ji,ny-3), gcosu(ji,ny-2), gcosu(ji,ny-1),  gcosu(ji,ny))
            CALL extra_1_east(pphiu(ji,ny-3), pphiu(ji,ny-2), pphiu(ji,ny-1), pphiu(ji,ny), gsinu(ji,ny-3), gsinu(ji,ny-2), gsinu(ji,ny-1),  gsinu(ji,ny))
            CALL extra_1_east(pphiv(ji,ny-3), pphiv(ji,ny-2), pphiv(ji,ny-1), pphiv(ji,ny), gcosv(ji,ny-3), gcosv(ji,ny-2), gcosv(ji,ny-1),  gcosv(ji,ny))
            CALL extra_1_east(pphiv(ji,ny-3), pphiv(ji,ny-2), pphiv(ji,ny-1), pphiv(ji,ny), gsinv(ji,ny-3), gsinv(ji,ny-2), gsinv(ji,ny-1),  gsinv(ji,ny))
            CALL extra_1_east(pphif(ji,ny-3), pphif(ji,ny-2), pphif(ji,ny-1), pphif(ji,ny), gcosf(ji,ny-3), gcosf(ji,ny-2), gcosf(ji,ny-1),  gcosf(ji,ny))
            CALL extra_1_east(pphif(ji,ny-3), pphif(ji,ny-2), pphif(ji,ny-1), pphif(ji,ny), gsinf(ji,ny-3), gsinf(ji,ny-2), gsinf(ji,ny-1),  gsinf(ji,ny))
         END DO         
         
         ! Souhernmost row:
         DO ji = 1, nx
            CALL extra_1_west(pphit(ji,4), pphit(ji,3), pphit(ji,2), pphit(ji,1), gcost(ji,4), gcost(ji,3), gcost(ji,2),  gcost(ji,1))
            CALL extra_1_west(pphit(ji,4), pphit(ji,3), pphit(ji,2), pphit(ji,1), gsint(ji,4), gsint(ji,3), gsint(ji,2),  gsint(ji,1)) 
            CALL extra_1_west(pphiu(ji,4), pphiu(ji,3), pphiu(ji,2), pphiu(ji,1), gcosu(ji,4), gcosu(ji,3), gcosu(ji,2),  gcosu(ji,1))
            CALL extra_1_west(pphiu(ji,4), pphiu(ji,3), pphiu(ji,2), pphiu(ji,1), gsinu(ji,4), gsinu(ji,3), gsinu(ji,2),  gsinu(ji,1)) 
            CALL extra_1_west(pphiv(ji,4), pphiv(ji,3), pphiv(ji,2), pphiv(ji,1), gcosv(ji,4), gcosv(ji,3), gcosv(ji,2),  gcosv(ji,1))
            CALL extra_1_west(pphiv(ji,4), pphiv(ji,3), pphiv(ji,2), pphiv(ji,1), gsinv(ji,4), gsinv(ji,3), gsinv(ji,2),  gsinv(ji,1)) 
            CALL extra_1_west(pphif(ji,4), pphif(ji,3), pphif(ji,2), pphif(ji,1), gcosf(ji,4), gcosf(ji,3), gcosf(ji,2),  gcosf(ji,1))
            CALL extra_1_west(pphif(ji,4), pphif(ji,3), pphif(ji,2), pphif(ji,1), gsinf(ji,4), gsinf(ji,3), gsinf(ji,2),  gsinf(ji,1)) 
         END DO
         
         ! Westernmost column:
         DO jj = 1, ny
            CALL extra_1_west(plamt(4,jj), plamt(3,jj), plamt(2,jj), plamt(1,jj), gcost(4,jj), gcost(3,jj), gcost(2,jj),  gcost(1,jj))
            CALL extra_1_west(plamt(4,jj), plamt(3,jj), plamt(2,jj), plamt(1,jj), gsint(4,jj), gsint(3,jj), gsint(2,jj),  gsint(1,jj))
            CALL extra_1_west(plamu(4,jj), plamu(3,jj), plamu(2,jj), plamu(1,jj), gcosu(4,jj), gcosu(3,jj), gcosu(2,jj),  gcosu(1,jj))
            CALL extra_1_west(plamu(4,jj), plamu(3,jj), plamu(2,jj), plamu(1,jj), gsinu(4,jj), gsinu(3,jj), gsinu(2,jj),  gsinu(1,jj))
            CALL extra_1_west(plamv(4,jj), plamv(3,jj), plamv(2,jj), plamv(1,jj), gcosv(4,jj), gcosv(3,jj), gcosv(2,jj),  gcosv(1,jj))
            CALL extra_1_west(plamv(4,jj), plamv(3,jj), plamv(2,jj), plamv(1,jj), gsinv(4,jj), gsinv(3,jj), gsinv(2,jj),  gsinv(1,jj))
            CALL extra_1_west(plamf(4,jj), plamf(3,jj), plamf(2,jj), plamf(1,jj), gcosf(4,jj), gcosf(3,jj), gcosf(2,jj),  gcosf(1,jj))
            CALL extra_1_west(plamf(4,jj), plamf(3,jj), plamf(2,jj), plamf(1,jj), gsinf(4,jj), gsinf(3,jj), gsinf(2,jj),  gsinf(1,jj))
         END DO

      END IF
         
      
   END SUBROUTINE angle2        !



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


   
   

   !!======================================================================
END MODULE mod_nemotools
