MODULE MOD_BILIN_2D
   !!-----------------------------------------------------------------
   !!
   !!         Method of interpolation : "bilinear 2D"
   !!         =======================================
   !!
   !!    L. Brodeau, February 2018
   !!    P. Mathiot, August 2010
   !!    L. BRODEAU, fall 2008
   !!    J.-M. MOLINES, P. MATHIOT , 2007
   !!
   !!   Naming conventions of variables and arrays in this module:
   !!   `_s` : "s" for "source domain" (global: `_src`)
   !!   `_t` : "t" for "target domain" (global: `_trg`)
   !!   `_e` : "e" for "extended in space" (frame extension of 2 points)
   !!
   !!-----------------------------------------------------------------
   USE mod_conf
   USE mod_manip, ONLY: EXTEND_ARRAY_2D_COOR, EXTEND_ARRAY_2D_DATA, DEGE_TO_DEGWE, FIND_NEAREST_POINT, DISTANCE
   USE mod_poly
   USE mod_grids, ONLY: mask_ignore_trg
   USE io_ezcdf,  ONLY: DUMP_FIELD ; !debug

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BILIN_2D

   INTEGER, PARAMETER  :: np_e = 4  !: Source grid spatial extension (4 => extra frame of 2 points)
   INTEGER, PARAMETER :: &
                                !! Only active if iverbose==2:
      &   idb = 0, & ! i-index of point to debug on target domain
      &   jdb = 0    ! j-index of point to debug on target domain

   !! Mapping for bilin:
   TYPE :: bln_map
      REAL(8)          :: ralfa
      REAL(8)          :: rbeta
      INTEGER          :: jip
      INTEGER          :: jjp
      INTEGER(1)       :: iqdrn
      INTEGER(2)       :: ipb ! ID of problem if any...
   END TYPE bln_map

   TYPE(bln_map), DIMENSION(:,:), ALLOCATABLE, PUBLIC, SAVE :: BILIN_MAP
   REAL(4),       DIMENSION(:,:), ALLOCATABLE, SAVE :: distance_to_np
   LOGICAL,                                    SAVE :: l_1st_call_bilin=.true.

   INTEGER,                                SAVE :: ni_s_e, nj_s_e ! shape of extended source arrays
   REAL(8), DIMENSION(:,:),   ALLOCATABLE, SAVE :: x_s_2d_e , y_s_2d_e, Z_s_2d_e ! extended coordinates


   CHARACTER(len=400), SAVE :: cf_wght_bilin
   LOGICAL,    PARAMETER    :: l_save_distance_to_np=.TRUE. !: for each point of target grid, shows the distance to the nearest point

   REAL(8), PARAMETER :: repsilon = 1.E-9

   !LOGICAL, SAVE :: l_last_y_row_missing = .FALSE.

   !! PUBLIC:

   LOGICAL,                                 PUBLIC, SAVE :: l_skip_bilin_mapping=.false.


CONTAINS


   SUBROUTINE BILIN_2D_INIT( kewper, pX1, pY1, pX2, pY2 )
      !!==============================================================================
      !!  Input :
      !!             kewper : east-west periodicity
      !!             pX1    : 2D source longitude array of shape (ni,nj)
      !!             pY1    : 2D source latitude  array of shape (ni,nj)
      !!             pX2    : 2D target longitude array (ni,nj) or (ni,1)
      !!             pY2    : 2D target latitude  array (ni,nj) or (nj,1)
      !!==============================================================================
      INTEGER,                 INTENT(in) :: kewper
      REAL(8), DIMENSION(:,:), INTENT(in) :: pX1, pY1 ! source 2D arrays of longitude and latitude
      REAL(8), DIMENSION(:,:), INTENT(in) :: pX2, pY2 ! target   "          "         "
      !!==============================================================================
      !!
      LOGICAL :: lefw
      !!==============================================================================
      !!
      WRITE(6,*) ''; WRITE(6,*) ''
      WRITE(6,*) '###################################################################'
      WRITE(6,*) '#                     BILIN INITIALIZATION'
      WRITE(6,*) '###################################################################'
      WRITE(6,*) ''

      IF( (SIZE(pX2,1)/=ni_trg).OR.(SIZE(pX2,2)/=nj_trg) ) &
         & CALL STOP_THIS( '[BILIN_2D_INIT()] => shape inconsistency on target domain!' )

      WRITE(6,'("  * Allocating array `BILIN_MAP`: ",i4," x ",i4)') ni_trg, nj_trg
      ALLOCATE ( BILIN_MAP(ni_trg,nj_trg) )
      IF (l_save_distance_to_np) THEN
         WRITE(6,'("   * Allocating array `distance_to_np`: ",i4," x ",i4)') ni_trg, nj_trg
         ALLOCATE ( distance_to_np(ni_trg,nj_trg) )
      END IF
      WRITE(6,*) '  * Allocations done...'
      WRITE(6,*) ''
      WRITE(cf_wght_bilin,'("sosie_mapping_",a,".nc")') TRIM(cpat)
      WRITE(6,*) '  * Mapping file is "',TRIM(cf_wght_bilin),'" !'

      ni_s_e = ni_src+np_e
      nj_s_e = nj_src+np_e

      WRITE(6,'("  * Allocating space-extended source coordinates arrays: ",i4.4,"x",i4.4)')  ni_s_e, nj_s_e
      ALLOCATE( x_s_2d_e(ni_s_e,nj_s_e) , y_s_2d_e(ni_s_e,nj_s_e) )
      WRITE(6,'("  * Allocating space-extended source field array: ",i4.4,"x",i4.4)')  ni_s_e, nj_s_e
      ALLOCATE( Z_s_2d_e(ni_s_e,nj_s_e) )
      !!
      WRITE(6,'("  * Filling space-extended source coordinates arrays")')

      PRINT *, 'lolo: i_orca_src =', i_orca_src
      CALL EXTEND_ARRAY_2D_COOR( kewper, pX1, pY1, x_s_2d_e, y_s_2d_e, is_orca_grid=i_orca_src )
      
      INQUIRE(FILE=cf_wght_bilin, EXIST=lefw )
      IF ( lefw ) THEN
         WRITE(6,*) '    => it was found in current directory !'
         WRITE(6,*) '      => still! Make sure that this is really the one you need...'
         CALL RD_MAPPING_AB(cf_wght_bilin, BILIN_MAP(:,:))
         WRITE(6,*) '    => mapping and weights read into ', TRIM(cf_wght_bilin), ' !'
         WRITE(6,*) '    => will skip routine MAPPING_BL !'
         l_skip_bilin_mapping = .TRUE.
         !!
      ELSE
         WRITE(6,*) '    => Not found in the current directory !'
         WRITE(6,*) '    => Going to create "', TRIM(cf_wght_bilin), '"'
         WRITE(6,*) '      => This might be time consuming if your grids are big,'
         WRITE(6,*) '         but it ONLY needs to be done once...'
         WRITE(6,*) '         Therefore, consider keeping this file for future'
         WRITE(6,*) '         interpolations using the same "source-target" setup...'
         l_skip_bilin_mapping = .FALSE.
         !!
         WRITE(6,*) '  * Calling `MAPPING_BL()` to fill mapping `BILIN_MAP` once for all!'
         !!LOLO: using extended arrays and not pX1, pY1:
         !LOLO:CALL MAPPING_BL( kewper, x_s_2d_e, y_s_2d_e, pX2, pY2, BILIN_MAP,  pmsk_dom_trg=mask_ignore_trg )
         CALL MAPPING_BL( kewper, x_s_2d_e, y_s_2d_e, pX2, pY2, BILIN_MAP )
         CALL SAVE_BL_2D_MP( pX2, pY2, BILIN_MAP )  ! Saving mapping into netCDF file if relevant...!
      END IF
      WRITE(6,*) ''
      !!
      WRITE(6,*) '  Initializations for BILIN done !'
      WRITE(6,*) ''
      WRITE(6,*) '###################################################################'
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
   END SUBROUTINE BILIN_2D_INIT



   SUBROUTINE BILIN_2D( k_ew_per, pX1, pY1, pZ1, pX2, pY2, pZ2,  mask_domain_trg )
      !!================================================================
      !!
      !! INPUT :     k_ew_per : east-west periodicity
      !!                        k_ew_per = -1  --> no periodicity
      !!                        k_ew_per >= 0  --> periodicity with overlap of k_ew_per points
      !!             pX1  : 2D source longitude array of shape (ni,nj)
      !!             pY1  : 2D source latitude  array of shape (ni,nj)
      !!             pZ1  : source field on source grid  "    "
      !!
      !!             pX2  : 2D target longitude array (ni,nj) or (ni,1)
      !!             pY2  : 2D target latitude  array (ni,nj) or (nj,1)
      !!
      !! OUTPUT :
      !!             pZ2    : field extrapolated from source to target grid
      !!
      !!
      !! OPTIONAL IN:
      !!      * mask_domain_trg: ignore (dont't treat) regions of the target domain where mask_domain_trg==0 !
      !!
      !!================================================================

      !! Input/Output arguments
      INTEGER,                   INTENT(in)  :: k_ew_per
      REAL(8),   DIMENSION(:,:), INTENT(in)  :: pX1, pY1
      REAL(wpl), DIMENSION(:,:), INTENT(in)  :: pZ1
      REAL(8),   DIMENSION(:,:), INTENT(in)  :: pX2, pY2
      REAL(wpl), DIMENSION(:,:), INTENT(out) :: pZ2
      !INTEGER,                 INTENT(in)  :: ithrd ! # OMP thread
      INTEGER(1), OPTIONAL, DIMENSION(:,:), INTENT(in) :: mask_domain_trg
      !! Local variables
      INTEGER :: nx2, ny2, iqd, iP, jP
      REAL(8) :: alpha, beta, rmeanv
      INTEGER :: ji, jj
      !!================================================================

      nx2 = SIZE(pZ2,1)
      ny2 = SIZE(pZ2,2)

      !PRINT *, 'Debug: BILIN_2D() => shape of pX1 =', SIZE(pX1,1), SIZE(pX1,2)

      IF( l_1st_call_bilin )  CALL BILIN_2D_INIT( k_ew_per, pX1, pY1, pX2, pY2 )
      

      CALL EXTEND_ARRAY_2D_DATA( k_ew_per, x_s_2d_e, y_s_2d_e,   REAL(pZ1,8), Z_s_2d_e,     is_orca_grid=i_orca_src )
      !! => from now on, `x_s_2d_e,y_s_2d_e,Z_s_2d_e` should be used instead of `pX1,pY1,pZ1` !
      
      !STOP'LOLO: mod_bilin_2d.f90' !debug
      
      pZ2(:,:) = rmissval ! Flagging non-interpolated output points

      !PRINT *, 'DEBUG bilin_2d, shape BILIN_MAP =', SIZE(BILIN_MAP,1), SIZE(BILIN_MAP,2), '#', 1
      !PRINT *, 'DEBUG bilin_2d, shape mask_ignore_trg =', SIZE(mask_ignore_trg,1), SIZE(mask_ignore_trg,2), '#', 1

      DO ji=1, nx2
         WHERE ( BILIN_MAP(ji,:)%jip < 1 ) mask_ignore_trg(ji,:) = 0
         WHERE ( BILIN_MAP(ji,:)%jjp < 1 ) mask_ignore_trg(ji,:) = 0

         BILIN_MAP(ji,:)%jip = MAX( BILIN_MAP(ji,:)%jip , 1 )  ! so no i or j <= 0
         BILIN_MAP(ji,:)%jjp = MAX( BILIN_MAP(ji,:)%jjp , 1 )  ! so no i or j <= 0

      END DO

      DO jj=1, ny2
         DO ji=1, nx2
            iP    = BILIN_MAP(ji,jj)%jip
            jP    = BILIN_MAP(ji,jj)%jjp
            iqd = BILIN_MAP(ji,jj)%iqdrn
            alpha = BILIN_MAP(ji,jj)%ralfa
            beta  = BILIN_MAP(ji,jj)%rbeta
            !!
            IF ( (ABS(degE_to_degWE(x_s_2d_e(iP,jP))-degE_to_degWE(pX2(ji,jj)))<1.E-5) &
               & .AND. (ABS(y_s_2d_e(iP,jP)-pY2(ji,jj))<1.E-5) ) THEN
               !! COPY:
               IF (iverbose>0) WRITE(6,*) ' *** BILIN_2D: "identical point" detected (crit: 1.E-5) => copying value, no interpolation!'
               pZ2(ji,jj) = Z_s_2d_e(iP,jP)
            ELSE
               !! INTERPOLATION:
               !pZ2(ji,jj) = INTERP_BL(k_ew_per, iP, jP, iqd, alpha, beta, Z_s_2d_e)
               !pZ2(ji,jj) = INTERP_BL(k_ew_per, ithrd, ji, jj, Z_s_2d_e)

               CALL INTERP_BL( k_ew_per, ji, jj, Z_s_2d_e, pZ2(ji,jj) )

            END IF
         END DO
      END DO

      pZ2(:,:) = pZ2(:,:)*REAL(mask_ignore_trg(:,:), 4) + REAL(1-mask_ignore_trg(:,:), 4)*(-9995.) ! masking problem points as in mask_ignore_trg


      !IF ( l_1st_call_bilin(ithrd) ) THEN
      !   !! Is the very last Y row fully masked! lolo and on a ORCA grid!!!
      !   IF ( i_orca_trg >= 4 ) THEN
      !      PRINT *, i_orca_trg, 'LOLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
      !      rmeanv = SUM(pZ2(:,ny2))/nx2
      !      l_last_y_row_missing(ithrd) = ( (rmeanv < rmissval + 0.1).AND.(rmeanv > rmissval - 0.1) )
      !   END IF
      !END IF

      !WRITE(6,*) ' l_last_y_row_missing =>', l_last_y_row_missing
      !IF ( i_orca_trg == 4 ) WRITE(6,*) ' Target grid is an ORCA grid with north-pole T-point folding!'
      !IF ( i_orca_trg == 6 ) WRITE(6,*) ' Target grid is an ORCA grid with north-pole F-point folding!'

      !! Correcting last missing band if relevant: LOLO: should use lbc_lnk no ????
      !IF ( l_last_y_row_missing(ithrd) ) THEN
      !   IF ( i_orca_trg == 4 ) THEN
      !      pZ2(2:nx2/2           ,ny2)   = pZ2(nx2:nx2-nx2/2-2:-1,ny2-2)
      !      pZ2(nx2:nx2-nx2/2-2:-1,ny2)   = pZ2(2:nx2/2           ,ny2-2)
      !   END IF
      !   IF ( i_orca_trg == 6 ) THEN
      !      pZ2(2:nx2/2             ,ny2) = pZ2(nx2-1:nx2-nx2/2+1:-1,ny2-1)
      !      pZ2(nx2-1:nx2-nx2/2+1:-1,ny2) = pZ2(2:nx2/2             ,ny2-1)
      !   END IF
      !END IF

      l_1st_call_bilin = .FALSE.

   END SUBROUTINE BILIN_2D

   SUBROUTINE SAVE_BL_2D_MP( pX2, pY2, pbln_map )
      !!=========================================================================
      REAL(8),       DIMENSION(:,:), INTENT(in) :: pX2, pY2
      TYPE(bln_map), DIMENSION(:,:), INTENT(in) :: pbln_map
      !!=========================================================================
      WRITE(6,*) ''
      IF ( l_skip_bilin_mapping ) THEN
         WRITE(6,*) '  ==> NOT writing bilin mapping file because it was found...'
      ELSE
         WRITE(6,*) '  ==> writing bilin mapping into "', TRIM(cf_wght_bilin),'" !'
         CALL SAVE_BL_MAPPING_NC( cf_wght_bilin, pX2, pY2, pbln_map, rmissval )
      END IF
      !!
      IF (l_save_distance_to_np) DEALLOCATE( distance_to_np )
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
   END SUBROUTINE SAVE_BL_2D_MP


   !FUNCTION INTERP_BL(k_ew_per, kiP, kjP, kqd, pa, pb, Z_in)
   SUBROUTINE INTERP_BL(k_ew_per, ilt, jlt, Z_in, pres )

      INTEGER,                 INTENT(in) :: k_ew_per
      INTEGER,                 INTENT(in) :: ilt, jlt   ! LOCAL (local omp domain) coordinates of treated point on target domain
      !INTEGER,                 INTENT(in) :: kiP, kjP, kqd
      !REAL(8),                 INTENT(in) :: pa, pb
      REAL(8), DIMENSION(:,:), INTENT(in) :: Z_in
      REAL(4),                 INTENT(out) :: pres

      INTEGER :: kiP, kjP, kqd
      REAL(8) :: pa, pb

      !REAL(4) :: INTERP_BL
      REAL(4) ::  wup, w1, w2, w3, w4
      INTEGER :: ki, kj, nxi, nyi, kiPm1, kiPp1
      INTEGER :: i1=0, j1=0, i2=0, j2=0, i3=0, j3=0, i4=0, j4=0
      INTEGER :: Nitl

      !! Choose the 4 interpolation points, according to sector and nearest point (kiP, kjP)

      !!   o<--o        x<--o         o<--x         o<--o
      !! 1 |   ^ NE   2 |   ^ SE    3 |   ^ SW    4 |   ^ NW
      !!   v   |        v   |         v   |         v   |
      !!   x-->o        o-->o         o-->o         o-->x

      !! Source domain is used in its whole...
      nxi = SIZE(Z_in,1)
      nyi = SIZE(Z_in,2)

      ki = ilt
      kj = jlt

      !!
      Nitl = SIZE(BILIN_MAP(:,:)%jip,1)
      IF( ki > Nitl ) THEN
         PRINT *, 'PROBLEM: ki > Nitl ! ki, Nitl, ilt ', ki, Nitl, ilt
         STOP
      END IF

      kiP = BILIN_MAP(ki,kj)%jip
      kjP = BILIN_MAP(ki,kj)%jjp
      kqd = BILIN_MAP(ki,kj)%iqdrn
      pa  = BILIN_MAP(ki,kj)%ralfa
      pb  = BILIN_MAP(ki,kj)%rbeta

      kiPm1 = kiP-1
      kiPp1 = kiP+1
      IF ( (kiPm1 ==   0  ).AND.(k_ew_per>=0) )  kiPm1 = nxi - k_ew_per
      IF ( (kiPp1 == nxi+1).AND.(k_ew_per>=0) )  kiPp1 = 1   + k_ew_per


      SELECT CASE (kqd)

      CASE (1)  ! nearest point is the bottom left corner point of local mesh
         i1=kiP   ; j1 = kjP  ! local mesh is located NE of nearest point
         i2=kiPp1 ; j2 = kjP
         i3=kiPp1 ; j3 = kjP+1
         i4=kiP   ; j4 = kjP+1

      CASE (2)  ! nearest point is the top left corner point of mesh
         i1=kiP   ; j1 = kjP    ! local mesh is located SE of nearest point
         i2=kiP   ; j2 = kjP-1
         i3=kiPp1 ; j3 = kjP-1
         i4=kiPp1 ; j4 = kjP

      CASE (3)  ! nearest point is the top righ corner point of mesh
         i1=kiP   ; j1 = kjP   ! local mesh is located SW of nearest point
         i2=kiPm1 ; j2 = kjP
         i3=kiPm1 ; j3 = kjP-1
         i4=kiP   ; j4 = kjP-1

      CASE (4)  ! nearest point is the bottom right corner point of mesh
         i1=kiP   ; j1 = kjP  ! local mesh is located NW of nearest point
         i2=kiP   ; j2 = kjP+1
         i3=kiPm1 ; j3 = kjP+1
         i4=kiPm1 ; j4 = kjP

      END SELECT

      !! compute sum weight above target point
      w1=REAL( (1. - pa)*(1. - pb) , 4)
      w2=REAL(       pa *(1. - pb) , 4)
      w3=REAL(       pa * pb       , 4)
      w4=REAL( (1. - pa)* pb       , 4)

      wup = w1 + w2 + w3 + w4

      !IF ( (i1==0).OR.(j1==0).OR.(i2==0).OR.(j2==0).OR.(i3==0).OR.(j3==0).OR.(i4==0).OR.(j4==0) ) THEN
      !   WRITE(6,*) ' WARNING: INTERP_BL => at least one of the i,j index is zero!'
      !END IF

      ! interpolate with non-masked  values, above target point

      !PRINT *, 'LOLO: j1,j2,j3,j4=', j1,j2,j3,j4

      IF ( wup == 0. ) THEN
         pres = -9998.
      ELSEIF ( (i1<1).OR.(i2<1).OR.(i3<1).OR.(i4<1) ) THEN
         pres = -9997.
      ELSEIF ( (j1<1).OR.(j2<1).OR.(j3<1).OR.(j4<1) ) THEN
         pres = -9996.
      ELSEIF ( (j1>nyi).OR.(j2>nyi).OR.(j3>nyi).OR.(j4>nyi) ) THEN
         pres = -9995.
      ELSEIF ( (i1>nxi).OR.(i2>nxi).OR.(i3>nxi).OR.(i4>nxi) ) THEN
         pres = -9994.
      ELSE
         pres = REAL(  ( Z_in(i1,j1)*w1 + Z_in(i2,j2)*w2 + Z_in(i3,j3)*w3 + Z_in(i4,j4)*w4 )/wup  , 4 )
      ENDIF

   END SUBROUTINE INTERP_BL




   SUBROUTINE MAPPING_BL(k_ew_per, plon_s, plat_s, plon_t, plat_t, pbln_map,   pmsk_dom_trg)
      !!----------------------------------------------------------------------------
      !!            ***  SUBROUTINE MAPPING_BL  ***
      !!
      !!   ** Purpose:  Write file of position, weight need by interpolation
      !!   *  Extract of CDFTOOLS cdfweight.f90 writen by Jean Marc Molines
      !!
      !! OPTIONAL:
      !!      * pmsk_dom_trg: ignore (dont't treat) regions of the target domain where pmsk_dom_trg==0 !
      !!----------------------------------------------------------------------------
      !!
      INTEGER,                 INTENT(in) :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in) :: plon_s, plat_s
      REAL(8), DIMENSION(:,:), INTENT(in) :: plon_t, plat_t
      TYPE(bln_map), DIMENSION(:,:), INTENT(out) :: pbln_map
      !!
      INTEGER(1), OPTIONAL, DIMENSION(:,:), INTENT(in) :: pmsk_dom_trg

      INTEGER :: &
         &     iqd, iqd0, iqd_old, &
         &     iP, jP,             &
         &     nx_s, ny_s, nx_t, ny_t, &
         &     ji, jj,   &
         &     iPm1, iPp1,  &
         &     jPm1, jPp1,  &
         &     iproblem

      REAL(8) ::  &
         &  xP, yP, &
         &  hPp, &            !: local maximum metrics
         &  lonP, latP,     & !: coordinates of the nearest point  (NP)
         &  lonN, latN, hN, & !: grid point North of NP, true heding from NP
         &  lonE, latE, hE, & !: grid point East of NP, true heding from NP
         &  lonS, latS, hS, & !: grid point South of NP, true heding from NP
         &  lonW, latW, hW, & !: grid point West of NP, true heding from NP
         &  hP                !: true heading of target point from NP

      REAL(8), DIMENSION(0:4) :: &
         &    loni, lati    !: the 4 grid points around target (1-4) + the target (0)

      !! To save in the netcdf file:
      INTEGER(4),    DIMENSION(:,:), ALLOCATABLE :: ki_nrst, kj_nrst
      INTEGER(1),    DIMENSION(:,:), ALLOCATABLE :: kmsk_ignr_trg
      !!
      REAL(8) :: zalfa, zbeta
      LOGICAL :: l_ok, lagain, lpdebug
      INTEGER :: icpt, ithrd
      !!
      !!DEBUG:
      !INTEGER :: jx
      CHARACTER(len=80) :: cf_tmp !debug
      LOGICAL :: l_skip_P
      !!----------------------------------------------------------------------------

      ithrd = 1 ! no OpenMP !
      !IF( PRESENT(ithread) ) ithrd = ithread

      nx_s = SIZE(plon_s,1) ! so these are extended arrays size!!!
      ny_s = SIZE(plon_s,2)
      !WRITE(6,'("Debug/MAPPING_BL: shape of plon_s for thread #",i1,": ",i4.4,"x",i4.4)') ithrd, nx_s, ny_s

      nx_t = SIZE(plon_t,1)
      ny_t = SIZE(plon_t,2)

      !WRITE(6,'("Debug/MAPPING_BL: size of target domain for thread #",i1,": ",i4.4,"x",i4.4)') 0, nx_t, ny_t


      ALLOCATE ( kmsk_ignr_trg(nx_t,ny_t), ki_nrst(nx_t,ny_t), kj_nrst(nx_t,ny_t) )

      pbln_map(:,:)%jip    = 0
      pbln_map(:,:)%jjp    = 0
      pbln_map(:,:)%iqdrn  = 0
      pbln_map(:,:)%ralfa  = 0.
      pbln_map(:,:)%rbeta  = 0.
      pbln_map(:,:)%ipb    = 0
      kmsk_ignr_trg(:,:)   = 1
      ki_nrst(:,:)         = 0
      kj_nrst(:,:)         = 0

      IF ( PRESENT(pmsk_dom_trg) ) kmsk_ignr_trg(:,:) = pmsk_dom_trg(:,:)

      !! DEBUG: checking input fields for each different thread:
      !WRITE(cf_tmp,'("in_mbl_lon_src_",i2.2,".nc")') 0
      !CALL DUMP_FIELD(REAL(plon_s(:,:),4), cf_tmp, 'lon')
      !WRITE(cf_tmp,'("in_mbl_lat_src_",i2.2,".nc")') 0
      !CALL DUMP_FIELD(REAL(plat_s(:,:),4), cf_tmp, 'lat')
      !!
      !WRITE(cf_tmp,'("in_mbl_lon_trg_",i2.2,".nc")') 0
      !CALL DUMP_FIELD(REAL(plon_t,4), cf_tmp, 'lon')
      !WRITE(cf_tmp,'("in_mbl_lat_trg_",i2.2,".nc")') 0
      !CALL DUMP_FIELD(REAL(plat_t,4), cf_tmp, 'lat')
      
      CALL FIND_NEAREST_POINT( plon_t, plat_t, plon_s, plat_s, ki_nrst, kj_nrst ) !,  &
      !LOLO:&                     pmsk_dom_trg=kmsk_ignr_trg )
      
      !!LOLO: ok since we work with extended arrays points i=1 and i=nx_s can be avoided:
      IF( k_ew_per>=0 ) THEN
         WHERE( ki_nrst(:,:)==1    ) ki_nrst(:,:) = nx_s - (np_e - 1 + k_ew_per)
         WHERE( ki_nrst(:,:)==nx_s ) ki_nrst(:,:) =   1  + (np_e - 1 + k_ew_per)         
         !IF( iP==1      ) iP = nx_s - (np_e - 1 + k_ew_per) !3 lilo still have to use ewper (here it's 0! in the example I'm running!)
         !IF( iP==2      ) iP = nx_s - (np_e - 2 + k_ew_per) !2 lilo
         !IF( iP==nx_s   ) iP =   1  + (np_e - 1 + k_ew_per) !3
         !IF( iP==nx_s-1 ) iP =   1  + (np_e - 2 + k_ew_per) !2
      END IF
      
      !debug:
      !WRITE(cf_tmp,'("in_mbl_nrst_ji",i2.2,".nc")') ithread
      !WRITE(cf_tmp,'("in_mbl_nrst_ji",i2.2,".nc")') 0
      !PRINT *, 'Debug/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      !CALL DUMP_FIELD(REAL(ki_nrst,4), cf_tmp, 'ji')
      !
      !WRITE(cf_tmp,'("in_mbl_nrst_jj",i2.2,".nc")') ithread
      !WRITE(cf_tmp,'("in_mbl_nrst_jj",i2.2,".nc")') 0
      !PRINT *, 'Debug/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      !CALL DUMP_FIELD(REAL(kj_nrst,4), cf_tmp, 'jj')
      !debug.

      !$OMP BARRIER
      !STOP'lilo'
      ! jusqu'ici, ok!

      !CALL DUMP_FIELD(REAL(kmsk_ignr_trg,4), 'kmsk_ignr_trg.nc', 'msk') ;  STOP'kmsk_ignr_trg.nc !!!'

      !! PART I
      DO jj = 1, ny_t
         DO ji = 1, nx_t

            lpdebug = ( (iverbose==2).AND.(ji==idb).AND.(jj==jdb) )

            IF ( kmsk_ignr_trg(ji,jj)==1 ) THEN
               !! => exclude regions that do not exist on source domain (kmsk_ignr_trg==0) and
               !! points for which the nearest point was not found (kmsk_ignr_trg==-1 or -2)

               !! Now deal with horizontal interpolation
               !! set longitude of input point in accordance with lon ( [lon0, 360+lon0 [ )
               xP = plon_t(ji,jj)
               yP = plat_t(ji,jj)

               iP = ki_nrst(ji,jj)
               jP = kj_nrst(ji,jj)
               
               IF(lpdebug) WRITE(6,*) ' *** #DEBUG: xP, yP =', xP, yP
               IF(lpdebug) WRITE(6,*) ' *** #DEBUG: iP, jP, nx_s, ny_s =', iP, jP, nx_s, ny_s

               IF ( (iP/=INT(rmissval)).AND.(jP/=INT(rmissval)).AND.(jP<ny_s) ) THEN

                  CALL GIVE_NGHBR_POINTS_SRC( k_ew_per, iP, jP, nx_s, ny_s, plon_s, plat_s,  &
                     &                        l_skip_P, iPm1, iPp1, jPm1, jPp1, lonP, latP,    &
                     &                        lonN, latN, lonE, latE, lonS, latS, lonW, latW )

                  IF ( .NOT. l_skip_P ) THEN

                     !! Restore target point longitude between 0 and 360
                     xP = MOD(xP,360._8)

                     !!LB:
                     !! All this HEADING stuf is aimed at finding in which source grid cell
                     !! the target point is comprised with regards to the nearest point
                     !! (there is actually 4 possible adjacent cells NE, SE, SW and NW)

                     !! Compute heading of target point and neighbours from the nearest point
                     hP = heading(lonP,  xP,latP,  yP)  ! target point
                     hN = heading(lonP,lonN,latP,latN)  ! 'north' on the grid
                     hE = heading(lonP,lonE,latP,latE)  ! 'east' on the grid
                     hS = heading(lonP,lonS,latP,latS)  ! 'south' on the grid
                     hW = heading(lonP,lonW,latP,latW)  ! 'west' on the grid
                     !!  => returns an angle !



                     !! First attemp, and generally the good one to find iqd !
                     !! ********************************************************
                     !! Determine the sector in wich the target point is
                     !! located: ( from 1, to 4 resp. adjacent cells NE, SE, SW,
                     !! NW of the grid) which mesh from the nearest point
                     !! point-of-view !
                     !!
                     !!   o--o        x--o         o--x         o--o
                     !! 1 |  | NE   2 |  | SE    3 |  | SW    4 |  | NW
                     !!   x--o        o--o         o--o         o--x
                     !!
                     iqd = 4
                     !!
                     ! To avoid problem with the GW meridian, pass to -180, 180 when working around GW
                     hPp = degE_to_degWE(hP)
                     !!
                     IF ( hN > hE ) hN = hN -360._8
                     IF ( hPp > hN .AND. hPp <= hE ) iqd=1
                     IF ( hP > hE  .AND. hP <= hS )  iqd=2
                     IF ( hP > hS  .AND. hP <= hW )  iqd=3
                     IF ( hP > hW  .AND. hPp <= hN)  iqd=4

                     iqd0    = iqd
                     iqd_old = iqd
                     !IF(lpdebug) WRITE(6,*) ' *** #DEBUG: first find of iqd =', iqd0

                     loni(0) = xP ;    lati(0) = yP      ! fill loni, lati for 0 = target point
                     loni(1) = lonP ;  lati(1) = latP    !                     1 = nearest point

                     IF (l_save_distance_to_np) distance_to_np(ji,jj) = DISTANCE(xP, lonP, yP, latP)

                     !! Problem is that sometimes, in the case of really twisted
                     !! meshes this method screws up, iqd is not what it
                     !! should be, so need the following loop on the value of iqd:
                     icpt = 0
                     lagain = .TRUE.
                     DO WHILE ( lagain )

                        SELECT CASE ( iqd ) ! point 2 3 4 are counter clockwise in the respective sector
                        CASE ( 1 )
                           loni(2) = lonE ; lati(2) = latE
                           loni(3) = MOD(plon_s(iPp1,jPp1), 360._8) ; lati(3) = plat_s(iPp1,jPp1)
                           loni(4) = lonN ; lati(4) = latN
                        CASE ( 2 )
                           loni(2) = lonS ; lati(2) = latS
                           loni(3) = MOD(plon_s(iPp1,jPm1), 360._8) ; lati(3) = plat_s(iPp1,jPm1)
                           loni(4) = lonE ; lati(4) = latE
                        CASE ( 3 )
                           loni(2) = lonW ; lati(2) = latW
                           loni(3) = MOD(plon_s(iPm1,jPm1), 360._8) ; lati(3) = plat_s(iPm1,jPm1)
                           loni(4) = lonS ; lati(4) = latS
                        CASE ( 4 )
                           loni(2) = lonN ; lati(2) = latN
                           loni(3) = MOD(plon_s(iPm1,jPp1), 360._8) ; lati(3) = plat_s(iPm1,jPp1)
                           loni(4) = lonW ; lati(4) = latW
                        END SELECT

                        WHERE ( loni <= 0.0 )  loni = loni + 360._8  ! P. Mathiot: Some bug with ERA40 grid

                        !! The tests!!!
                        l_ok = L_InPoly ( loni(1:4), lati(1:4), xP, yP )    ! $$
                        IF(lpdebug) WRITE(6,*) ' *** #DEBUG: l_ok =', l_ok

                        IF ( (.NOT. l_ok).AND.(yP < 88.) ) THEN
                           !! Mhhh... Seems like the "heading()" approach
                           !! screwed up... i.e the point xP,yP is not into the
                           !! mesh corresponding to current iquadran, trying all
                           !! other adjacent meshes to find if it belongs to one
                           !! of them...
                           icpt  = icpt + 1
                           iqd = icpt
                        ELSE
                           !! Everything okay! Point is inside the the mesh
                           !! corresponding to current iquadran ! :D
                           lagain = .FALSE.
                           IF ( (icpt>0).AND.(lpdebug) ) WRITE(6,*) ' --- iquadran corrected thanks to iterative test (old,new) =>', iqd_old, iqd
                        END IF
                        IF ( icpt == 5 ) THEN
                           lagain = .FALSE. ! simply give up
                           iqd = iqd0 ! Giving what first method gave
                        END IF
                     END DO !DO WHILE ( lagain )

                     !! resolve a non linear system of equation for alpha and beta
                     !! ( the non dimensional coordinates of target point)
                     CALL LOCAL_COORD(loni, lati, zalfa, zbeta, iproblem)
                     pbln_map(ji,jj)%ipb = iproblem

                     !LOLO: mark this in ID_problem => case when iqd = iqd0 ( IF ( icpt == 5 ) ) above!!!
                     IF ( icpt == 5 ) pbln_map(ji,jj)%ipb = 2 ! IDing the screw-up from above...

                     IF (lpdebug) THEN
                        WRITE(6,*) ' *** #DEBUG :'
                        WRITE(6,*) 'Nearest point :',lonP,  latP,  hP, hPp
                        WRITE(6,*) 'North point :',  lonN , latN , hN
                        WRITE(6,*) 'East  point :',  lonE , latE , hE
                        WRITE(6,*) 'South point :',  lonS , latS , hS
                        WRITE(6,*) 'West  point :',  lonW , latW , hW
                        WRITE(6,*) 'iqd =',iqd
                        WRITE(6,*) ''
                        WRITE(6,*) ' Nearest 4 points :'
                        WRITE(6,*) 'Pt 1 :',loni(1), lati(1)
                        WRITE(6,*) 'Pt 2 :',loni(2), lati(2)
                        WRITE(6,*) 'Pt 3 :',loni(3), lati(3)
                        WRITE(6,*) 'Pt 4 :',loni(4), lati(4)
                        WRITE(6,*) ''
                     END IF

                     !! Saving into arrays to be written at the end:
                     pbln_map(ji,jj)%iqdrn = iqd
                     pbln_map(ji,jj)%ralfa = zalfa
                     pbln_map(ji,jj)%rbeta = zbeta
                     pbln_map(ji,jj)%ipb   = 0

                  END IF ! IF ((iPm1 < 1).OR.(jPm1 < 1).OR.(iPp1 > nx_s))
               END IF ! IF ( (iP/=INT(rmissval)).AND.(jP/=INT(rmissval)).AND.(jP<ny_s) )
            END IF ! IF ( kmsk_ignr_trg(ji,jj)==1 )
            
         ENDDO !DO ji = 1, nx_t
      ENDDO !DO jj = 1, ny_t

      !lolo: UGLY!
      pbln_map(:,:)%jip = ki_nrst(:,:)
      pbln_map(:,:)%jjp = kj_nrst(:,:)

      !CALL DUMP_FIELD(REAL(pbln_map(:,:)%ralfa,4), 'alpha_trg.nc', 'alpha'); STOP'alpha_trg.nc'
      
      WHERE ( (ki_nrst == INT(rmissval)).OR.(kj_nrst == INT(rmissval)) )
         pbln_map(:,:)%ralfa  = rmissval
         pbln_map(:,:)%rbeta  = rmissval
      END WHERE

      ! Awkwardly fixing problematic points but remembering them in ID_problem

      ! Negative values that are actually 0
      WHERE ( ((pbln_map(:,:)%ralfa < 0.).AND.(pbln_map(:,:)%ralfa > -repsilon)) ) pbln_map(:,:)%ralfa = 0.0
      WHERE ( ((pbln_map(:,:)%rbeta < 0.).AND.(pbln_map(:,:)%rbeta > -repsilon)) ) pbln_map(:,:)%rbeta = 0.0

      WHERE ( (pbln_map(:,:)%ralfa > rmissval).AND.(pbln_map(:,:)%ralfa < 0.) )
         pbln_map(:,:)%ralfa = 0.5
         pbln_map(:,:)%ipb = 4
      END WHERE
      WHERE ( pbln_map(:,:)%ralfa > 1. )
         pbln_map(:,:)%ralfa = 0.5
         pbln_map(:,:)%ipb =  5
      END WHERE

      WHERE ( (pbln_map(:,:)%rbeta > rmissval).AND.(pbln_map(:,:)%rbeta < 0.) )
         pbln_map(:,:)%rbeta = 0.5
         pbln_map(:,:)%ipb = 6
      END WHERE
      WHERE ( pbln_map(:,:)%rbeta > 1. )
         pbln_map(:,:)%rbeta = 0.5
         pbln_map(:,:)%ipb = 7
      END WHERE

      ! iquadran was not found:
      WHERE ( pbln_map(:,:)%iqdrn < 1 )
         pbln_map(:,:)%iqdrn = 1 ! maybe bad... but at least reported in ID_problem ...
         pbln_map(:,:)%ipb = 44
      END WHERE

      WHERE (kmsk_ignr_trg <= -1) pbln_map(:,:)%ipb = -1 ! Nearest point was not found by "FIND_NEAREST"
      WHERE (kmsk_ignr_trg ==  0) pbln_map(:,:)%ipb = -2 ! No idea if possible... #lolo
      WHERE (kmsk_ignr_trg <  -2) pbln_map(:,:)%ipb = -3 ! No idea if possible... #lolo


      !debug:
      !WRITE(cf_tmp,'("in_mbl_alfa_thread",i2.2,".nc")') 0
      !PRINT *, 'Debug/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      !CALL DUMP_FIELD(REAL(pbln_map(:,:)%ralfa,4), cf_tmp, 'alfa')
      !!
      !WRITE(cf_tmp,'("in_mbl_beta_thread",i2.2,".nc")') 0
      !PRINT *, 'Debug/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      !CALL DUMP_FIELD(REAL(pbln_map(:,:)%rbeta,4), cf_tmp, 'beta')
      !debug.
      !$OMP BARRIER
      
      !BILIN_MAP(:,:) = pbln_map(:,:)

      DEALLOCATE ( kmsk_ignr_trg, ki_nrst, kj_nrst )

   END SUBROUTINE MAPPING_BL




   SUBROUTINE LOCAL_COORD(pxlam, pxphi, pa, pb, ipb)

      !!----------------------------------------------------------
      !!           ***  SUBROUTINE  local_coord    ***
      !!
      !!  ** Purpose : Compute the local coordinate in a grid cell
      !!
      !!  ** Method : from N. Daget Web page :
      !!       http://aton.cerfacs.fr/~daget/TECHREPORT/TR_CMGC_06_18_html/node8.html
      !!
      !! * history:
      !!      Original : J.M. Molines ( May 2007)
      !!----------------------------------------------------------

      !! * Arguments
      REAL(8), DIMENSION(0:4), INTENT(in)  :: pxlam, pxphi
      REAL(8)                , INTENT(out) :: pa, pb
      INTEGER                , INTENT(out) :: ipb  !: 0 if everything went fine, > 0 otherwize!

      !! * Local variables
      REAL(8) :: &
         &      zalpha=0. , zbeta=0., zresmax=0.1, zres, &
         &      zdeta, zdalp, zdbet, zdlam, zdphi, z1, z2, z3

      REAL(8), DIMENSION(2,2):: za
      REAL(8), DIMENSION(0:4):: zxlam

      INTEGER :: niter=0

      INTEGER, PARAMETER :: itermax=200 !: maximum of iteration

      ipb = 0

      zxlam(:) = pxlam(:)       !: save input longitude in workinh array

      ! when near the 0 deg line and we must work in the frame -180 180
      IF ((ABS(zxlam(1)-zxlam(4))>=180.).OR.(ABS(zxlam(1)-zxlam(2))) >= 180.  &
         &                            .OR.(ABS(zxlam(1)-zxlam(3))  >= 180. )) &
         &  zxlam = degE_to_degWE(zxlam)

      zres=1000. ; zdlam=0.5 ;  zdphi=0.5 ;  zalpha=0. ;  zbeta=0. ;  niter=0

      DO WHILE ( (zres > zresmax).AND.(niter < itermax) )

         z1 = (zxlam(2) - zxlam(1))
         z2 = (zxlam(1) - zxlam(4))
         z3 = (zxlam(3) - zxlam(2))

         za(1,1) =  z1 + (z2 + z3 )*zbeta
         za(1,2) = -z2 + (z2 + z3 )*zalpha

         za(2,1) = pxphi(2) - pxphi(1) + (pxphi(1) - pxphi(4) + pxphi(3) - pxphi(2))*zbeta
         za(2,2) = pxphi(4) - pxphi(1) + (pxphi(1) - pxphi(4) + pxphi(3) - pxphi(2))*zalpha

         ! Determinant
         zdeta = det(za(1,1), za(1,2), za(2,1), za(2,2) )

         !! Solution of
         !! | zdlam |        | zdalp |
         !! |       | =  za .|       |
         !! | zdphi |        | zdbet |

         zdeta = ( SIGN(1._8,zdeta)*MAX(ABS(zdeta), repsilon) )  ! just to avoid FPE division by zero sometimes...

         zdalp = DET(  zdlam , za(1,2), zdphi  , za(2,2) ) / zdeta
         zdbet = DET( za(1,1), zdlam  , za(2,1), zdphi   ) / zdeta

         !! Compute residual ( loop criteria)
         zres = sqrt(zdalp*zdalp + zdbet*zdbet )

         !! Compute alpha and beta from 1rst guess :
         zalpha = zalpha + zdalp
         zbeta  = zbeta  + zdbet

         !! Compute corresponding lon/lat for this alpha, beta
         zdlam = zxlam(0) - ((1.-zalpha)*(1-zbeta)*zxlam(1) + zalpha*(1-zbeta)*zxlam(2)  &
            &                    +  zalpha*zbeta*zxlam(3) + (1-zalpha)*zbeta*zxlam(4))
         zdphi = pxphi(0)  - ((1.-zalpha)*(1-zbeta)*pxphi(1)  + zalpha*(1-zbeta)*pxphi(2)   &
            &                    +  zalpha*zbeta*pxphi(3)  + (1-zalpha)*zbeta*pxphi(4))

         niter = niter + 1  ! increment iteration counter

      END DO   ! loop until zres small enough (or itermax reach )

      IF ( niter >= itermax )   THEN
         zalpha = 0.5
         zbeta  = 0.5
         ipb    = 11
      END IF

      pa = zalpha
      pb = zbeta

      !! Problem if the 4 latitudes surrounding 'lati' are equal!
      IF ( (pxphi(1)==pxphi(2)).AND.(pxphi(2)==pxphi(3)).AND.(pxphi(3)==pxphi(4)) ) THEN
         pa  = 0.5
         pb  = 0.5
         ipb = 12
      END IF
      
   END SUBROUTINE LOCAL_COORD


   FUNCTION DET(p1, p2, p3, p4)
      !!----------------------------------------------------------
      !!          ***  FUNCTION DET   ***
      !!
      !!    ** Purpose : compute determinant
      !!
      !! * history:
      !!     J.M. Molines, may 2007
      !!----------------------------------------------------------
      REAL(8),INTENT(in) :: p1, p2, p3, p4
      REAL(8)            :: DET
      DET = p1*p4 - p2*p3
   END FUNCTION DET


   FUNCTION heading(plona, plonb, plata, platb)
      !!--------------------------------------------------------------
      !!           ***   FUNCTION HEADING   ***
      !!
      !!   ** Purpose: Compute true heading between point a and b
      !!
      !!   ** Method : suppose that the 2 points are not too far away from each other
      !!               so that heading can be computed with loxodromy
      !!
      !!      "heading" is the angle (in degrees) that going from point A to point B makes:
      !!               * 100% northward =>   0 deg.
      !!               * 100% eastward  =>  90 deg.
      !!               * 100% southward => 180 deg.
      !!                 ...
      !!
      !!  * history
      !!         J.M. Molines, may 2007
      !!--------------------------------------------------------------
      !!  Arguments
      REAL(8), INTENT(in) :: plata, plona, platb, plonb
      REAL(8) :: heading

      !!  Local variables
      REAL(8)  :: zpi, zconv, angled, xa, xb, ya, yb, xb_xa, rr

      zpi=ACOS(-1._8)
      zconv=zpi/180.  ! for degree to radian conversion

      !! There is a problem if the Greenwich meridian pass between a and b
      IF (iverbose>1) WRITE(6,*) ' * HEADIN() => Plonb  Plona ' , plonb, plona
      xa=plona*zconv
      xb=plonb*zconv

      rr = MAX(ABS(tan(zpi/4.-zconv*plata/2.)), repsilon)  !lolo just to avoid FPE sometimes
      ya = -LOG(rr)

      rr = MAX(ABS(tan(zpi/4.-zconv*platb/2.)), repsilon)  !lolo just to avoid FPE sometimes
      yb = -LOG(rr)

      IF (iverbose>1) WRITE(6,*) ' * HEADIN() =>  xa_xb , modulo 2pi', xb-xa, MOD((xb-xa),2*zpi)
      xb_xa=MOD((xb-xa),2*zpi)

      IF ( xb_xa >=  zpi ) xb_xa = xb_xa -2*zpi
      IF ( xb_xa <= -zpi ) xb_xa = xb_xa +2*zpi
      IF (iverbose>1)  WRITE(6,*) ' * HEADIN() => yb-ya, xb_xa ', yb-ya , xb_xa

      angled = ATAN2(xb_xa, yb - ya)

      heading=angled*180./zpi
      IF (heading < 0) heading = heading + 360._8

   END FUNCTION heading


   
   SUBROUTINE SAVE_BL_MAPPING_NC(cf_out, plon, plat, pbln_map, rflag )
      !!=========================================================================
      USE netcdf
      USE io_ezcdf, ONLY : sherr      
      !!
      CHARACTER(len=*),              INTENT(in) :: cf_out
      REAL(8),       DIMENSION(:,:), INTENT(in) :: plon, plat
      TYPE(bln_map), DIMENSION(:,:), INTENT(in) :: pbln_map
      REAL(4),                       INTENT(in) :: rflag
      !!=========================================================================
      INTEGER                                   :: id_f, id_x, id_y, id_lo, id_la
      CHARACTER(LEN=400), PARAMETER   ::     &
         &    cabout = 'Created with SOSIE interpolation environement => https://github.com/brodeau/sosie/'
      INTEGER          :: nx, ny, id_v1, id_v2, id_v3, id_v4, id_v5, id_v6, id_dnp
      CHARACTER(len=80), PARAMETER :: crtn = 'SAVE_BL_MAPPING_NC'
      !!=========================================================================
      nx = SIZE(pbln_map,1)
      ny = SIZE(pbln_map,2)
      !!
      !!           CREATE NETCDF OUTPUT FILE :
      CALL sherr( NF90_CREATE(cf_out, NF90_NETCDF4, id_f),  crtn,cf_out,'dummy')

      CALL sherr( NF90_DEF_DIM(id_f, 'x',  nx, id_x), crtn,cf_out,'dummy')
      CALL sherr( NF90_DEF_DIM(id_f, 'y',  ny, id_y), crtn,cf_out,'dummy')

      CALL sherr( NF90_DEF_VAR(id_f, 'lon',  NF90_DOUBLE, (/id_x,id_y/), id_lo, deflate_level=5),  crtn,cf_out,'lon')
      CALL sherr( NF90_DEF_VAR(id_f, 'lat',  NF90_DOUBLE, (/id_x,id_y/), id_la, deflate_level=5),  crtn,cf_out,'lat')
      !!
      CALL sherr( NF90_DEF_VAR(id_f, 'iP',   NF90_INT,    (/id_x,id_y/), id_v1, deflate_level=5),  crtn,cf_out,'iP'   )
      CALL sherr( NF90_DEF_VAR(id_f, 'jP',   NF90_INT,    (/id_x,id_y/), id_v2, deflate_level=5),  crtn,cf_out,'jP'   )
      CALL sherr( NF90_DEF_VAR(id_f, 'iqd',  NF90_INT,    (/id_x,id_y/), id_v3, deflate_level=5),  crtn,cf_out,'iqd'  )
      CALL sherr( NF90_DEF_VAR(id_f, 'alfa', NF90_DOUBLE, (/id_x,id_y/), id_v4, deflate_level=5),  crtn,cf_out,'alfa' )
      CALL sherr( NF90_DEF_VAR(id_f, 'beta', NF90_DOUBLE, (/id_x,id_y/), id_v5, deflate_level=5),  crtn,cf_out,'beta' )
      CALL sherr( NF90_DEF_VAR(id_f, 'ipb',  NF90_INT,    (/id_x,id_y/), id_v6, deflate_level=5),  crtn,cf_out,'ipb'  )

      IF ( l_save_distance_to_np ) THEN
         CALL sherr( NF90_DEF_VAR(id_f, 'dist_np', NF90_REAL, (/id_x,id_y/), id_dnp, deflate_level=5), crtn,cf_out,'dist_np' )
         CALL sherr( NF90_PUT_ATT(id_f, id_dnp, 'long_name', 'Distance to nearest point'),             crtn,cf_out,'dist_np' )
         CALL sherr( NF90_PUT_ATT(id_f, id_dnp, 'units'    , 'km'                       ),             crtn,cf_out,'dist_np' )
      END IF

      IF ( rflag /= 0. ) THEN
         CALL sherr( NF90_PUT_ATT(id_f, id_v1, '_FillValue',  INT(rflag)  ), crtn,cf_out,'iP   (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v2, '_FillValue',  INT(rflag)  ), crtn,cf_out,'jP   (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v3, '_FillValue',  INT(rflag)  ), crtn,cf_out,'iqd  (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v4, '_FillValue', REAL(rflag,8)), crtn,cf_out,'alfa (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v5, '_FillValue', REAL(rflag,8)), crtn,cf_out,'beta (masking)' )
      END IF
      !lolo

      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'Info', 'File containing mapping/weight information for bilinear interpolation with SOSIE.'), &
         &      crtn,cf_out,'dummy')
      CALL sherr( NF90_PUT_ATT(id_f, NF90_GLOBAL, 'About', trim(cabout)),  crtn,cf_out,'dummy')

      CALL sherr( NF90_ENDDEF(id_f),  crtn,cf_out,'dummy') ! END OF DEFINITION

      CALL sherr( NF90_PUT_VAR(id_f, id_lo, plon),     crtn,cf_out,'lon')
      CALL sherr( NF90_PUT_VAR(id_f, id_la, plat),     crtn,cf_out,'lat')

      CALL sherr( NF90_PUT_VAR(id_f, id_v1, pbln_map(:,:)%jip   ),  crtn,cf_out,'iP')
      CALL sherr( NF90_PUT_VAR(id_f, id_v2, pbln_map(:,:)%jjp   ),  crtn,cf_out,'jP')
      CALL sherr( NF90_PUT_VAR(id_f, id_v3, pbln_map(:,:)%iqdrn ),  crtn,cf_out,'iqd')
      CALL sherr( NF90_PUT_VAR(id_f, id_v4, pbln_map(:,:)%ralfa ),  crtn,cf_out,'alfa')
      CALL sherr( NF90_PUT_VAR(id_f, id_v5, pbln_map(:,:)%rbeta ),  crtn,cf_out,'beta')
      CALL sherr( NF90_PUT_VAR(id_f, id_v6, pbln_map(:,:)%ipb   ),  crtn,cf_out,'ipb')

      IF ( l_save_distance_to_np ) CALL sherr( NF90_PUT_VAR(id_f, id_dnp, distance_to_np), crtn,cf_out,'lon' )

      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_out,'dummy')

   END SUBROUTINE SAVE_BL_MAPPING_NC


   SUBROUTINE  RD_MAPPING_AB(cf_in, pbln_map)
      !!
      USE netcdf
      USE io_ezcdf, ONLY : sherr
      !!
      CHARACTER(len=*),              INTENT(in)  :: cf_in
      TYPE(bln_map), DIMENSION(:,:), INTENT(out) :: pbln_map
      !!
      INTEGER :: id_f
      INTEGER :: id_v1, id_v2, id_v3, id_v4, id_v5, id_v6
      !!
      CHARACTER(len=80), PARAMETER :: crtn = 'RD_MAPPING_AB'
      !!
      CALL sherr( NF90_OPEN(cf_in, NF90_NOWRITE, id_f),  crtn,cf_in,'dummy' )
      !!
      CALL sherr( NF90_INQ_VARID(id_f, 'iP',   id_v1),  crtn,cf_in,'iP'   )
      CALL sherr( NF90_INQ_VARID(id_f, 'jP',   id_v2),  crtn,cf_in,'jP'   )
      CALL sherr( NF90_INQ_VARID(id_f, 'iqd',  id_v3),  crtn,cf_in,'iqd'  )
      CALL sherr( NF90_INQ_VARID(id_f, 'alfa', id_v4),  crtn,cf_in,'alfa' )
      CALL sherr( NF90_INQ_VARID(id_f, 'beta', id_v5),  crtn,cf_in,'beta' )
      CALL sherr( NF90_INQ_VARID(id_f, 'ipb',  id_v6),  crtn,cf_in,'ipb'  )
      !!
      CALL sherr( NF90_GET_VAR(id_f, id_v1, pbln_map(:,:)%jip   ), crtn,cf_in,'iP'   )
      CALL sherr( NF90_GET_VAR(id_f, id_v2, pbln_map(:,:)%jjp   ), crtn,cf_in,'jP'   )
      CALL sherr( NF90_GET_VAR(id_f, id_v3, pbln_map(:,:)%iqdrn ), crtn,cf_in,'iqd'  )
      CALL sherr( NF90_GET_VAR(id_f, id_v4, pbln_map(:,:)%ralfa ), crtn,cf_in,'alfa' )
      CALL sherr( NF90_GET_VAR(id_f, id_v5, pbln_map(:,:)%rbeta ), crtn,cf_in,'beta' )
      CALL sherr( NF90_GET_VAR(id_f, id_v6, pbln_map(:,:)%ipb   ), crtn,cf_in,'ipb'  )
      !!
      CALL sherr( NF90_CLOSE(id_f),  crtn,cf_in,'dummy' )
      !!
   END SUBROUTINE RD_MAPPING_AB





   SUBROUTINE GIVE_NGHBR_POINTS_SRC( kewp, iX, jX, nxs, nys, xlon_src, xlat_src,  &
      &                              lskip, iXm1, iXp1, jXm1, jXp1, plonX, platX, &
      &                              plonN, platN, plonE, platE, plonS, platS, plonW, platW )
      !!==============================================================================
      !!
      !! INPUT:
      !!   * kewp : easr-west periodicity
      !!   * iX   : i-index of nearest point on source grid
      !!   * jX   : j-index of nearest point on source grid
      !!   * nxs, nys : shape of input array GLOBAL ???
      !!   * xlon_src : 2D array of source longitude
      !!   * xlat_src : 2D array of source latitude
      !!
      !! OUTPUT:
      !!   * lskip : if TRUE, skip this point...
      !!   * plonX, platX : coordinates of nearest point on source grid
      !!   * plonN, platN, plonE, platE, plonS, platS, plonW, platW : coordinates of the 4 points surrounding the nearest point on source grid
      !!
      !!==========================================================================================================
      INTEGER,                 INTENT(in)  :: kewp, iX, jX, nxs, nys
      REAL(8), DIMENSION(:,:), INTENT(in)  :: xlon_src, xlat_src
      LOGICAL,                 INTENT(out) :: lskip
      INTEGER,                 INTENT(out) :: iXm1, iXp1, jXm1, jXp1
      REAL(8),                 INTENT(out) :: plonX, platX, plonN, platN, plonE, platE, plonS, platS, plonW, platW
      !!
      !INTEGER :: nxs, nys
      !!==========================================================================================================
      !nxs = SIZE(xlon_src,1)
      !nys = SIZE(xlon_src,2)
      !
      ! jX<ny1 < last upper row was an extrapolation!
      iXm1 = iX-1
      iXp1 = iX+1
      jXm1 = jX-1
      jXp1 = jX+1

      plonX = 0.
      platX = 0.
      plonN = 0.
      platN = 0.
      plonE = 0.
      platE = 0.
      plonS = 0.
      platS = 0.
      plonW = 0.
      platW = 0.

      !! LOLO: remove from now, it's problematic in OMP with x-decomposition !!!!
      !!      => all EWP shit should be done before OMP???
      !IF ( iXm1 == 0 ) THEN
      !   !! We are in the extended case !!!
      !   IF ( kewp>=0 ) iXm1 = nxs - kewp
      !END IF
      !IF ( iXp1 == nxs+1 ) THEN
      !   IF ( kewp>=0 ) iXp1 = 1   + kewp
      !END IF
      !IF(lpdebug) WRITE(6,*) ' *** #DEBUG: iXm1, iXp1 =', iXm1, iXp1
      !IF(lpdebug) WRITE(6,*) ' *** #DEBUG: jXm1, jXp1 =', jXm1, jXp1

      lskip = ((iXm1 < 1).OR.(jXm1 < 1).OR.(iXp1 > nxs).OR.(jXp1 > nys))

      IF ( .NOT. lskip) THEN

         plonX = MOD(xlon_src(iX,jX)  , 360._8)
         plonN = MOD(xlon_src(iX,jXp1), 360._8)
         plonE = MOD(xlon_src(iXp1,jX), 360._8)
         plonS = MOD(xlon_src(iX,jXm1), 360._8)
         plonW = MOD(xlon_src(iXm1,jX), 360._8)

         platX = xlat_src(iX,jX)   ! nearest point
         platN = xlat_src(iX,jXp1) ! N (grid)
         platE = xlat_src(iXp1,jX) ! E (grid)
         platS = xlat_src(iX,jXm1) ! S (grid)
         platW = xlat_src(iXm1,jX) ! W (grid)

      ELSE
         IF(iverbose>0) WRITE(6,*) 'WARNING: mod_bilin_2d.f90 => bound problem => ',xlon_src(iX,jX),xlon_src(iX,jX),nxs,nys,iX,jX
         IF(iverbose>0) WRITE(6,*) '          iXm1, iXp1, nxs =', iXm1, iXp1, nxs
         IF(iverbose>0) WRITE(6,*) '          jXm1, jXp1, nys =', jXm1, jXp1, nys
         !IF(iverbose>0) WRITE(6,*) '         => ignoring current nearest point for i,j =', ji, jj, '(of target domain)'
         IF(iverbose>0) WRITE(6,*) ''
      END IF

   END SUBROUTINE GIVE_NGHBR_POINTS_SRC


END MODULE MOD_BILIN_2D
