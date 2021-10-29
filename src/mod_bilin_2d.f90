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
   !!-----------------------------------------------------------------

   USE mod_conf
   USE mod_manip, ONLY: EXT_NORTH_TO_90_REGG, DEGE_TO_DEGWE, FIND_NEAREST_POINT, DISTANCE
   USE mod_poly
   USE mod_grids, ONLY: mask_ignore_trg

   USE io_ezcdf, ONLY : DUMP_FIELD !lolodbg
   
   IMPLICIT NONE

   PRIVATE

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

   TYPE(bln_map), DIMENSION(:,:), ALLOCATABLE, PUBLIC, SAVE :: bilin_map
   REAL(4),       DIMENSION(:,:), ALLOCATABLE, SAVE :: distance_to_np
   LOGICAL,       DIMENSION(:),   ALLOCATABLE, SAVE :: l_1st_call_bilin
   !--------


   CHARACTER(len=400), SAVE :: cf_wght_bilin
   LOGICAL,    PARAMETER    :: l_save_distance_to_np=.TRUE. !: for each point of target grid, shows the distance to the nearest point

   REAL(8), PARAMETER :: repsilon = 1.E-9

   !LOGICAL, SAVE :: l_last_y_row_missing = .FALSE.

   !! PUBLIC:

   LOGICAL,                                 PUBLIC, SAVE :: l_skip_bilin_mapping

   PUBLIC :: BILIN_2D_INIT, BILIN_2D_WRITE_MAPPING, BILIN_2D, MAPPING_BL, INTERP_BL

CONTAINS


   SUBROUTINE BILIN_2D_INIT()
      !!==============================================================================
      !!
      USE io_ezcdf, ONLY : TEST_XYZ  ! , DUMP_2D_FIELD
      !!
      !REAL(8),    DIMENSION(:,:),           INTENT(in) :: px_src, py_src
      !REAL(8),    DIMENSION(:,:),           INTENT(in) :: px_trg, py_trg
      !INTEGER(1), DIMENSION(:,:), OPTIONAL, INTENT(in) :: mask_domain_trg
      !!
      LOGICAL :: lefw
      !!==============================================================================
      !!
      WRITE(6,*) ''; WRITE(6,*) ''
      WRITE(6,*) '###################################################################'
      WRITE(6,*) '#                  BILINEAR 2D INITIALIZATION'
      WRITE(6,*) '###################################################################'
      WRITE(6,*) ''
      WRITE(6,'("   * Allocating array bilin_map: ",i5," x ",i5)') ni_trg, nj_trg
      ALLOCATE ( bilin_map(ni_trg,nj_trg) )

      IF (l_save_distance_to_np) ALLOCATE ( distance_to_np(ni_trg,nj_trg) )

      ALLOCATE ( l_1st_call_bilin(Nthrd) )
      l_1st_call_bilin(:)     = .TRUE.

      WRITE(6,*) '  * Allocations done...'
      WRITE(6,*) ''
      WRITE(cf_wght_bilin,'("sosie_mapping_",a,".nc")') TRIM(cpat)

      WRITE(6,*) '  * Mapping file is "',TRIM(cf_wght_bilin),'" !'

      INQUIRE(FILE=cf_wght_bilin, EXIST=lefw )
      IF ( lefw ) THEN
         WRITE(6,*) '    => it was found in current directory !'
         WRITE(6,*) '      => still! Make sure that this is really the one you need...'
         CALL RD_MAPPING_AB(cf_wght_bilin, bilin_map(:,:))
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
      END IF

      WRITE(6,*) '  Initializations for BILIN_2D done !'
      WRITE(6,*) ''
      WRITE(6,*) '###################################################################'
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
   END SUBROUTINE BILIN_2D_INIT


   SUBROUTINE BILIN_2D_WRITE_MAPPING( )
      !!
      WRITE(6,*) ''
      IF ( l_skip_bilin_mapping ) THEN
         WRITE(6,*) '  ==> NOT writing bilin mapping file because it was found...'
      ELSE
         WRITE(6,*) '  ==> writing bilin mapping into "', TRIM(cf_wght_bilin),'" !'
         CALL P2D_MAPPING_AB( cf_wght_bilin, lon_trg, lat_trg, bilin_map, rmissval )
      END IF
      !!
      IF (l_save_distance_to_np) DEALLOCATE( distance_to_np )
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
   END SUBROUTINE BILIN_2D_WRITE_MAPPING



   SUBROUTINE BILIN_2D( k_ew_per, pX1, pY1, pZ1, pX2, pY2, pZ2, ithrd,  mask_domain_trg )
      !!================================================================
      !!
      !! INPUT :     k_ew_per : east-west periodicity
      !!                        k_ew_per = -1  --> no periodicity
      !!                        k_ew_per >= 0  --> periodicity with overlap of k_ew_per points
      !!             pX1   : 2D source longitude array of shape (ni,nj)
      !!             pY1   : 2D source latitude  array of shape (ni,nj)
      !!             pZ1   : source field on source grid  "    "
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
      !!
      !! Input/Output arguments
      INTEGER,                 INTENT(in)  :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pX1, pY1
      REAL(4), DIMENSION(:,:), INTENT(in)  :: pZ1
      REAL(8), DIMENSION(:,:), INTENT(in)  :: pX2, pY2
      REAL(4), DIMENSION(:,:), INTENT(out) :: pZ2
      INTEGER,                 INTENT(in)  :: ithrd ! # OMP thread
      INTEGER(1), OPTIONAL, DIMENSION(:,:), INTENT(in) :: mask_domain_trg
      !! Local variables
      INTEGER :: nx2, ny2, iqd, iP, jP
      REAL(8) :: alpha, beta, rmeanv
      INTEGER :: ji, jj, iom1, iom2

      nx2 = SIZE(pZ2,1)
      ny2 = SIZE(pZ2,2)

      iom1 = io1(ithrd)
      iom2 = io2(ithrd)
      

      PRINT *, 'LOLOdbg: BILIN_2D() => shape of pX1 =', SIZE(pX1,1), SIZE(pX1,2)

      !IF ( l_1st_call_bilin(ithrd) ) l_last_y_row_missing = .FALSE.

      pZ2(:,:) = rmissval ! Flagging non-interpolated output points

      PRINT *, 'LOLODBG bilin_2d, shape bilin_map =', SIZE(bilin_map,1), SIZE(bilin_map,2), '#', ithrd
      PRINT *, 'LOLODBG bilin_2d, shape mask_ignore_trg =', SIZE(mask_ignore_trg,1), SIZE(mask_ignore_trg,2), '#', ithrd
      PRINT *, '  io1(ithrd), io2(ithrd)', iom1, iom2, '#', ithrd

      DO ji=iom1, iom2
         WHERE ( bilin_map(ji,:)%jip < 1 ) mask_ignore_trg(ji,:) = 0
         WHERE ( bilin_map(ji,:)%jjp < 1 ) mask_ignore_trg(ji,:) = 0

         bilin_map(ji,:)%jip = MAX( bilin_map(ji,:)%jip , 1 )  ! so no i or j <= 0
         bilin_map(ji,:)%jjp = MAX( bilin_map(ji,:)%jjp , 1 )  ! so no i or j <= 0
         
      END DO
      !WHERE ( bilin_map(iom1:iom2,:)%jip < 1 ) mask_ignore_trg = 0
      !WHERE ( bilin_map(iom1:iom2,:)%jjp < 1 ) mask_ignore_trg = 0

      !!WHERE ( (IMETRICS(:,:,3 < 1) ) mask_ignore_trg = 0 ; ! iqd => problem in interp ORCA2->ORCA1 linked to iqd < 1

      !bilin_map(iom1:iom2,:)%jip = MAX( bilin_map(iom1:iom2,:)%jip , 1 )  ! so no i or j <= 0
      !bilin_map(iom1:iom2,:)%jjp = MAX( bilin_map(iom1:iom2,:)%jjp , 1 )  ! so no i or j <= 0

      DO jj=1, ny2
         DO ji=1, nx2
            iP    = bilin_map(ji+iom1-1,jj)%jip
            jP    = bilin_map(ji+iom1-1,jj)%jjp
            iqd = bilin_map(ji+iom1-1,jj)%iqdrn
            alpha = bilin_map(ji+iom1-1,jj)%ralfa
            beta  = bilin_map(ji+iom1-1,jj)%rbeta
            !!
            IF ( (ABS(degE_to_degWE(pX1(iP,jP))-degE_to_degWE(pX2(ji,jj)))<1.E-5) .AND. (ABS(pY1(iP,jP)-pY2(ji,jj))<1.E-5) ) THEN
               !! COPY:
               IF (iverbose>0) WRITE(6,*) ' *** BILIN_2D: "identical point" detected (crit: 1.E-5) => copying value, no interpolation!'
               pZ2(ji,jj) = pZ1(iP,jP)
            ELSE
               !! INTERPOLATION:
               !pZ2(ji,jj) = INTERP_BL(k_ew_per, iP, jP, iqd, alpha, beta, pZ1)
               !pZ2(ji,jj) = INTERP_BL(k_ew_per, ithrd, ji, jj, pZ1)

               CALL INTERP_BL( k_ew_per, ithrd, ji, jj, pZ1, pZ2(ji,jj) )
               
            END IF
         END DO
      END DO

      pZ2(:,:) = pZ2(:,:)*REAL(mask_ignore_trg(iom1:iom2,:), 4) + REAL(1-mask_ignore_trg(iom1:iom2,:), 4)*(-9995.) ! masking problem points as in mask_ignore_trg


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

      l_1st_call_bilin(ithrd) = .FALSE.

   END SUBROUTINE BILIN_2D


   !FUNCTION INTERP_BL(k_ew_per, kiP, kjP, kqd, pa, pb, Z_in)
   SUBROUTINE INTERP_BL(k_ew_per, ithrd, ilt, jlt, Z_in, pres )
   
      INTEGER,                 INTENT(in) :: k_ew_per
      INTEGER,                 INTENT(in) :: ithrd, ilt, jlt   ! LOCAL (local omp domain) coordinates of treated point on target domain
      !INTEGER,                 INTENT(in) :: kiP, kjP, kqd
      !REAL(8),                 INTENT(in) :: pa, pb
      REAL(4), DIMENSION(:,:), INTENT(in) :: Z_in
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

      !PRINT *, 'LOLOINTERP_BL: ithrd, shape(io1), i1, i2 =>', INT(ithrd,1), INT(SIZE(io1),1), INT(io1(ithrd),2), INT(io2(ithrd),2)
      
      !ki = ilt + io1(ithrd) - 1
      ki = io1(ithrd) - 1 + ilt
      kj = jlt

      !!
      Nitl = SIZE(bilin_map(:,:)%jip,1)
      IF( ki > Nitl ) THEN
         PRINT *, 'PROBLEM: ki > Nitl ! ithrd, ki, Nitl, ilt, io1(ithrd) ', INT(ithrd,1), ki, Nitl, ilt, io1(ithrd)
         STOP
      END IF
      
      kiP = bilin_map(ki,kj)%jip
      kjP = bilin_map(ki,kj)%jjp
      kqd = bilin_map(ki,kj)%iqdrn
      pa  = bilin_map(ki,kj)%ralfa
      pb  = bilin_map(ki,kj)%rbeta
      
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
         pres = ( Z_in(i1,j1)*w1 + Z_in(i2,j2)*w2 + Z_in(i3,j3)*w3 + Z_in(i4,j4)*w4 )/wup
      ENDIF

   END SUBROUTINE INTERP_BL




   SUBROUTINE MAPPING_BL(k_ew_per, plon_src, plat_src, plon_trg, plat_trg,  ithread, pmsk_dom_trg)

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
      REAL(8), DIMENSION(:,:), INTENT(in) :: plon_src, plat_src
      REAL(8), DIMENSION(:,:), INTENT(in) :: plon_trg, plat_trg
      !!
      INTEGER,    OPTIONAL,                 INTENT(in) :: ithread
      INTEGER(1), OPTIONAL, DIMENSION(:,:), INTENT(in) :: pmsk_dom_trg

      INTEGER :: &
         &     iqd, iqd0, iqd_old, &
         &     iP, jP,             &
         &     nxi, nyi, nxo, nyo, &
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
      TYPE(bln_map), DIMENSION(:,:), ALLOCATABLE :: zbln_map  ! local array on thread !
      INTEGER(4),    DIMENSION(:,:), ALLOCATABLE :: ki_nrst, kj_nrst
      INTEGER(1),    DIMENSION(:,:), ALLOCATABLE :: kmsk_ignr_trg
      !!
      REAL(8) :: zalfa, zbeta
      LOGICAL :: l_ok, lagain, lpdebug
      INTEGER :: icpt, ithrd, iom1, iom2
      !!
      !!DEBUG:
      !INTEGER :: jx
      CHARACTER(len=80) :: cf_tmp !lolodbg
      !!----------------------------------------------------------------------------

      ithrd = 1 ! no OpenMP !
      IF( PRESENT(ithread) ) ithrd = ithread
      iom1 = io1(ithrd)
      iom2 = io2(ithrd)

      nxi = SIZE(plon_src,1)
      nyi = SIZE(plon_src,2)
      WRITE(6,'("LOLOdbg/MAPPING_BL: shape of plon_src for thread #",i1,": ",i4.4,"x",i4.4)') ithrd, nxi, nyi
      

      

      nxo = SIZE(plon_trg,1)
      nyo = SIZE(plon_trg,2)

      IF( (ithrd==1).AND.((iom1/=1).OR.(iom2/=nxo)) ) THEN
         WRITE(6,*) 'ERROR in "MAPPING_BL", wrong OMP partitioning...'
         STOP
      END IF

      WRITE(6,'("LOLOdbg/MAPPING_BL: size of target domain for thread #",i1,": ",i4.4,"x",i4.4)') ithrd, nxo, nyo


      ALLOCATE ( zbln_map(nxo,nyo), kmsk_ignr_trg(nxo,nyo), ki_nrst(nxo,nyo), kj_nrst(nxo,nyo) )

      zbln_map(:,:)%jip    = 0
      zbln_map(:,:)%jjp    = 0
      zbln_map(:,:)%iqdrn  = 0
      zbln_map(:,:)%ralfa  = 0.
      zbln_map(:,:)%rbeta  = 0.
      zbln_map(:,:)%ipb    = 0
      kmsk_ignr_trg(:,:) = 1
      ki_nrst(:,:)      = 0
      kj_nrst(:,:)      = 0

      IF ( PRESENT(pmsk_dom_trg) ) kmsk_ignr_trg(:,:) = pmsk_dom_trg(:,:)

      !! DEBUG: checking input fields for each different thread:
      WRITE(cf_tmp,'("in_mbl_lon_src_",i2.2,".nc")') ithrd
      CALL DUMP_FIELD(REAL(plon_src(:,:),4), cf_tmp, 'lon')
      WRITE(cf_tmp,'("in_mbl_lat_src_",i2.2,".nc")') ithrd
      CALL DUMP_FIELD(REAL(plat_src(:,:),4), cf_tmp, 'lat')
      !!
      WRITE(cf_tmp,'("in_mbl_lon_trg_",i2.2,".nc")') ithrd
      CALL DUMP_FIELD(REAL(plon_trg,4), cf_tmp, 'lon')
      WRITE(cf_tmp,'("in_mbl_lat_trg_",i2.2,".nc")') ithrd
      CALL DUMP_FIELD(REAL(plat_trg,4), cf_tmp, 'lat')      
      !STOP

      CALL FIND_NEAREST_POINT( plon_trg, plat_trg, plon_src, plat_src, ki_nrst, kj_nrst,  &
         &                     ithread=ithrd, pmsk_dom_trg=kmsk_ignr_trg )

      !lolodbg:
      WRITE(cf_tmp,'("in_mbl_nrst_ji",i2.2,".nc")') ithread
      PRINT *, 'LOLOdbg/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      CALL DUMP_FIELD(REAL(ki_nrst,4), cf_tmp, 'ji')
      !
      WRITE(cf_tmp,'("in_mbl_nrst_jj",i2.2,".nc")') ithread
      PRINT *, 'LOLOdbg/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      CALL DUMP_FIELD(REAL(kj_nrst,4), cf_tmp, 'jj')
      !lolodbg.

      !$OMP BARRIER
      !STOP'lilo'

      DO jj = 1, nyo
         DO ji = 1, nxo

            lpdebug = ( (iverbose==2).AND.(ji==idb).AND.(jj==jdb) )

            IF ( kmsk_ignr_trg(ji,jj)==1 ) THEN
               !! => exclude regions that do not exist on source domain (kmsk_ignr_trg==0) and
               !! points for which the nearest point was not found (kmsk_ignr_trg==-1 or -2)

               !! Now deal with horizontal interpolation
               !! set longitude of input point in accordance with lon ( [lon0, 360+lon0 [ )
               xP = plon_trg(ji,jj)
               yP = plat_trg(ji,jj)

               iP = ki_nrst(ji,jj)
               jP = kj_nrst(ji,jj)

               IF(lpdebug) WRITE(6,*) ' *** #DEBUG: xP, yP =', xP, yP
               IF(lpdebug) WRITE(6,*) ' *** #DEBUG: iP, jP, nxi, nyi =', iP, jP, nxi, nyi

               IF ( (iP/=INT(rmissval)).AND.(jP/=INT(rmissval)).AND.(jP<nyi) ) THEN
                  ! jP<ny1 < last upper row was an extrapolation!
                  iPm1 = iP-1
                  iPp1 = iP+1
                  jPm1 = jP-1
                  jPp1 = jP+1

                  IF ( iPm1 == 0 ) THEN
                     !! We are in the extended case !!!
                     IF ( k_ew_per>=0 ) iPm1 = nxi - k_ew_per
                  END IF

                  IF ( iPp1 == nxi+1 ) THEN
                     IF ( k_ew_per>=0 ) iPp1 = 1   + k_ew_per
                  END IF

                  IF(lpdebug) WRITE(6,*) ' *** #DEBUG: iPm1, iPp1 =', iPm1, iPp1
                  IF(lpdebug) WRITE(6,*) ' *** #DEBUG: jPm1, jPp1 =', jPm1, jPp1

                  IF ((iPm1 < 1).OR.(jPm1 < 1).OR.(iPp1 > nxi)) THEN
                     IF(iverbose>0) WRITE(6,*) 'WARNING: mod_bilin_2d.f90 => bound problem => ',xP,yP,nxi,nyi,iP,jP
                     IF(iverbose>0) WRITE(6,*) '          iPm1, iPp1, nxi =', iPm1, iPp1, nxi
                     IF(iverbose>0) WRITE(6,*) '          jPm1, jPp1, nyi =', jPm1, jPp1, nyi
                     IF(iverbose>0) WRITE(6,*) '         => ignoring current nearest point for i,j =', ji, jj, '(of target domain)'
                     IF(iverbose>0) WRITE(6,*) ''
                  ELSE

                     lonP = MOD(plon_src(iP,jP)  , 360._8) ; latP = plat_src(iP,jP)   ! nearest point
                     lonN = MOD(plon_src(iP,jPp1), 360._8) ; latN = plat_src(iP,jPp1) ! N (grid)
                     lonE = MOD(plon_src(iPp1,jP), 360._8) ; latE = plat_src(iPp1,jP) ! E (grid)
                     lonS = MOD(plon_src(iP,jPm1), 360._8) ; latS = plat_src(iP,jPm1) ! S (grid)
                     lonW = MOD(plon_src(iPm1,jP), 360._8) ; latW = plat_src(iPm1,jP) ! W (grid)

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

                     IF (l_save_distance_to_np) distance_to_np(ji+iom1-1,jj) = DISTANCE(xP, lonP, yP, latP)

                     !! Problem is that sometimes, in the case of really twisted
                     !! meshes this method screws up, iqd is not what it
                     !! should be, so need the following loop on the value of iqd:
                     icpt = 0
                     lagain = .TRUE.
                     DO WHILE ( lagain )

                        SELECT CASE ( iqd ) ! point 2 3 4 are counter clockwise in the respective sector
                        CASE ( 1 )
                           loni(2) = lonE ; lati(2) = latE
                           loni(3) = MOD(plon_src(iPp1,jPp1), 360._8) ; lati(3) = plat_src(iPp1,jPp1)
                           loni(4) = lonN ; lati(4) = latN
                        CASE ( 2 )
                           loni(2) = lonS ; lati(2) = latS
                           loni(3) = MOD(plon_src(iPp1,jPm1), 360._8) ; lati(3) = plat_src(iPp1,jPm1)
                           loni(4) = lonE ; lati(4) = latE
                        CASE ( 3 )
                           loni(2) = lonW ; lati(2) = latW
                           loni(3) = MOD(plon_src(iPm1,jPm1), 360._8) ; lati(3) = plat_src(iPm1,jPm1)
                           loni(4) = lonS ; lati(4) = latS
                        CASE ( 4 )
                           loni(2) = lonN ; lati(2) = latN
                           loni(3) = MOD(plon_src(iPm1,jPp1), 360._8) ; lati(3) = plat_src(iPm1,jPp1)
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
                     zbln_map(ji,jj)%ipb = iproblem

                     !LOLO: mark this in ID_problem => case when iqd = iqd0 ( IF ( icpt == 5 ) ) above!!!
                     IF ( icpt == 5 ) zbln_map(ji,jj)%ipb = 2 ! IDing the screw-up from above...

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
                     zbln_map(ji,jj)%iqdrn = iqd
                     zbln_map(ji,jj)%ralfa = zalfa
                     zbln_map(ji,jj)%rbeta = zbeta
                     zbln_map(ji,jj)%ipb   = 0

                  END IF ! IF ((iPm1 < 1).OR.(jPm1 < 1).OR.(iPp1 > nxi))
               END IF ! IF ( (iP/=INT(rmissval)).AND.(jP/=INT(rmissval)).AND.(jP<nyi) )
            END IF ! IF ( kmsk_ignr_trg(ji,jj)==1 )
         ENDDO
      ENDDO

      !lolo: UGLY!
      zbln_map(:,:)%jip = ki_nrst(:,:)
      zbln_map(:,:)%jjp = kj_nrst(:,:)

      WHERE ( (ki_nrst == INT(rmissval)).OR.(kj_nrst == INT(rmissval)) )
         zbln_map(:,:)%ralfa  = rmissval
         zbln_map(:,:)%rbeta  = rmissval
      END WHERE

      ! Awkwardly fixing problematic points but remembering them in ID_problem

      ! Negative values that are actually 0
      WHERE ( ((zbln_map(:,:)%ralfa < 0.).AND.(zbln_map(:,:)%ralfa > -repsilon)) ) zbln_map(:,:)%ralfa = 0.0
      WHERE ( ((zbln_map(:,:)%rbeta < 0.).AND.(zbln_map(:,:)%rbeta > -repsilon)) ) zbln_map(:,:)%rbeta = 0.0

      WHERE ( (zbln_map(:,:)%ralfa > rmissval).AND.(zbln_map(:,:)%ralfa < 0.) )
         zbln_map(:,:)%ralfa = 0.5
         zbln_map(:,:)%ipb = 4
      END WHERE
      WHERE ( zbln_map(:,:)%ralfa > 1. )
         zbln_map(:,:)%ralfa = 0.5
         zbln_map(:,:)%ipb =  5
      END WHERE

      WHERE ( (zbln_map(:,:)%rbeta > rmissval).AND.(zbln_map(:,:)%rbeta < 0.) )
         zbln_map(:,:)%rbeta = 0.5
         zbln_map(:,:)%ipb = 6
      END WHERE
      WHERE ( zbln_map(:,:)%rbeta > 1. )
         zbln_map(:,:)%rbeta = 0.5
         zbln_map(:,:)%ipb = 7
      END WHERE

      ! iquadran was not found:
      WHERE ( zbln_map(:,:)%iqdrn < 1 )
         zbln_map(:,:)%iqdrn = 1 ! maybe bad... but at least reported in ID_problem ...
         zbln_map(:,:)%ipb = 44
      END WHERE

      WHERE (kmsk_ignr_trg <= -1) zbln_map(:,:)%ipb = -1 ! Nearest point was not found by "FIND_NEAREST"
      WHERE (kmsk_ignr_trg ==  0) zbln_map(:,:)%ipb = -2 ! No idea if possible... #lolo
      WHERE (kmsk_ignr_trg <  -2) zbln_map(:,:)%ipb = -3 ! No idea if possible... #lolo


      !lolodbg:
      WRITE(cf_tmp,'("in_mbl_alfa_thread",i2.2,".nc")') ithread
      PRINT *, 'LOLOdbg/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      CALL DUMP_FIELD(REAL(zbln_map(:,:)%ralfa,4), cf_tmp, 'alfa')
      !!
      WRITE(cf_tmp,'("in_mbl_beta_thread",i2.2,".nc")') ithread
      PRINT *, 'LOLOdbg/MAPPING_BL: saving '//TRIM(cf_tmp)//' !!!'
      CALL DUMP_FIELD(REAL(zbln_map(:,:)%rbeta,4), cf_tmp, 'beta')
      !lolodbg.


      !$OMP BARRIER
      STOP'mod_bilin_2d.f90:MAPPING_BL()=>lilo2'

      
      bilin_map(iom1:iom2,:) = zbln_map(:,:)

      DEALLOCATE ( zbln_map, kmsk_ignr_trg, ki_nrst, kj_nrst )

   END SUBROUTINE MAPPING_BL




   SUBROUTINE LOCAL_COORD(xlam, xphi, pa, pb, ipb)

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
      REAL(8), DIMENSION(0:4), INTENT(in)  :: xlam, xphi
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

      zxlam = xlam       !: save input longitude in workinh array

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

         za(2,1) = xphi(2) - xphi(1) + (xphi(1) - xphi(4) + xphi(3) - xphi(2))*zbeta
         za(2,2) = xphi(4) - xphi(1) + (xphi(1) - xphi(4) + xphi(3) - xphi(2))*zalpha

         ! Determinant
         zdeta = det(za(1,1), za(1,2), za(2,1), za(2,2) )

         !! Solution of
         !! | zdlam |        | zdalp |
         !! |       | =  za .|       |
         !! | zdphi |        | zdbet |

         zdeta = ( SIGN(1._8,zdeta)*MAX(ABS(zdeta), repsilon) )  ! just to avoid FPE division by zero sometimes...

         zdalp = det(zdlam,  za(1,2) , zdphi, za(2,2)  ) / zdeta
         zdbet = det(za(1,1)  , zdlam, za(2,1)   ,zdphi) / zdeta

         !! Compute residual ( loop criteria)
         zres = sqrt(zdalp*zdalp + zdbet*zdbet )

         !! Compute alpha and beta from 1rst guess :
         zalpha = zalpha + zdalp
         zbeta  = zbeta  + zdbet

         !! Compute corresponding lon/lat for this alpha, beta
         zdlam = zxlam(0) - ((1.-zalpha)*(1-zbeta)*zxlam(1) + zalpha*(1-zbeta)*zxlam(2)  &
            &                    +  zalpha*zbeta*zxlam(3) + (1-zalpha)*zbeta*zxlam(4))
         zdphi = xphi(0)  - ((1.-zalpha)*(1-zbeta)*xphi(1)  + zalpha*(1-zbeta)*xphi(2)   &
            &                    +  zalpha*zbeta*xphi(3)  + (1-zalpha)*zbeta*xphi(4))

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
      IF ( (xphi(1)==xphi(2)).AND.(xphi(2)==xphi(3)).AND.(xphi(3)==xphi(4)) ) THEN
         pa  = 0.5
         pb  = 0.5
         ipb = 12
      END IF

   END SUBROUTINE LOCAL_COORD


   FUNCTION det(p1, p2, p3, p4)
      !!----------------------------------------------------------
      !!          ***  FUNCTION DET   ***
      !!
      !!    ** Purpose : compute determinant
      !!
      !! * history:
      !!     J.M. Molines may 2007
      !!----------------------------------------------------------
      REAL(8),INTENT(in) :: p1, p2, p3, p4
      REAL(8) :: det
      det = p1*p4 - p2*p3
   END FUNCTION det


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



   SUBROUTINE P2D_MAPPING_AB(cf_out, plon, plat, pbln_map, rflag )

      USE netcdf
      USE io_ezcdf, ONLY : sherr

      CHARACTER(len=*),              INTENT(in) :: cf_out
      REAL(8),       DIMENSION(:,:), INTENT(in) :: plon, plat
      TYPE(bln_map), DIMENSION(:,:), INTENT(in) :: pbln_map
      REAL(4),                       INTENT(in) :: rflag
      !!
      INTEGER                                   :: id_f, id_x, id_y, id_lo, id_la
      CHARACTER(LEN=400), PARAMETER   ::     &
         &    cabout = 'Created with SOSIE interpolation environement => https://github.com/brodeau/sosie/'
      !!
      INTEGER          :: nx, ny, id_v1, id_v2, id_v3, id_v4, id_v5, id_v6, id_dnp
      CHARACTER(len=80), PARAMETER :: crtn = 'P2D_MAPPING_AB'

      nx = SIZE(pbln_map,1)
      ny = SIZE(pbln_map,2)


      !!           CREATE NETCDF OUTPUT FILE :
      !!           ---------------------------
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

   END SUBROUTINE P2D_MAPPING_AB


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

END MODULE MOD_BILIN_2D
