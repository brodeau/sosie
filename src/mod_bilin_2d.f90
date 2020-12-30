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
   USE mod_manip
   USE mod_poly

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=400), SAVE :: cf_wght_bilin
   LOGICAL,            SAVE :: l_skip_bilin_mapping
   LOGICAL,    PARAMETER    :: l_save_distance_to_np=.TRUE. !: for each point of target grid, shows the distance to the nearest point

   REAL(8), PARAMETER :: repsilon = 1.E-9

   LOGICAL, SAVE :: l_last_y_row_missing

   INTEGER, PARAMETER :: &
      &                     iverbose = 0  , &  ! 0 to 2...
      &                     idb      = 0  , &  ! i-index of point to debug on target domain (iverbose = 2 !)
      &                     jdb      = 0       ! j-index     "         "         "          "


   PUBLIC :: BILIN_2D_INIT, BILIN_2D_WRITE_MAPPING, BILIN_2D, MAPPING_BL, INTERP_BL


CONTAINS



   SUBROUTINE BILIN_2D_INIT( )
      !!
      LOGICAL :: lefw
      !!
      WRITE(6,*) ''; WRITE(6,*) ''
      WRITE(6,*) '###################################################################'
      WRITE(6,*) '#                  BILINEAR 2D INITIALIZATION'
      WRITE(6,*) '###################################################################'
      WRITE(6,*) ''
      WRITE(6,'("   * Allocating array bilin_map: ",i5," x ",i5)') ni_trg, nj_trg
      ALLOCATE ( bilin_map(ni_trg,nj_trg) )
      IF (l_save_distance_to_np) ALLOCATE ( distance_to_np(ni_trg,nj_trg) )
      WRITE(6,*) '  * Allocation done...'
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
         CALL P2D_MAPPING_AB( cf_wght_bilin, lon_trg, lat_trg, bilin_map, rflg )
      END IF
      !!
      IF (l_save_distance_to_np) DEALLOCATE( distance_to_np )
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
   END SUBROUTINE BILIN_2D_WRITE_MAPPING





   SUBROUTINE BILIN_2D( k_ew_per, X10, Y10, Z1, X20, Y20, Z2, ithrd,  mask_domain_trg )
      !!================================================================
      !!
      !! INPUT :     k_ew_per : east-west periodicity
      !!                        k_ew_per = -1  --> no periodicity
      !!                        k_ew_per >= 0  --> periodicity with overlap of k_ew_per points
      !!             X10   : 2D source longitude array (ni,nj) or (ni,1)
      !!             Y10   : 2D source latitude  array (ni,nj) or (nj,1)
      !!             Z1    : source field on source grid
      !!
      !!             X20   : 2D target longitude array (ni,nj) or (ni,1)
      !!             Y20   : 2D target latitude  array (ni,nj) or (nj,1)
      !!
      !! OUTPUT :
      !!             Z2    : field extrapolated from source to target grid
      !!
      !!
      !! OPTIONAL IN:
      !!      * mask_domain_trg: ignore (dont't treat) regions of the target domain where mask_domain_trg==0 !
      !!
      !!================================================================
      USE io_ezcdf, ONLY : TEST_XYZ  ! , DUMP_2D_FIELD
      !!
      !! Input/Output arguments
      INTEGER,                 INTENT(in)  :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X10, Y10
      REAL(4), DIMENSION(:,:), INTENT(in)  :: Z1
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X20, Y20
      REAL(4), DIMENSION(:,:), INTENT(out) :: Z2
      INTEGER,                 INTENT(in)  :: ithrd ! # OMP thread
      INTEGER(1), OPTIONAL, DIMENSION(:,:), INTENT(in) :: mask_domain_trg
      !! Local variables
      INTEGER :: nx1, ny1, ny1w, nx2, ny2, iqd, iP, jP

      REAL(8) :: alpha, beta, rmeanv, ymx
      LOGICAL :: l_add_extra_j
      INTEGER :: ji, jj

      REAL(8),    DIMENSION(:,:), ALLOCATABLE :: X1, Y1, X2, Y2, X1w, Y1w
      REAL(4),    DIMENSION(:,:), ALLOCATABLE :: Z1w
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: mask_ignore_trg !, msk_res

      CHARACTER(len=2)   :: ctype

      ctype = TEST_XYZ(X10,Y10,Z1)
      nx1 = SIZE(Z1,1)
      ny1 = SIZE(Z1,2)
      ALLOCATE ( X1(nx1,ny1) , Y1(nx1,ny1) )


      IF ( ctype == '1d' ) THEN
         DO jj=1, ny1
            X1(:,jj) = X10(:,1)
         END DO
         DO ji=1, nx1
            Y1(ji,:) = Y10(:,1)
         END DO
      ELSE
         X1 = X10 ; Y1 = Y10
      END IF


      !! Working arrays for source domain:
      l_add_extra_j = .FALSE.
      ny1w = ny1

      IF ( l_reg_src .AND. (ny1 > 20) .AND. (nx1 > 20) ) THEN !lolo, ensure it's a map, not a something fishy...
         ymx = Y1(nx1/2,ny1)
         IF ( (ymx < 90.) .AND. ( 2.*ymx - Y1(nx1/2,ny1-1) >= 90. ) ) THEN
            ny1w = ny1 + 2
            l_add_extra_j = .TRUE.
            !!
            IF ( l_first_call_interp_routine(ithrd) ) THEN
               WRITE(6,*) ''
               WRITE(6,*) '  ------ W A R N I N G ! ! ! ------'
               WRITE(6,*) ' *** your source grid is regular and seems to include the north pole.'
               WRITE(6,*) '     => yet the highest latitude in the latitude array is ', ymx
               WRITE(6,*) '     => will generate and use two extra upper J row where lat=90 and lat=90+dy on all 2D source arrays !!! '
            END IF
         END IF
      END IF

      ALLOCATE ( X1w(nx1,ny1w) , Y1w(nx1,ny1w) , Z1w(nx1,ny1w) )


      !! Is source grid supposed to include North pole?
      !! - if yes, and suppose that the highest latitude "lat_max" on the grid is
      !!   something like 89. or 89.5 then we need to add an extra upper-row for the latitude lat=90 !
      !!  => so for each lonwe interpolate X(lon,90) = X(lon+180,lat_max)
      IF ( l_add_extra_j ) THEN
         CALL EXT_NORTH_TO_90_REGG( X1, Y1, Z1,  X1w, Y1w, Z1w )
      ELSE
         X1w(:,:) = X1(:,:)
         Y1w(:,:) = Y1(:,:)
         Z1w(:,:) = Z1(:,:)
      END IF

      DEALLOCATE ( X1, Y1 )

      ctype = '00'
      ctype = TEST_XYZ(X20, Y20, Z2)
      nx2 = SIZE(Z2,1)
      ny2 = SIZE(Z2,2)

      ALLOCATE ( X2(nx2,ny2) , Y2(nx2,ny2) , mask_ignore_trg(nx2,ny2) ) !, msk_res(nx2,ny2) )

      mask_ignore_trg(:,:) = 1
      IF ( PRESENT(mask_domain_trg) ) mask_ignore_trg(:,:) = mask_domain_trg(:,:)

      IF ( ctype == '1d' ) THEN
         DO jj=1, ny2
            X2(:,jj) = X20(:,1)
         END DO
         DO ji=1, nx2
            Y2(ji,:) = Y20(:,1)
         END DO
      ELSE
         X2 = X20
         Y2 = Y20
      END IF

      IF ( l_first_call_interp_routine(ithrd) ) THEN

         l_last_y_row_missing = .FALSE.

         IF( .NOT. l_skip_bilin_mapping ) CALL MAPPING_BL(k_ew_per, X1w, Y1w, X2, Y2,  ithread=ithrd, pmsk_dom_trg=mask_ignore_trg)

      END IF

      Z2(:,:) = rflg ! Flagging non-interpolated output points

      mask_ignore_trg(:,:) = 1
      WHERE ( bilin_map(io1(ithrd):io2(ithrd),:)%jip < 1 ) mask_ignore_trg = 0
      WHERE ( bilin_map(io1(ithrd):io2(ithrd),:)%jjp < 1 ) mask_ignore_trg = 0

      !WHERE ( (IMETRICS(:,:,3 < 1) ) mask_ignore_trg = 0 ; ! iqd => problem in interp ORCA2->ORCA1 linked to iqd < 1


      bilin_map(io1(ithrd):io2(ithrd),:)%jip = MAX( bilin_map(io1(ithrd):io2(ithrd),:)%jip , 1 )  ! so no i or j <= 0
      bilin_map(io1(ithrd):io2(ithrd),:)%jjp = MAX( bilin_map(io1(ithrd):io2(ithrd),:)%jjp , 1 )  ! so no i or j <= 0

      DO jj=1, ny2
         DO ji=1, nx2
            iP    = bilin_map(ji+io1(ithrd)-1,jj)%jip
            jP    = bilin_map(ji+io1(ithrd)-1,jj)%jjp
            iqd = bilin_map(ji+io1(ithrd)-1,jj)%iqdrn
            alpha = bilin_map(ji+io1(ithrd)-1,jj)%ralfa
            beta  = bilin_map(ji+io1(ithrd)-1,jj)%rbeta
            !!
            IF ( (ABS(degE_to_degWE(X1w(iP,jP))-degE_to_degWE(X2(ji,jj)))<1.E-5) .AND. (ABS(Y1w(iP,jP)-Y2(ji,jj))<1.E-5) ) THEN
               !! COPY:
               IF (iverbose>0) WRITE(6,*) ' *** BILIN_2D: "identical point" detected (crit: 1.E-5) => copying value, no interpolation!'
               Z2(ji,jj) = Z1w(iP,jP)
            ELSE
               !! INTERPOLATION:
               Z2(ji,jj) = INTERP_BL(k_ew_per, iP, jP, iqd, alpha, beta, Z1w)
            END IF
         END DO
      END DO

      Z2 = Z2*REAL(mask_ignore_trg, 4) + REAL(1-mask_ignore_trg, 4)*(-9995.) ! masking problem points as in mask_ignore_trg


      IF ( l_first_call_interp_routine(ithrd) ) THEN
         !! Is the very last Y row fully masked! lolo and on a ORCA grid!!!
         IF ( i_orca_trg >= 4 ) THEN
            rmeanv = SUM(Z2(:,ny2))/nx2
            l_last_y_row_missing = ( (rmeanv < rflg + 0.1).AND.(rmeanv > rflg - 0.1) )
         END IF
      END IF

      !WRITE(6,*) ' l_last_y_row_missing =>', l_last_y_row_missing
      !IF ( i_orca_trg == 4 ) WRITE(6,*) ' Target grid is an ORCA grid with north-pole T-point folding!'
      !IF ( i_orca_trg == 6 ) WRITE(6,*) ' Target grid is an ORCA grid with north-pole F-point folding!'

      !! Correcting last missing band if relevant: LOLO: should use lbc_lnk no ????
      IF ( l_last_y_row_missing ) THEN
         IF ( i_orca_trg == 4 ) THEN
            Z2(2:nx2/2           ,ny2)   = Z2(nx2:nx2-nx2/2-2:-1,ny2-2)
            Z2(nx2:nx2-nx2/2-2:-1,ny2)   = Z2(2:nx2/2           ,ny2-2)
         END IF
         IF ( i_orca_trg == 6 ) THEN
            Z2(2:nx2/2             ,ny2) = Z2(nx2-1:nx2-nx2/2+1:-1,ny2-1)
            Z2(nx2-1:nx2-nx2/2+1:-1,ny2) = Z2(2:nx2/2             ,ny2-1)
         END IF
      END IF

      DEALLOCATE ( X1w, Y1w, Z1w, X2, Y2, mask_ignore_trg )

      l_first_call_interp_routine(ithrd) = .FALSE.

   END SUBROUTINE BILIN_2D


   FUNCTION INTERP_BL(k_ew_per, jiP, jjP, iqd, xa, xb, Z_in)

      IMPLICIT none

      INTEGER,                 INTENT(in) :: k_ew_per
      INTEGER,                 INTENT(in) :: jiP, jjP, iqd
      REAL(8),                 INTENT(in) :: xa, xb
      REAL(4), DIMENSION(:,:), INTENT(in) :: Z_in

      REAL(4) :: INTERP_BL
      REAL(4) ::  wup, w1, w2, w3, w4
      INTEGER  :: nxi, nyi, jiPm1, jiPp1, &
         &        i1=0, j1=0, i2=0, j2=0, i3=0, j3=0, i4=0, j4=0

      !! Choose the 4 interpolation points, according to sector and nearest point (jiP, jjP)

      !!   o<--o        x<--o         o<--x         o<--o
      !! 1 |   ^ NE   2 |   ^ SE    3 |   ^ SW    4 |   ^ NW
      !!   v   |        v   |         v   |         v   |
      !!   x-->o        o-->o         o-->o         o-->x

      nxi = SIZE(Z_in,1)
      nyi = SIZE(Z_in,2)

      jiPm1 = jiP-1
      jiPp1 = jiP+1
      IF ( (jiPm1 ==   0  ).AND.(k_ew_per>=0) )  jiPm1 = nxi - k_ew_per
      IF ( (jiPp1 == nxi+1).AND.(k_ew_per>=0) )  jiPp1 = 1   + k_ew_per


      SELECT CASE (iqd)

      CASE (1)  ! nearest point is the bottom left corner point of local mesh
         i1=jiP   ; j1 = jjP  ! local mesh is located NE of nearest point
         i2=jiPp1 ; j2 = jjP
         i3=jiPp1 ; j3 = jjP+1
         i4=jiP   ; j4 = jjP+1

      CASE (2)  ! nearest point is the top left corner point of mesh
         i1=jiP   ; j1 = jjP    ! local mesh is located SE of nearest point
         i2=jiP   ; j2 = jjP-1
         i3=jiPp1 ; j3 = jjP-1
         i4=jiPp1 ; j4 = jjP

      CASE (3)  ! nearest point is the top righ corner point of mesh
         i1=jiP   ; j1 = jjP   ! local mesh is located SW of nearest point
         i2=jiPm1 ; j2 = jjP
         i3=jiPm1 ; j3 = jjP-1
         i4=jiP   ; j4 = jjP-1

      CASE (4)  ! nearest point is the bottom right corner point of mesh
         i1=jiP   ; j1 = jjP  ! local mesh is located NW of nearest point
         i2=jiP   ; j2 = jjP+1
         i3=jiPm1 ; j3 = jjP+1
         i4=jiPm1 ; j4 = jjP

      END SELECT

      !! compute sum weight above target point
      w1=REAL( (1. - xa)*(1. - xb) , 4)
      w2=REAL(       xa *(1. - xb) , 4)
      w3=REAL(       xa * xb       , 4)
      w4=REAL( (1. - xa)* xb       , 4)

      wup = w1 + w2 + w3 + w4

      !IF ( (i1==0).OR.(j1==0).OR.(i2==0).OR.(j2==0).OR.(i3==0).OR.(j3==0).OR.(i4==0).OR.(j4==0) ) THEN
      !   WRITE(6,*) ' WARNING: INTERP_BL => at least one of the i,j index is zero!'
      !END IF

      ! interpolate with non-masked  values, above target point

      IF ( wup == 0. ) THEN
         INTERP_BL = -9998.
      ELSEIF ( (i1<1).OR.(i2<1).OR.(i3<1).OR.(i4<1) ) THEN
         INTERP_BL = -9997.
      ELSEIF ( (j1<1).OR.(j2<1).OR.(j3<1).OR.(j4<1) ) THEN
         INTERP_BL = -9996.
      ELSEIF ( (j1>nyi).OR.(j2>nyi).OR.(j3>nyi).OR.(j4>nyi) ) THEN
         INTERP_BL = -9995.
      ELSEIF ( (i1>nxi).OR.(i2>nxi).OR.(i3>nxi).OR.(i4>nxi) ) THEN
         INTERP_BL = -9994.
      ELSE
         INTERP_BL = ( Z_in(i1,j1)*w1 + Z_in(i2,j2)*w2 + Z_in(i3,j3)*w3 + Z_in(i4,j4)*w4 )/wup
      ENDIF

   END FUNCTION INTERP_BL




   SUBROUTINE MAPPING_BL(k_ew_per, pX, pY, plon_trg, plat_trg,  ithread, pmsk_dom_trg)

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
      USE io_ezcdf
      USE mod_poly, ONLY : L_InPoly
      !!
      INTEGER,                 INTENT(in) :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in) :: pX, pY
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

      INTEGER(4), DIMENSION(:,:), ALLOCATABLE :: i_nrst_src, j_nrst_src

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
      TYPE(bln_map), DIMENSION(:,:),   ALLOCATABLE :: pbln_map  ! local array on thread !
      !!
      INTEGER(1), DIMENSION(:,:),   ALLOCATABLE :: mask_ignore_trg

      REAL(8) :: alpha, beta
      LOGICAL :: l_ok, lagain, lpdebug
      INTEGER :: icpt, ithrd
      
      ithrd = 1 ! no OpenMP !
      IF( PRESENT(ithread) ) ithrd = ithread
      
      nxi = size(pX,1)
      nyi = size(pX,2)

      nxo = size(plon_trg,1)
      nyo = SIZE(plon_trg,2)

      ALLOCATE ( pbln_map(nxo,nyo), mask_ignore_trg(nxo,nyo), i_nrst_src(nxo, nyo), j_nrst_src(nxo, nyo) )
      pbln_map(:,:)%jip    = 0
      pbln_map(:,:)%jjp    = 0
      pbln_map(:,:)%iqdrn  = 0
      pbln_map(:,:)%ralfa  = 0.
      pbln_map(:,:)%rbeta  = 0.
      pbln_map(:,:)%ipb    = 0
      mask_ignore_trg(:,:) = 1
      i_nrst_src(:,:)      = 0
      j_nrst_src(:,:)      = 0

      IF ( PRESENT(pmsk_dom_trg) ) mask_ignore_trg(:,:) = pmsk_dom_trg(:,:)

      CALL FIND_NEAREST_POINT( plon_trg, plat_trg, pX, pY, i_nrst_src, j_nrst_src,  ithread=ithrd, pmsk_dom_trg=mask_ignore_trg )

      DO jj = 1, nyo
         DO ji = 1, nxo

            lpdebug = ( (iverbose==2).AND.(ji==idb).AND.(jj==jdb) )

            IF(lpdebug) WRITE(6,*) ' *** LOLO debug:', idb, jdb

            IF ( mask_ignore_trg(ji,jj)==1 ) THEN
               !! => exclude regions that do not exist on source domain (mask_ignore_trg==0) and
               !! points for which the nearest point was not found (mask_ignore_trg==-1 or -2)

               !! Now deal with horizontal interpolation
               !! set longitude of input point in accordance with lon ( [lon0, 360+lon0 [ )
               xP = plon_trg(ji,jj)
               yP = plat_trg(ji,jj)

               iP = i_nrst_src(ji,jj)
               jP = j_nrst_src(ji,jj)

               IF(lpdebug) WRITE(6,*) ' *** LOLO debug: xP, yP =', xP, yP
               IF(lpdebug) WRITE(6,*) ' *** LOLO debug: iP, jP, nxi, nyi =', iP, jP, nxi, nyi

               IF ( (iP/=INT(rflg)).AND.(jP/=INT(rflg)).AND.(jP<nyi) ) THEN   ! jP<ny1 < last upper row was an extrapolation!

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

                  IF(lpdebug) WRITE(6,*) ' *** LOLO debug: iPm1, iPp1 =', iPm1, iPp1
                  IF(lpdebug) WRITE(6,*) ' *** LOLO debug: jPm1, jPp1 =', jPm1, jPp1

                  IF ((iPm1 < 1).OR.(jPm1 < 1).OR.(iPp1 > nxi)) THEN
                     IF(iverbose>0) WRITE(6,*) 'WARNING: mod_bilin_2d.f90 => bound problem => ',xP,yP,nxi,nyi,iP,jP
                     IF(iverbose>0) WRITE(6,*) '          iPm1, iPp1, nxi =', iPm1, iPp1, nxi
                     IF(iverbose>0) WRITE(6,*) '          jPm1, jPp1, nyi =', jPm1, jPp1, nyi
                     IF(iverbose>0) WRITE(6,*) '         => ignoring current nearest point for i,j =', ji, jj, '(of target domain)'
                     IF(iverbose>0) WRITE(6,*) ''
                  ELSE

                     lonP = MOD(pX(iP,jP)  , 360._8) ; latP = pY(iP,jP)   ! nearest point
                     lonN = MOD(pX(iP,jPp1), 360._8) ; latN = pY(iP,jPp1) ! N (grid)
                     lonE = MOD(pX(iPp1,jP), 360._8) ; latE = pY(iPp1,jP) ! E (grid)
                     lonS = MOD(pX(iP,jPm1), 360._8) ; latS = pY(iP,jPm1) ! S (grid)
                     lonW = MOD(pX(iPm1,jP), 360._8) ; latW = pY(iPm1,jP) ! W (grid)

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
                     !IF(lpdebug) WRITE(6,*) ' *** LOLO debug: first find of iqd =', iqd0

                     loni(0) = xP ;    lati(0) = yP      ! fill loni, lati for 0 = target point
                     loni(1) = lonP ;  lati(1) = latP    !                     1 = nearest point

                     IF (l_save_distance_to_np) distance_to_np(ji+io1(ithrd)-1,jj) = DISTANCE(xP, lonP, yP, latP)

                     !! Problem is that sometimes, in the case of really twisted
                     !! meshes this method screws up, iqd is not what it
                     !! shoule be, so need the following loop on the value of iqd:
                     icpt = 0
                     lagain = .TRUE.
                     DO WHILE ( lagain )

                        SELECT CASE ( iqd ) ! point 2 3 4 are counter clockwise in the respective sector
                        CASE ( 1 )
                           loni(2) = lonE ; lati(2) = latE
                           loni(3) = MOD(pX(iPp1,jPp1), 360._8) ; lati(3) = pY(iPp1,jPp1)
                           loni(4) = lonN ; lati(4) = latN
                        CASE ( 2 )
                           loni(2) = lonS ; lati(2) = latS
                           loni(3) = MOD(pX(iPp1,jPm1), 360._8) ; lati(3) = pY(iPp1,jPm1)
                           loni(4) = lonE ; lati(4) = latE
                        CASE ( 3 )
                           loni(2) = lonW ; lati(2) = latW
                           loni(3) = MOD(pX(iPm1,jPm1), 360._8) ; lati(3) = pY(iPm1,jPm1)
                           loni(4) = lonS ; lati(4) = latS
                        CASE ( 4 )
                           loni(2) = lonN ; lati(2) = latN
                           loni(3) = MOD(pX(iPm1,jPp1), 360._8) ; lati(3) = pY(iPm1,jPp1)
                           loni(4) = lonW ; lati(4) = latW
                        END SELECT

                        WHERE ( loni <= 0.0 )  loni = loni + 360._8  ! P. Mathiot: Some bug with ERA40 grid

                        !! The tests!!!
                        l_ok = L_InPoly ( loni(1:4), lati(1:4), xp, yp )    ! $$
                        IF(lpdebug) WRITE(6,*) ' *** LOLO debug: l_ok =', l_ok

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
                     CALL LOCAL_COORD(loni, lati, alpha, beta, iproblem)
                     pbln_map(ji,jj)%ipb = iproblem

                     !LOLO: mark this in ID_problem => case when iqd = iqd0 ( IF ( icpt == 5 ) ) above!!!
                     IF ( icpt == 5 ) pbln_map(ji,jj)%ipb = 2 ! IDing the screw-up from above...

                     IF (lpdebug) THEN
                        WRITE(6,*) ' *** LOLO debug :'
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
                     pbln_map(ji,jj)%ralfa = alpha
                     pbln_map(ji,jj)%rbeta = beta
                     pbln_map(ji,jj)%ipb   = 0

                  END IF ! IF ((iPm1 < 1).OR.(jPm1 < 1).OR.(iPp1 > nxi))
               END IF ! IF ( (iP/=INT(rflg)).AND.(jP/=INT(rflg)).AND.(jP<nyi) )
            END IF ! IF ( mask_ignore_trg(ji,jj)==1 )
         ENDDO
      ENDDO

      !lolo: UGLY!
      pbln_map(:,:)%jip = i_nrst_src(:,:)
      pbln_map(:,:)%jjp = j_nrst_src(:,:)

      WHERE ( (i_nrst_src == INT(rflg)).OR.(j_nrst_src == INT(rflg)) )
         pbln_map(:,:)%ralfa  = rflg
         pbln_map(:,:)%rbeta  = rflg
      END WHERE

      !! Awkwardly fixing problematic points but remembering them in ID_problem

      !! Negative values that are actually 0
      WHERE ( ((pbln_map(:,:)%ralfa < 0.).AND.(pbln_map(:,:)%ralfa > -repsilon)) ) pbln_map(:,:)%ralfa = 0.0
      WHERE ( ((pbln_map(:,:)%rbeta < 0.).AND.(pbln_map(:,:)%rbeta > -repsilon)) ) pbln_map(:,:)%rbeta = 0.0

      WHERE ( (pbln_map(:,:)%ralfa > rflg).AND.(pbln_map(:,:)%ralfa < 0.) )
         pbln_map(:,:)%ralfa = 0.5
         pbln_map(:,:)%ipb = 4
      END WHERE
      WHERE ( pbln_map(:,:)%ralfa > 1. )
         pbln_map(:,:)%ralfa = 0.5
         pbln_map(:,:)%ipb =  5
      END WHERE

      WHERE ( (pbln_map(:,:)%rbeta > rflg).AND.(pbln_map(:,:)%rbeta < 0.) )
         pbln_map(:,:)%rbeta = 0.5
         pbln_map(:,:)%ipb = 6
      END WHERE
      WHERE ( pbln_map(:,:)%rbeta > 1. )
         pbln_map(:,:)%rbeta = 0.5
         pbln_map(:,:)%ipb = 7
      END WHERE

      !! iquadran was not found:
      WHERE ( pbln_map(:,:)%iqdrn < 1 )
         pbln_map(:,:)%iqdrn = 1 ! maybe bad... but at least reported in ID_problem ...
         pbln_map(:,:)%ipb = 44
      END WHERE

      WHERE (mask_ignore_trg <= -1) pbln_map(:,:)%ipb = -1 ! Nearest point was not found by "FIND_NEAREST"
      WHERE (mask_ignore_trg ==  0) pbln_map(:,:)%ipb = -2 ! No idea if possible... #lolo
      WHERE (mask_ignore_trg <  -2) pbln_map(:,:)%ipb = -3 ! No idea if possible... #lolo


      IF( ithrd > 0 ) THEN
         !! OMP i-decomposition:
         bilin_map(io1(ithrd):io2(ithrd),:) = pbln_map(:,:)
      ELSE
         bilin_map(:,:) = pbln_map(:,:)
      END IF

      DEALLOCATE ( pbln_map, mask_ignore_trg, i_nrst_src, j_nrst_src )

   END SUBROUTINE MAPPING_BL





   SUBROUTINE LOCAL_COORD(xlam, xphi, xa, xb, ipb)

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

      IMPLICIT NONE

      !! * Arguments
      REAL(8), DIMENSION(0:4), INTENT(in)  :: xlam, xphi
      REAL(8)                , INTENT(out) :: xa, xb
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

      xa = zalpha
      xb = zbeta

      !! Problem if the 4 latitudes surrounding 'lati' are equal!
      IF ( (xphi(1)==xphi(2)).AND.(xphi(2)==xphi(3)).AND.(xphi(3)==xphi(4)) ) THEN
         xa  = 0.5
         xb  = 0.5
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
      IMPLICIT NONE
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
      IMPLICIT NONE
      !!  Arguments
      REAL(8), INTENT(in) :: plata, plona, platb, plonb
      REAL(8) :: heading

      !!  Local variables
      REAL(8)  :: zpi, zconv, angled, xa, xb, ya, yb, xb_xa, rr

      zpi=ACOS(-1._8)
      zconv=zpi/180.  ! for degree to radian conversion

      !! There is a problem if the Greenwich meridian pass between a and b
      IF (iverbose>0) PRINT *,' Plonb  Plona ' , plonb, plona
      xa=plona*zconv
      xb=plonb*zconv

      rr = MAX(ABS(tan(zpi/4.-zconv*plata/2.)), repsilon)  !lolo just to avoid FPE sometimes
      ya = -LOG(rr)

      rr = MAX(ABS(tan(zpi/4.-zconv*platb/2.)), repsilon)  !lolo just to avoid FPE sometimes
      yb = -LOG(rr)

      IF (iverbose>0) WRITE(6,*) ' xa_xb , modulo 2pi', xb-xa, MOD((xb-xa),2*zpi)
      xb_xa=MOD((xb-xa),2*zpi)

      IF ( xb_xa >=  zpi ) xb_xa = xb_xa -2*zpi
      IF ( xb_xa <= -zpi ) xb_xa = xb_xa +2*zpi
      IF (iverbose>0)  print *, 'yb -ya, xb_xa ',yb -ya , xb_xa

      angled = ATAN2(xb_xa, yb - ya)

      heading=angled*180./zpi
      IF (heading < 0) heading = heading + 360._8

   END FUNCTION heading



   SUBROUTINE P2D_MAPPING_AB(cf_out, plon, plat, pbln_map, vflag )

      USE netcdf
      USE io_ezcdf, ONLY : sherr

      CHARACTER(len=*),              INTENT(in) :: cf_out
      REAL(8),       DIMENSION(:,:), INTENT(in) :: plon, plat
      TYPE(bln_map), DIMENSION(:,:), INTENT(in) :: pbln_map
      REAL(8),                       INTENT(in) :: vflag
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

      IF ( vflag /= 0. ) THEN
         CALL sherr( NF90_PUT_ATT(id_f, id_v1, '_FillValue', INT(vflag)), crtn,cf_out,'iP   (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v2, '_FillValue', INT(vflag)), crtn,cf_out,'jP   (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v3, '_FillValue', INT(vflag)), crtn,cf_out,'iqd  (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v4, '_FillValue',     vflag ), crtn,cf_out,'alfa (masking)' )
         CALL sherr( NF90_PUT_ATT(id_f, id_v5, '_FillValue',     vflag ), crtn,cf_out,'beta (masking)' )
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
