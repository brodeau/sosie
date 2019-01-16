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

   !! To be read into the netcdf file only at "l_first_call_interp_routine = .TRUE."
   REAL(8),    DIMENSION(:,:,:), ALLOCATABLE, SAVE :: RAB       !: alpha, beta
   INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: IMETRICS  !: iP, jP, iqdrn at each point
   INTEGER(2), DIMENSION(:,:),   ALLOCATABLE, SAVE :: IPB       !: problem ID


   PRIVATE

   LOGICAL, PARAMETER :: ldebug  = .FALSE.

   REAL(8), PARAMETER :: repsilon = 1.E-9

   INTEGER  :: iqdrn, iqdrn0, iqdrn_old !: grid sector from 1 to 4 (clockwise, 1=NE) in wich target point
   INTEGER  :: iP, jP   !: location of nearest point on input grid

   LOGICAL, SAVE :: l_last_y_row_missing

   PUBLIC :: BILIN_2D, MAPPING_BL, INTERP_BL


CONTAINS


   SUBROUTINE BILIN_2D(k_ew_per, X10, Y10, Z1, X20, Y20, Z2, cnpat,  mask_domain_trg)

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
      !!             cnpat : name of current configuration pattern
      !!                      -> to recognise the mapping/weight file
      !!
      !! OUTPUT :
      !!             Z2    : field extrapolated from source to target grid
      !!
      !!
      !! OPTIONAL IN:
      !!      * mask_domain_trg: ignore (dont't treat) regions of the target domain where mask_domain_trg==0 !
      !!
      !!================================================================

      USE io_ezcdf, ONLY : RD_MAPPING_AB, P2D_MAPPING_AB, TEST_XYZ  ! , DUMP_2D_FIELD

      !! Input/Output arguments
      INTEGER,                 INTENT(in)  :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X10, Y10
      REAL(4), DIMENSION(:,:), INTENT(in)  :: Z1
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X20, Y20
      REAL(4), DIMENSION(:,:), INTENT(out) :: Z2
      CHARACTER(len=*),        INTENT(in)  :: cnpat
      INTEGER(1), OPTIONAL ,DIMENSION(:,:), INTENT(in) :: mask_domain_trg

      !! Local variables
      INTEGER :: nx1, ny1, ny1w, nx2, ny2

      REAL(8) :: alpha, beta, rmeanv, ymx
      LOGICAL :: l_add_extra_j, lefw
      INTEGER :: icpt, ji, jj

      REAL(8),    DIMENSION(:,:), ALLOCATABLE :: X1, Y1, X2, Y2, X1w, Y1w
      REAL(4),    DIMENSION(:,:), ALLOCATABLE :: Z1w
      INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: mask_ignore_trg !, msk_res

      CHARACTER(len=2)   :: ctype
      CHARACTER(len=400) :: cf_wght

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
            IF ( l_first_call_interp_routine ) THEN
               PRINT *, ''
               PRINT *, '  ------ W A R N I N G ! ! ! ------'
               PRINT *, ' *** your source grid is regular and seems to include the north pole.'
               PRINT *, '     => yet the highest latitude in the latitude array is ', ymx
               PRINT *, '     => will generate and use two extra upper J row where lat=90 and lat=90+dy on all 2D source arrays !!! '
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
         X2 = X20 ; Y2 = Y20
      END IF

      WRITE(cf_wght,'("sosie_mapping_",a,".nc")') TRIM(cnpat)


      IF ( l_first_call_interp_routine ) THEN

         l_last_y_row_missing = .FALSE.

         !! Testing if the file containing weights exists or if we need to create it
         !! (2nd option might be pretty time-consuming!!!
         PRINT*,'';PRINT*,'********************************************************'
         INQUIRE(FILE=cf_wght, EXIST=lefw )
         IF ( lefw ) THEN
            PRINT *, 'Mapping file ', TRIM(cf_wght), ' was found!'
            PRINT *, 'Still! Insure that this is really the one you need!!!'
            PRINT *, 'No need to build it, skipping routine MAPPING_BL !'
         ELSE
            PRINT *, 'No mapping file found in the current directory!'
            PRINT *, 'We are going to build it: ', TRIM(cf_wght)
            PRINT *, 'This is very time consuming, but only needs to be done once...'
            PRINT *, 'Therefore, you should keep this file for any future interpolation'
            PRINT *, 'using the same "source-target" setup'
            CALL MAPPING_BL(k_ew_per, X1w, Y1w, X2, Y2, cf_wght,  mask_domain_trg=mask_ignore_trg)
         END IF
         PRINT *, ''; PRINT *, 'MAPPING_BL OK';
         PRINT*,'********************************************************';PRINT*,'';PRINT*,''

         !! We read the mapping metrics in the netcdf file (regardless of
         !! whether the mapping file was just created or not) => maybe not that
         !! smart but ensure that we saved the right stuff in the netcdf mapping
         !! file...
         ALLOCATE ( IMETRICS(nx2,ny2,3), RAB(nx2,ny2,2), IPB(nx2,ny2) )
         CALL RD_MAPPING_AB(cf_wght, IMETRICS, RAB, IPB)
         PRINT *, ''; PRINT *, 'Mapping and weights read into ', TRIM(cf_wght); PRINT *, ''
      END IF


      Z2(:,:) = rflg ! Flagging non-interpolated output points
      icpt = 0

      mask_ignore_trg(:,:) = 1
      WHERE ( (IMETRICS(:,:,1) < 1) ) mask_ignore_trg = 0
      WHERE ( (IMETRICS(:,:,2) < 1) ) mask_ignore_trg = 0

      !WHERE ( (IMETRICS(:,:,3 < 1) ) mask_ignore_trg = 0 ; ! iqdrn => problem in interp ORCA2->ORCA1 linked to iqdrn < 1 !!! LOLO

      IMETRICS(:,:,1:2) = MAX( IMETRICS(:,:,1:2) , 1 )  ! so no i or j <= 0

      DO jj=1, ny2
         DO ji=1, nx2
            iP    = IMETRICS(ji,jj,1)
            jP    = IMETRICS(ji,jj,2)
            iqdrn = IMETRICS(ji,jj,3)
            alpha = RAB(ji,jj,1)
            beta  = RAB(ji,jj,2)
            !!
            IF ( (ABS(degE_to_degWE(X1w(iP,jP))-degE_to_degWE(X2(ji,jj)))<1.E-5) .AND. (ABS(Y1w(iP,jP)-Y2(ji,jj))<1.E-5) ) THEN
               !! COPY:
               IF (ldebug) WRITE(6,*)' *** BILIN_2D: "identical point" detected (crit: 1.E-5) => copying value, no interpolation!'
               Z2(ji,jj) = Z1w(iP,jP)
            ELSE
               !! INTERPOLATION:
               Z2(ji,jj) = INTERP_BL(k_ew_per, iP, jP, iqdrn, alpha, beta, Z1w)
            END IF
         END DO
      END DO

      Z2 = Z2*REAL(mask_ignore_trg, 4) + REAL(1-mask_ignore_trg, 4)*-9995. ! masking problem points as in mask_ignore_trg


      IF ( l_first_call_interp_routine ) THEN
         !! Is the very last Y row fully masked! lolo and on a ORCA grid!!!
         IF ( i_orca_trg >= 4 ) THEN
            rmeanv = SUM(Z2(:,ny2))/nx2
            l_last_y_row_missing = ( (rmeanv < rflg + 0.1).AND.(rmeanv > rflg - 0.1) )
         END IF
      END IF

      !PRINT *, ' l_last_y_row_missing =>', l_last_y_row_missing
      !IF ( i_orca_trg == 4 ) PRINT *, ' Target grid is an ORCA grid with north-pole T-point folding!'
      !IF ( i_orca_trg == 6 ) PRINT *, ' Target grid is an ORCA grid with north-pole F-point folding!'

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

      l_first_call_interp_routine = .FALSE.

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
      !   PRINT *, ' WARNING: INTERP_BL => at least one of the i,j index is zero!'
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




   SUBROUTINE MAPPING_BL(k_ew_per, X1, Y1, lon_trg, lat_trg, cf_w,  mask_domain_trg)

      !!----------------------------------------------------------------------------
      !!            ***  SUBROUTINE MAPPING_BL  ***
      !!
      !!   ** Purpose:  Write file of position, weight need by interpolation
      !!   *  Extract of CDFTOOLS cdfweight.f90 writen by Jean Marc Molines
      !!
      !! OPTIONAL:
      !!      * mask_domain_trg: ignore (dont't treat) regions of the target domain where mask_domain_trg==0 !
      !!----------------------------------------------------------------------------

      USE io_ezcdf
      USE mod_poly, ONLY : L_InPoly

      IMPLICIT NONE

      INTEGER,                 INTENT(in) :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in) :: X1, Y1
      REAL(8), DIMENSION(:,:), INTENT(in) :: lon_trg, lat_trg
      CHARACTER(len=*)       , INTENT(in) :: cf_w ! file containing mapping pattern
      INTEGER(1), OPTIONAL ,DIMENSION(:,:), INTENT(in) :: mask_domain_trg

      LOGICAL, PARAMETER :: l_save_distance_to_np=.TRUE. !: for each point of target grid, shows the distance to the nearest point found...
      INTEGER :: &
         &     nxi, nyi, nxo, nyo, &
         &     ji, jj,   &
         &     iPm1, iPp1,  &
         &     jPm1, jPp1,  &
         &     iproblem

      INTEGER(4), DIMENSION(:,:), ALLOCATABLE :: i_nrst_in, j_nrst_in

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
      REAL(8),    DIMENSION(:,:,:), ALLOCATABLE :: ZAB       !: alpha, beta
      INTEGER(4), DIMENSION(:,:,:), ALLOCATABLE :: MTRCS  !: iP, jP, iqdrn at each point
      INTEGER(2), DIMENSION(:,:),   ALLOCATABLE :: ID_problem
      INTEGER(1), DIMENSION(:,:),   ALLOCATABLE :: mask_ignore_trg
      REAL(4),    DIMENSION(:,:),   ALLOCATABLE :: distance_to_np

      REAL(8) :: alpha, beta
      LOGICAL :: l_ok, lagain
      INTEGER :: icpt, idb, jdb

      nxi = size(X1,1)
      nyi = size(X1,2)

      nxo = size(lon_trg,1)
      nyo = size(lon_trg,2)

      ALLOCATE ( ZAB(nxo,nyo,2), MTRCS(nxo,nyo,3), ID_problem(nxo,nyo), mask_ignore_trg(nxo,nyo), &
         &       i_nrst_in(nxo, nyo), j_nrst_in(nxo, nyo) )
      ZAB(:,:,:)      = 0.0
      MTRCS(:,:,:)    = 0
      ID_problem(:,:) = 0
      mask_ignore_trg(:,:) = 1

      IF (l_save_distance_to_np) ALLOCATE ( distance_to_np(nxo, nyo) )

      IF ( PRESENT(mask_domain_trg) ) mask_ignore_trg(:,:) = mask_domain_trg(:,:)

      CALL FIND_NEAREST_POINT( lon_trg, lat_trg, X1, Y1, i_nrst_in, j_nrst_in,   mask_domain_trg=mask_ignore_trg )


      idb = 0 ! i-index of point to debug on target domain
      jdb = 0 ! j-index of point to debug on target domain

      DO jj = 1, nyo
         DO ji = 1, nxo

            IF((ji==idb).AND.(jj==jdb)) PRINT *, ' *** LOLO debug:', idb, jdb

            IF ( mask_ignore_trg(ji,jj)==1 ) THEN
               !! => exclude regions that do not exist on source domain (mask_ignore_trg==0) and
               !! points for which the nearest point was not found (mask_ignore_trg==-1 or -2)

               !! Now deal with horizontal interpolation
               !! set longitude of input point in accordance with lon ( [lon0, 360+lon0 [ )
               xP = lon_trg(ji,jj)
               yP = lat_trg(ji,jj)

               iP = i_nrst_in(ji,jj)
               jP = j_nrst_in(ji,jj)

               IF((ji==idb).AND.(jj==jdb)) PRINT *, ' *** LOLO debug: xP, yP =', xP, yP
               IF((ji==idb).AND.(jj==jdb)) PRINT *, ' *** LOLO debug: iP, jP, nxi, nyi =', iP, jP, nxi, nyi

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

                  IF((ji==idb).AND.(jj==jdb)) PRINT *, ' *** LOLO debug: iPm1, iPp1 =', iPm1, iPp1
                  IF((ji==idb).AND.(jj==jdb)) PRINT *, ' *** LOLO debug: jPm1, jPp1 =', jPm1, jPp1

                  IF ((iPm1 < 1).OR.(jPm1 < 1).OR.(iPp1 > nxi)) THEN
                     PRINT *, 'WARNING: mod_bilin_2d.f90 => bound problem => ',xP,yP,nxi,nyi,iP,jP
                     PRINT *, '          iPm1, iPp1, nxi =', iPm1, iPp1, nxi
                     PRINT *, '          jPm1, jPp1, nyi =', jPm1, jPp1, nyi
                     PRINT *, '         => ignoring current nearest point for i,j =', ji, jj, '(of target domain)'
                     PRINT *, ''
                  ELSE

                     lonP = MOD(X1(iP,jP)  , 360._8) ; latP = Y1(iP,jP)   ! nearest point
                     lonN = MOD(X1(iP,jPp1), 360._8) ; latN = Y1(iP,jPp1) ! N (grid)
                     lonE = MOD(X1(iPp1,jP), 360._8) ; latE = Y1(iPp1,jP) ! E (grid)
                     lonS = MOD(X1(iP,jPm1), 360._8) ; latS = Y1(iP,jPm1) ! S (grid)
                     lonW = MOD(X1(iPm1,jP), 360._8) ; latW = Y1(iPm1,jP) ! W (grid)

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



                     !! First attemp, and generally the good one to find iqdrn !
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
                     iqdrn = 4
                     !!
                     ! To avoid problem with the GW meridian, pass to -180, 180 when working around GW
                     hPp = degE_to_degWE(hP)
                     !!
                     IF ( hN > hE ) hN = hN -360._8
                     IF ( hPp > hN .AND. hPp <= hE ) iqdrn=1
                     IF ( hP > hE  .AND. hP <= hS )  iqdrn=2
                     IF ( hP > hS  .AND. hP <= hW )  iqdrn=3
                     IF ( hP > hW  .AND. hPp <= hN)  iqdrn=4

                     !iqdrn0    = iqdrn
                     !iqdrn_old = iqdrn
                     !IF((ji==idb).AND.(jj==jdb)) PRINT *, ' *** LOLO debug: first find of iqdrn =', iqdrn0

                     loni(0) = xP ;    lati(0) = yP      ! fill loni, lati for 0 = target point
                     loni(1) = lonP ;  lati(1) = latP    !                     1 = nearest point

                     IF (l_save_distance_to_np) distance_to_np(ji,jj) = DISTANCE(xP, lonP, yP, latP)

                     !! Problem is that sometimes, in the case of really twisted
                     !! meshes this method screws up, iqdrn is not what it
                     !! shoule be, so need the following loop on the value of iqdrn:
                     icpt = 0
                     lagain = .TRUE.
                     DO WHILE ( lagain )

                        SELECT CASE ( iqdrn ) ! point 2 3 4 are counter clockwise in the respective sector
                        CASE ( 1 )
                           loni(2) = lonE ; lati(2) = latE
                           loni(3) = MOD(X1(iPp1,jPp1), 360._8) ; lati(3) = Y1(iPp1,jPp1)
                           loni(4) = lonN ; lati(4) = latN
                        CASE ( 2 )
                           loni(2) = lonS ; lati(2) = latS
                           loni(3) = MOD(X1(iPp1,jPm1), 360._8) ; lati(3) = Y1(iPp1,jPm1)
                           loni(4) = lonE ; lati(4) = latE
                        CASE ( 3 )
                           loni(2) = lonW ; lati(2) = latW
                           loni(3) = MOD(X1(iPm1,jPm1), 360._8) ; lati(3) = Y1(iPm1,jPm1)
                           loni(4) = lonS ; lati(4) = latS
                        CASE ( 4 )
                           loni(2) = lonN ; lati(2) = latN
                           loni(3) = MOD(X1(iPm1,jPp1), 360._8) ; lati(3) = Y1(iPm1,jPp1)
                           loni(4) = lonW ; lati(4) = latW
                        END SELECT

                        WHERE ( loni <= 0.0 )  loni = loni + 360._8  ! P. Mathiot: Some bug with ERA40 grid

                        !! The tests!!!
                        l_ok = L_InPoly ( loni(1:4), lati(1:4), xp, yp )    ! $$
                        IF((ji==idb).AND.(jj==jdb)) PRINT *, ' *** LOLO debug: l_ok =', l_ok

                        IF ( (.NOT. l_ok).AND.(yP < 88.) ) THEN
                           !! Mhhh... Seems like the "heading()" approach
                           !! screwed up... i.e the point xP,yP is not into the
                           !! mesh corresponding to current iquadran, trying all
                           !! other adjacent meshes to find if it belongs to one
                           !! of them...
                           icpt  = icpt + 1
                           iqdrn = icpt
                        ELSE
                           !! Everything okay! Point is inside the the mesh
                           !! corresponding to current iquadran ! :D
                           lagain = .FALSE.
                           IF ( (icpt>0).AND.(ji==idb).AND.(jj==jdb) ) PRINT *, ' --- iquadran corrected thanks to iterative test (old,new) =>', iqdrn_old, iqdrn
                        END IF
                        IF ( icpt == 5 ) THEN
                           lagain = .FALSE. ! simply give up
                           iqdrn = iqdrn0 ! Giving what first method gave
                        END IF
                     END DO !DO WHILE ( lagain )

                     !! resolve a non linear system of equation for alpha and beta
                     !! ( the non dimensional coordinates of target point)
                     CALL LOCAL_COORD(loni, lati, alpha, beta, iproblem)
                     ID_problem(ji,jj) = iproblem

                     !LOLO: mark this in ID_problem => case when iqdrn = iqdrn0 ( IF ( icpt == 5 ) ) above!!!
                     !lolo IF ( icpt == 5 ) ID_problem(ji,jj) = 2 ! IDing the screw-up from above...

                     !IF (ldebug) THEN
                     IF ((ji==idb).AND.(jj==jdb)) THEN
                        PRINT *, ' *** LOLO debug :'
                        PRINT *, 'Nearest point :',lonP,  latP,  hP, hPp
                        PRINT *, 'North point :',  lonN , latN , hN
                        PRINT *, 'East  point :',  lonE , latE , hE
                        PRINT *, 'South point :',  lonS , latS , hS
                        PRINT *, 'West  point :',  lonW , latW , hW
                        PRINT *, 'iqdrn =',iqdrn
                        PRINT *, ''
                        PRINT *, ' Nearest 4 points :'
                        PRINT *, 'Pt 1 :',loni(1), lati(1)
                        PRINT *, 'Pt 2 :',loni(2), lati(2)
                        PRINT *, 'Pt 3 :',loni(3), lati(3)
                        PRINT *, 'Pt 4 :',loni(4), lati(4)
                        PRINT *, ''
                     END IF

                     !! Saving into arrays to be written at the end:
                     !MTRCS(ji,jj,:) = (/ iP, jP, iqdrn /)
                     MTRCS(ji,jj,3) = iqdrn
                     ZAB(ji,jj,:)   = (/ alpha, beta /)

                  END IF ! IF ((iPm1 < 1).OR.(jPm1 < 1).OR.(iPp1 > nxi).OR.(jPp1 > nyi))

               END IF

            END IF
         ENDDO
      ENDDO

      !lolo: UGLY!
      MTRCS(:,:,1) = i_nrst_in(:,:)
      MTRCS(:,:,2) = j_nrst_in(:,:)

      WHERE ( (i_nrst_in == INT(rflg)).OR.(j_nrst_in == INT(rflg)) )
         ZAB(:,:,1) = rflg
         ZAB(:,:,2) = rflg
      END WHERE

      !! Awkwardly fixing problematic points but remembering them in ID_problem

      !! Negative values that are actually 0
      WHERE ( ((ZAB(:,:,1) < 0.).AND.(ZAB(:,:,1) > -repsilon)) ) ZAB(:,:,1) = 0.0
      WHERE ( ((ZAB(:,:,2) < 0.).AND.(ZAB(:,:,2) > -repsilon)) ) ZAB(:,:,2) = 0.0

      WHERE ( (ZAB(:,:,1) > rflg).AND.(ZAB(:,:,1) < 0.) )
         ZAB(:,:,1) = 0.5
         ID_problem(:,:) = 4
      END WHERE
      WHERE ( ZAB(:,:,1) > 1. )
         ZAB(:,:,1) = 0.5
         ID_problem(:,:) =  5
      END WHERE

      WHERE ( (ZAB(:,:,2) > rflg).AND.(ZAB(:,:,2) < 0.) )
         ZAB(:,:,2) = 0.5
         ID_problem(:,:) = 6
      END WHERE
      WHERE ( ZAB(:,:,2) > 1. )
         ZAB(:,:,2) = 0.5
         ID_problem(:,:) = 7
      END WHERE

      !! iquadran was not found:
      WHERE ( MTRCS(:,:,3) < 1 )
         MTRCS(:,:,3) = 1 ! maybe bad... but at least reported in ID_problem ...
         ID_problem(:,:) = 1
      END WHERE

      WHERE (mask_ignore_trg <= -1) ID_problem = -1 ! Nearest point was not found by "FIND_NEAREST"
      WHERE (mask_ignore_trg ==  0) ID_problem = -2 ! No idea if possible... #lolo
      WHERE (mask_ignore_trg <  -2) ID_problem = -3 ! No idea if possible... #lolo

      !! Print metrics and weight into a netcdf file 'cf_w':
      CALL P2D_MAPPING_AB(cf_w, lon_trg, lat_trg, MTRCS, ZAB, rflg, ID_problem,  d2np=distance_to_np)

      DEALLOCATE ( i_nrst_in, j_nrst_in, MTRCS, ZAB, ID_problem, mask_ignore_trg )
      IF ( l_save_distance_to_np ) DEALLOCATE ( distance_to_np )

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
      IF ( ldebug) print *,' Plonb  Plona ' , plonb, plona
      xa=plona*zconv
      xb=plonb*zconv

      rr = MAX(ABS(tan(zpi/4.-zconv*plata/2.)), repsilon)  !lolo just to avoid FPE sometimes
      ya = -LOG(rr)

      rr = MAX(ABS(tan(zpi/4.-zconv*platb/2.)), repsilon)  !lolo just to avoid FPE sometimes
      yb = -LOG(rr)

      IF (ldebug) PRINT *,' xa_xb , modulo 2pi', xb-xa, MOD((xb-xa),2*zpi)
      xb_xa=MOD((xb-xa),2*zpi)

      IF ( xb_xa >=  zpi ) xb_xa = xb_xa -2*zpi
      IF ( xb_xa <= -zpi ) xb_xa = xb_xa +2*zpi
      IF (ldebug)  print *, 'yb -ya, xb_xa ',yb -ya , xb_xa

      angled = ATAN2(xb_xa, yb - ya)

      heading=angled*180./zpi
      IF (heading < 0) heading = heading + 360._8

   END FUNCTION heading

END MODULE MOD_BILIN_2D
