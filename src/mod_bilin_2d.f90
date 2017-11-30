MODULE MOD_BILIN_2D

   !!-----------------------------------------------------------------
   !!
   !!         Method of interpolation : "bilinear 2D"
   !!         =======================================
   !!
   !!    L. Brodeau, September 2017
   !!    P. Mathiot, August 2010
   !!    L. BRODEAU, fall 2008
   !!    J.-M. MOLINES, P. MATHIOT , 2007
   !!
   !!-----------------------------------------------------------------

   USE mod_conf
   USE mod_manip

   IMPLICIT NONE

   !! To be read into the netcdf file only at "l_first_call_interp_routine = .TRUE."
   REAL(8), DIMENSION(:,:,:), ALLOCATABLE, SAVE :: RAB       !: alpha, beta
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: IMETRICS  !: iP, jP, iquadran at each point


   PRIVATE

   LOGICAL, PARAMETER :: ldebug  = .FALSE.

   REAL(8), PARAMETER :: repsilon = 1.E-9

   INTEGER  :: iquadran !: grid sector from 1 to 4 (clockwise, 1=NE) in wich target point
   INTEGER  :: iP, jP   !: location of nearest point on input grid

   LOGICAL, SAVE :: l_last_y_row_missing

   PUBLIC :: BILIN_2D


CONTAINS


   SUBROUTINE BILIN_2D(k_ew_per, X10, Y10, Z1, X20, Y20, Z2, cnpat)

      !!================================================================
      !!
      !! INPUT :     k_ew_per : east-west periodicity
      !!                        k_ew_per = -1  --> no periodicity
      !!                        k_ew_per >= 0  --> periodicity with overlap of k_ew_per points
      !!             X10   : 2D source longitude array (ni*nj) or (ni*1)
      !!             Y10   : 2D source latitude  array (ni*nj) or (nj*1)
      !!             Z1    : source field on source grid
      !!
      !!             X20   : 2D target longitude array (ni*nj) or (ni*1)
      !!             Y20   : 2D target latitude  array (ni*nj) or (nj*1)
      !!
      !!             cnpat : name of current configuration pattern
      !!                      -> to recognise the mapping/weight file
      !!
      !! OUTPUT :
      !!             Z2    : field extrapolated from source to target grid
      !!
      !!================================================================

      USE io_ezcdf ; !lolo?

      !! Input/Output arguments
      INTEGER,                 INTENT(in)  :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X10, Y10
      REAL(4), DIMENSION(:,:), INTENT(in)  :: Z1
      REAL(8), DIMENSION(:,:), INTENT(in)  :: X20, Y20
      REAL(4), DIMENSION(:,:), INTENT(out) :: Z2
      CHARACTER(len=*),      INTENT(in)  :: cnpat


      !! Local variables
      INTEGER :: nx1, ny1, nx2, ny2

      INTEGER, PARAMETER :: n_extd = 4    ! total source grid extension (2 points south, 2 points north)

      REAL(8) :: alpha, beta, rmeanv
      LOGICAL :: lefw
      INTEGER :: cpt, ji, jj, ni1, nj1

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: &
         &    X1, Y1, X2, Y2,   &
         &    lon_in , lat_in, Z_in

      CHARACTER(len=2)   :: ctype

      CHARACTER(LEN=400) :: cf_wght

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

      ctype = '00'
      ctype = TEST_XYZ(X20, Y20, Z2)
      nx2 = size(Z2,1) ; ny2 = size(Z2,2)
      ALLOCATE ( X2(nx2,ny2) , Y2(nx2,ny2) )


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


      !!                       S T A R T

      WRITE(cf_wght,'("sosie_mapping_",a,".nc")') trim(cnpat)


      !! Extending the source 2D domain with a frame of 2 points:
      !!    We extend initial 2D array with a frame, adding n_extd points in each
      !!    dimension This is really needed specially for preserving good east-west
      !!    perdiodicity...

      ni1 = nx1  ;   nj1 = ny1 + n_extd

      ALLOCATE ( Z_in(ni1,nj1), lon_in(ni1,nj1), lat_in(ni1,nj1) )
      CALL FILL_EXTRA_NORTH_SOUTH(X1, Y1, REAL(Z1,8), lon_in, lat_in, Z_in)

      DEALLOCATE (X1, Y1)

      IF (ldebug) THEN
         CALL PRTMASK(REAL(lon_in,4), 'lon_in_extended.nc', 'lon')
         CALL PRTMASK(REAL(lat_in,4), 'lat_in_extended.nc', 'lat')
         CALL PRTMASK(REAL(Z_in,4),     'Z_in_extended.nc',  'Z' )
      END IF
      !STOP 'LOLO! after FILL_EXTRA_NORTH_SOUTH in mod_bilin_2d.f90' ! for only testing what FILL_EXTRA_NORTH_SOUTH actually does...


      IF ( l_first_call_interp_routine ) THEN

         l_last_y_row_missing = .FALSE.

         !! Testing if the file containing weights exists or if we need to create it
         !! (the latter is very time consuming!!!
         PRINT*,'';PRINT*,'********************************************************'
         INQUIRE(FILE=cf_wght, EXIST=lefw )
         IF ( lefw ) THEN
            PRINT *, 'Mapping file ', TRIM(cf_wght), ' was found!'
            PRINT *, 'Still! Insure that this is really the one you need!!!'
            PRINT *, 'No need to build it, skipping routine MAPPING !'
         ELSE
            PRINT *, 'No mapping file found in the current directory!'
            PRINT *, 'We are going to build it: ', trim(cf_wght)
            PRINT *, 'This is very time consuming, but only needs to be done once...'
            PRINT *, 'Therefore, you should keep this file for any future interpolation'
            PRINT *, 'using the same "source-target" setup'

            CALL MAPPING(k_ew_per, lon_in, lat_in, X2, Y2, cf_wght)

         END IF
         PRINT *, ''; PRINT *, 'MAPPING OK';
         PRINT*,'********************************************************';PRINT*,'';PRINT*,''
      END IF

      Z2 = -777.0      ! Flagging non-interpolated output points

      cpt = 0

      IF ( l_first_call_interp_routine ) THEN

         ALLOCATE ( IMETRICS(nx2,ny2,3), RAB(nx2,ny2,2) )
         
         CALL RD_MAPPING_AB(cf_wght, IMETRICS, RAB)
         PRINT *, ''; PRINT *, 'Mapping and weights read into ', trim(cf_wght); PRINT *, ''

      END IF

      DO jj=1, ny2
         DO ji=1, nx2

            iP       = IMETRICS(ji,jj,1)
            jP       = IMETRICS(ji,jj,2)
            iquadran = IMETRICS(ji,jj,3)

            alpha    = RAB(ji,jj,1)
            beta     = RAB(ji,jj,2)

            IF ( (iP == INT(rflg)).OR.(jP == INT(rflg)) ) THEN
               Z2(ji,jj) = rflg ; ! masking
            ELSE
               !! INTERPOLATION !
               Z2(ji,jj) = REAL(INTERP_BL(k_ew_per, iP, jP, alpha, beta, Z_in), 4)
               !!
            END IF

         END DO
      END DO

      IF ( l_first_call_interp_routine ) THEN
         !! Is the very last Y row fully masked! lolo and on a ORCA grid!!!
         IF ( i_orca_out >= 4 ) THEN
            rmeanv = SUM(Z2(:,ny2))/nx2
            l_last_y_row_missing = ( (rmeanv < rflg + 0.1).AND.(rmeanv > rflg - 0.1) )
         END IF
      END IF

      !PRINT *, ' l_last_y_row_missing =>', l_last_y_row_missing
      !IF ( i_orca_out == 4 ) PRINT *, ' Target grid is an ORCA grid with north-pole T-point folding!'
      !IF ( i_orca_out == 6 ) PRINT *, ' Target grid is an ORCA grid with north-pole F-point folding!'

      !! Correcting last missing band if relevant:
      IF ( l_last_y_row_missing ) THEN
         IF ( i_orca_out == 4 ) THEN
            Z2(2:nx2/2           ,ny2)   = Z2(nx2:nx2-nx2/2-2:-1,ny2-2)
            Z2(nx2:nx2-nx2/2-2:-1,ny2)   = Z2(2:nx2/2           ,ny2-2)
         END IF
         IF ( i_orca_out == 6 ) THEN
            Z2(2:nx2/2             ,ny2) = Z2(nx2-1:nx2-nx2/2+1:-1,ny2-1)
            Z2(nx2-1:nx2-nx2/2+1:-1,ny2) = Z2(2:nx2/2             ,ny2-1)
         END IF
      END IF

      !! Deallocation :
      DEALLOCATE ( Z_in, lon_in, lat_in, X2, Y2)

      l_first_call_interp_routine = .FALSE.

      !!DEALLOCATE ( IMETRICS, RAB )

   END SUBROUTINE BILIN_2D






   FUNCTION INTERP_BL(k_ew_per, iP, jP, xa, xb, Z_in)

      IMPLICIT none

      INTEGER,                 INTENT(in) :: k_ew_per
      INTEGER,                 INTENT(in) :: iP, jP
      REAL(8),                 INTENT(in) :: xa, xb
      REAL(8), DIMENSION(:,:), INTENT(in) :: Z_in

      REAL(8) :: INTERP_BL, wup, w1, w2, w3, w4
      INTEGER  :: nxi, iPm1, iPp1, &
         &        i1=0, j1=0, i2=0, j2=0, i3=0, j3=0, i4=0, j4=0

      !! Choose the 4 interpolation points, according to sector and nearest point (iP, jP)

      !!   o<--o        x<--o         o<--x         o<--o
      !! 1 |   ^ NE   2 |   ^ SE    3 |   ^ SW    4 |   ^ NW
      !!   v   |        v   |         v   |         v   |
      !!   x-->o        o-->o         o-->o         o-->x

      nxi = SIZE(Z_in,1)

      iPm1 = iP-1
      iPp1 = iP+1
      IF ( (iPm1 ==   0  ).AND.(k_ew_per>=0) )  iPm1 = nxi - k_ew_per
      IF ( (iPp1 == nxi+1).AND.(k_ew_per>=0) )  iPp1 = 1   + k_ew_per


      SELECT CASE (iquadran)

      CASE (1)  ! nearest point is the bottom left corner point of local mesh
         i1=iP   ; j1 = jP  ! local mesh is located NE of nearest point
         i2=iPp1 ; j2 = jP
         i3=iPp1 ; j3 = jP+1
         i4=iP   ; j4 = jP+1

      CASE (2)  ! nearest point is the top left corner point of mesh
         i1=iP   ; j1 = jP    ! local mesh is located SE of nearest point
         i2=iP   ; j2 = jP-1
         i3=iPp1 ; j3 = jP-1
         i4=iPp1 ; j4 = jP

      CASE (3)  ! nearest point is the top righ corner point of mesh
         i1=iP   ; j1 = jP   ! local mesh is located SW of nearest point
         i2=iPm1 ; j2 = jP
         i3=iPm1 ; j3 = jP-1
         i4=iP   ; j4 = jP-1

      CASE (4)  ! nearest point is the bottom right corner point of mesh
         i1=iP   ; j1 = jP  ! local mesh is located NW of nearest point
         i2=iP   ; j2 = jP+1
         i3=iPm1 ; j3 = jP+1
         i4=iPm1 ; j4 = jP

      END SELECT

      !! compute sum weight above target point
      w1=(1 - xa)*(1 - xb)
      w2=     xa *(1 - xb)
      w3=     xa * xb
      w4=(1 - xa)* xb

      wup = w1 + w2 + w3 + w4

      IF ( (i1==0).OR.(j1==0).OR.(i2==0).OR.(j2==0).OR.(i3==0).OR.(j3==0).OR.(i4==0).OR.(j4==0) ) THEN
         PRINT *, ' WARNING: INTERP_BL => at least one of the i,j index is zero!'
      END IF
      
      ! interpolate with non-masked  values, above target point

      IF ( wup == 0. ) THEN
         INTERP_BL = -666.
      ELSEIF ( (i1==0).OR.(j1==0).OR.(i2==0).OR.(j2==0).OR.(i3==0).OR.(j3==0).OR.(i4==0).OR.(j4==0) ) THEN
         INTERP_BL = -777.
      ELSE
         INTERP_BL = ( Z_in(i1,j1)*w1 + Z_in(i2,j2)*w2 + Z_in(i3,j3)*w3 + Z_in(i4,j4)*w4 )/wup
      ENDIF
      
   END FUNCTION INTERP_BL




   SUBROUTINE MAPPING(k_ew_per, lon_in, lat_in, lon_out, lat_out, cf_w)

      !!----------------------------------------------------------------------------
      !!            ***  SUBROUTINE MAPPING  ***
      !!
      !!   ** Purpose:  Write file of position, weight need by interpolation
      !!   *  Extract of CDFTOOLS cdfweight.f90 writen by Jean Marc Molines
      !!
      !!
      !!----------------------------------------------------------------------------

      USE io_ezcdf

      IMPLICIT NONE

      INTEGER,                 INTENT(in) :: k_ew_per
      REAL(8), DIMENSION(:,:), INTENT(in) :: lon_in, lat_in
      REAL(8), DIMENSION(:,:), INTENT(in) :: lon_out, lat_out
      CHARACTER(len=400)     , INTENT(in) :: cf_w ! file containing mapping pattern


      INTEGER :: &
         &     nxi, nyi, nxo, nyo, &
         &     ji, jj,   &
         &     iP,  jP,  &
         &     iPm1, iPp1,  &
         &     iproblem

      INTEGER, DIMENSION(:,:), ALLOCATABLE :: i_nrst_in, j_nrst_in

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
      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: ZAB       !: alpha, beta
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: MTRCS  !: iP, jP, iquadran at each point
      INTEGER, DIMENSION(:,:),   ALLOCATABLE :: mask_metrics

      REAL(8) :: alpha, beta


      nxi = size(lon_in,1)
      nyi = size(lon_in,2)

      nxo = size(lon_out,1)
      nyo = size(lon_out,2)


      ALLOCATE ( ZAB(nxo,nyo,2), MTRCS(nxo,nyo,3) )

      ALLOCATE ( mask_metrics(nxo,nyo) )
      mask_metrics(:,:) = 0

      ALLOCATE ( i_nrst_in(ni_out, nj_out), j_nrst_in(ni_out, nj_out) )

      
      CALL FIND_NEAREST_POINT( lon_out, lat_out, lon_in, lat_in, i_nrst_in, j_nrst_in )
      
      DO jj = 1, nyo
         DO ji = 1, nxo

            !IF ( (ji == 1).AND.(MOD(jj,10)==0) ) PRINT *, ' *** j index =>',jj,'/',nyo

            !! Now deal with horizontal interpolation
            !! set longitude of input point in accordance with lon ( [lon0, 360+lon0 [ )
            xP = lon_out(ji,jj)
            yP = lat_out(ji,jj)

            iP = i_nrst_in(ji,jj)
            jP = j_nrst_in(ji,jj)

            IF ( (ip /= INT(rflg)).AND.(jp /= INT(rflg)) ) THEN

               iPm1 = iP-1
               iPp1 = iP+1

               IF ( iPm1 == 0 ) THEN
                  !! We are in the extended case !!!
                  IF ( k_ew_per>=0 ) iPm1 = nxi - k_ew_per
               END IF

               IF ( iPp1 == nxi+1 ) THEN
                  IF ( k_ew_per>=0 ) iPp1 = 1   + k_ew_per
               END IF

               IF ((iPm1 < 1).OR.(jP-1 < 1).OR.(iPp1 > nxi).OR.(jP+1 > nyi)) THEN !lulu
                  PRINT *, 'WARNING: mod_bilin_2d.f90 => bound problem => ',xP,yP,nxi,nyi,iP,jP
                  PRINT *, '         => ignoring current nearest point for i,j =', ji, jj, '(of target domain)'
                  PRINT *, ''
               ELSE

                  lonP = MOD(lon_in(iP,jP)  , 360._8) ; latP = lat_in(iP,jP)   ! nearest point
                  lonN = MOD(lon_in(iP,jP+1), 360._8) ; latN = lat_in(iP,jP+1) ! N (grid)
                  lonE = MOD(lon_in(iPp1,jP), 360._8) ; latE = lat_in(iPp1,jP) ! E (grid)
                  lonS = MOD(lon_in(iP,jP-1), 360._8) ; latS = lat_in(iP,jP-1) ! S (grid)
                  lonW = MOD(lon_in(iPm1,jP), 360._8) ; latW = lat_in(iPm1,jP) ! W (grid)

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

                  !! determine the sector in wich the target point is located:
                  !!  ( from 1, to 4 resp. adjacent cells NE, SE, SW, NW  of the grid)
                  !!  which mesh from the nearest point point-of-view !
                  !!
                  !!   o--o        x--o         o--x         o--o
                  !! 1 |  | NE   2 |  | SE    3 |  | SW    4 |  | NW
                  !!   x--o        o--o         o--o         o--x

                  iquadran = 4

                  ! to avoid problem with the GW meridian, pass to -180, 180 when working around GW
                  IF ( hP > 180. ) THEN
                     hPp = hP - 360._8
                  ELSE
                     hPp = hP
                  ENDIF

                  IF ( hN > hE ) hN = hN -360._8
                  IF ( hPp > hN .AND. hPp <= hE ) iquadran=1
                  IF ( hP > hE  .AND. hP <= hS )  iquadran=2
                  IF ( hP > hS  .AND. hP <= hW )  iquadran=3
                  IF ( hP > hW  .AND. hPp <= hN)  iquadran=4

                  loni(0) = xP ;    lati(0) = yP      ! fill loni, lati for 0 = target point
                  loni(1) = lonP ;  lati(1) = latP    !                     1 = nearest point

                  SELECT CASE ( iquadran ) ! point 2 3 4 are counter clockwise in the respective sector
                  CASE ( 1 )
                     loni(2) = lonE ; lati(2) = latE
                     loni(3) = MOD(lon_in(iPp1,jP+1), 360._8) ; lati(3) = lat_in(iPp1,jP+1)
                     loni(4) = lonN ; lati(4) = latN
                  CASE ( 2 )
                     loni(2) = lonS ; lati(2) = latS
                     loni(3) = MOD(lon_in(iPp1,jP-1), 360._8) ; lati(3) = lat_in(iPp1,jP-1)
                     loni(4) = lonE ; lati(4) = latE
                  CASE ( 3 )
                     loni(2) = lonW ; lati(2) = latW
                     loni(3) = MOD(lon_in(iPm1,jP-1), 360._8) ; lati(3) = lat_in(iPm1,jP-1)
                     loni(4) = lonS ; lati(4) = latS
                  CASE ( 4 )
                     loni(2) = lonN ; lati(2) = latN
                     loni(3) = MOD(lon_in(iPm1,jP+1), 360._8) ; lati(3) = lat_in(iPm1,jP+1)
                     loni(4) = lonW ; lati(4) = latW
                  END SELECT

                  WHERE ( loni <= 0.0 )  loni = loni + 360._8  ! P. Mathiot: Some bug with ERA40 grid

                  !! resolve a non linear system of equation for alpha and beta
                  !! ( the non dimensional coordinates of target point)
                  CALL LOCAL_COORD(loni, lati, alpha, beta, iproblem)
                  mask_metrics(ji,jj) = iproblem

                  !ELSE   ! point is outside the domaine, put dummy values
                  !   alpha = rflg
                  !   beta  = rflg
                  !   mask_metrics(ji,jj) = 1
                  !ENDIF


                  IF (ldebug) THEN
                     PRINT *, 'Nearest point :',lonP,  latP,  hP, hPp
                     PRINT *, 'North point :',  lonN , latN , hN
                     PRINT *, 'East  point :',  lonE , latE , hE
                     PRINT *, 'South point :',  lonS , latS , hS
                     PRINT *, 'West  point :',  lonW , latW , hW
                     PRINT *, 'iquadran =',iquadran
                     PRINT *, ''
                     PRINT *, ' Nearest 4 points :'
                     PRINT *, 'Pt 1 :',loni(1), lati(1)
                     PRINT *, 'Pt 2 :',loni(2), lati(2)
                     PRINT *, 'Pt 3 :',loni(3), lati(3)
                     PRINT *, 'Pt 4 :',loni(4), lati(4)
                     PRINT *, ''
                  END IF

                  !! Saving into arrays to be written at the end:
                  !MTRCS(ji,jj,:) = (/ iP, jP, iquadran /)
                  MTRCS(ji,jj,3) = iquadran
                  ZAB(ji,jj,:)   = (/ alpha, beta /)

               END IF !lulu

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
      !lolo.



      !! Awkwardly fixing problematic points but remembering them in mask_metrics

      !! Negative values that are actually 0
      WHERE ( ((ZAB(:,:,1) < 0.).AND.(ZAB(:,:,1) > -repsilon)) ) ZAB(:,:,1) = 0.0
      WHERE ( ((ZAB(:,:,2) < 0.).AND.(ZAB(:,:,2) > -repsilon)) ) ZAB(:,:,2) = 0.0

      WHERE ( (ZAB(:,:,1) < 0.).AND.(ZAB(:,:,1) > rflg) )
         ZAB(:,:,1) = 0.5
         mask_metrics(:,:) = -6
      END WHERE
      WHERE ( ZAB(:,:,1) > 1. )
         ZAB(:,:,1) = 0.5
         mask_metrics(:,:) =  6
      END WHERE

      WHERE ( (ZAB(:,:,2) < 0.).AND.(ZAB(:,:,2) > rflg) )
         ZAB(:,:,2) = 0.5
         mask_metrics(:,:) = mask_metrics(:,:) - 10
      END WHERE
      WHERE ( ZAB(:,:,2) > 1. )
         ZAB(:,:,2) = 0.5
         mask_metrics(:,:) = mask_metrics(:,:) + 10
      END WHERE


      !! Print metrics and weight into a netcdf file 'cf_w':
      CALL P2D_MAPPING_AB(cf_w, lon_out, lat_out, MTRCS, ZAB, rflg, mask_metrics)

      DEALLOCATE ( MTRCS, ZAB, mask_metrics )

   END SUBROUTINE MAPPING





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

      IF ((ABS(zxlam(1)-zxlam(4))>=180.).OR.(ABS(zxlam(1)-zxlam(2))) >= 180. &
         &                            .OR.(ABS(zxlam(1)-zxlam(3))  >= 180. )) THEN
         ! then we are near the 0 deg line and we must work in the frame -180 180
         WHERE ( zxlam >= 180. ) zxlam=zxlam -360._8
      ENDIF

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
         !IF ( maxval(xphi) < 89. ) THEN
         !   PRINT *, ''
         !   PRINT *, 'LOCAL_COORD of mod_bilin_2d.f90: problematic region found!!!!'
         !   PRINT *, ' --> relax, this may only occur on your masked regions...'
         !   PRINT *, ''
         !   PRINT *, ' *** input longitude:', xlam ; PRINT *, ''
         !   PRINT *, ' *** input latitude:' , xphi ; PRINT *, ''
         !END IF
         zalpha = 0.5 ; zbeta  = 0.5 ; ipb = 5
      END IF

      xa = zalpha
      xb = zbeta

      !! Problem if the 4 latitudes surrounding 'lati' are equal!
      IF ( (xphi(1)==xphi(2)).and.(xphi(2)==xphi(3)).and.(xphi(3)==xphi(4)) ) THEN
         xa  = 0.5
         xb  = 0.5
         ipb = 3
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
      !!            so that heading can be computed with loxodromy
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

      IF ( xb_xa >= zpi ) xb_xa = xb_xa -2*zpi
      IF ( xb_xa <= - zpi ) xb_xa = xb_xa +2*zpi
      IF (ldebug)  print *, 'yb -ya, xb_xa ',yb -ya , xb_xa

      angled = ATAN2(xb_xa, yb - ya)

      heading=angled*180./zpi
      IF (heading < 0) heading = heading + 360._8

   END FUNCTION heading

END MODULE MOD_BILIN_2D
