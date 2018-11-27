!! Taken from CDFTOOLS !
!!   => https://github.com/meom-group/CDFTOOLS
!!
!!
MODULE MOD_POLY
   !!======================================================================
   !!                     ***  MODULE  MOD_POLY  ***
   !! Determine if a given point is within a polygon or not. This module is
   !! inherited from de finite element mesh generator program (TRIGRID)
   !!=====================================================================
   !! History :  2.1  : 03/2006  : J.M. Molines : Port from trigrid
   !!            3.0  : 12/2010  : J.M. Molines : Doctor norm + Licence
   !!            4.0  : 03/2017  : J.M. Molines
   !!
   !!          Use algorithms developped in the late 80's for a finite element
   !!          mesh generator (TRIGRID) by R. Walters, C. Werner et Al.
   !!          Some original comments are maintained for references.
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   routines      : description
   !!   ReadPoly      : Read polygon file
   !!   PrepPoly      : Compute elements of the sides of polygon
   !!   InPoly        : Check if a point in in or out of a polygon
   !!   PointSlope    : Internal routine which computes slopes and coeff, of the side
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   PRIVATE

   INTEGER(4), PUBLIC , PARAMETER :: jpvert = 4  !: Number of vertex per polygon

   REAL(4), DIMENSION(jpvert+1) :: vertx, verty  ! 2dim. array of polygons and their X,Y  coordinates
   REAL(4)                      :: rmaxx, rmaxy  ! max x,y of polygon coordinates
   REAL(4)                      :: rminx, rminy  ! min x,y of polygon coordinates
   REAL(8), DIMENSION(jpvert)   :: slope         ! slope of the sides of polygone
   REAL(8), DIMENSION(jpvert)   :: ra, rb, rc    ! equation of side of polygon

   !PUBLIC  :: ReadPoly
   PUBLIC  :: L_InPoly
   PRIVATE :: PointSlope

   !!----------------------------------------------------------------------
   !! CDFTOOLS_4.0 , MEOM 2017
   !! $Id$
   !! Copyright (c) 2017, J.-M. Molines
   !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
   !! @class system
   !!----------------------------------------------------------------------

CONTAINS

   FUNCTION L_InPoly ( vplon, vplat, pxpoint, pypoint )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE InPoly  ***
      !!
      !! ** Purpose :  To see if a point is inside or outside of the specified
      !!        polygon.
      !!
      !! ** Method  :  Use the equation of the side of the polygon to determine
      !!               if a point is in (L_InPoly = true) or out (L_InPoly = false)
      !!
      !! References : Trigrid
      !!----------------------------------------------------------------------
      REAL(8), DIMENSION(4), INTENT(in) :: vplon, vplat
      REAL(8),               INTENT(in) :: pxpoint, pypoint  ! Position to check
      LOGICAL      :: L_InPoly
      !!
      INTEGER(4) :: ji, jj
      INTEGER(4) :: icross, inumvert
      REAL(8)    :: zxpt, zx, zy, zevenodd
      !!----------------------------------------------------------------------
      rmaxx = MAXVAL(vplon)
      rmaxy = MAXVAL(vplat)
      rminx = MINVAL(vplon)
      rminy = MINVAL(vplat)

      IF ( (rminx < 0.).OR.(pxpoint < 0.) ) STOP 'ERROR (L_InPoly of mod_poly.f90): we only expect positive longitudes!'

      IF ( (rminx < 10.).AND.(rmaxx > 350.) ) THEN
         !! Need to reorganize in frame -180. -- +180.:
         zx  =           SIGN(1.d0,180.-pxpoint )*MIN(pxpoint ,ABS(pxpoint -360.))
         vertx(1:4) = REAL( SIGN(1.d0,180.-vplon   )*MIN(vplon ,  ABS(  vplon -360.)) , 4 )
         rmaxx = MAXVAL(vertx(1:4))
         rminx = MINVAL(vertx(1:4))
      ELSE
         zx         = pxpoint
         vertx(1:4) = REAL( vplon , 4)
      END IF


      zy         = pypoint
      verty(1:4) = REAL( vplat(:) , 4 )




      ! Automatically close the polygon
      vertx(jpvert+1)=vertx(1)
      verty(jpvert+1)=verty(1)
      ! add dummy 0.001 to integer vertex coordinates... to avoid singular problem
      DO jj=1, jpvert+1
         IF ( (vertx(jj) - INT( vertx(jj) ) ) == 0 ) vertx(jj) = vertx(jj)+0.001
         IF ( (verty(jj) - INT( verty(jj) ) ) == 0 ) verty(jj) = verty(jj)+0.001
      END DO

      !PRINT *, '' ; PRINT *, ''
      !PRINT *, ' vertx = ', vertx
      !PRINT *, ''
      !PRINT *, ' verty = ', verty
      !PRINT *, ''

      inumvert = 4 !lolo
      !!
      DO ji = 1, inumvert-1
         CALL PointSlope ( slope(ji), vertx(ji), vertx(ji+1),     & !
            &                         verty(ji), verty(ji+1),       &
            &                     ra(ji), rb(ji), rc(ji) )
      END DO
      !       - ( ji = 1, inumvert-1 )
      CALL PointSlope ( slope(inumvert), vertx(inumvert), vertx(1), & !
         &                               verty(inumvert), verty(1), &
         &                        ra(inumvert), rb(inumvert), rc(inumvert) )

      !PRINT *, ' slope=', slope

      !
      !     - get the number of cross with the polygon boundary
      icross = 0
      !     - see if point falls in the max and min range of polygon
      IF ( zx <=  rmaxx ) THEN
         IF ( zx >=  rminx ) THEN
            IF ( zy <=  rmaxy ) THEN
               IF ( zy >=  rminy ) THEN
                  !               - step through the polygon boundaries
                  DO ji = 1, inumvert
                     !                 - see if slope = 9999 and if point is on same y axis
                     IF ( slope(ji) ==  9999 ) THEN
                        IF ( zx >=  vertx(ji) ) THEN
                           IF ( ji ==  inumvert ) THEN
                              IF (        ( (zy <=  verty(inumvert) ) .AND.    &
                                 &        (zy >   verty(1)        ) ) .OR.   &
                                 &      ( (zy >=  verty(inumvert) ) .AND.    &
                                 &        (zy <   verty(1)        ) ) ) THEN
                                 !                         - it has crossed the polygon boundary
                                 icross = icross + 1
                                 !                              if (zy == 398) print *, zx, zy, icross ,'A', ji
                              ENDIF ! ( zy test )
                           ELSEIF (       ( (zy <=  verty(ji)   ) .AND.        &
                              &           (zy >   verty(ji+1) ) ) .OR.       &
                              &         ( (zy >=  verty(ji)   ) .AND.        &
                              &           (zy <   verty(ji+1) ) ) ) THEN
                              !                       - it has crossed the polygon boundary
                              icross = icross + 1
                              !                              if (zy == 398) print *, zx, zy, icross,'B', ji
                           ENDIF !    ( ji = inumvert )
                        ENDIF   !    ( zx >= vertx(ji) )
                        !                   - see if normal slope (+ or -), and if point is not
                        !                    - higher or lower than y endpoints of the vertices
                     ELSEIF ( slope(ji) .NE. 0 ) THEN
                        zxpt = ( rc(ji) + zy ) / ra(ji)
                        IF ( ji ==  inumvert ) THEN
                           IF (            ( (zxpt <=  vertx(inumvert) ) .AND.   &
                              &            (zxpt >   vertx(1)        ) ) .OR.  &
                              &          ( (zxpt >=  vertx(inumvert) ) .AND.   &
                              &            (zxpt <   vertx(1)        ) ) ) THEN
                              IF ( zx >=  zxpt) THEN
                                 !                          - it has crossed the polygon boundary
                                 icross = icross + 1
                                 !                              if (zy == 398) print *, zx, zy, icross,'C', ji
                              ENDIF ! ( zx >= zxpt )
                           ENDIF !  ( zxpt test )
                        ELSEIF (            ( (zxpt <=  vertx(ji)   ) .AND.     &
                           &                (zxpt >   vertx(ji+1) ) ) .OR.    &
                           &              ( (zxpt >=  vertx(ji)   ) .AND.     &
                           &                (zxpt <   vertx(ji+1) ) ) ) THEN
                           IF ( zx >=   zxpt ) THEN
                              !                       - it has crossed the polygon boundary
                              icross = icross + 1
                              !                              if (zy == 398) print *, zx, zy, icross,'D', ji, slope(ji), zxpt
                           ENDIF ! ( zx >= zxpt )
                        ENDIF ! ( ji = inumvert )
                     ENDIF !  ( zxpt test )
                  END DO !  ( ji = 1, inumvert )
                  !          - decide how many times scanline crossed poly bounds
                  zevenodd = AMOD ( ( icross * 1.0 ), 2.0 )
                  IF ( zevenodd .NE. 0 ) THEN
                     !            - point is in polygon
                     L_InPoly = .TRUE.
                  ELSE
                     L_InPoly = .FALSE.
                  ENDIF
                  !            - ( zevenodd ne 0 )
               ELSE
                  L_InPoly = .FALSE.
               ENDIF
               !          - ( zy >= rminy )
            ELSE
               L_InPoly = .FALSE.
            ENDIF
            !           - ( zy <= rmaxy )
         ELSE
            L_InPoly = .FALSE.
         ENDIF
         !         - ( zx >= rminx )
      ELSE
         L_InPoly = .FALSE.
      ENDIF
      !       - ( zx <= rmaxx )

   END FUNCTION L_InPoly


   SUBROUTINE PointSlope ( pslup, pvertxa, pvertxb, pvertya, pvertyb, pax, pby, pcnstnt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE PointSlope  ***
      !!
      !! ** Purpose :   To get the slope and general equations of lines.
      !!
      !! ** Method  :  GIVEN: vertxa, vertxb, vertya, vertyb = endpoints of line section
      !!                            to operate on.
      !!             RETURNS: slup = slope of the line section
      !!              ax, by, cnstnt = general eqation of the line section.
      !!
      !! References :    trigrid
      !!----------------------------------------------------------------------
      REAL(8), INTENT(out) :: pslup
      REAL(4),  INTENT(in) :: pvertxa, pvertxb, pvertya, pvertyb
      REAL(8), INTENT(out) :: pax, pby, pcnstnt

      REAL(8) :: zvertxa, zvertxb, zvertya, zvertyb
      REAL(8) :: zrise, zrun
      !!----------------------------------------------------------------------

      zvertxa = pvertxa  ; zvertxb = pvertxb
      zvertya = pvertya  ; zvertyb = pvertyb

      zrise = zvertyb - zvertya
      zrun  = zvertxb - zvertxa

      IF ( zrun ==  0 ) THEN
         pslup = 9999
      ELSE
         pslup = zrise / zrun
      ENDIF

      IF ( ABS(pslup) <=  0.001 ) THEN
         pslup = 0.0
      ENDIF

      IF ( pslup ==  0 ) THEN
         pax = pslup
         pby = 1
         pcnstnt = zvertya
      ELSEIF ( pslup ==  9999 ) THEN
         pax = 1
         pby = 0
         pcnstnt = zvertxa
      ELSE
         pax = pslup
         pby = -1
         pcnstnt = ( pslup * zvertxa - zvertya )
      ENDIF

   END SUBROUTINE PointSlope

END MODULE MOD_POLY
