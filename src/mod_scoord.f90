MODULE mod_scoord

   IMPLICIT NONE

CONTAINS

   ! -- compute depth from s coodinates parameters --
   SUBROUTINE depth_from_scoord(zs_par, zz_bathy, zn_nx, zn_ny, zn_nz, zz_depth)

      USE mod_conf, ONLY : scoord_params
      IMPLICIT NONE

      TYPE( scoord_params ),INTENT(in)                 :: zs_par
      INTEGER,INTENT(in)                               :: zn_nx, zn_ny, zn_nz
      REAL,DIMENSION(zn_nx,zn_ny),INTENT(in)           :: zz_bathy
      REAL(4),DIMENSION(zn_nx,zn_ny,zn_nz),INTENT(out) :: zz_depth

      !! local
      REAL(8),DIMENSION(zn_nz) :: Cs_r,sc_r
      REAL(8) :: hc, hwater, cff_r, cff1_r, cff2_r, hinv
      INTEGER :: ji, jj, jk

      REAL(8),DIMENSION(zn_nx,zn_ny) :: Zt_avg1 ! free surface (should be an input too)

      Zt_avg1 = 0.

      SELECT CASE (zs_par%Vtransform )

         ! Song and Haidvogel
      CASE (1)
         PRINT *, ' This transformation has not been coded' ; STOP
         ! UCLA
      CASE (2)
         PRINT *, ' Vtransform = 2'
         !!
         SELECT CASE (zs_par%Vstretching )

         CASE (1)
            PRINT *, ' This stretching has not been coded' ; STOP

         CASE (2)
            PRINT *, ' This stretching has not been coded' ; STOP

         CASE (3)
            PRINT *, ' This stretching has not been coded' ; STOP

         CASE (4)
            PRINT *, ' Vstretching = 4'
            CALL compute_scoord_2_4(zs_par,zn_nz,Cs_r,sc_r)

            hc=zs_par%Tcline

            DO jj=1,zn_ny
               DO ji=1,zn_nx
                  DO jk=1,zs_par%Nlevels

                     cff_r  = hc*sc_r(jk)
                     cff1_r = Cs_r(jk)

                     hwater=zz_bathy(ji,jj)
                     hinv=1.0/(hc+hwater)
                     cff2_r=(cff_r+cff1_r*hwater)*hinv

                     zz_depth(ji,jj,jk) = REAL( Zt_avg1(ji,jj)+(Zt_avg1(ji,jj)+hwater)*cff2_r , 4 )
                  END DO
               END DO
            END DO

         END SELECT ! Vstretching

      END SELECT ! Vtransform


   END SUBROUTINE depth_from_scoord

   ! -- compute s coodinates functions --
   SUBROUTINE compute_scoord_2_4(zs_par,zn_nz,Cs_r,sc_r)

      USE mod_conf, ONLY : scoord_params
      IMPLICIT NONE

      TYPE( scoord_params ),INTENT(in) :: zs_par
      INTEGER, INTENT(in) :: zn_nz
      REAL(8),DIMENSION(zn_nz),INTENT(out) :: Cs_r, sc_r

      !! local
      REAL :: Nlev, theta_s, theta_b, zsc_r
      INTEGER :: k, N
      REAL(8) :: ds, Csur, Cbot

      !! copy in local variables
      Nlev = REAL(zs_par%Nlevels) ; theta_s = REAL(zs_par%theta_s) ; theta_b = REAL(zs_par%theta_b)
      N = zs_par%Nlevels

      ds=1.0/Nlev

      DO k=1,N
         zsc_r=ds*(REAL(k-Nlev)-0.5)
         sc_r(k)=zsc_r
         IF (theta_s.gt.0.0) THEN
            Csur=(1.0-COSH(theta_s*zsc_r))/    &
               &           (COSH(theta_s)-1.0)
         ELSE
            Csur=-zsc_r**2
         END IF

         IF (theta_b.gt.0.0) THEN
            Cbot=(EXP(theta_b*Csur)-1.0)/      &
               &           (1.0-EXP(-theta_b))
            Cs_r(k)=Cbot
         ELSE
            Cs_r(k)=Csur
         END IF
      END DO

   END SUBROUTINE compute_scoord_2_4

END MODULE mod_scoord
