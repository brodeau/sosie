MODULE MOD_BDROWN

   !USE io_ezcdf !LOLO
   !USE mod_manip

   IMPLICIT none

   PRIVATE

   PUBLIC :: bdrown, smoother

   LOGICAL, PARAMETER :: ldebug = .FALSE.

   REAL, PARAMETER :: ris2 = 1.0/SQRT(2.0)

CONTAINS

   SUBROUTINE BDROWN(k_ew, X, mask,   nb_inc, nb_smooth)

      !!#############################################################################
      !!
      !!  PURPOSE : fill continental areas of field X (defined by mask=0)
      !!  -------   using nearest surrounding sea points (defined by mask=1)
      !!            field X is absoluletly unchanged on mask=1 points
      !!
      !!  k_ew :  east-west periodicity on the input file/grid
      !!          k_ew = -1  --> no periodicity
      !!          k_ew >= 0  --> periodicity with overlap of k_ew points
      !!
      !!  X    :  treated array                             (2D array)
      !!  mask :  land-sea mask    INTEGER !!!!             (2D array)
      !!
      !! Optional:
      !!  * nb_inc : how far in terms of number of grid-point we extrapolate sea value into land
      !!                => default: nb_inc = 400
      !!                     (will normally stop before 400 iterations, when all land points have been treated!!!)
      !!
      !!  * nb_smooth : number of times the smoother is applied on masked region (mask=0)
      !!                => default: nb_smooth = 10
      !!
      !!
      !!                       Author : Laurent BRODEAU, 2014
      !!
      !!#############################################################################

      !! Arguments :
      INTEGER,                       INTENT(in)    :: k_ew
      REAL(4),    DIMENSION(:,:),    INTENT(inout) :: X
      INTEGER(1), DIMENSION(:,:),    INTENT(in)    :: mask

      INTEGER,    OPTIONAL,          INTENT(in)    :: nb_inc, nb_smooth

      !! Local :
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:) :: maskv, mask_coast, mtmp
      REAL(4),    ALLOCATABLE, DIMENSION(:,:) :: dold, xtmp

      INTEGER :: &
         &      ninc_max,      &
         &      nsmooth_max,          &
         &      ni, nj,        &
         &      jinc,          &
         &      ji, jj, jci,   &
         &      jim, jip

      INTEGER, DIMENSION(2) :: ivi, vim_per, vip_per

      INTEGER, PARAMETER :: jinc_debg = 2


      X = X * mask  ! we rather have 0s on continents than some fucked up high values...

      ninc_max = 200   ! will stop before when all land points have been treated!!!
      IF ( present(nb_inc) ) ninc_max = nb_inc

      nsmooth_max = 10
      IF ( present(nb_smooth) ) nsmooth_max = nb_smooth


      IF ( (size(X,1) /= size(mask,1)).OR.(size(X,2) /= size(mask,2)) ) THEN
         PRINT *, 'ERROR, mod_bdrown.F90 => BDROWN : size of data and mask do not match!!!'; STOP
      END IF

      ni = SIZE(X,1)
      nj = SIZE(X,2)





      !! Backing up original mask into mask2(:,:)
      ALLOCATE ( maskv(ni,nj), dold(ni,nj), xtmp(ni,nj), mask_coast(ni,nj), mtmp(ni,nj) )



      ivi = (/ 1 , ni /)


      IF (k_ew >= 0) THEN
         vim_per = (/ ni-k_ew ,  ni-1  /)
         vip_per = (/    2    , 1+k_ew /)
      END IF

      jinc = 0
      maskv(:,:) = mask(:,:)



      DO jinc = 1, ninc_max

         !! Quiting if no land point left:
         IF ( .NOT. (ANY(maskv == 0))  ) THEN
            IF ( ldebug ) PRINT *, 'BDROWN: No land points left! Leaving incursion loop at jinc =', jinc
            EXIT
         END IF

         dold(:,:) = X(:,:)




         !! Building mask of the coast-line (belonging to land points)

         mask_coast(:,:) = 0

         mask_coast(2:ni-1,2:nj-1) = (maskv(3:ni,2:nj-1) + maskv(2:ni-1,3:nj) + maskv(1:ni-2,2:nj-1) + maskv(2:ni-1,1:nj-2)) &
            &                     *(-(maskv(2:ni-1,2:nj-1)-1))


         !! West and East boundaries with periodicity
         !! ------------------------------------------

         IF (k_ew >= 0) THEN

            DO jci = 1, 2
               jim = vim_per(jci)  ! ji-1
               ji  = ivi(jci)      ! first ji = 1, then ji = ni
               jip = vip_per(jci)  ! ji+1
               mask_coast(ji,2:nj-1) = (maskv(jip,2:nj-1) + maskv(ji,3:nj) + maskv(jim,2:nj-1) + maskv(ji,1:nj-2)) &
                  &                     *(-(maskv(ji,2:nj-1)-1))
            END DO

         ELSE

            !! West LBC:
            mask_coast(1,2:nj-1)  = (maskv(2,2:nj-1) + maskv(1,3:nj)     + maskv( 1,1:nj-2))*(-(maskv( 1,2:nj-1) -1))
            !! East LBC:
            mask_coast(ni,2:nj-1) = (maskv(ni,3:nj) + maskv(ni-1,2:nj-1) + maskv(ni,1:nj-2))*(-(maskv(ni,2:nj-1) -1))

         END IF



         ! -------
         ! jj=1
         ! -------
         mask_coast(2:ni-1,1) = (maskv(3:ni,1) + maskv(2:ni-1,2) + maskv(1:ni-2,1)) &
            &                     *(-(maskv(2:ni-1,1)-1))
         !!
         !! ji=1, jj=1
         IF (k_ew >= 0) THEN
            mask_coast(1,1) = (maskv(2,1) + maskv(1,2) + maskv(ni-k_ew,1))*(-(maskv(1,1)-1))
         ELSE
            mask_coast(1,1) = (maskv(2,1) + maskv(1,2)                   )*(-(maskv(1,1)-1))
         END IF
         ! ji=ni, jj=1
         IF (k_ew >= 0) THEN
            mask_coast(ni,1) = (maskv(1+k_ew,1) + maskv(ni,2) + maskv(ni,1))*(-(maskv(ni,1)-1))
         ELSE
            mask_coast(ni,1) = (                  maskv(ni,2) + maskv(ni,1))*(-(maskv(ni,1)-1))
         END IF

         ! jj=nj
         ! -------
         mask_coast(2:ni-1,nj) = (maskv(3:ni,nj) + maskv(1:ni-2,nj) + maskv(2:ni-1,nj-1)) &
            &                     *(-(maskv(2:ni-1,nj)-1))
         !! ji=1, jj=nj
         IF (k_ew >= 0) THEN
            mask_coast(1,nj)  = (maskv(2,nj)            + maskv(ni-k_ew,nj)  + maskv(1,nj-1)   )*(-(maskv(1,nj) -1))
         ELSE
            mask_coast(1,nj)  = (maskv(2,nj)                                 + maskv(1,nj-1)   )*(-(maskv(1,nj) -1))
         END IF
         !! ji=ni, jj=nj
         IF (k_ew >= 0) THEN
            mask_coast(ni,nj) = (maskv(k_ew+1,nj) + maskv(ni-1,nj)     + maskv(ni-1,nj-1))    *(-(maskv(ni,nj) -1))
         ELSE
            mask_coast(ni,nj) = (          maskv(ni-1,nj)     + maskv(ni-1,nj-1))             *(-(maskv(ni,nj) -1))
         END IF


         !! mask_coast is fine now
         mtmp(:,:) = mask_coast(:,:)
         mask_coast(:,:) = 0

         WHERE ( mtmp(:,:) > 0 )
            mask_coast = 1
         END WHERE


         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(REAL(maskv,4), 'maskv.nc', 'lsm')
         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(REAL(mask_coast,4), 'mask_coast.nc', 'lsm') !; STOP 'mod_bdrown.F90 => boo!'
         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(X, 'data_X_before.nc', 'lsm')
         !STOP



         !! mask_coast done, time to fill the coastline points with values from the nearest se points
         !! -----------------------------------------------------------------------------------------

         !! Center of the domain:
         DO jj = 2, nj-1
            DO ji = 2, ni-1
               IF ( mask_coast(ji,jj) == 1 ) THEN
                  X(ji,jj) = 1./(maskv(ji+1,jj)+maskv(ji,jj+1)+maskv(ji-1,jj)+maskv(ji,jj-1) + &
                     & ris2*(maskv(ji+1,jj+1)+maskv(ji-1,jj+1)+maskv(ji-1,jj-1)+maskv(ji+1,jj-1)))*( &
                     & maskv(ji+1,jj)*dold(ji+1,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                     & maskv(ji-1,jj)*dold(ji-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                     & ris2*maskv(ji+1,jj+1)*dold(ji+1,jj+1) + ris2*maskv(ji-1,jj+1)*dold(ji-1,jj+1) + &
                     & ris2*maskv(ji-1,jj-1)*dold(ji-1,jj-1) + ris2*maskv(ji+1,jj-1)*dold(ji+1,jj-1)  )
               END IF
            END DO
         END DO



         DO jci = 1, 2


            ji  = ivi(jci)      ! first ji = 1, then ji = ni


            IF (k_ew >= 0) THEN

               !! West and East boundaries with periodicity


               jim = vim_per(jci)  ! ji-1
               jip = vip_per(jci)  ! ji+1

               DO jj = 2, nj-1
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(jip,jj)+maskv(ji,jj+1)+maskv(jim,jj)+maskv(ji,jj-1) + &
                        & ris2*(maskv(jip,jj+1)+maskv(jim,jj+1)+maskv(jim,jj-1)+maskv(jip,jj-1)))*( &
                        & maskv(jip,jj)*dold(jip,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & maskv(jim,jj)*dold(jim,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & ris2*maskv(jip,jj+1)*dold(jip,jj+1) + ris2*maskv(jim,jj+1)*dold(jim,jj+1) + &
                        & ris2*maskv(jim,jj-1)*dold(jim,jj-1) + ris2*maskv(jip,jj-1)*dold(jip,jj-1)  )
                  END IF
               END DO

            ELSE

               !! West & East LBCs when not east-west periodicity, extrapolating lineraly
               IF ( ji == 1 ) THEN
                  DO jj = 2, nj-1
                     IF ( mask_coast(ji,jj) == 1 ) THEN
                        X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj+1)+maskv(ji,jj-1) + &
                           & ris2*maskv(2,jj+1)+ris2*maskv(2,jj-1))*( &
                           & maskv(2,jj)*dold(2,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                           & maskv(ji,jj-1)*dold(ji,jj-1) + &
                           & ris2*maskv(2,jj+1)*dold(2,jj+1) + &
                           & ris2*maskv(2,jj-1)*dold(2,jj-1)  )
                     END IF
                  END DO
               END IF
               IF ( ji == ni ) THEN
                  DO jj = 2, nj-1
                     IF ( mask_coast(ji,jj) == 1 ) THEN
                        X(ji,jj) = 1./( maskv(ji,jj+1)+maskv(ni-1,jj)+maskv(ji,jj-1) + &
                           & ris2*maskv(ni-1,jj+1)+ris2*maskv(ni-1,jj-1))*( &
                           & maskv(ji,jj+1)*dold(ji,jj+1) + &
                           & maskv(ni-1,jj)*dold(ni-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                           & ris2*maskv(ni-1,jj+1)*dold(ni-1,jj+1) + &
                           & ris2*maskv(ni-1,jj-1)*dold(ni-1,jj-1) )
                     END IF
                  END DO
               END IF
            END IF
         END DO





         !! Center of Top row:
         jj = nj
         DO ji = 2, ni-1
            IF ( mask_coast(ji,jj) == 1 ) THEN
               X(ji,jj) = 1./( maskv(ji+1,jj)+maskv(ji-1,jj)+maskv(ji,jj-1) + &
                  & ris2*maskv(ji-1,jj-1)+ris2*maskv(ji+1,jj-1) )*( &
                  & maskv(ji+1,jj)*dold(ji+1,jj) + &
                  & maskv(ji-1,jj)*dold(ji-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                  & ris2*maskv(ji-1,jj-1)*dold(ji-1,jj-1) + ris2*maskv(ji+1,jj-1)*dold(ji+1,jj-1)  )
            END IF
         END DO

         !! West and East corner of top row:
         DO jci = 1, 2

            ji  = ivi(jci)      ! first ji = 1, then ji = ni

            IF (k_ew >= 0) THEN
               jim = vim_per(jci)  ! ji-1
               jip = vip_per(jci)  ! ji+1
               IF ( mask_coast(ji,jj) == 1 ) THEN
                  X(ji,jj) = 1./(maskv(jip,jj)+maskv(jim,jj)+maskv(ji,jj-1) + &
                     & ris2*maskv(jim,jj-1)+ris2*maskv(jip,jj-1))*( &
                     & maskv(jip,jj)*dold(jip,jj) + &
                     & maskv(jim,jj)*dold(jim,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                     & ris2*maskv(jim,jj-1)*dold(jim,jj-1) + ris2*maskv(jip,jj-1)*dold(jip,jj-1)  )
               END IF

               ! No E-W periodicity:
            ELSE
               IF ( ji == 1 ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj-1) + &
                        & ris2*maskv(2,jj-1))*( &
                        & maskv(2,jj)*dold(2,jj) + &
                        & maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & ris2*maskv(2,jj-1)*dold(2,jj-1)  )
                  END IF
               END IF
               IF ( ji == ni ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(ni-1,jj)+maskv(ji,jj-1) + &
                        & ris2*maskv(ni-1,jj-1))*( &
                        & maskv(ni-1,jj)*dold(ni-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & ris2*maskv(ni-1,jj-1)*dold(ni-1,jj-1)  )
                  END IF
               END IF

            END IF
         END DO




         !! Center of Bottom row:
         jj = 1
         DO ji = 2, ni-1
            IF ( mask_coast(ji,jj) == 1 ) THEN
               X(ji,jj) = 1./(maskv(ji+1,jj)+maskv(ji,jj+1)+maskv(ji-1,jj) + &
                  & ris2*maskv(ji+1,jj+1)+ris2*maskv(ji-1,jj+1) )*( &
                  & maskv(ji+1,jj)*dold(ji+1,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                  & maskv(ji-1,jj)*dold(ji-1,jj) + &
                  & ris2*maskv(ji+1,jj+1)*dold(ji+1,jj+1) + ris2*maskv(ji-1,jj+1)*dold(ji-1,jj+1) )
            END IF
         END DO

         !! West and East corner of bottom row:
         DO jci = 1, 2

            ji  = ivi(jci)      ! first ji = 1, then ji = ni

            IF (k_ew >= 0) THEN
               jim = vim_per(jci)  ! ji-1
               jip = vip_per(jci)  ! ji+1
               IF ( mask_coast(ji,jj) == 1 ) THEN
                  X(ji,jj) = 1./(maskv(jip,jj)+maskv(ji,jj+1)+maskv(jim,jj) + &
                     & ris2*maskv(jip,jj+1)+ris2*maskv(jim,jj+1) )*( &
                     & maskv(jip,jj)*dold(jip,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                     & maskv(jim,jj)*dold(jim,jj) + &
                     & ris2*maskv(jip,jj+1)*dold(jip,jj+1) + ris2*maskv(jim,jj+1)*dold(jim,jj+1) )
               END IF

               !! No E-W periodicity:
            ELSE
               IF ( ji == 1 ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj+1) + &
                        & ris2*maskv(2,jj+1) )*( &
                        & maskv(2,jj)*dold(2,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & ris2*maskv(2,jj+1)*dold(2,jj+1) )
                  END IF
               END IF
               IF ( ji == ni ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(ji,jj+1)+maskv(ni-1,jj) + &
                        & ris2*maskv(ni-1,jj+1) )*( &
                        & maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & maskv(ni-1,jj)*dold(ni-1,jj) + &
                        & ris2*maskv(ni-1,jj+1)*dold(ni-1,jj+1) )
                  END IF
               END IF
            END IF
         END DO





         !! Loosing land for the next iteration:
         !! -----------------------------------
         maskv = maskv + mask_coast
         !! -----------------------------------


         !IF ( jinc == jinc_debg) THEN
         !   CALL DUMP_2D_FIELD(X, 'data_X_after.nc', 'lsm')
         !   STOP 'after lolo'
         !END IF


      END DO

      !! Time to smooth over land! (what's been drowned):
      mtmp = 1 - mask ! 1 over continents, 0 over seas!
      CALL SMOOTHER(k_ew, X,  nb_smooth=nsmooth_max, msk=mtmp)
      !! *** l_exclude_mask_points=.true. would be stupid here,
      !!       it's actually good if sea values are used and are
      !!       propagating inland in the present CASE

      !CALL DUMP_2D_FIELD(X, 'drowned_final.nc', 'lsm') ;     !     STOP 'lolo'

      DEALLOCATE ( maskv, mtmp, xtmp, dold, mask_coast )

      IF ( ldebug ) PRINT *, 'BDROWN: jinc =', jinc

   END SUBROUTINE BDROWN




   SUBROUTINE SMOOTHER(k_ew, X,  nb_smooth, msk, l_exclude_mask_points)
      !!#############################################################################
      !!
      !!  PURPOSE : Smooth a fied with a nearest-points box-car typ of smoothing
      !!
      !!  k_ew :  east-west periodicity on the input file/grid
      !!          k_ew = -1  --> no periodicity
      !!          k_ew >= 0  --> periodicity with overlap of k_ew points
      !!
      !!  X    :  treated array                             (2D array)
      !!
      !! Optional:
      !!  * nb_smooth  : number of times the smoother is applied on masked region (mask=0)
      !!                => default: nb_smooth = 10
      !!  * msk : mask array that defines where the smoothing should be applied
      !!                 => where msk==1: smoothing applies
      !!                 => where msk==0: original values of X will be preserved
      !!
      !!  * l_exclude_mask_points: if true, the smoothing process will not use any value
      !!                           from points that belong to regions where msk==0
      !!
      !!#############################################################################

      !! Arguments :
      INTEGER,                    INTENT(in)           :: k_ew
      REAL(4),    DIMENSION(:,:), INTENT(inout)        :: X
      INTEGER,    OPTIONAL,                 INTENT(in) :: nb_smooth
      INTEGER(1), OPTIONAL, DIMENSION(:,:), INTENT(in) :: msk
      LOGICAL,    OPTIONAL                , INTENT(in) :: l_exclude_mask_points

      REAL(4),    ALLOCATABLE, DIMENSION(:,:) :: xorig, xtmp, rdnm
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:) :: nbp

      INTEGER, DIMENSION(2) :: ivi, vim_per, vip_per

      LOGICAL :: l_mask, l_emp

      INTEGER :: &
         &      nsmooth_max,          &
         &      ni, nj,        &
         &      ji, jci,   &
         &      jim, jip, js

      REAL(4), PARAMETER :: &
         &  w0 = 0.35   ! weight given to the point i,j in the boxcar process


      nsmooth_max = 10
      IF ( PRESENT(nb_smooth) ) nsmooth_max = nb_smooth

      l_mask = PRESENT(msk)

      l_emp  = PRESENT(l_exclude_mask_points)

      ni = SIZE(X,1)
      nj = SIZE(X,2)

      IF ( (l_emp).AND.(.NOT. l_mask) ) THEN
         PRINT *, 'PROBLEM in SMOOTH (mod_bdrown.f90): you need to provide a "msk"'
         PRINT *, '                                    if you set l_exclude_mask_points=.true.!'
         STOP
      END IF


      ALLOCATE ( xtmp(ni,nj) )

      IF (l_emp) THEN
         ALLOCATE ( rdnm(ni,nj) , nbp(ni,nj) )
         rdnm(:,:) = 0.25
         nbp(:,:)  = 0
      END IF

      IF ( l_mask ) THEN
         ALLOCATE ( xorig(ni,nj) )
         xorig(:,:) = X(:,:)
      END IF

      ivi = (/ 1 , ni /)

      IF (k_ew >= 0) THEN
         vim_per = (/ ni-k_ew ,  ni-1  /)
         vip_per = (/    2    , 1+k_ew /)
      END IF


      DO js = 1, nsmooth_max

         xtmp(:,:) = X(:,:)

         IF ( l_emp ) xtmp(:,:) = xtmp(:,:)*REAL(msk(:,:),4)

         !! Center of the domain:
         !! ---------------------
         IF ( l_emp ) THEN
            !! -- SMOOTHER is excluding masked points...'
            nbp(2:ni-1,2:nj-1) =   msk(3:ni,2:nj-1) + msk(2:ni-1,3:nj) + msk(1:ni-2,2:nj-1) + msk(2:ni-1,1:nj-2) & ! Number of non-masked surrounded points
               &                 + msk(3:ni,3:nj)   + msk(3:ni,1:nj-2) + msk(1:ni-2,1:nj-2) + msk(1:ni-2,3:nj)

            rdnm(2:ni-1,2:nj-1) = 1. / MAX( REAL( &
               &          msk(3:ni,2:nj-1) + msk(2:ni-1,3:nj) + msk(1:ni-2,2:nj-1) + msk(2:ni-1,1:nj-2)   &
               & + ris2*( msk(3:ni,3:nj)   + msk(3:ni,1:nj-2) + msk(1:ni-2,1:nj-2) + msk(1:ni-2,3:nj) )  ,4), 0.01 )

            X(2:ni-1,2:nj-1) = w0   *xtmp(2:ni-1,2:nj-1) &
               & + (1.-w0)*(         xtmp(3:ni,2:nj-1) + xtmp(2:ni-1,3:nj) + xtmp(1:ni-2,2:nj-1) + xtmp(2:ni-1,1:nj-2)    &
               &            + ris2*( xtmp(3:ni,3:nj)   + xtmp(3:ni,1:nj-2) + xtmp(1:ni-2,1:nj-2) + xtmp(1:ni-2,3:nj) )  ) &
               &                                  * rdnm(2:ni-1,2:nj-1)

            WHERE( nbp(2:ni-1,2:nj-1) == 0 ) X(2:ni-1,2:nj-1) = xorig(2:ni-1,2:nj-1)  ! We do nothing for pixels that have zero neighbor!!!

         ELSE
            !IF ( l_mask ) PRINT *, ' -- SMOOTH is NOT excluding masked points! (despite presence of "msk")'
            X(2:ni-1,2:nj-1) = w0   *xtmp(2:ni-1,2:nj-1) &
               & + (1.-w0)*(         xtmp(3:ni,2:nj-1) + xtmp(2:ni-1,3:nj) + xtmp(1:ni-2,2:nj-1) + xtmp(2:ni-1,1:nj-2)   &
               &            + ris2*( xtmp(3:ni,3:nj)   + xtmp(3:ni,1:nj-2) + xtmp(1:ni-2,1:nj-2) + xtmp(1:ni-2,3:nj) )  ) &
               &                       / ( 4. * (1. + ris2) )

         END IF


         !! we can use east-west periodicity:
         !! ---------------------------------
         IF (k_ew >= 0) THEN
            DO jci = 1, 2
               jim = vim_per(jci)  ! ji-1
               ji  = ivi(jci)      ! first ji = 1, then ji = ni
               jip = vip_per(jci)  ! ji+1


               IF ( l_emp ) THEN
                  !! -- SMOOTHER is excluding masked points...'
                  nbp(ji,2:nj-1) =   msk(jip,2:nj-1) + msk(ji,3:nj)    + msk(jim,2:nj-1) + msk(ji,1:nj-2) & ! Number of non-masked surrounded points
                     &             + msk(jip,3:nj)   + msk(jip,1:nj-2) + msk(jim,1:nj-2) + msk(jim,3:nj)

                  rdnm(ji,2:nj-1) = 1./MAX( REAL(msk(jip,2:nj-1) + msk(ji,3:nj)    + msk(jim,2:nj-1) + msk(ji,1:nj-2)   &
                     &                  + ris2*( msk(jip,3:nj)   + msk(jip,1:nj-2) + msk(jim,1:nj-2) + msk(jim,3:nj) ),4),0.01)

                  X(ji,2:nj-1) = w0*xtmp(ji,2:nj-1) &
                     & + (1.-w0)*(         xtmp(jip,2:nj-1) + xtmp(ji,3:nj)    + xtmp(jim,2:nj-1) + xtmp(ji,1:nj-2)    &
                     &            + ris2*( xtmp(jip,3:nj)   + xtmp(jip,1:nj-2) + xtmp(jim,1:nj-2) + xtmp(jim,3:nj) ) ) &
                     &                   * rdnm(ji,2:nj-1)

                  WHERE( nbp(ji,2:nj-1) == 0 ) X(ji,2:nj-1) = xorig(ji,2:nj-1)  ! We do nothing for pixels that have zero neighbor!!!

               ELSE
                  X(ji,2:nj-1) = w0*xtmp(ji,2:nj-1) &
                     &       + (1.-w0)*( xtmp(jip,2:nj-1) + xtmp(ji,3:nj) + xtmp(jim,2:nj-1) + xtmp(ji,1:nj-2)  &
                     &            + ris2*( xtmp(jip,3:nj)   + xtmp(jip,1:nj-2) + xtmp(jim,1:nj-2) + xtmp(jim,3:nj) ) ) &
                     &                   / ( 4. * (1. + ris2) )

               END IF


            END DO
         END IF

         !! Smoothing is applied only where msk==1, values of X remain unchanged elsewhere:
         IF ( l_mask ) X(:,2:nj-1) = msk(:,2:nj-1)*X(:,2:nj-1) - (msk(:,2:nj-1) - 1)*xorig(:,2:nj-1)

      END DO

      DEALLOCATE ( xtmp )
      IF (l_mask) DEALLOCATE (  xorig )
      IF (l_emp)  DEALLOCATE ( rdnm , nbp )

   END SUBROUTINE SMOOTHER


END MODULE MOD_BDROWN
