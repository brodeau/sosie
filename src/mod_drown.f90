MODULE MOD_DROWN

   !USE io_ezcdf
   !USE mod_manip
   
   IMPLICIT none

   PRIVATE

   PUBLIC :: drown, smooth

   LOGICAL, PARAMETER :: ldebug = .FALSE.

CONTAINS
   
   SUBROUTINE DROWN(k_ew, X, mask,   nb_inc, nb_smooth)

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
      !!                => default: nb_smooth = 5
      !!
      !!
      !!                       Author : Laurent BRODEAU, 2014
      !!
      !!#############################################################################

      !! Arguments :
      INTEGER,                       INTENT(in)    :: k_ew
      REAL(4),    DIMENSION(:,:),    INTENT(inout) :: X
      INTEGER(2), DIMENSION(:,:),    INTENT(in)    :: mask

      INTEGER,    OPTIONAL,          INTENT(in)    :: nb_inc, nb_smooth

      !! Local :
      INTEGER(2), ALLOCATABLE, DIMENSION(:,:) :: maskv, mask_coast, mtmp
      REAL(4),    ALLOCATABLE, DIMENSION(:,:) :: dold, xtmp

      INTEGER :: &
         &      ninc_max,      &
         &      nsmooth_max,          &
         &      ni, nj,        &
         &      jinc,          &
         &      ji, jj, jci,   &
         &      jim, jip

      REAL(4), PARAMETER :: rr = 0.707

      INTEGER, DIMENSION(2) :: ivi, vim_per, vip_per

      INTEGER, PARAMETER :: jinc_debg = 2


      ninc_max = 200   ! will stop before when all land points have been treated!!!
      IF ( present(nb_inc) ) ninc_max = nb_inc

      nsmooth_max = 5
      IF ( present(nb_smooth) ) nsmooth_max = nb_smooth


      IF ( (size(X,1) /= size(mask,1)).OR.(size(X,2) /= size(mask,2)) ) THEN
         PRINT *, 'ERROR, mod_drown.F90 => DROWN : size of data and mask do not match!!!'; STOP
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
            IF ( ldebug ) PRINT *, 'DROWN: No land points left! Leaving incursion loop at jinc =', jinc
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
         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(REAL(mask_coast,4), 'mask_coast.nc', 'lsm') !; STOP 'mod_drown.F90 => boo!'
         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(X, 'data_X_before.nc', 'lsm')
         !STOP



         !! mask_coast done, time to fill the coastline points with values from the nearest se points
         !! -----------------------------------------------------------------------------------------

         !! Center of the domain:
         DO jj = 2, nj-1
            DO ji = 2, ni-1
               IF ( mask_coast(ji,jj) == 1 ) THEN
                  X(ji,jj) = 1./(maskv(ji+1,jj)+maskv(ji,jj+1)+maskv(ji-1,jj)+maskv(ji,jj-1) + &
                     & rr*(maskv(ji+1,jj+1)+maskv(ji-1,jj+1)+maskv(ji-1,jj-1)+maskv(ji+1,jj-1)))*( &
                     & maskv(ji+1,jj)*dold(ji+1,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                     & maskv(ji-1,jj)*dold(ji-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                     & rr*maskv(ji+1,jj+1)*dold(ji+1,jj+1) + rr*maskv(ji-1,jj+1)*dold(ji-1,jj+1) + &
                     & rr*maskv(ji-1,jj-1)*dold(ji-1,jj-1) + rr*maskv(ji+1,jj-1)*dold(ji+1,jj-1)  )
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
                        & rr*(maskv(jip,jj+1)+maskv(jim,jj+1)+maskv(jim,jj-1)+maskv(jip,jj-1)))*( &
                        & maskv(jip,jj)*dold(jip,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & maskv(jim,jj)*dold(jim,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & rr*maskv(jip,jj+1)*dold(jip,jj+1) + rr*maskv(jim,jj+1)*dold(jim,jj+1) + &
                        & rr*maskv(jim,jj-1)*dold(jim,jj-1) + rr*maskv(jip,jj-1)*dold(jip,jj-1)  )
                  END IF
               END DO

            ELSE

               !! West & East LBCs when not east-west periodicity, extrapolating lineraly
               IF ( ji == 1 ) THEN
                  DO jj = 2, nj-1
                     IF ( mask_coast(ji,jj) == 1 ) THEN
                        X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj+1)+maskv(ji,jj-1) + &
                           & rr*maskv(2,jj+1)+rr*maskv(2,jj-1))*( &
                           & maskv(2,jj)*dold(2,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                           & maskv(ji,jj-1)*dold(ji,jj-1) + &
                           & rr*maskv(2,jj+1)*dold(2,jj+1) + &
                           & rr*maskv(2,jj-1)*dold(2,jj-1)  )
                     END IF
                  END DO
               END IF
               IF ( ji == ni ) THEN
                  DO jj = 2, nj-1
                     IF ( mask_coast(ji,jj) == 1 ) THEN
                        X(ji,jj) = 1./( maskv(ji,jj+1)+maskv(ni-1,jj)+maskv(ji,jj-1) + &
                           & rr*maskv(ni-1,jj+1)+rr*maskv(ni-1,jj-1))*( &
                           & maskv(ji,jj+1)*dold(ji,jj+1) + &
                           & maskv(ni-1,jj)*dold(ni-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                           & rr*maskv(ni-1,jj+1)*dold(ni-1,jj+1) + &
                           & rr*maskv(ni-1,jj-1)*dold(ni-1,jj-1) )
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
                  & rr*maskv(ji-1,jj-1)+rr*maskv(ji+1,jj-1) )*( &
                  & maskv(ji+1,jj)*dold(ji+1,jj) + &
                  & maskv(ji-1,jj)*dold(ji-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                  & rr*maskv(ji-1,jj-1)*dold(ji-1,jj-1) + rr*maskv(ji+1,jj-1)*dold(ji+1,jj-1)  )
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
                     & rr*maskv(jim,jj-1)+rr*maskv(jip,jj-1))*( &
                     & maskv(jip,jj)*dold(jip,jj) + &
                     & maskv(jim,jj)*dold(jim,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                     & rr*maskv(jim,jj-1)*dold(jim,jj-1) + rr*maskv(jip,jj-1)*dold(jip,jj-1)  )
               END IF

               ! No E-W periodicity:
            ELSE
               IF ( ji == 1 ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj-1) + &
                        & rr*maskv(2,jj-1))*( &
                        & maskv(2,jj)*dold(2,jj) + &
                        & maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & rr*maskv(2,jj-1)*dold(2,jj-1)  )
                  END IF
               END IF
               IF ( ji == ni ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(ni-1,jj)+maskv(ji,jj-1) + &
                        & rr*maskv(ni-1,jj-1))*( &
                        & maskv(ni-1,jj)*dold(ni-1,jj) + maskv(ji,jj-1)*dold(ji,jj-1) + &
                        & rr*maskv(ni-1,jj-1)*dold(ni-1,jj-1)  )
                  END IF
               END IF

            END IF
         END DO




         !! Center of Bottom row:
         jj = 1
         DO ji = 2, ni-1
            IF ( mask_coast(ji,jj) == 1 ) THEN
               X(ji,jj) = 1./(maskv(ji+1,jj)+maskv(ji,jj+1)+maskv(ji-1,jj) + &
                  & rr*maskv(ji+1,jj+1)+rr*maskv(ji-1,jj+1) )*( &
                  & maskv(ji+1,jj)*dold(ji+1,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                  & maskv(ji-1,jj)*dold(ji-1,jj) + &
                  & rr*maskv(ji+1,jj+1)*dold(ji+1,jj+1) + rr*maskv(ji-1,jj+1)*dold(ji-1,jj+1) )
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
                     & rr*maskv(jip,jj+1)+rr*maskv(jim,jj+1) )*( &
                     & maskv(jip,jj)*dold(jip,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                     & maskv(jim,jj)*dold(jim,jj) + &
                     & rr*maskv(jip,jj+1)*dold(jip,jj+1) + rr*maskv(jim,jj+1)*dold(jim,jj+1) )
               END IF

               !! No E-W periodicity:
            ELSE
               IF ( ji == 1 ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(2,jj)+maskv(ji,jj+1) + &
                        & rr*maskv(2,jj+1) )*( &
                        & maskv(2,jj)*dold(2,jj) + maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & rr*maskv(2,jj+1)*dold(2,jj+1) )
                  END IF
               END IF
               IF ( ji == ni ) THEN
                  IF ( mask_coast(ji,jj) == 1 ) THEN
                     X(ji,jj) = 1./(maskv(ji,jj+1)+maskv(ni-1,jj) + &
                        & rr*maskv(ni-1,jj+1) )*( &
                        & maskv(ji,jj+1)*dold(ji,jj+1) + &
                        & maskv(ni-1,jj)*dold(ni-1,jj) + &
                        & rr*maskv(ni-1,jj+1)*dold(ni-1,jj+1) )
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

      !! Time to smooth what's been drowned:
      mtmp = 1 - mask ! 1 over continents, 0 over seas!
      CALL SMOOTH(k_ew, X,  nb_smooth=nsmooth_max, mask_apply=mtmp)

      !CALL DUMP_2D_FIELD(X, 'drowned_final.nc', 'lsm') ;     !     STOP 'lolo'

      DEALLOCATE ( maskv, mtmp, xtmp, dold, mask_coast )

      IF ( ldebug ) PRINT *, 'DROWN: jinc =', jinc

   END SUBROUTINE DROWN

   


   SUBROUTINE SMOOTH(k_ew, X,  nb_smooth, mask_apply)
      
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
      !!                => default: nb_smooth = 5
      !!  * mask_apply : mask array that defines where the smoothing should be applied
      !!                 => where mask_apply==1: smoothing applies
      !!                 => where mask_apply==0: original values of X will be preserved
      !!
      !!#############################################################################

      !! Arguments :
      INTEGER,                    INTENT(in)           :: k_ew
      REAL(4),    DIMENSION(:,:), INTENT(inout)        :: X     
      INTEGER,    OPTIONAL,                 INTENT(in) :: nb_smooth
      INTEGER(2), OPTIONAL, DIMENSION(:,:), INTENT(in) :: mask_apply

      REAL(4),    ALLOCATABLE, DIMENSION(:,:) :: xorig, xtmp

      INTEGER, DIMENSION(2) :: ivi, vim_per, vip_per

      LOGICAL :: l_mask
      
      INTEGER :: &
         &      nsmooth_max,          &
         &      ni, nj,        &
         &      ji, jci,   &
         &      jim, jip, js
      
      nsmooth_max = 5
      IF ( PRESENT(nb_smooth) ) nsmooth_max = nb_smooth

      l_mask = .FALSE.
      IF ( PRESENT(mask_apply) ) l_mask = .TRUE.
      
      ni = SIZE(X,1)
      nj = SIZE(X,2)


      ALLOCATE ( xtmp(ni,nj) )
      IF ( l_mask ) ALLOCATE (  xorig(ni,nj) )

      ivi = (/ 1 , ni /)
      
      IF (k_ew >= 0) THEN
         vim_per = (/ ni-k_ew ,  ni-1  /)
         vip_per = (/    2    , 1+k_ew /)
      END IF

      IF ( l_mask ) xorig(:,:) = X(:,:)  
      
      DO js = 1, nsmooth_max
         
         xtmp(:,:) = X(:,:)

         !! Center of the domain:
         X(2:ni-1,2:nj-1) = 0.35*xtmp(2:ni-1,2:nj-1) &
            &           + 0.65*0.25*(  xtmp(3:ni,2:nj-1)   + xtmp(2:ni-1,3:nj)     &
            &                        + xtmp(1:ni-2,2:nj-1) + xtmp(2:ni-1,1:nj-2) )
         
         !! we can use east-west periodicity:
         IF (k_ew >= 0) THEN
            DO jci = 1, 2
               jim = vim_per(jci)  ! ji-1
               ji  = ivi(jci)      ! first ji = 1, then ji = ni
               jip = vip_per(jci)  ! ji+1
               
               X(ji,2:nj-1) = 0.35*xtmp(ji,2:nj-1) &
                  &       + 0.65*0.25*(xtmp(jip,2:nj-1) + xtmp(ji,3:nj) + xtmp(jim,2:nj-1) + xtmp(ji,1:nj-2) )
               
            END DO
         END IF

         !! Smoothing is applied only where mask_apply==1, values of X remain unchanged elsewhere:
         IF ( l_mask ) X(:,2:nj-1) = mask_apply(:,2:nj-1)*X(:,2:nj-1) - (mask_apply(:,2:nj-1) - 1)*xorig(:,2:nj-1)
         
      END DO
      
      DEALLOCATE ( xtmp )
      IF ( l_mask ) DEALLOCATE (  xorig )
      
   END SUBROUTINE SMOOTH
   
   
END MODULE MOD_DROWN
