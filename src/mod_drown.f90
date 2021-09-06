MODULE MOD_DROWN

   !USE io_ezcdf !LOLO
   !USE mod_manip

   IMPLICIT none

   PRIVATE

   PUBLIC :: drown, smoother

   LOGICAL, PARAMETER :: ldebug = .FALSE.

   REAL, PARAMETER :: ris2 = 1.0/SQRT(2.0)

CONTAINS

   SUBROUTINE DROWN(k_ew, X, mask,   nb_inc)

      !!#############################################################################
      !!
      !!  PURPOSE : fill continental areas of field X (defined by mask=0)
      !!  -------   using weighted box average (gaussian weight) of surronded
      !!            sea point (defined by mask=1)
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
      !!
      !!                       Author : Laurent BRODEAU, 2014
      !!                                Pierre Mathiot,  2020 : update the drowning algo and add openMP instruction
      !!
      !!#############################################################################

      !! Arguments :
      INTEGER,                       INTENT(in)    :: k_ew
      REAL(4),    DIMENSION(:,:),    INTENT(inout) :: X
      INTEGER(1), DIMENSION(:,:),    INTENT(in)    :: mask

      INTEGER,    OPTIONAL,          INTENT(in)    :: nb_inc

      !! Local :
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:) :: maskv, mask_coast
      REAL(4),    ALLOCATABLE, DIMENSION(:,:) :: xtmp

      INTEGER :: &
         &      ninc_max,      &
         &      ni, nj, ns,    &
         &      jinc,          &
         &      ji, jj, jim,   &
         &      jjm, jic, jjc, ji0, jip1, jim1, jjp1, jjm1

      REAL(4) :: datmsk, summsk, zweight

      INTEGER, PARAMETER :: jinc_debg = 2

      !REAL(8) :: tSTART, tEND, omp_get_wtime 
      !tSTART = omp_get_wtime() 

      X = X * mask  ! we rather have 0s on continents than some fucked up high values...

      ninc_max = 200   ! will stop before when all land points have been treated!!!
      IF ( present(nb_inc) ) ninc_max = nb_inc


      IF ( (size(X,1) /= size(mask,1)).OR.(size(X,2) /= size(mask,2)) ) THEN
         PRINT *, 'ERROR, mod_drown.F90 => DROWN : size of data and mask do not match!!!'; STOP
      END IF

      ni = SIZE(X,1)
      nj = SIZE(X,2)

      !! Backing up original mask into mask2(:,:)
      ALLOCATE ( maskv(ni,nj), xtmp(ni,nj), mask_coast(ni,nj) )

      jinc = 0
      maskv(:,:) = mask(:,:)

      DO jinc = 1, ninc_max

         !! Quiting if no land point left:
         IF ( .NOT. (ANY(maskv == 0))  ) THEN
            IF ( ldebug ) PRINT *, 'DROWN: No land points left! Leaving incursion loop at jinc =', jinc
            EXIT
         END IF

         !! Building mask of the coast-line (belonging to land points)
         mask_coast(:,:) = 0
         !$OMP PARALLEL DEFAULT(NONE) SHARED(k_ew,maskv,mask_coast,ni,nj) PRIVATE(ji,jj,jjp1,jjm1,ji0,jip1,jim1)
         !$OMP DO SCHEDULE(STATIC)
         DO jj = 1, nj
            DO ji = 1, ni
               !
               jjp1 = jj+1 ; jjm1 = jj-1 ;
               IF ( jjp1 > nj ) jjp1 = nj
               IF ( jjm1 <  1 ) jjm1 = 1

               ji0 = ji ; jip1 = ji+1 ; jim1 = ji-1 ;
               IF  ( k_ew >= 0 ) THEN
                  IF ( ji0  > ni - k_ew ) ji0  = ji0  - ni + 2 * k_ew
                  IF ( ji0  <  1 + k_ew ) ji0  = ji0  + ni - 2 * k_ew
                  IF ( jip1 > ni - k_ew ) jip1 = jip1 - ni + 2 * k_ew ! W boundary
                  IF ( jip1 <  1 + k_ew ) jip1 = jip1 + ni - 2 * k_ew ! W boundary
                  IF ( jim1 > ni - k_ew ) jim1 = jim1 - ni + 2 * k_ew ! W boundary
                  IF ( jim1 <  1 + k_ew ) jim1 = jim1 + ni - 2 * k_ew ! E boundary
               ELSE
                  IF ( jip1 > ni ) jip1 = ni
                  IF ( jim1 <  1 ) jim1 = 1
               END IF
               !
               ! mask_coast
               mask_coast(ji,jj) = MIN( maskv(jim1,jj) + maskv(jip1,jj) + maskv(ji0,jjm1) + maskv(ji0,jjp1) , 1 ) * ( 1 - maskv(ji0,jj) )
               !
            END DO
         END DO
         !$OMP END DO
         !$OMP END PARALLEL
         !
         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(REAL(maskv,4), 'maskv.nc', 'lsm')
         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(REAL(mask_coast,4), 'mask_coast.nc', 'lsm') !; STOP 'mod_drown.F90 => boo!'
         !IF ( jinc == jinc_debg) CALL DUMP_2D_FIELD(X, 'data_X_before.nc', 'lsm')
         !STOP


         !! mask_coast done, time to fill the coastline points with values from the nearest se points
         !! -----------------------------------------------------------------------------------------
         ns=MIN(jinc,50)
         xtmp = X
         
         !$OMP PARALLEL DEFAULT(NONE) SHARED(ni,nj,mask_coast,ns,k_ew,maskv,xtmp,X,jinc ) PRIVATE(zweight,summsk,datmsk,ji,jj,jjm,jim,jic,jjc)
         !$OMP DO SCHEDULE(DYNAMIC)
         DO jj = 1, nj
            DO ji = 1, ni
               !
               ! update only coastal point
               IF ( mask_coast(ji,jj) == 1 ) THEN
                  !
                  summsk = 0.0
                  datmsk = 0.0
                  !
                  ! compute weighted average in a box center on ji,jj 
                  DO jjm=-ns,ns
                     DO jim=-ns,ns
                        !
                        ! box index definition
                        jic = ji + jim ; jjc = jj + jjm ;
                        IF  ( k_ew >= 0 ) THEN
                           IF ( jic > ni - k_ew ) jic = jic - ni + 2 * k_ew ! W boundary
                           IF ( jic <  1 + k_ew ) jic = jic + ni - 2 * k_ew ! E boundary
                        END IF
                        !
                        IF (jic >= 1 .AND. jic <= ni .AND. jjc >= 1 .AND. jjc <= nj) THEN
                           !
                           ! compute gaussian weight
                           zweight = EXP(-1.*(jjm**2+jim**2)/jinc**2)*maskv(jic,jjc)
                           summsk=summsk + zweight
                           datmsk=datmsk + zweight*X(jic,jjc)
                        END IF
                     END DO
                  END DO
                  !
                  ! compute mean over the defined box with gaussian weight (only where data valid)
                  xtmp(ji,jj) = datmsk / summsk
               END IF
            END DO
         END DO
         !$OMP END DO
         !$OMP END PARALLEL
         !
         ! update X
         X = xtmp
         !
         ! Loosing land for the next iteration:
         maskv = maskv + mask_coast
         !
         IF ( ldebug ) PRINT *, 'DROWN: jinc =', jinc
         !
      END DO

      !CALL PRTMASK(X, 'drowned_final.nc', 'lsm') ;          STOP 'lolo'


      DEALLOCATE ( maskv, xtmp, mask_coast )

      !tEND = omp_get_wtime() 
      !PRINT *, '~Drowning algo took~', tEND - tSTART, '~seconds~'

   END SUBROUTINE DROWN




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
         PRINT *, 'PROBLEM in SMOOTH (mod_drown.f90): you need to provide a "msk"'
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


END MODULE MOD_DROWN
