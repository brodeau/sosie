MODULE MOD_DROWN3D

   USE io_ezcdf !LOLO
   !USE mod_manip
   USE, INTRINSIC :: ieee_arithmetic

   IMPLICIT none

   PRIVATE

   PUBLIC :: drown3d, persistence_topbot

   LOGICAL, PARAMETER :: ldebug = .FALSE.

   REAL, PARAMETER :: ris2 = 1.0/SQRT(2.0)

CONTAINS

   SUBROUTINE DROWN3D(k_ew, X, mask, cmethod)

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
      INTEGER,                         INTENT(in)    :: k_ew
      REAL(4),    DIMENSION(:,:,:),    INTENT(inout) :: X
      INTEGER(1), DIMENSION(:,:,:),    INTENT(in)    :: mask

      CHARACTER(LEN=4), INTENT(in) :: cmethod

      !! Local :
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:,:) :: maskv
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:)   :: mask_coast
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:,:) :: mask_tofill, mask_dta
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:,:) :: zmask, zmask_tofill
      REAL(4),    ALLOCATABLE, DIMENSION(:,:)   :: xtmp,zdebug
      REAL(4),    ALLOCATABLE, DIMENSION(:,:,:) :: zX, zdata
      REAL(4),                 DIMENSION(3)     :: zw, zweight

      INTEGER :: &
         &      ninc_max,          &
         &      ni, nj, ns, nk,    &
         &      jinc,              &
         &      ji, jj, jim, jk,   &
         &      jjm, jic, jjc, ji0, jip1, jim1, jjp1, jjm1, jkm1, jkp1

      INTEGER :: ntodo=1, ntodo_b=0, ntodo_bb=0, coast_dir

      REAL(4) :: datmsk, summsk

      INTEGER, PARAMETER :: jinc_debg = 2

      ! flood filling
      INTEGER :: iib, ijb, ii, ij, iip1, iim1, ijp1, ijm1, ip, idist, iptest, np
      INTEGER(1), ALLOCATABLE, DIMENSION(:,:) :: ioptm
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: ipile
      !REAL(8) :: tSTART, tEND, omp_get_wtime 
      !tSTART = omp_get_wtime() 

      ntodo=1 ; ntodo_b=0 ; ntodo_bb=0

      X = X * mask  ! we rather have 0s on continents than some fucked up high values...

      IF ( (size(X,1) /= size(mask,1)) .OR. (size(X,2) /= size(mask,2)) .OR. (size(X,3) /= size(mask,3)) ) THEN
         PRINT *, 'ERROR, mod_drown.F90 => DROWN : size of data and mask do not match!!!'; STOP
      END IF

      ni = SIZE(X,1)
      nj = SIZE(X,2)
      nk = SIZE(X,3)

      ninc_max = MAX(ni,nj)   ! will stop before when all land points have been treated!!!

      ! set weight for the upper (1) and lower slice (3)
      zw(1) = 0.001
      zw(2) = 1.000
      zw(3) = 0.001

      !! Backing up original mask into mask2(:,:)
      ALLOCATE ( mask_tofill(ni,nj,nk), zmask(ni,nj,nk), mask_dta(ni,nj,nk) )
      ALLOCATE ( zdata(ni,nj,nk))
      
      zdata=mask ; CALL DUMP_FIELD(zdata, 'mask_after_zint.nc', 'mask_after_zint') !#LB

      ! set to NaN value still at on 'land'
      mask_dta=mask
      DO jk = 1,nk
         DO jj = 1, nj
            DO ji = 1, ni
               IF ( ISNAN(X(ji,jj,jk)) ) THEN
                  mask_dta(ji,jj,jk) = 0
               END IF
            END DO
         END DO
      END DO
      WHERE (mask_dta == 0) X = 12345.

      IF (TRIM(cmethod) == 'HOLE') mask_dta = mask_dta * mask

      !! Building mask defining where to fill
      CALL MSK_TOFILL(mask, mask_dta, mask_tofill, k_ew, cmethod)
      zdata=mask_tofill ; CALL DUMP_FIELD(zdata, 'mask_tofill.nc', 'mask_tofill') !#LB

      DO jinc = 1, ninc_max

         IF ( TRIM(cmethod) == 'LAND' .AND. (jinc == 500) ) CALL DUMP_FIELD(X, 'data_jinc500.nc', 'data_jinc500')

         IF (ntodo == ntodo_b .AND. ntodo /= ntodo_bb) THEN
            coast_dir=0
            PRINT *, 'Fill vertically ',jinc
            ! is applying persistence_topbot here good enough ? => will decrease amount of memory and initialisation used later.

         ELSE IF (ntodo /= ntodo_b) THEN
            coast_dir=1
         ELSE
            PRINT *, 'Still closed pool (probably issue with missing north fold treatment in msk_coat and other algo), proceed to filling land'

            ! set to NaN value still at on 'land'
            WHERE ( X == 12345. )
               X = ieee_value(1., ieee_quiet_nan)
            END WHERE

            EXIT
         END IF

         ntodo_bb= ntodo_b
         ntodo_b = ntodo

         ns=MIN(jinc,40)

         zmask(:,:,:) = mask_dta(:,:,:)
         zdata(:,:,:) = X(:,:,:)
         
         !$OMP PARALLEL DEFAULT(PRIVATE) &
            !$OMP SHARED(nk, nj, ni, ns, zmask, mask_tofill, zdata, k_ew, X, mask_dta, coast_dir, jinc, zw, cmethod)

         IF ( TRIM( cmethod ) == 'HOLE' ) THEN
            ALLOCATE ( ioptm(-ns:ns,-ns:ns  ) )
            ALLOCATE ( ipile((2*ns+1)**2/2,5) )
         END IF

         ALLOCATE ( zmask_tofill(ni,nj,3)  )
         ALLOCATE ( maskv(ni,nj,3), xtmp(ni,nj), mask_coast(ni,nj) )
         ALLOCATE ( zX(ni,nj,3))

         !$OMP DO SCHEDULE(STATIC,1)
         DO jk = 1,nk
            ! define upper and lower level
            jkp1=MIN(nk,jk+1)
            jkm1=MAX( 1,jk-1)

            ! define upper, current and lower level mask
            maskv(:,:,1) = zmask(:,:,jkm1)
            maskv(:,:,2) = zmask(:,:,jk  )
            maskv(:,:,3) = zmask(:,:,jkp1)

            zmask_tofill(:,:,1) = mask_tofill(:,:,jkm1)
            zmask_tofill(:,:,2) = mask_tofill(:,:,jk  )
            zmask_tofill(:,:,3) = mask_tofill(:,:,jkp1)

            zX(:,:,1) = zdata(:,:,jkm1)
            zX(:,:,2) = zdata(:,:,jk  )
            zX(:,:,3) = zdata(:,:,jkp1)

            !! Building mask of the coast-line (belonging to land points)
            CALL MSK_COAST(mask_coast, zmask_tofill, maskv, k_ew, coast_dir)

            !! mask_coast done, time to fill the coastline points with values from the nearest se points
            !! -----------------------------------------------------------------------------------------
            xtmp(:,:) = zdata(:,:,jk)
           
            DO jj = 1, nj
               DO ji = 1, ni

                  ! update only coastal point
                  IF ( mask_coast(ji,jj) == 1 ) THEN

                     ! initialise variables
                     summsk = 0.0
                     datmsk = 0.0
                     !
                     IF ( TRIM(cmethod) == 'HOLE' ) THEN

                     ioptm(: ,: ) = 1
                     ioptm(0 ,0 ) = 0 
                     ipile(:,:)=HUGE(1)

                     ! initialised pile
                     ipile(1,:) = [ji,jj,0,0,0]
                     np = 1

                     ! loop until the pile size is 0 or if the pool is larger than the critical size
                     DO WHILE ( np /= 0 );

                        ! load closest cell to seed
                        ip = MINLOC(ipile(:,5),DIM=1)

                        ! next point
                        ii=ipile(ip,1); ij=ipile(ip,2); iib=ipile(ip,3); ijb=ipile(ip,4); idist=ipile(ip,5)

                        ! check neighbour cells and update pile ( assume E-W periodicity )
                        IF ( k_ew >= 0 ) THEN
                           IF ( ii == ni+1 ) ii=1 + k_ew
                           IF ( ii == 0    ) ii=ni-k_ew
                           iip1=ii+1; IF ( iip1 == ni+1) iip1=1 + k_ew
                           iim1=ii-1; IF ( iim1 == 0   ) iim1=ni-k_ew
                        ELSE
                           IF ( ii == ni+1 ) ii=ni
                           IF ( ii == 0    ) ii=1
                           iip1=ii+1; IF ( iip1 == ni+1 ) iip1=ni
                           iim1=ii-1; IF ( iim1 == 0    ) iim1=1
                        END IF
                        ijp1=MIN(ij+1,nj)  ! north fold not treated
                        ijm1=MAX(ij-1,1)
            
                        ! compute gaussian weight
                        zweight(:) = EXP(-1.*idist**2/jinc**2)*maskv(ii,ij,:)*zw(:)
                        summsk=summsk + SUM(zweight)
                        datmsk=datmsk + SUM(zweight(:)*zX(ii,ij,:)*maskv(ii,ij,:))

                        ! update pile size
                        ipile(ip,:) = [0,0,0,0,HUGE(1)]; np=np-1
                     
                        ! check neighbour cells and update pile
                        IF (zmask_tofill(ii, ijp1,2) == 1 .AND. ioptm(iib, MIN( ns,ijb+1)) == 1 .AND. idist <= ns - 1) THEN
                           np=np+1; ip = FINDLOC(ipile(:,5), HUGE(1),DIM=1); ipile(ip,:)=[ii  ,ijp1,iib ,ijb+1, idist+1]
                           ioptm (iib,ijb+1) = 0 
                        END IF
                        IF (zmask_tofill(ii, ijm1,2) == 1 .AND. ioptm(iib, MAX(-ns,ijb-1)) == 1 .AND. idist <= ns - 1) THEN
                           np=np+1; ip = FINDLOC(ipile(:,5), HUGE(1),DIM=1); ipile(ip,:)=[ii  ,ijm1,iib ,ijb-1, idist+1]
                           ioptm (iib,ijb-1) = 0 
                        END IF
                        IF (zmask_tofill(iip1, ij,2) == 1 .AND. ioptm(MIN( ns,iib+1), ijb) == 1 .AND. idist <= ns - 1) THEN
                           np=np+1; ip = FINDLOC(ipile(:,5), HUGE(1),DIM=1); ipile(ip,:)=[iip1,ij ,iib+1,ijb, idist+1]
                           ioptm (iib+1,ijb) = 0 
                        END IF
                        IF (zmask_tofill(iim1, ij,2) == 1 .AND. ioptm(MAX(-ns,iib-1), ijb) == 1 .AND. idist <= ns - 1) THEN
                           np=np+1; ip = FINDLOC(ipile(:,5), HUGE(1),DIM=1); ipile(ip,:)=[iim1,ij ,iib-1,ijb, idist+1]
                           ioptm (iib-1,ijb) = 0 
                        END IF
                     END DO

                     ELSE IF ( (TRIM( cmethod ) == 'LAND') .OR. (TRIM( cmethod ) == 'BOUN') ) THEN
                     ! compute weighted average in a box center on ji,jj 
                     DO jjm=-ns,ns
                        DO jim=-ns,ns
                           !
                           ! box index definition
                           jic = ji + jim ; jjc = jj + jjm ;
                           IF  ( k_ew >= 0 ) THEN
                              IF ( jic > ni - k_ew ) jic = jic - ni + k_ew ! W boundary
                              IF ( jic <  1 + k_ew ) jic = jic + ni - k_ew ! E boundary
                           END IF
                           !
                           IF (jic >= 1 .AND. jic <= ni .AND. jjc >= 1 .AND. jjc <= nj) THEN
                              !
                              ! compute gaussian weight
                              zweight(:) = EXP(-1.*(jjm**2+jim**2)/jinc**2)*maskv(jic,jjc,:)*zw(:)
                              summsk=summsk + SUM(zweight)
                              datmsk=datmsk + SUM(zweight(:)*zX(jic,jjc,:)*maskv(jic,jjc,:))
                           END IF
                        END DO
                     END DO
                    ELSE
                       PRINT *, TRIM(cmethod),' is unknown'
                       STOP
                    END IF
                     !
                     ! compute mean over the defined box with gaussian weight (only where data valid)
                     xtmp(ji,jj) = datmsk / summsk
                  END IF
               END DO
            END DO
            !
            ! update X
            X(:,:,jk) = xtmp(:,:)
            !
            ! Loosing land for the next iteration:
            mask_dta(:,:,jk) = zmask(:,:,jk) + mask_coast(:,:)
            !
         END DO 
         !$OMP END DO

         IF ( TRIM(cmethod ) == 'HOLE' ) DEALLOCATE(ioptm, ipile)
         DEALLOCATE ( maskv, xtmp, mask_coast, zmask_tofill )
         DEALLOCATE ( zX )
 
         !$OMP END PARALLEL
         !
         IF ( ldebug ) PRINT *, 'DROWN: jinc =', jinc
         ntodo = SUM(INT(mask_tofill)-INT(mask_dta))
         PRINT *, 'DROWN: jinc =', jinc, ntodo

         !! Quiting if no land point left:
         IF ( SUM(INT(mask_tofill)-INT(mask_dta)) == 0  ) THEN
            IF ( ldebug ) PRINT *, 'DROWN: No land points left! Leaving incursion loop at jinc =', jinc
            EXIT
         END IF
         !
      END DO

      ! set to NaN value still at on 'land'
      WHERE ( mask_tofill == 0 )
         X = ieee_value(1., ieee_quiet_nan)
      END WHERE

   !tEND = omp_get_wtime() 
   !PRINT *, '~Drowning algo took~', tEND - tSTART, '~seconds~'

   END SUBROUTINE DROWN3D

   SUBROUTINE MSK_TOFILL(mask, mask_dta, mask_tofill, k_ew, cmethod)
      INTEGER(1), DIMENSION(:,:,:), INTENT(in)   :: mask, mask_dta
      INTEGER(1), DIMENSION(:,:,:), INTENT(out)  :: mask_tofill
      INTEGER(4),                   INTENT(in)   :: k_ew

      CHARACTER(LEN=4), INTENT(in) :: cmethod

      INTEGER(1), ALLOCATABLE, DIMENSION(:,:) :: zmask
      INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: mbkt, mikt

      INTEGER(4) :: ni, nj, nk
      INTEGER(4) :: ji, jj, jk, jm
      INTEGER(4) :: ji0,jip1,jim1,jjm1,jjp1

      INTEGER(4) :: ztmp

      ni = SIZE(mask,1)
      nj = SIZE(mask,2)
      nk = SIZE(mask,3)
      
      IF (TRIM(cmethod) == 'LAND') THEN
         mask_tofill(:,:,:) = 1
      ELSE IF (TRIM(cmethod) == 'HOLE') THEN
         mask_tofill(:,:,:) = mask(:,:,:)
      ELSE IF (TRIM(cmethod) == 'BOUN') THEN

         mask_tofill = mask_dta
         
         ALLOCATE(zmask(ni,nj))

         ! open 2 cell on the side
         DO jk=1,nk
            DO jm = 1,2
               zmask(:,:) = mask_tofill(:,:,jk)
               DO ji=1,ni
                  DO jj=1,nj
                     ! manage periodicity
                     jjp1 = jj+1 ; jjm1 = jj-1 ;
                     IF ( jjp1 > nj ) jjp1 = nj
                     IF ( jjm1 <  1 ) jjm1 = 1
         
                     ji0 = ji ; jip1 = ji+1 ; jim1 = ji-1 ;
                     IF  ( k_ew >= 0 ) THEN
                        IF ( ji0  > ni - k_ew ) ji0  = ji0  - ni + k_ew
                        IF ( ji0  <  1 + k_ew ) ji0  = ji0  + ni - k_ew
                        IF ( jip1 > ni - k_ew ) jip1 = jip1 - ni + k_ew ! W boundary
                        IF ( jip1 <  1 + k_ew ) jip1 = jip1 + ni - k_ew ! E boundary
                        IF ( jim1 > ni - k_ew ) jim1 = jim1 - ni + k_ew ! W boundary
                        IF ( jim1 <  1 + k_ew ) jim1 = jim1 + ni - k_ew ! E boundary
                     ELSE
                        IF ( jip1 > ni ) jip1 = ni
                        IF ( jim1 <  1 ) jim1 = 1
                     END IF

                     ! find coastal point (include corner)
                     IF ( zmask( ji0, jj) == 0 ) THEN
                        ztmp =        zmask(jim1,jjm1) + zmask(jim1,jj) + zmask(jim1,jjp1)
                        ztmp = ztmp + zmask(ji0 ,jjm1) +                  zmask(ji0 ,jjp1)
                        ztmp = ztmp + zmask(jip1,jjm1) + zmask(jip1,jj) + zmask(jip1,jjp1)
                        mask_tofill(ji0,jj,jk) = MIN(ztmp,1)
                     END IF

                  END DO
               END DO
            END DO
         END DO

         DEALLOCATE(zmask)
 
      END IF

   END SUBROUTINE

   SUBROUTINE MSK_COAST(mask_coast, mask_tofill, mask, k_ew, idir)

      INTEGER(1), DIMENSION(:,:,:),    INTENT(in)     :: mask
      INTEGER(1), DIMENSION(:,:,:),    INTENT(in)     :: mask_tofill
      INTEGER(1), DIMENSION(:,:)  ,    INTENT(out)    :: mask_coast
      INTEGER(4) :: k_ew, idir

      INTEGER :: &
         &      ni, nj,        &
         &      ji, jj,        &
         &      ji0, jip1, jim1, jjp1, jjm1

      INTEGER(4) :: zsum

      ni = SIZE(mask,1)
      nj = SIZE(mask,2)

      !! Building mask of the coast-line (belonging to land points)
      mask_coast(:,:) = 0
      DO jj = 1, nj
         DO ji = 1, ni
            !
            jjp1 = jj+1 ; jjm1 = jj-1 ;
            IF ( jjp1 > nj ) jjp1 = nj
            IF ( jjm1 <  1 ) jjm1 = 1

            ji0 = ji ; jip1 = ji+1 ; jim1 = ji-1 ;
            IF  ( k_ew >= 0 ) THEN
               IF ( ji0  > ni - k_ew ) ji0  = ji0  - ni + k_ew
               IF ( ji0  <  1 + k_ew ) ji0  = ji0  + ni - k_ew
               IF ( jip1 > ni - k_ew ) jip1 = jip1 - ni + k_ew ! W boundary
               IF ( jip1 <  1 + k_ew ) jip1 = jip1 + ni - k_ew ! W boundary
               IF ( jim1 > ni - k_ew ) jim1 = jim1 - ni + k_ew ! W boundary
               IF ( jim1 <  1 + k_ew ) jim1 = jim1 + ni - k_ew ! E boundary
            ELSE
               IF ( jip1 > ni ) jip1 = ni
               IF ( jim1 <  1 ) jim1 = 1
            END IF
            !
            ! mask_coast: inside mask_tofill (value = 1) and not in mask (value = 0) and along coast of mask
            zsum = idir * ( mask(jim1,jj,2) + mask(jip1,jj,2) + mask(ji0,jjm1,2) + mask(ji0,jjp1,2) ) + (1-idir) * ( mask(ji0,jj,3) + mask(ji0,jj,1) )
            mask_coast(ji,jj) = MIN( zsum , 1 ) * ( 1 - mask(ji0,jj,2)) * mask_tofill(ji0,jj,2)
            !
         END DO
      END DO
      !
   END SUBROUTINE MSK_COAST

   SUBROUTINE PERSISTENCE_TOPBOT(rdata)
      REAL(4), DIMENSION(:,:,:), INTENT(inout):: rdata

      INTEGER(4), ALLOCATABLE, DIMENSION(:,:) :: mbkt, mikt

      INTEGER(4) :: ni, nj, nk
      INTEGER(4) :: ji, jj, jk

      PRINT *, 'compute persistence top/bot'

      ni = SIZE(rdata,1)
      nj = SIZE(rdata,2)
      nk = SIZE(rdata,3)
      
      ALLOCATE(mbkt(ni,nj), mikt(ni,nj))

      ! find bottom level
      mbkt(:,:) = nk
      DO jj=1,nj
         DO ji=1,ni
            DO jk=1,nk-1
               IF ( (.NOT. ISNAN(rdata(ji,jj,jk))) .AND. ISNAN(rdata(ji,jj,jk+1)) ) THEN
                  mbkt(ji,jj) = jk
               END IF
            END DO
         END DO
      END DO

     ! find top level
      mikt(:,:) = 1
      DO jj=1,nj
         DO ji=1,ni
            DO jk=2,nk
               IF ( (.NOT. ISNAN(rdata(ji,jj,jk))) .AND. ISNAN(rdata(ji,jj,jk-1)) ) THEN
                  mikt(ji,jj) = jk
               END IF
            END DO
         END DO
      END DO

      ! persistence for the 2 cells above the top of below the bottom
      DO jj=1,nj
         DO ji=1,ni
            ! use 2 because akima need to access this data 
            DO jk =1,2
               rdata(ji,jj,MAX(1 ,mikt(ji,jj)-jk)) = rdata(ji,jj,MAX(1 ,mikt(ji,jj)-jk+1))
               rdata(ji,jj,MIN(nk,mbkt(ji,jj)+jk)) = rdata(ji,jj,MIN(nk,mbkt(ji,jj)+jk-1))
            END DO
         END DO
      END DO

      DEALLOCATE(mikt, mbkt)

   END SUBROUTINE


END MODULE MOD_DROWN3D
