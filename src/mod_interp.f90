MODULE MOD_INTERP

   USE mod_conf       !* important parameters, namelist and misc routines
   USE mod_manip      !* misc. manipulation of 2D arrays
   USE mod_drown      !* extrapolation over masked surfaces
   USE mod_akima_2d   !* Akima method algorithm
   USE mod_bilin_2d   !* Bi-linear method (for handling irregular source grids)
   USE mod_akima_1d   !* 1D Akima method for vertical interpolation

   USE mod_nemotools, ONLY: lbc_lnk
   !USE io_ezcdf ; !LOLOdebug

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: INTERP_2D, INTERP_3D

CONTAINS

   SUBROUTINE INTERP_2D()

      !USE io_ezcdf ; !LOLOdebug
      !! ================
      !! 2D INTERPOLATION
      !! ================

      !! lon-aranging or lat-flipping field
      IF ( nlat_inc_in == -1 ) CALL FLIP_UD(data_in)
      IF ( nlon_inc_in == -1 ) CALL LONG_REORG_2D(i_chg_lon, data_in)

      mask_in = mask_in_b    ! re-filling the mask with trusted values...

      IF ( l_drown_src ) THEN
         !CALL DUMP_2D_FIELD(data_in, '1_before_drown.nc', cv_in)
         !CALL DUMP_2D_FIELD(REAL(mask_in(:,:,1),4), '1_mask_before_drown.nc', 'mask')
         !! Extrapolate sea values over land :
         CALL DROWN(ewper, data_in, mask_in(:,:,1), nb_inc=idrown)
         !CALL DUMP_2D_FIELD(data_in, '2_after_drown.nc', cv_in)
      ELSE
         PRINT *, '-------------------'
         PRINT *, 'DROWN NOT CALLED!!!'
         PRINT *, '-------------------'
      END IF

      IF ( ismooth > 0 ) THEN
         !! First, may apply a smoothing on "data_in" in case target grid is much coarser than the source grid!
         WRITE(6,'("     --- ",a,": smoothing ",i2," times!")') TRIM(cv_in), ismooth
         PRINT *, ' Smoothing '//TRIM(cv_in)//'!', ismooth, ' times'
         CALL SMOOTH(ewper, data_in,  nb_smooth=ismooth, mask_apply=mask_in(:,:,1))
      END IF


      !! Call interpolation procedure :
      !! ------------------------------

      SELECT CASE(cmethod)

      CASE('akima')
         CALL akima_2d(ewper, lon_in, lat_in, data_in, lon_out, lat_out, data_out)

      CASE('bilin')
         CALL bilin_2d(ewper, lon_in, lat_in, data_in, lon_out, lat_out, data_out, cpat,  mask_domain_trg=IGNORE)

      CASE('no_xy')
         WRITE(6,*) 'ERROR (mod_interp.f90): method "no_xy" makes no sense for 2D interp!'
         STOP

      CASE DEFAULT
         PRINT *, 'Interpolation method ', cmethod, ' is unknown!!!' ; STOP
      END SELECT


      !! If target grid extends too much in latitude compared to source grid, need to
      !! extrapolate a bit at bottom and top of the domain :
      IF (lregout) CALL extrp_hl(data_out)

      !! Applying bound corrections
      WHERE ( data_out > vmax )  data_out = rmaskvalue
      WHERE ( data_out < vmin )  data_out = rmaskvalue

      !lolo:
      !! If overshoot of latitudes between target and source domain (target has higher values than source):
      !! => apply a drown because the relevant areas were masked (even if lmout==false)!
      IF (jj_ex_btm > 0) THEN
         CALL DROWN(ewper_out, data_out, mask_out(:,:,1),  nb_inc=idrown, nb_smooth=5)
      END IF
      !lolo.

      !! Masking result if requested
      IF ( lmout ) THEN
         WHERE (mask_out(:,:,1) == 0)  data_out = rmaskvalue
      ELSE
         rmaskvalue = 0.
      ENDIF

      !! If target grid is an ORCA grid, calling "lbc_lnk":
      IF ( i_orca_out > 0 ) CALL lbc_lnk( i_orca_out, data_out, c_orca_out, 1.0_8 )

   END SUBROUTINE INTERP_2D



   SUBROUTINE INTERP_3D()

      !! ================
      !! 3D INTERPOLATION
      !! ================

      INTEGER :: ji, jj, jk, jklast=0
      REAL(8) :: zmax_in, zmax_out

      CHARACTER(len=128) :: cfdbg !DEBUG

      !! Interpolation of each of the nk_in levels onto the target grid
      !! --------------------------------------------------------------


      !! Need to do some sort of vertical drown downward to propagate the bottom
      !! value (last water pomit mask==1) down into the sea-bed...
      !DO jj = 1, nj_in
      !   DO ji = 1, ni_in
      !      IF ( mask_in(ji,jj,1) == 1 ) THEN  ! ocean at first level!
      !         IF ((jj==678).AND.(ji==1108)) PRINT *, ' LOLO, before deepening', data3d_in(ji,jj,:) !LOLOdebug
      !         jk = 1
      !         DO WHILE ( jk <= nk_in )
      !            IF ( (jklast == 0).AND.(jk < nk_in) ) THEN
      !               IF ( (mask_in(ji,jj,jk) == 1).AND.(mask_in(ji,jj,jk+1) == 0) ) THEN
      !                  jklast = jk ! last point with water before seabed...
      !                  !IF ((jj==678).AND.(ji==1108)) PRINT *, 'LOLO: bottom point at jk=', jklast
      !               END IF
      !            END IF
      !            IF ( jk > jklast ) data3d_in(ji,jj,jk) = data3d_in(ji,jj,jklast) ! persistence!
      !            jk = jk + 1
      !         END DO
      !      END IF
      !      IF ((jj==678).AND.(ji==1108)) PRINT *, ' LOLO, after deepening', data3d_in(ji,jj,:) !LOLOdebug
      !   END DO
      !END DO
      !jklast = 0
      !!LOLOdebug:
      !DO jk = 1, nk_in
      !   WRITE(cfdbg,'("data_in_drowned_h+b_lev",i2.2,".nc")') jk ; PRINT *, ' *** saving ', TRIM(cfdbg)
      !   CALL DUMP_2D_FIELD(data3d_in(:,:,jk), TRIM(cfdbg), 'data')
      !END DO
      !!LOLOdebug.



      DO jk = 1, nk_in

         PRINT *, '### Preparing source field at level : ', jk

         IF ( cmethod /= 'no_xy' ) THEN !LOLO: WHY????
            IF ( nlat_inc_in == -1 ) CALL FLIP_UD(data3d_in(:,:,jk))
            IF ( nlon_inc_in == -1 ) CALL LONG_REORG_2D(i_chg_lon, data3d_in(:,:,jk))
         END IF

         IF ( l_drown_src ) THEN
            !! Extrapolate sea values over land :
            WRITE(6,'("     --- ",a,": Extrapolating source data over land at level #",i3.3)') TRIM(cv_in), jk
            CALL DROWN(ewper, data3d_in(:,:,jk), mask_in(:,:,jk), nb_inc=idrown, nb_smooth=5 )
         ELSE
            PRINT *, '-------------------'
            PRINT *, 'DROWN NOT CALLED!!!'
            PRINT *, '-------------------'
         END IF

         IF ( ismooth > 0 ) THEN
            IF ( TRIM(cmethod) == 'no_xy' ) THEN
               PRINT *, 'ERROR: makes no sense to perform "no_xy" vertical interpolation and to have ismooth > 0 !'
               STOP
            END IF
            !! First, may apply a smoothing on "data_in" in case target grid is much coarser than the source grid!
            WRITE(6,'("     --- ",a,": Smoothing level #",i3.3," ",i2," times!")') TRIM(cv_in), jk, ismooth
            CALL SMOOTH(ewper, data3d_in(:,:,jk),  nb_smooth=ismooth, mask_apply=mask_in(:,:,jk))
         END IF

      END DO !DO jk = 1, nk_in



      !! Prevent last level to cause problem when not a single water point (last level is only rock!)
      IF ( SUM(mask_in(:,:,nk_in)) == 0) data3d_in(:,:,nk_in) = data3d_in(:,:,nk_in-1) ! persistence!

      !LOLOdebug:
      !DO jk = nk_in/2, nk_in
      !   WRITE(cfdbg,'("data_to_be_interpolated_lev",i2.2,".nc")') jk ; PRINT *, ' *** saving ', TRIM(cfdbg)
      !   CALL DUMP_2D_FIELD(data3d_in(:,:,jk), TRIM(cfdbg), cv_in)
      !END DO
      !LOLOdebug.


      PRINT *, ''
      PRINT *, ' 3D field prepared at all levels, ready to be interpolated...'
      PRINT *, ''


      !! Now! 3D input field ready to be interpolated...
      DO jk = 1, nk_in

         IF (TRIM(cmethod) /= 'no_xy' ) PRINT *, '  *** interpolating at level ', jk

         SELECT CASE(TRIM(cmethod))

         CASE('akima')
            CALL akima_2d(ewper, lon_in,  lat_in,  data3d_in(:,:,jk), &
               &              lon_out, lat_out, data3d_tmp(:,:,jk))
            IF ( trim(ctype_z_in) == 'z' ) THEN
               !! we don't need horizontal interpolation, all levels are flat
               depth_in_trgt2d(:,:,jk) = depth_in(1,1,jk)
            ELSE
               !! input is sigma, layers are non-flat
               CALL akima_2d(ewper, lon_in,  lat_in, depth_in(:,:,jk),       &
                  &              lon_out, lat_out,   depth_in_trgt2d(:,:,jk) )
            ENDIF

         CASE('bilin')
            CALL bilin_2d(ewper, lon_in,  lat_in,  data3d_in(:,:,jk), &
               &              lon_out, lat_out, data3d_tmp(:,:,jk), cpat)

            IF ( trim(ctype_z_in) == 'z' ) THEN
               !! we don't need horizontal interpolation, all levels are flat
               depth_in_trgt2d(:,:,jk) = depth_in(1,1,jk)
               IF ( l_identical_levels ) depth_in_trgt2d(:,:,jk) = depth_out(1,1,jk) !lolo
            ELSE
               !! input is sigma, layers are non-flat
               CALL bilin_2d(ewper, lon_in,  lat_in, depth_in(:,:,jk), &
                  &              lon_out, lat_out,   depth_in_trgt2d(:,:,jk), cpat)
            ENDIF

         CASE('no_xy')
            PRINT *, '  *** copying at level ', jk
            data3d_tmp(:,:,jk) = data3d_in(:,:,jk)
            IF ( TRIM(ctype_z_in) == 'z' ) THEN
               !! we don't need horizontal interpolation, all levels are flat
               depth_in_trgt2d(:,:,jk) = depth_in(1,1,jk)
            ELSE
               !! input is sigma, layers are non-flat
               depth_in_trgt2d(:,:,jk) = depth_in(:,:,jk)
            ENDIF

         CASE DEFAULT
            PRINT *, 'Interpolation method "', trim(cmethod), '" is unknown!!!'; STOP
         END SELECT

         IF (lregout) CALL extrp_hl(data3d_tmp(:,:,jk))

      END DO !DO jk = 1, nk_in

      PRINT *, ''

      data3d_out(:,:,:) = rmaskvalue ! Masking everything



      !! Time for vertical interpolation


      !IF ( .NOT. (TRIM(cf_x_out)  == 'spheric') ) THEN  !LOLO WTF was this???
      IF ( .NOT. l_identical_levels ) THEN  ! ( l_identical_levels always false for S-coordinates...)

         depth_in  = ABS(depth_in)
         depth_out = ABS(depth_out)
         depth_in_trgt2d  = ABS(depth_in_trgt2d)

         !! Need to perform a vertical interpolation from data3d_tmp to data3d_out :


         PRINT *, ''
         IF ( (TRIM(ctype_z_in) == 'z').AND.(TRIM(ctype_z_out) == 'z') ) THEN
            
            !! Go for the vectorial routine...
            PRINT *, ' *** CALLING AKIMA_1D_3D for vertical interpolation !!!'
            CALL AKIMA_1D( depth_in_trgt2d(1,1,:), data3d_tmp, depth_out(1,1,:), data3d_out(:,:,:), rmaskvalue )
            
         ELSE


            DO jj = 1, nj_out
               DO ji = 1, ni_out

                  !! RD dev notes : we need to make sure that the depth vector for both in and out
                  !! are from smallest to largest value so that persistance works
                  IF ( trim(ctype_z_in) == 'sigma' ) THEN
                     CALL FLIP_UD(depth_in_trgt2d(ji,jj,:))
                     CALL FLIP_UD(data3d_tmp(ji,jj,:))
                  ENDIF

                  IF ( trim(ctype_z_out) == 'sigma' ) THEN
                     CALL FLIP_UD(depth_out(ji,jj,:))
                  ENDIF

                  !! RD dev notes : we compare the depth from source depth vector and target depth vector
                  !! at the same horizontal location : compare depth_in_trgt2d and depth_out
                  zmax_in  = MAXVAL(depth_in_trgt2d(ji,jj,:))
                  zmax_out = MAXVAL(depth_out(ji,jj,:))

                  IF ( zmax_out > zmax_in ) THEN
                     !! Must find the last target level less deep than zmax_in
                     jklast = 1
                     DO WHILE ( jklast < nk_out )
                        IF ( depth_out(ji,jj,jklast+1) > zmax_in ) EXIT
                        jklast = jklast + 1
                     END DO
                  ELSE
                     jklast = nk_out
                  END IF

                  IF ( (mask_out(ji,jj,1) == 1) .OR. (.NOT. lmout) ) THEN
                     IF ( lmout ) THEN  ! adapting nlev if masking target
                        nlev = 1
                        !! RD while loop causes seg fault in debug
                        DO jk=1,nk_out
                           IF ( mask_out(ji,jj,jk) == 1 ) nlev = nlev + 1
                        ENDDO
                        nlev = nlev - 1
                     ELSE
                        nlev = jklast
                     END IF
                     !!
                     IF ( (MOD(ji,100)==0).AND.(MOD(jj,100)==0) ) PRINT *, '  ... calling AKIMA_1D for ji_out,jj_out =', ji,jj
                     !!
                     CALL AKIMA_1D( depth_in_trgt2d(ji,jj,:),data3d_tmp(ji,jj,:),  &
                        &           depth_out(ji,jj,1:nlev), data3d_out(ji,jj,1:nlev))
                     !!
                     !! Apply persistance at the bottom if target depth goes deeper that source depth
                     !! RD dev notes : I think the indices were off. If jklast is the last target level
                     !! that can be properly computed then we want to apply persistance to jklast + 1 to nk_out
                     !! btw, interp from z to sigma works slightly better without this on my test case
                     IF ( (jklast > 0).AND.(jklast < nk_out) ) THEN
                        DO jk = jklast+1, nk_out
                           data3d_out(ji,jj,jk) = data3d_out(ji,jj,jklast)
                        END DO
                     END IF

                     !! RD dev notes : when interpolating to sigma, need to reverse again arrays
                     IF ( trim(ctype_z_out) == 'sigma' ) THEN
                        CALL FLIP_UD(depth_out(ji,jj,:))
                        CALL FLIP_UD(data3d_out(ji,jj,:))
                     ENDIF

                  END IF
               END DO
            END DO
         END IF  !IF ( (TRIM(ctype_z_in) == 'z').AND.(TRIM(ctype_z_out) == 'z') ) THEN

      ELSE

         WRITE(6,*) '  *** skipping vertical interpolation as source and target vertical levels are identical!'
         data3d_out = data3d_tmp ! target levels are same than source levels

      END IF  !IF ( .NOT. l_identical_levels )      ! => no vertical interpolation required...


      !! avoid working with 3D arrays as a whole : produce SEGV on some machines (small)
      !! RD : replaced out of bounds values by vmin/vmax not rmaskvalue
      !! as it induced spval instead of zero on BGC fields
      DO jk=1,nk_out
         WHERE ((data3d_out(:,:,jk) > vmax).and.(data3d_out(:,:,jk) /= rmaskvalue)) &
            &   data3d_out(:,:,jk) = vmax
         WHERE ((data3d_out(:,:,jk) < vmin).and.(data3d_out(:,:,jk) /= rmaskvalue)) &
            &   data3d_out(:,:,jk) = vmin
         
         !! If target grid is an ORCA grid, calling "lbc_lnk":
         IF ( i_orca_out > 0 ) CALL lbc_lnk( i_orca_out, data3d_out(:,:,jk), c_orca_out, 1.0_8 )
         
         IF ( lmout ) data3d_out(:,:,jk) = data3d_out(:,:,jk)*mask_out(:,:,jk) + (1. - mask_out(:,:,jk))*rmaskvalue
         
      END DO


      
      PRINT *, ''
      PRINT *, ' 3D interpolation done!'
      PRINT *, ''

   END SUBROUTINE INTERP_3D





   SUBROUTINE extrp_hl(X2d)

      !! Extrapolating data at top/north bottom/south with persistence!

      REAL(4), DIMENSION(ni_out, nj_out), INTENT(inout) :: X2d

      INTEGER :: jj

      IF (jj_ex_top > 0) THEN
         DO jj=jj_ex_top,(nlat_inc_out+1)/2*nj_out+(1 - nlat_inc_out)/2,nlat_inc_out
            X2d(:,jj) = X2d(:,jj_ex_top-nlat_inc_out)
         END DO
      END IF

      IF (jj_ex_btm > 0) THEN
         DO jj=(nlat_inc_out + 1)/2+(1 - nlat_inc_out)/2*nj_out,jj_ex_btm,nlat_inc_out
            X2d(:,jj) = X2d(:,jj_ex_btm+nlat_inc_out)
            !PRINT *, ' LOLO: jj_out=', jj, 'values of jj_out=',jj_ex_btm+nlat_inc_out
         END DO
      END IF

   END SUBROUTINE extrp_hl



END MODULE MOD_INTERP
