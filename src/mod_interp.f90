MODULE MOD_INTERP

   USE mod_conf       !* important parameters, namelist and misc routines
   USE mod_manip      !* misc. manipulation of 2D arrays
   USE mod_drown      !* extrapolation over masked surfaces
   USE mod_akima_2d   !* Akima method algorithm
   USE mod_bilin_2d   !* Bi-linear method (for handling irregular source grids)
   USE mod_akima_1d   !* 1D Akima method for vertical interpolation

   !USE io_ezcdf ; !LOLOdebug

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: INTERP_2D, INTERP_3D

CONTAINS

   SUBROUTINE INTERP_2D()

      !! ================
      !! 2D INTERPOLATION
      !! ================

      !! lon-aranging or lat-flipping field
      IF ( nlat_inc_in == -1 ) CALL FLIP_UD_2D(data_in)
      IF ( nlon_inc_in == -1 ) CALL LONG_REORG_2D(i_chg_lon, data_in)

      mask_in = mask_in_b    ! re-filling the mask with trusted values...

      IF ( ldrown ) THEN
         !! Extrapolate sea values over land :
         IF ( lmout ) THEN
            CALL DROWN(ewper, data_in, mask_in(:,:,1),  nb_inc=100, nb_smooth=0)
         ELSE            
            CALL DROWN(ewper, data_in, mask_in(:,:,1))
         END IF
      ELSE
         PRINT *, '-------------------'
         PRINT *, 'DROWN NOT CALLED!!!'
         PRINT *, '-------------------'
      END IF

      !! Call interpolation procedure :
      !! ------------------------------

      !! First, may apply a smoothing on "data_in" in case target grid is much coarser than the source grid!
      !! => CALL SMOOTH
      
      SELECT CASE(cmethod)

      CASE('akima')
         CALL akima_2d(ewper, lon_in, lat_in, data_in, lon_out, lat_out, data_out)

      CASE('bilin')
         CALL bilin_2d(ewper, lon_in, lat_in, data_in, lon_out, lat_out, data_out, cpat)

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

      !!   Masking result
      IF ( lmout ) THEN
         WHERE (mask_out(:,:,1) == 0)  data_out = rmaskvalue
      ELSE
         rmaskvalue = 0.
      ENDIF

   END SUBROUTINE INTERP_2D



   SUBROUTINE INTERP_3D()

      !! ================
      !! 3D INTERPOLATION
      !! ================

      INTEGER :: ji, jj, jk, jk_last=0
      REAL(8) :: zmax_in, zmax_out

      !! Interpolation of each of the nk_in levels onto the target grid
      !! --------------------------------------------------------------
      DO jk = 1, nk_in

         PRINT *, ' Level : ', jk


         IF ( nlat_inc_in == -1 ) CALL FLIP_UD_2D(data3d_in(:,:,jk))
         IF ( nlon_inc_in == -1 ) CALL LONG_REORG_2D(i_chg_lon, data3d_in(:,:,jk))

         IF ( ldrown ) THEN
            !! Extrapolate sea values over land :
            !WRITE(6,*) '*** Extrapolating source data over continents with DROWN on level jk =', jk
            IF ( lmout ) THEN
               CALL DROWN(ewper, data3d_in(:,:,jk), mask_in(:,:,jk),  nb_inc=100, nb_smooth=0)
            ELSE            
               CALL DROWN(ewper, data3d_in(:,:,jk), mask_in(:,:,jk))
            END IF
            !LOLOdebug:
            !WRITE(cfdbg,'("data_in_drowned_lev",i2.2,".nc")') jk
            !IF ( jk == 17 ) THEN
            !   !CALL PRTMASK(data3d_in(:,:,jk), trim(cfdbg), 'data')
            !   CALL PRTMASK(real(mask_in(:,:,jk),4), 'lolo.nc', 'data')
            !   STOP
            !END IF
            !LOLOdebug.
         ELSE
            PRINT *, '-------------------'
            PRINT *, 'DROWN NOT CALLED!!!'
            PRINT *, '-------------------'
         END IF

         !! First, may apply a smoothing on "data3d_in(:,:,jk)" in case target grid is much coarser than the source grid!
         !! => CALL SMOOTH

         SELECT CASE(cmethod)

         CASE('akima')
            CALL akima_2d(ewper, lon_in,  lat_in,  data3d_in(:,:,jk), &
               &              lon_out, lat_out, data3d_tmp(:,:,jk))

         CASE('bilin')
            CALL bilin_2d(ewper, lon_in,  lat_in,  data3d_in(:,:,jk), &
               &              lon_out, lat_out, data3d_tmp(:,:,jk), cpat)

         CASE('no_xy')
            data3d_tmp(:,:,jk) = data3d_in(:,:,jk)

         CASE DEFAULT
            PRINT *, 'Interpolation method "', trim(cmethod), '" is unknown!!!'; STOP
         END SELECT

         IF (lregout) CALL extrp_hl(data3d_tmp(:,:,jk))


      END DO

      !! Masking everything :
      data3d_out(:,:,:) = rmaskvalue

      IF ( .NOT. (trim(cf_x_out)  == 'spheric') ) THEN

         depth_in  = ABS(depth_in)
         depth_out = ABS(depth_out)

         zmax_in  = MAXVAL(depth_in)
         zmax_out = MAXVAL(depth_out)

         IF ( zmax_out > zmax_in ) THEN
            !! Must find the last target level less deep than zmax_in
            jk_last = 1
            DO WHILE ( jk_last < nk_out )
               IF ( depth_out(jk_last+1) > zmax_in ) EXIT
               jk_last = jk_last + 1
            END DO
         END IF

         !! Need to perform a vertical interpolation from data3d_tmp to data3d_out :
         DO jj = 1, nj_out
            DO ji = 1, ni_out

               nlev = nk_out
               IF ( (mask_out(ji,jj,1) == 1) .OR. (.NOT. lmout) ) THEN
                  IF ( lmout ) THEN  ! adapting nlev if masking target
                     nlev = 1
                     !! RD while loop causes seg fault in debug
                     DO jk=1,nk_out
                        IF ( mask_out(ji,jj,jk) == 1 ) nlev = nlev + 1
                     ENDDO
                     nlev = nlev - 1
                  END IF

                  CALL AKIMA_1D( REAL(depth_in(:)      ,4), data3d_tmp(ji,jj,:),    &
                     &           REAL(depth_out(1:nlev),4), data3d_out(ji,jj,1:nlev)  )

               END IF
            END DO
         END DO

         !! Assuring persistance at the bottom if target depth goes deeper that source depth
         IF ( (jk_last > 0).AND.(jk_last < nk_out) ) THEN
            DO jk = jk_last-1, nk_out
               data3d_out(:,:,jk) = data3d_out(:,:,jk_last-2)
            END DO
         END IF

      ELSE
         data3d_out = data3d_tmp ! target levels are same than source levels
      END IF                     ! => no vertical interpolation required...

      !! avoid working with 3D arrays as a whole : produce SEGV on some machines (small)
      !! RD : replaced out of bounds values by vmin/vmax not rmaskvalue
      !! as it induced spval instead of zero on BGC fields
      DO jk=1,nk_out
         WHERE ((data3d_out(:,:,jk) > vmax).and.(data3d_out(:,:,jk) /= rmaskvalue)) &
            &   data3d_out(:,:,jk) = vmax
         WHERE ((data3d_out(:,:,jk) < vmin).and.(data3d_out(:,:,jk) /= rmaskvalue)) &
            &   data3d_out(:,:,jk) = vmin
         IF ( lmout ) data3d_out(:,:,jk) = data3d_out(:,:,jk)*mask_out(:,:,jk) + (1. - mask_out(:,:,jk))*rmaskvalue

      END DO

   END SUBROUTINE INTERP_3D





   SUBROUTINE extrp_hl(X2d)

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
         END DO
      END IF

   END SUBROUTINE extrp_hl



END MODULE MOD_INTERP
