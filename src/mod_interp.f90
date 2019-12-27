MODULE MOD_INTERP

   USE mod_conf       !* important parameters, namelist and misc routines
   USE mod_manip      !* misc. manipulation of 2D arrays
   USE mod_drown      !* extrapolation over masked surfaces
   USE mod_akima_2d   !* Akima method algorithm
   USE mod_bilin_2d   !* Bi-linear method (for handling irregular source grids)
   USE mod_akima_1d   !* 1D Akima method for vertical interpolation
   USE mod_grids

   USE mod_nemotools, ONLY: lbc_lnk
   USE io_ezcdf,      ONLY: DUMP_FIELD, TEST_XYZ ; !LOLOdebug

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: INTERP_2D, INTERP_3D

CONTAINS

   SUBROUTINE INTERP_2D( jt, Nt)

      !! ================
      !! 2D INTERPOLATION
      !! ================

      INTEGER, INTENT(in) :: jt, Nt

      INTEGER :: i1,j1, i2,j2, jtr, ji, jj
      INTEGER :: ni_src_x, nj_src_x
      CHARACTER(len=2) :: ctype
      INTEGER, PARAMETER :: n_extd = 4    ! source grid extension
      REAL(8), DIMENSION(:,:), ALLOCATABLE       :: X1, Y1
      REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: X2, Y2
      REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: data_src_x
      REAL(8), DIMENSION(:,:), ALLOCATABLE, SAVE :: X1_x, Y1_x
      REAL(8), DIMENSION(4),                SAVE :: xy_range_src

      !! lon-aranging or lat-flipping field
      IF ( nlat_icr_src == -1 ) CALL FLIP_UD(data_src)
      IF ( nlon_icr_src == -1 ) CALL LONG_REORG_2D(i_chg_lon, data_src)

      mask_src = mask_src_b    ! re-filling the mask with trusted values...

      IF ( l_drown_src ) THEN
         !! Extrapolate sea values over land :
          IF ( idrown%l_msk_chg ) CALL CREATE_LSM( 'source', cf_lsm_src, cv_lsm_src, mask_src(:,:,1),  xfield=data_src )
         CALL DROWN(ewper_src, data_src, mask_src(:,:,1), nb_inc=idrown%np_penetr, nb_smooth=idrown%nt_smooth)

         !loloIF ( l_save_drwn ) data_src_drowned(:,:,1) = data_src(:,:)

      ELSE
         PRINT *, '-------------------'
         PRINT *, 'DROWN NOT CALLED!!!'
         PRINT *, '-------------------'
      END IF

      IF ( ismooth > 0 ) THEN
         !! First, may apply a smoothing on "data_src" in case target grid is much coarser than the source grid!
         WRITE(6,'("     --- ",a,": smoothing ",i4," times!")') TRIM(cv_src), ismooth
         PRINT *, ' Smoothing '//TRIM(cv_src)//'!', ismooth, ' times'
         CALL SMOOTHER(ewper_src, data_src,  nb_smooth=ismooth, msk=mask_src(:,:,1), l_exclude_mask_points=.true.)
      END IF


      !! East-West and South-North border extension of source domain:
      !! ------------------------------------------------------------

      !!LOLO: some stuffs here can be done only at first call! #FIXME
      
      ctype = TEST_XYZ(lon_src, lat_src, data_src)
         
      !! Source extended domain:
      ni_src_x = ni_src + n_extd
      nj_src_x = nj_src + n_extd
      
      IF ( jt == 1 ) THEN
         PRINT *, '  *** allocating data_src_x, X1_x, Y1_x:', ni_src_x,nj_src_x, '(',ni_src, nj_src,')'
         ALLOCATE ( data_src_x(ni_src_x,nj_src_x), X1_x(ni_src_x,nj_src_x), Y1_x(ni_src_x,nj_src_x) )
         PRINT *, '      => allocation done!'
      END IF
      !PRINT *, '  *** allocating data_src_x:', ni_src_x,nj_src_x, '(',ni_src, nj_src,')'      
      !ALLOCATE ( data_src_x(ni_src_x,nj_src_x) )

      

      IF ( jt == 1 ) THEN
         !! Going to work with 2D longitude,latitude arrays (source domain) => X1, Y1 (regardless regularity of source grid)
         IF ( ctype == '1d' ) THEN
            FORALL (jj = 1:nj_src) X1_x(3:ni_src_x-2,jj+2) = lon_src(:,1)
            FORALL (ji = 1:ni_src) Y1_x(ji+2,3:nj_src_x-2) = lat_src(:,1)
         ELSE
            X1_x(3:ni_src_x-2,3:nj_src_x-2) = lon_src(:,:)
            Y1_x(3:ni_src_x-2,3:nj_src_x-2) = lat_src(:,:)
         END IF
      END IF
      
      data_src_x(3:ni_src_x-2,3:nj_src_x-2) = REAL( data_src(:,:) , 8)      
      !CALL DUMP_FIELD(REAL(X1_x,4), 'X1_ext.nc', 'var')
      !CALL DUMP_FIELD(REAL(Y1_x,4), 'Y1_ext.nc', 'var')
      !CALL DUMP_FIELD(REAL(data_src_x,4), 'data_src_x.nc', 'var')
      

      !LOLO: dodgy to call the same array as input and output :
      CALL FILL_EXTRA_BANDS(jt, ewper_src, X1_x(3:ni_src_x-2,3:nj_src_x-2), Y1_x(3:ni_src_x-2,3:nj_src_x-2), data_src_x(3:ni_src_x-2,3:nj_src_x-2), &
         &                            X1_x, Y1_x, data_src_x,  is_orca_grid=i_orca_src)
      
      !CALL DUMP_FIELD(REAL(data_src_x,4), 'data_src_ext.nc', 'var') !,   xlon=X1_x, xlat=Y1_x)


      IF ( jt == 1 ) THEN
         !!LOLO #FIXME put all these squaring-coordinates stuff somewhere else, like in mod_grid..., also known from a module,
         !! Make target lat and lon rectangle and not vector!
         ctype = '00'
         ctype = TEST_XYZ(lon_trg, lat_trg, data_trg)
         ALLOCATE ( X2(ni_trg,nj_trg) , Y2(ni_trg,nj_trg) )
         IF ( ctype == '1d' ) THEN
            FORALL (jj=1:nj_trg) X2(:,jj) = lon_trg(:,1)
            FORALL (ji=1:ni_trg) Y2(ji,:) = lat_trg(:,1)
         ELSE
            X2 = lon_trg ; Y2 = lat_trg
         END IF
      END IF





      !! Call interpolation procedure :
      !! ------------------------------

      SELECT CASE(cmethod)

      CASE('akima')


         IF ( jt == 1 ) THEN

            xy_range_src(:) = (/ MINVAL(X1_x),MAXVAL(X1_x) , MINVAL(Y1_x),MAXVAL(Y1_x) /)  !! Checking if the target grid does not overlap source grid :

            IF ( Nthrd > 1 ) THEN
               PRINT *, ' *** Allocating "ixy_mapping" => ', ni_trg, nj_trg, 2
            END IF
            ALLOCATE ( ixy_mapping(ni_trg, nj_trg, 2) )
            ixy_mapping(:,:,:) = 0

            !$OMP PARALLEL DO
            DO jtr = 1, Nthrd
!$             PRINT *, ' Running "FIND_SRC_CELL" on OMP thread #', INT(jtr,1)
               CALL FIND_SRC_CELL( X1_x, Y1_x, xy_range_src, X2(i_b_l(jtr):i_b_r(jtr),:), Y2(i_b_l(jtr):i_b_r(jtr),:), ixy_mapping(i_b_l(jtr):i_b_r(jtr),:,:) )
               !CALL FIND_SRC_CELL( X1_x, Y1_x, xy_range_src, X2(:,j_b_l(jtr):j_b_r(jtr)), Y2(:,j_b_l(jtr):j_b_r(jtr)), ixy_mapping(:,j_b_l(jtr):j_b_r(jtr),:) )
            END DO
            !$OMP END PARALLEL DO
            
         END IF ! IF ( jt == 1 )
         
         
!$       PRINT *, ''
         
         !$OMP PARALLEL DO
         DO jtr = 1, Nthrd
!$          PRINT *, ' Running "AKIMA_2D" on OMP thread #', INT(jtr,1)
            !! ewper_src useless now that extension is done above???? right?
            CALL AKIMA_2D( ewper_src, X1_x, Y1_x, data_src_x, &
               &           X2(i_b_l(jtr):i_b_r(jtr),:), Y2(i_b_l(jtr):i_b_r(jtr),:), data_trg(i_b_l(jtr):i_b_r(jtr),:), &
               &           ixy_mapping(i_b_l(jtr):i_b_r(jtr),:,:), jtr )
         END DO
         !$OMP END PARALLEL DO

         
      CASE('bilin')

         IF ( jt == 1 ) THEN
            ALLOCATE ( l_first_call_interp_routine(Nthrd) )
            l_first_call_interp_routine(:) = .TRUE.
            ALLOCATE ( IMETRICS(ni_trg,nj_trg,3), RAB(ni_trg,nj_trg,2), IPB(ni_trg,nj_trg), DIST_NP(ni_trg,nj_trg) )
            !ALLOCATE ( IMETRICS(ni_trg,nj_trg,3), RAB(ni_trg,nj_trg,2), IPB(ni_trg,nj_trg) )
         END IF
         
         !$OMP PARALLEL DO
         DO jtr = 1, Nthrd
!$          PRINT *, ' Running "bilin_2d" on OMP thread #', INT(jtr,1)
            CALL bilin_2d( ewper_src, X1_x, Y1_x, REAL(data_src_x,4), &
               &                      X2(i_b_l(jtr):i_b_r(jtr),:), Y2(i_b_l(jtr):i_b_r(jtr),:), data_trg(i_b_l(jtr):i_b_r(jtr),:), &
               &           cpat, jtr,  mask_domain_trg=IGNORE)
         END DO
         !$OMP END PARALLEL DO


      CASE('no_xy')
         WRITE(6,*) 'ERROR (mod_interp.f90): method "no_xy" makes no sense for 2D interp!'
         STOP

      CASE DEFAULT
         PRINT *, 'Interpolation method ', cmethod, ' is unknown!!!' ; STOP
      END SELECT

      !DEALLOCATE ( data_src_x )
      !DEALLOCATE ( X1_x, Y1_x )   !! lolo a la fin !!!


      !! If target grid extends too much in latitude compared to source grid, need to
      !! extrapolate a bit at bottom and top of the domain :
      IF (l_reg_trg) CALL extrp_hl(data_trg)


      !! Applying bound corrections
      !! LOLO should update mask_trg accordingly !!!?
      WHERE ( data_trg > vmax )  data_trg = rmiss_val
      WHERE ( data_trg < vmin )  data_trg = rmiss_val

      !lolo:
      !! If overshoot of latitudes between target and source domain (target has higher values than source):
      !! => apply a drown because the relevant areas were masked (even if lmout==false)!
      IF (jj_ex_btm > 0) THEN
         CALL DROWN(ewper_trg, data_trg, mask_trg(:,:,1), nb_inc=idrown%np_penetr, nb_smooth=idrown%nt_smooth)
      END IF
      !lolo.

      IF ( ibx_xtra_sm(0) > 0 ) THEN
         WRITE(6,'("     --- ",a,": post-interp extra smoothing ",i4," times!")') TRIM(cv_out), ibx_xtra_sm(0)
         i1=ibx_xtra_sm(1)
         j1=ibx_xtra_sm(2)
         i2=ibx_xtra_sm(3)
         j2=ibx_xtra_sm(4)
         PRINT *, '          => on box:', i1,j1, i2,j2
         !CALL SMOOTHER(ewper_trg, data_trg(i1:i2,j1:j2),  nb_smooth=ibx_xtra_sm(0), msk=mask_trg(i1:i2,j1:j2,1), l_exclude_mask_points=.TRUE.)
         CALL SMOOTHER(ewper_trg, data_trg(i1:i2,j1:j2),  nb_smooth=ibx_xtra_sm(0))
         PRINT *, ''
      END IF

      IF ( ismooth_out > 0 ) THEN
         WRITE(6,'("     --- ",a,": post-interp smoothing ",i4," times!")') TRIM(cv_out), ismooth_out
         CALL SMOOTHER(ewper_trg, data_trg,  nb_smooth=ismooth_out, msk=mask_trg(:,:,1), l_exclude_mask_points=.true.)
         PRINT *, ''
      END IF



      !! Masking result if requested
      IF ( lmout ) THEN
         WHERE (mask_trg(:,:,1) == 0)  data_trg = rmiss_val
      ELSE
         rmiss_val = 0.
      ENDIF

      !! If target grid is an ORCA grid, calling "lbc_lnk":
      IF ( i_orca_trg > 0 ) CALL lbc_lnk( i_orca_trg, data_trg, c_orca_trg, 1.0_8 )

      !IF ( jt == Nt ) DEALLOCATE ( X2, Y2, ixy_mapping )
            
   END SUBROUTINE INTERP_2D



   SUBROUTINE INTERP_3D()

      !! ================
      !! 3D INTERPOLATION
      !! ================

      INTEGER :: i1,j1, i2,j2
      INTEGER :: ji, jj, jk, jklast=0
      REAL(8) :: zmax_src, zmax_trg

      !! Interpolation of each of the nk_src levels onto the target grid
      !! --------------------------------------------------------------

      !CALL DUMP_2D_FIELD(data3d_src(:,nj_src/2,:), '00_Slice_before_anything.tmp', 's')


      DO jk = 1, nk_src
         PRINT *, '### Preparing source field at level : ', jk

         IF ( cmethod /= 'no_xy' ) THEN !LOLO: WHY????
            IF ( nlat_icr_src == -1 ) CALL FLIP_UD(data3d_src(:,:,jk))
            IF ( nlon_icr_src == -1 ) CALL LONG_REORG_2D(i_chg_lon, data3d_src(:,:,jk))
         END IF

         IF ( l_drown_src ) THEN
            !! Extrapolate sea values over land :
            WRITE(6,'("     --- ",a,": Extrapolating source data over land at level #",i3.3)') TRIM(cv_src), jk
            PRINT *, 'LOLO: calling DROWN with: ', idrown%np_penetr, idrown%nt_smooth
            CALL DROWN(ewper_src, data3d_src(:,:,jk), mask_src(:,:,jk), nb_inc=idrown%np_penetr, nb_smooth=idrown%nt_smooth )
            !lolo IF ( l_save_drwn ) data_src_drowned(:,:,jk) = data3d_src(:,:,jk)

            !CALL DUMP_2D_FIELD(data3d_src(:,nj_src/2,:), '01_Slice_in_just_after_horiz_drown.tmp', 's')

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
            !! First, may apply a smoothing on "data_src" in case target grid is much coarser than the source grid!
            WRITE(6,'("     --- ",a,": Smoothing level #",i3.3," ",i2," times!")') TRIM(cv_src), jk, ismooth
            CALL SMOOTHER(ewper_src, data3d_src(:,:,jk),  nb_smooth=ismooth, msk=mask_src(:,:,jk), l_exclude_mask_points=.true.)
         END IF

      END DO !DO jk = 1, nk_src



      !! Need to do some sort of vertical drown downward to propagate the bottom
      !! value (last water pomit mask==1) down into the sea-bed...
      !! => calling drown vertical slice by vertical slices (zonal vertical slices)
      PRINT *, ''
      PRINT *, '### Extrapolating bottom value of source field downward into sea-bed!'
      DO jj = 1, nj_src
         !CALL DUMP_2D_FIELD(data3d_src(:,jj,:), 'slice_before', 's')
         CALL DROWN( -1, data3d_src(:,jj,:), mask_src(:,jj,:), nb_inc=100, nb_smooth=5 )
         !CALL DUMP_2D_FIELD(data3d_src(:,jj,:), 'slice_after', 's')
      END DO
      PRINT *, '   => Done!'
      PRINT *, ''

      !CALL DUMP_2D_FIELD(data3d_src(:,nj_src/2,:), '02_Slice_in_just_after_vert_drown.tmp', 's')




      !! Prevent last level to f-word shit up when not a single water point (last level is only rock!)
      IF ( SUM(mask_src(:,:,nk_src)) == 0) data3d_src(:,:,nk_src) = data3d_src(:,:,nk_src-1) ! persistence!

      !LOLOdebug:
      !DO jk = nk_src/2, nk_src
      !   WRITE(cfdbg,'("data_to_be_srcterpolated_lev",i2.2,".nc")') jk ; PRINT *, ' *** saving ', TRIM(cfdbg)
      !   CALL DUMP_2D_FIELD(data3d_src(:,:,jk), TRIM(cfdbg), cv_src)
      !END DO
      !LOLOdebug.


      PRINT *, ''
      PRINT *, ' 3D field prepared at all levels, ready to be interpolated...'
      PRINT *, ''



      !CALL DUMP_2D_FIELD(data3d_src(:,nj_src/2,:), '03_Slice_in_before_interp.tmp', 's')

      !! Now! 3D input field ready to be interpolated...
      DO jk = 1, nk_src

         IF (TRIM(cmethod) /= 'no_xy' ) PRINT *, '  *** interpolating at level ', jk

         SELECT CASE(TRIM(cmethod))

         CASE('akima')
            PRINT *, 'mod_interp.f90 not yet!'
            STOP
            !CALL akima_2d(ewper_src, lon_src,  lat_src,  data3d_src(:,:,jk), &
            !   &              lon_trg, lat_trg, data3d_tmp(:,:,jk), 1)
            IF ( trim(ctype_z_src) == 'z' ) THEN
               !! we don't need horizontal interpolation, all levels are flat
               depth_src_trgt2d(:,:,jk) = depth_src(1,1,jk)
            ELSE
               !! input is sigma, layers are non-flat
               !CALL akima_2d(ewper_src, lon_src,  lat_src, depth_src(:,:,jk),       &
               !   &              lon_trg, lat_trg,   depth_src_trgt2d(:,:,jk), 1 )
            ENDIF

            !CASE('bilin')
            !   CALL bilin_2d(ewper_src, lon_src,  lat_src,  data3d_src(:,:,jk), &
            !      &              lon_trg, lat_trg, data3d_tmp(:,:,jk), cpat, 1)

            IF ( trim(ctype_z_src) == 'z' ) THEN
               !! we don't need horizontal interpolation, all levels are flat
               depth_src_trgt2d(:,:,jk) = depth_src(1,1,jk)
               IF ( l_identical_levels ) depth_src_trgt2d(:,:,jk) = depth_trg(1,1,jk) !lolo
            ELSE
               !! input is sigma, layers are non-flat
               !CALL bilin_2d(ewper_src, lon_src,  lat_src, depth_src(:,:,jk), &
               !   &              lon_trg, lat_trg,   depth_src_trgt2d(:,:,jk), cpat, 1)
            ENDIF

         CASE('no_xy')
            PRINT *, '  *** copying at level ', jk
            data3d_tmp(:,:,jk) = data3d_src(:,:,jk)
            IF ( TRIM(ctype_z_src) == 'z' ) THEN
               !! we don't need horizontal interpolation, all levels are flat
               depth_src_trgt2d(:,:,jk) = depth_src(1,1,jk)
            ELSE
               !! input is sigma, layers are non-flat
               depth_src_trgt2d(:,:,jk) = depth_src(:,:,jk)
            ENDIF

         CASE DEFAULT
            PRINT *, 'Interpolation method "', trim(cmethod), '" is unknown!!!'; STOP
         END SELECT

         IF (l_reg_trg) CALL extrp_hl(data3d_tmp(:,:,jk))

      END DO !DO jk = 1, nk_src

      PRINT *, ''

      data3d_trg(:,:,:) = rmiss_val ! Masking everything



      !! Time for vertical interpolation


      !IF ( .NOT. (TRIM(cf_x_trg)  == 'spheric') ) THEN  !LOLO WTF was this???
      IF ( .NOT. l_identical_levels ) THEN  ! ( l_identical_levels always false for S-coordinates...)

         depth_src  = ABS(depth_src)
         depth_trg = ABS(depth_trg)
         depth_src_trgt2d  = ABS(depth_src_trgt2d)

         !! Need to perform a vertical interpolation from data3d_tmp to data3d_trg :


         PRINT *, ''
         IF ( (TRIM(ctype_z_src) == 'z').AND.(TRIM(ctype_z_trg) == 'z') ) THEN

            !! Go for the vectorial routine...
            PRINT *, ' *** CALLING AKIMA_1D_3D for vertical interpolation !!!'
            WRITE(6,'("      => ",i4.4," x ",i4.4," x ",i3.3," to ",i4.4," x ",i4.4," x ",i3.3)') &
               &       SIZE(data3d_tmp,1), SIZE(data3d_tmp,2), SIZE(depth_src_trgt2d(1,1,:),1), &
               &       ni_trg, nj_trg, SIZE(depth_trg(1,1,:),1)

            CALL AKIMA_1D( depth_src_trgt2d(1,1,:), data3d_tmp, depth_trg(1,1,:), data3d_trg(:,:,:), rmiss_val )

            PRINT *, '  --- AKIMA_1D_3D done!'
            PRINT *, ''

         ELSE


            DO jj = 1, nj_trg
               DO ji = 1, ni_trg

                  !! RD dev notes : we need to make sure that the depth vector for both in and out
                  !! are from smallest to largest value so that persistance works
                  IF ( trim(ctype_z_src) == 'sigma' ) THEN
                     CALL FLIP_UD(depth_src_trgt2d(ji,jj,:))
                     CALL FLIP_UD(data3d_tmp(ji,jj,:))
                  ENDIF

                  IF ( trim(ctype_z_trg) == 'sigma' ) THEN
                     CALL FLIP_UD(depth_trg(ji,jj,:))
                  ENDIF

                  !! RD dev notes : we compare the depth from source depth vector and target depth vector
                  !! at the same horizontal location : compare depth_src_trgt2d and depth_trg
                  zmax_src  = MAXVAL(depth_src_trgt2d(ji,jj,:))
                  zmax_trg = MAXVAL(depth_trg(ji,jj,:))

                  IF ( zmax_trg > zmax_src ) THEN
                     !! Must find the last target level less deep than zmax_src
                     jklast = 1
                     DO WHILE ( jklast < nk_trg )
                        IF ( depth_trg(ji,jj,jklast+1) > zmax_src ) EXIT
                        jklast = jklast + 1
                     END DO
                  ELSE
                     jklast = nk_trg
                  END IF

                  IF ( (mask_trg(ji,jj,1) == 1) .OR. (.NOT. lmout) ) THEN
                     IF ( lmout ) THEN  ! adapting nlev if masking target
                        nlev = 1
                        !! RD while loop causes seg fault in debug
                        DO jk=1,nk_trg
                           IF ( mask_trg(ji,jj,jk) == 1 ) nlev = nlev + 1
                        ENDDO
                        nlev = nlev - 1
                     ELSE
                        nlev = jklast
                     END IF
                     !!
                     IF ( (MOD(ji,100)==0).AND.(MOD(jj,100)==0) ) PRINT *, '  ... calling AKIMA_1D for ji_trg,jj_trg =', ji,jj
                     !!
                     CALL AKIMA_1D( depth_src_trgt2d(ji,jj,:),data3d_tmp(ji,jj,:),  &
                        &           depth_trg(ji,jj,1:nlev), data3d_trg(ji,jj,1:nlev))
                     !!
                     !! Apply persistance at the bottom if target depth goes deeper that source depth
                     !! RD dev notes : I think the indices were off. If jklast is the last target level
                     !! that can be properly computed then we want to apply persistance to jklast + 1 to nk_trg
                     !! btw, interp from z to sigma works slightly better without this on my test case
                     IF ( (jklast > 0).AND.(jklast < nk_trg) ) THEN
                        DO jk = jklast+1, nk_trg
                           data3d_trg(ji,jj,jk) = data3d_trg(ji,jj,jklast)
                        END DO
                     END IF

                     !! RD dev notes : when interpolating to sigma, need to reverse again arrays
                     IF ( trim(ctype_z_trg) == 'sigma' ) THEN
                        CALL FLIP_UD(depth_trg(ji,jj,:))
                        CALL FLIP_UD(data3d_trg(ji,jj,:))
                     ENDIF

                  END IF
               END DO
            END DO
         END IF  !IF ( (TRIM(ctype_z_src) == 'z').AND.(TRIM(ctype_z_trg) == 'z') ) THEN

      ELSE

         WRITE(6,*) '  *** skipping vertical interpolation as source and target vertical levels are identical!'
         data3d_trg = data3d_tmp ! target levels are same than source levels

      END IF  !IF ( .NOT. l_identical_levels )      ! => no vertical interpolation required...


      !! avoid working with 3D arrays as a whole : produce SEGV on some machines (small)
      !! RD : replaced out of bounds values by vmin/vmax not rmiss_val
      !! as it induced spval instead of zero on BGC fields
      DO jk=1,nk_trg
         !! LOLO: make the following consistent with vmin/vmax "bound correction" in 2D case????
         WHERE ((data3d_trg(:,:,jk) > vmax).and.(data3d_trg(:,:,jk) /= rmiss_val)) &
            &   data3d_trg(:,:,jk) = vmax
         WHERE ((data3d_trg(:,:,jk) < vmin).and.(data3d_trg(:,:,jk) /= rmiss_val)) &
            &   data3d_trg(:,:,jk) = vmin


         IF ( ibx_xtra_sm(0) > 0 ) THEN
            WRITE(6,'("     --- ",a,": post-interp extra smoothing ",i4," times at level ",i3.3,"!")') TRIM(cv_out), ibx_xtra_sm(0), jk
            i1=ibx_xtra_sm(1)
            j1=ibx_xtra_sm(2)
            i2=ibx_xtra_sm(3)
            j2=ibx_xtra_sm(4)
            PRINT *, '          => on box:', i1,j1, i2,j2
            CALL SMOOTHER(ewper_trg, data3d_trg(i1:i2,j1:j2,jk),  nb_smooth=ibx_xtra_sm(0))
            PRINT *, ''
         END IF

         IF ( ismooth_out > 0 ) THEN
            WRITE(6,'("     --- ",a,": post-interp smoothing ",i4," times at level ",i3.3,"!")') TRIM(cv_out), ismooth_out, jk
            CALL SMOOTHER(ewper_trg, data3d_trg(:,:,jk),  nb_smooth=ismooth_out, msk=mask_trg(:,:,jk), l_exclude_mask_points=.TRUE.)
            PRINT *, ''
         END IF

         !! If target grid is an ORCA grid, calling "lbc_lnk":
         IF ( i_orca_trg > 0 ) CALL lbc_lnk( i_orca_trg, data3d_trg(:,:,jk), c_orca_trg, 1.0_8 )

         IF ( lmout ) data3d_trg(:,:,jk) = data3d_trg(:,:,jk)*mask_trg(:,:,jk) + (1. - mask_trg(:,:,jk))*rmiss_val

      END DO



      PRINT *, ''
      PRINT *, ' 3D interpolation done!'
      PRINT *, ''

   END SUBROUTINE INTERP_3D





   SUBROUTINE extrp_hl(X2d)

      !! Extrapolating data at top/north bottom/south with persistence!

      REAL(4), DIMENSION(ni_trg, nj_trg), INTENT(inout) :: X2d

      INTEGER :: jj

      IF (jj_ex_top > 0) THEN
         DO jj=jj_ex_top,(nlat_icr_trg+1)/2*nj_trg+(1 - nlat_icr_trg)/2,nlat_icr_trg
            X2d(:,jj) = X2d(:,jj_ex_top-nlat_icr_trg)
         END DO
      END IF

      IF (jj_ex_btm > 0) THEN
         DO jj=(nlat_icr_trg + 1)/2+(1 - nlat_icr_trg)/2*nj_trg,jj_ex_btm,nlat_icr_trg
            X2d(:,jj) = X2d(:,jj_ex_btm+nlat_icr_trg)
            !PRINT *, ' LOLO: jj_trg=', jj, 'values of jj_trg=',jj_ex_btm+nlat_icr_trg
         END DO
      END IF

   END SUBROUTINE extrp_hl



END MODULE MOD_INTERP
