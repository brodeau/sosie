MODULE MOD_GRIDS

   USE mod_conf
   USE mod_scoord

   IMPLICIT none

   PRIVATE

   PUBLIC :: SRC_DOMAIN, &
      &      TRG_DOMAIN, &
      &      IS_ORCA_NORTH_FOLD, &
      &      TERMINATE

   INTEGER :: &
      &     ji, jj, jt0, jz0, &
      &     n1, n2, n3, nrec,   &
      &     if0, iv0

   REAL(wpl) :: rmv, dx, dy
   LOGICAL :: lmval


CONTAINS




   SUBROUTINE SRC_DOMAIN()

      USE io_ezcdf  !: => only to get time array !LB???

      INTEGER :: jt

      !! Determine input dimensions from input file :
      CALL know_dim_in()

      !! Allocate input arrays with input dimensions :
      ALLOCATE ( data_in(ni_in,nj_in), mask_in(ni_in,nj_in,nk_in), data_in_b(ni_in,nj_in),    &
         &     mask_in_b(ni_in,nj_in,nk_in), vt0(Ntr0), vt(Ntr) )

      vt(:) = 0.
      
      IF ( lregin ) THEN
         ALLOCATE ( lon_in(ni_in,1),     lat_in(nj_in,1)     )
      ELSE
         ALLOCATE ( lon_in(ni_in,nj_in), lat_in(ni_in,nj_in) )
      END IF

      IF ( l_int_3d ) ALLOCATE ( data3d_in(ni_in,nj_in,nk_in), depth_in(ni_in,nj_in,nk_in) )
      !! if input grid is in terrain-following, need array for input bathy
      IF ( l_int_3d .AND. trim(ctype_z_in) == 'sigma' ) ALLOCATE ( bathy_in(ni_in,nj_in) )
      !! Filling time array :
      IF ( ltime  )  THEN
         IF ( lct ) THEN       ! time is being controlled
            DO jt = 1, Ntr
               vt(jt) = t0 + t_stp*REAL(jt)
            END DO
            nb_att_t = 1
            vatt_info_t(:)%cname = 'null'
            vatt_info_t(1)%cname = 'units'
            vatt_info_t(1)%itype = 2 ! char
            vatt_info_t(1)%val_char = 'unknown'
            vatt_info_t(1)%ilength = LEN('unknown')
            !!
         ELSE    ! we get time vector and its attributes into input file
            CALL GETVAR_1D(cf_in, cv_t_in, vt0) ;  vt(:) = vt0(j_start:j_stop)
            CALL GETVAR_ATTRIBUTES(cf_in, cv_t_in,  nb_att_t, vatt_info_t)
         END IF
      END IF

      CALL get_src_conf()

      max_lon_in = maxval(lon_in);   max_lat_in = maxval(lat_in)
      min_lon_in = minval(lon_in);   min_lat_in = minval(lat_in)

      PRINT *, ''
      gt_orca_in = IS_ORCA_NORTH_FOLD( lat_in , cname_long=trim(cv_lon_in) )
      i_orca_in  = gt_orca_in%ifld_nord
      IF ( i_orca_in == 4 ) PRINT *, ' Input grid is an ORCA grid with north-pole T-point folding!'
      IF ( i_orca_in == 6 ) PRINT *, ' Input grid is an ORCA grid with north-pole F-point folding!'
      PRINT *, ''

   END SUBROUTINE SRC_DOMAIN




   SUBROUTINE TRG_DOMAIN()

      !lolo: remettre: USE io_ezcdf, ONLY: getvar_attributes
      USE io_ezcdf

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: xdum

      IF ( TRIM(cmethod) == 'no_xy' ) THEN
         lregout = lregin
         !! Geting them from source file:
         cv_lon_out = cv_lon_in
         cv_lat_out = cv_lat_in
         CALL GETVAR_ATTRIBUTES(cf_x_in, cv_lon_in, nb_att_lon, vatt_info_lon)
         CALL GETVAR_ATTRIBUTES(cf_x_in, cv_lat_in, nb_att_lat, vatt_info_lat)
      END IF

      CALL know_dim_out()

      !! Allocate output arrays with output dimensions :
      ALLOCATE ( mask_out(ni_out,nj_out,nk_out), data_out(ni_out,nj_out), IGNORE(ni_out,nj_out) )

      IF ( .NOT. lmout ) mask_out(:,:,:) = 1 ; !lolo

      IF ( lregout ) THEN
         ALLOCATE ( lon_out(ni_out, 1) , lat_out(nj_out, 1), lon_out_b(ni_out, 1) )
      ELSE
         ALLOCATE ( lon_out(ni_out, nj_out) , lat_out(ni_out, nj_out), lon_out_b(ni_out, nj_out) )
      END IF

      IF (l_int_3d) THEN
         ALLOCATE ( depth_in_tmp(ni_out, nj_out, nk_in), depth_out(ni_out,nj_out,nk_out) )
         ALLOCATE ( data3d_tmp(ni_out, nj_out, nk_in), data3d_out(ni_out,nj_out,nk_out) )
      END IF

      !! If output grid is terrain-following, then allocate for bathy_out
      IF ( l_int_3d .AND. trim(ctype_z_out) == 'sigma' ) ALLOCATE ( bathy_out(ni_out,nj_out) )
      IF ( l_int_3d .AND. trim(ctype_z_out) == 'sigma' ) ALLOCATE ( Cs_rho(nk_out), Sc_rho(nk_out) )

      IF ( l_int_3d .AND. trim(ctype_z_out) == 'sigma' ) CALL compute_scoord_2_4(ssig_out,nk_out,Cs_rho,Sc_rho)

      jj_ex_top = 0 ; jj_ex_btm = 0

      IF ( TRIM(cmethod) == 'no_xy' ) THEN

         CALL get_trg_conf()

         !! The very few things we have to do if no 2D horizontal interpolation needed:
         max_lat_out  = max_lat_in
         nlat_inc_out = nlat_inc_in

         CALL rd_vgrid(nk_out, cf_z_out, cv_z_out, depth_out(1,1,:))
         WRITE(6,*) ''; WRITE(6,*) 'Output Depths ='; PRINT *, depth_out(1,1,:) ; WRITE(6,*) ''
         CALL GETVAR_ATTRIBUTES(cf_z_out, cv_z_out,  nb_att_z, vatt_info_z)
         DO ji=1,ni_out
            DO jj=1,nj_out
               depth_out(ji,jj,:) = depth_out(1,1,:)
            ENDDO
         ENDDO

         lon_out   = lon_in
         lon_out_b = lon_in
         lat_out   = lat_in

      ELSE

         !! Normal case => 2D horizontal interpolation will be performed !

         CALL get_trg_conf()

         max_lat_out = maxval(lat_out) ;   min_lat_out = minval(lat_out) ;
         PRINT*,''
         WRITE(6,*) 'Latitude max on  input grid =', max_lat_in
         WRITE(6,*) 'Latitude max on output grid =', max_lat_out
         PRINT*,''
         WRITE(6,*) 'Latitude min on  input grid =', min_lat_in
         WRITE(6,*) 'Latitude min on output grid =', min_lat_out
         PRINT*,''

         !! Building IGNORE mask:
         IGNORE = 1
         IF ( (.NOT.l_glob_lon_wize).OR.(.NOT.l_glob_lat_wize) ) ALLOCATE ( xdum(ni_out,nj_out) )
         IF ( .NOT.l_glob_lon_wize ) THEN
            WRITE(*,'("  => going to disregard points of target domain with lon < ",f7.2," and lon > ",f7.2)'), lon_min_1,lon_max_1
            IF ( lregout ) THEN
               DO jj=1,nj_out
                  xdum(:,jj) = lon_out(:,1)
               END DO
            ELSE
               xdum = lon_out
            END IF
            xdum = SIGN(1._8,180._8-xdum)*MIN(xdum,ABS(xdum-360._8)) ! like lon_out but between -180 and +180 !
            !CALL DUMP_2D_FIELD(REAL(xdum,4), 'lon_out_180-180.nc', 'lon') ; !#lolo
            WHERE ( xdum < lon_min_1 ) IGNORE=0
            WHERE ( xdum > lon_max_1 ) IGNORE=0
         END IF
         IF ( .NOT.l_glob_lat_wize ) THEN
            WRITE(*,'("  => going to disregard points of target domain with lat < ",f7.2," and lat > ",f7.2)'), min_lat_in, max_lat_in
            IF ( lregout ) THEN
               DO ji=1,ni_out
                  xdum(ji,:) = lat_out(:,1)
               END DO
            ELSE
               xdum = lat_out
            END IF
            WHERE ( xdum < min_lat_in ) IGNORE=0
            WHERE ( xdum > max_lat_in ) IGNORE=0
            PRINT *, ''            
         END IF
         IF ( (.NOT.l_glob_lon_wize).OR.(.NOT.l_glob_lat_wize) ) DEALLOCATE ( xdum )
         

         !! Is target latitude increasing with j : 1 = yes | -1 = no
         nlat_inc_out = 1

         IF ( lregout ) THEN
            IF (lat_out(1,1) > lat_out(nj_out,1)) nlat_inc_out = -1
         ELSE
            IF (lat_out(1,1) > lat_out(1,nj_out)) nlat_inc_out = -1
         END IF

         !! Find first jj that exits max_lat_in --> a mettre ailleurs!
         IF ( lregout ) THEN
            IF  ( max_lat_out >  max_lat_in ) THEN
               jj_ex_top = (nlat_inc_out + 1)/2 + (1 - nlat_inc_out)/2*nj_out
               DO WHILE ( lat_out(jj_ex_top,1) < max_lat_in )
                  jj_ex_top = jj_ex_top + nlat_inc_out
               END DO
            END IF
            !! Find first ji that exits max_lat_in --> a mettre ailleurs!
            IF  ( min_lat_out <  min_lat_in ) THEN
               jj_ex_btm = (nlat_inc_out + 1)/2*nj_out + (1 - nlat_inc_out)/2
               DO WHILE ( lat_out(jj_ex_btm,1) > min_lat_in )
                  jj_ex_btm = jj_ex_btm - nlat_inc_out
               END DO
            END IF
         END IF

         IF (jj_ex_top > 0) THEN
            jj_ex_top = jj_ex_top - nlat_inc_out
            WRITE(6,*) 'jj_ex_top =', jj_ex_top, lat_out(jj_ex_top,1)
         END IF

         IF (jj_ex_btm > 0) THEN
            jj_ex_btm = jj_ex_btm + nlat_inc_out
            WRITE(6,*) 'jj_ex_btm =', jj_ex_btm, lat_out(jj_ex_btm,1)
            !IF (jj_ex_btm > 0) jj_ex_btm = jj_ex_btm + 3*nlat_inc_out !lolo Works but 3 points is too much!!!
         END IF

      END IF

      !! LOLO: masking target mask
      IF (jj_ex_btm > 0) THEN
         DO jj=(nlat_inc_out + 1)/2+(1 - nlat_inc_out)/2*nj_out,jj_ex_btm,nlat_inc_out
            mask_out(:,jj,:) = 0
         END DO
         !debug: CALL DUMP_2D_FIELD(REAL(mask_out(:,:,1),4), 'mask_out.nc', 'lsm')
      END IF

      ! Type of target grid (only matters if ORCA grid...)
      gt_orca_out = IS_ORCA_NORTH_FOLD( lon_out , cname_long=trim(cv_lon_out) )
      i_orca_out = gt_orca_out%ifld_nord
      c_orca_out = gt_orca_out%cgrd_type
      IF ( i_orca_out == 4 ) PRINT *, ' Target grid is an ORCA '//c_orca_out//' grid with north-pole T-point folding!'
      IF ( i_orca_out == 6 ) PRINT *, ' Target grid is an ORCA '//c_orca_out//' grid with north-pole F-point folding!'
      PRINT *, ''

   END SUBROUTINE TRG_DOMAIN





   SUBROUTINE TERMINATE()

      !data_in = 0.0; mask_in = 0.0; data_in_b = 0.0; lon_in = 0.0; lat_in = 0.0
      !mask_in_b = 0.0 ; vt0 = 0.0; vt = 0.0

      DEALLOCATE ( data_in, mask_in, data_in_b, lon_in, lat_in, mask_in_b, vt0, vt )

      !mask_out = 0.0; data_out = 0.0; lat_out = 0.0; lon_out = 0.0; lon_out_b = 0.0

      DEALLOCATE ( mask_out, data_out, lat_out, lon_out, lon_out_b )

      IF (l_int_3d) THEN
         !data3d_in = 0.0 ; depth_in = 0.0
         DEALLOCATE ( data3d_in, depth_in )
         !data3d_tmp = 0.0;  depth_out = 0.0;  data3d_out = 0.0
         DEALLOCATE ( data3d_tmp, depth_out, data3d_out )
         DEALLOCATE ( depth_in_tmp )
         IF (trim(ctype_z_in) == 'sigma' )  DEALLOCATE ( bathy_in )
         IF (trim(ctype_z_out) == 'sigma' ) DEALLOCATE ( bathy_out )
      END IF

   END SUBROUTINE TERMINATE






   !! LOCAL SUBROUTINES
   !! ~~~~~~~~~~~~~~~~~

   SUBROUTINE get_src_conf()

      USE io_ezcdf
      USE mod_manip
      USE mod_scoord
      !! Local :
      INTEGER :: ji, jj, jk
      REAL    :: rval_thrshld, lon_min_2, lon_max_2
      LOGICAL :: l_loc1, l_loc2
      REAL(wpl), DIMENSION(:,:,:), ALLOCATABLE :: z3d_tmp

      !! Getting grid on source domain:
      CALL rd_grid(-1, lregin, cf_x_in, cv_lon_in, cv_lat_in, lon_in, lat_in)

      IF ( l_int_3d ) THEN
         IF ( trim(ctype_z_in) == 'sigma' ) THEN
            !! read input bathymetry
            CALL GETVAR_2D(if0,iv0,cf_bathy_in, cv_bathy_in, 0, 0, 0, bathy_in(:,:))
            !! compute 3D depth_in for input variable from bathy and sigma parameters
            CALL depth_from_scoord(ssig_in, bathy_in, ni_in, nj_in, nk_in, depth_in)
         ELSEIF ( trim(ctype_z_in) == 'z' ) THEN
            !! in z case, the depth vector is copied at each grid-point
            CALL rd_vgrid(nk_in, cf_z_in, cv_z_in, depth_in(1,1,:))
            WRITE(6,*) ''; WRITE(6,*) 'Input Depths ='; PRINT *, depth_in(1,1,:) ; WRITE(6,*) ''
            DO ji=1,ni_in
               DO jj=1,nj_in
                  depth_in(ji,jj,:) = depth_in(1,1,:)
               ENDDO
            ENDDO
         ELSE
            PRINT*,''; PRINT *, 'Not a valid input vertical coordinate' ; PRINT*,''
         ENDIF

         IF ( trim(ctype_z_in) == 'z' ) THEN
            PRINT*,''; WRITE(6,*) 'Input has z coordinates and depth vector is =', depth_in(1,1,:); PRINT*,''
         ELSEIF ( trim(ctype_z_in) == 'sigma' ) THEN
            PRINT*,''; WRITE(6,*) 'Input has sigma coordinates and depth range is ', MINVAL(depth_in), &
               &                           ' to ', MAXVAL(depth_in) ; PRINT*,''
         ELSE
            PRINT*,''; WRITE(6,*) 'You should not see this' ; STOP
         ENDIF
      END IF

      lon_min_1 = MINVAL(lon_in)
      lon_max_1 = MAXVAL(lon_in)
      PRINT *, ' *** Minimum longitude on source domain before reorg. : ', REAL(lon_min_1,4)
      PRINT *, ' *** Maximum longitude on source domain before reorg. : ', REAL(lon_max_1,4)

      IF ( lregin ) THEN
         !! Fixing input 1D longitude:
         CALL FIX_LONG_1D(ni_in, lon_in(:,1), nlon_inc_in, i_chg_lon)
         !! Fixing input 1D latitude:
         CALL FIX_LATD_1D(nj_in, lat_in(:,1), nlat_inc_in)
      ELSE
         WHERE ( lon_in < 0. )  lon_in = lon_in + 360.
      END IF

      lon_min_2 = MINVAL(lon_in)
      lon_max_2 = MAXVAL(lon_in)
      PRINT *, ' *** Minimum longitude on source domain now: ', lon_min_2
      PRINT *, ' *** Maximum longitude on source domain now: ', lon_max_2

      ! lolo: IMPROVE! This is disgusting:
      l_loc1 = (lon_min_1 <  0.).AND.(lon_min_1 > -175.).AND.(lon_max_1 >  0. ).AND.(lon_max_1 <  175.)
      l_loc2 = (lon_min_2 >= 0.).AND.(lon_min_2 <   2.5).AND.(lon_max_2 >357.5).AND.(lon_max_2 <= 360.)
      IF ( (.NOT. l_loc1).AND.(l_loc2) ) THEN
         l_glob_lon_wize = .TRUE.
         PRINT *, 'Looks like global setup (longitude-wise at least...)'
      ELSE
         PRINT *, 'Looks like regional setup (longitude-wise at least...)'
         l_glob_lon_wize = .FALSE.
         WRITE(*,'("  => going to disregard points of target domain with lon < ",f7.2," and lon > ",f7.2)'), lon_min_1,lon_max_1
      END IF
      PRINT *, ''
      
      l_glob_lat_wize = .TRUE.
      IF ( MAXVAL(lat_in) < 88. ) l_glob_lat_wize =.FALSE.
            
      !! Getting land-sea mask on source domain
      mask_in(:,:,:) = 1 ! by default everything is considered sea (helps for smoothing when no LSM)

      IF ( ldrown ) THEN

         IF ( TRIM(cf_lsm_in) == 'missing_value' ) THEN

            WRITE(6,*) 'Opening land-sea mask from missing_value on input data!'
            CALL CHECK_4_MISS(cf_in, cv_in, lmval, rmv, ca_missval)
            IF ( .NOT. lmval ) THEN
               PRINT *, 'ERROR (get_src_conf of mod_grids.f90) : '//TRIM(cv_in)//' has no missing value attribute!'
               PRINT *, '      (in '//TRIM(cf_in)//')'
               STOP
            END IF
            ALLOCATE ( z3d_tmp(ni_in,nj_in,nk_in) )
            !! Read data field (at time 1 if time exists) :
            IF ( ltime  ) jt0 = 1
            IF ( l_int_3d ) THEN
               CALL GETVAR_3D(if0, iv0, cf_in, cv_in, Ntr,      jt0, z3d_tmp)
            ELSE
               IF ( l3d ) jz0 = jplev
               CALL GETVAR_2D(if0, iv0, cf_in, cv_in, Ntr, jz0, jt0, z3d_tmp(:,:,1))
            END IF
            mask_in = 1
            WHERE ( z3d_tmp == rmv ) mask_in = 0
            DEALLOCATE ( z3d_tmp )

         ELSEIF ( (TRIM(cf_lsm_in) == 'val+').OR.(TRIM(cf_lsm_in) == 'val-').OR.(TRIM(cf_lsm_in) == 'value') ) THEN
            READ(cv_lsm_in,*) rval_thrshld
            IF (TRIM(cf_lsm_in) == 'val+')  WRITE(6,*) ' Land-sea mask is defined from values >=', rval_thrshld
            IF (TRIM(cf_lsm_in) == 'val-')  WRITE(6,*) ' Land-sea mask is defined from values <=', rval_thrshld
            IF (TRIM(cf_lsm_in) == 'value') WRITE(6,*) ' Land-sea mask is defined from values ==', rval_thrshld
            ALLOCATE ( z3d_tmp(ni_in,nj_in,nk_in) )
            !! Read data field (at time 1 if time exists) :
            IF ( ltime  ) jt0 = 1
            IF ( l_int_3d ) THEN
               CALL GETVAR_3D(if0, iv0, cf_in, cv_in, Ntr,      jt0, z3d_tmp)
            ELSE
               IF ( l3d ) jz0 = jplev
               CALL GETVAR_2D(if0, iv0, cf_in, cv_in, Ntr, jz0, jt0, z3d_tmp(:,:,1))
            END IF
            mask_in = 1
            IF (TRIM(cf_lsm_in) == 'val+') THEN
               WHERE ( z3d_tmp >= rval_thrshld ) mask_in = 0
            END IF
            IF (TRIM(cf_lsm_in) == 'val-') THEN
               WHERE ( z3d_tmp <= rval_thrshld ) mask_in = 0
            END IF
            IF (TRIM(cf_lsm_in) == 'value') THEN
               WHERE ( z3d_tmp == rval_thrshld ) mask_in = 0
            END IF
            DEALLOCATE ( z3d_tmp )

         ELSEIF ((TRIM(cf_lsm_in)=='nan').OR.(TRIM(cf_lsm_in)=='NaN')) THEN
            !! NaN values are considered mask!
            WRITE(6,*) ' Land-sea mask is defined from values larger than', vmax
            ALLOCATE ( z3d_tmp(ni_in,nj_in,nk_in) )
            !! Read data field (at time 1 if time exists) :
            IF ( ltime  ) jt0 = 1
            IF ( l_int_3d ) THEN
               CALL GETVAR_3D(if0, iv0, cf_in, cv_in, Ntr,      jt0, z3d_tmp)
            ELSE
               IF ( l3d ) jz0 = jplev
               CALL GETVAR_2D(if0, iv0, cf_in, cv_in, Ntr, jz0, jt0, z3d_tmp(:,:,1))
            END IF
            mask_in = 1
            DO jk = 1, nk_in
               DO jj =  1, nj_in
                  DO ji =  1, ni_in
                     IF ( ISNAN(z3d_tmp(ji,jj,jk)) ) mask_in(ji,jj,jk) = 0
                  END DO
               END DO
            END DO
            !debug: CALL DUMP_2D_FIELD(REAL(mask_in(:,:,1),4), 'mask_in.nc', 'mask')
            DEALLOCATE ( z3d_tmp )

         ELSE

            !! We are reading source land-sea mask from a file:
            !! -----------------------------------------------
            !! checking dimension:
            CALL DIMS(cf_lsm_in, cv_lsm_in, n1, n2, n3, nrec)

            IF ( l3d .AND. ( jplev > 1 ) ) THEN
               IF ( n3 == nk_in ) THEN
                  WRITE(6,*) 'Opening 3D land-sea mask on source grid for level', jplev
                  PRINT *, trim(cv_lsm_in)
                  !! if terrain-following, open the 2d mask, not sure interp one single level works
                  IF (trim(ctype_z_in) == 'sigma' ) THEN
                     CALL GETMASK_2D(cf_lsm_in, cv_lsm_in, mask_in(:,:,1))
                  ELSEIF (trim(ctype_z_in) == 'z' ) THEN
                     CALL GETMASK_2D(cf_lsm_in, cv_lsm_in, mask_in(:,:,1), jlev=jplev)
                  ELSE
                     STOP
                  ENDIF
               ELSE
                  WRITE(6,*) 'PROBLEM! You want to interpolate level', jplev
                  WRITE(6,*) 'but your source land-sea mask is not 3D!'
                  WRITE(6,*) '=> set ldrown to false in the namelist'
                  WRITE(6,*) 'If you want to "drown" a level other than the surface,'
                  WRITE(6,*) 'please provide a 3D input land-sea mask'
                  STOP
               END IF
            ELSEIF ( l_int_3d ) THEN
               !! RD: dims reads n3 = -1 on lsm, needs to force n3 to ROMS Nlevels
               !! RD: maybe there is a more elegant way to do that
               IF (trim(ctype_z_in) == 'sigma' ) n3 = ssig_in%Nlevels
               IF ( n3 == nk_in ) THEN
                  !! if terrain-following, read 2D mask and copy it on all levels
                  IF (trim(ctype_z_in) == 'sigma' ) THEN
                     WRITE(6,*) 'Opening 2D land-sea mask file on source grid: ', trim(cf_lsm_in)
                     CALL GETMASK_2D(cf_lsm_in, cv_lsm_in, mask_in(:,:,1))
                     DO jz0=2,nk_in
                        mask_in(:,:,jz0) = mask_in(:,:,1)
                     ENDDO
                  ELSEIF (trim(ctype_z_in) == 'z' ) THEN
                     WRITE(6,*) 'Opening 3D land-sea mask file on source grid, ', trim(cv_lsm_in)
                     CALL GETMASK_3D(cf_lsm_in, cv_lsm_in, mask_in)
                  ELSE
                     STOP
                  ENDIF
               ELSE
                  WRITE(6,*) 'We need to open the 3D source land-sea mask,'
                  WRITE(6,*) 'but the vertical dimension of it does not match!'
                  STOP
               END IF
            ELSE
               WRITE(6,*) 'Opening land-sea mask file on source grid, ', trim(cv_lsm_in)
               CALL GETMASK_2D(cf_lsm_in, cv_lsm_in, mask_in(:,:,1))
            END IF
         END IF

      END IF ! IF ( ldrown )

      !! Need to modify the mask if lon or lat have been modified :
      IF ( nlat_inc_in == -1 ) CALL FLIP_UD_3D(mask_in)
      IF ( nlon_inc_in == -1 ) CALL LONG_REORG_3D(i_chg_lon, mask_in)




   END SUBROUTINE get_src_conf


   SUBROUTINE get_trg_conf()

      USE io_ezcdf
      USE mod_scoord

      REAL(wpl), DIMENSION(:,:,:), ALLOCATABLE :: z3d_tmp


      IF ( (lregout).AND.(TRIM(cf_x_out) == 'spheric') ) THEN

         !! Building target grid:
         READ(cv_lon_out,*) dx ; READ(cv_lat_out,*) dy
         cv_lon_out = 'lon'           ; cv_lat_out = 'lat'
         WRITE(6,*) '  * dx, dy =', dx, dy
         WRITE(6,*) '  * ni_out, nj_out =', ni_out, nj_out ;  PRINT*,''
         DO ji = 1, ni_out
            lon_out(ji,1) = dx/2.0 + dx*REAL(ji - 1 , 8)
         END DO
         DO jj = 1, nj_out
            lat_out(jj,1) = -90 + dy/2.0 + dy*REAL(jj - 1 , 8)
         END DO

         WRITE(6,*) ''; WRITE(6,*) 'Target Longitude array (deg.E):'; PRINT *, lon_out; WRITE(6,*) ''
         WRITE(6,*) 'Target Latitude array (deg.N):';  PRINT *, lat_out; WRITE(6,*) ''; WRITE(6,*) ''


      ELSE
         !! Getting target grid from netcdf file:
         CALL rd_grid(ivect, lregout, cf_x_out, cv_lon_out, cv_lat_out, lon_out, lat_out)
         !!
         IF ( lregout ) THEN
            WRITE(6,*) ''; WRITE(6,*) 'Target Longitude array (deg.E):'; PRINT *, lon_out; WRITE(6,*) ''
            WRITE(6,*) 'Target Latitude array (deg.N):';  PRINT *, lat_out; WRITE(6,*) ''; WRITE(6,*) ''
         END IF
      END IF

      !! Netcdf attributes for longitude and latitude:
      IF ( TRIM(cf_x_out) == 'spheric' ) THEN
         nb_att_lon = 1
         vatt_info_lon(:)%cname = 'null'
         vatt_info_lon(1)%cname = 'units'
         vatt_info_lon(1)%itype = 2 ! char
         vatt_info_lon(1)%val_char = 'degrees_east'
         vatt_info_lon(1)%ilength = LEN('degrees_east')
         nb_att_lat = 1
         vatt_info_lat(:)%cname = 'null'
         vatt_info_lat(1)%cname = 'units'
         vatt_info_lat(1)%itype = 2 ! char
         vatt_info_lat(1)%val_char = 'degrees_west'
         vatt_info_lat(1)%ilength = LEN('degrees_west')
      ELSE
         !! Geting them from target file:
         CALL GETVAR_ATTRIBUTES(cf_x_out, cv_lon_out, nb_att_lon, vatt_info_lon)
         CALL GETVAR_ATTRIBUTES(cf_x_out, cv_lat_out, nb_att_lat, vatt_info_lat)
      END IF

      lon_out_b = lon_out
      WHERE ( lon_out < 0. ) lon_out = lon_out + 360.

      IF ( l_int_3d ) THEN
         IF ( trim(cf_x_out)  == 'spheric') THEN
            cf_z_out = cf_z_in ;  cv_z_out = cv_z_in         !Important
         END IF

         IF ( trim(ctype_z_out) == 'sigma' ) THEN
            !! read bathy for output grid
            CALL GETVAR_2D(if0,iv0,cf_bathy_out, cv_bathy_out, 0, 0, 0, bathy_out(:,:))
            !! compute target depth on output grid from bathy_out and ssig_out params
            CALL depth_from_scoord(ssig_out, bathy_out, ni_out, nj_out, ssig_out%Nlevels, depth_out)
            CALL GETVAR_ATTRIBUTES(cf_bathy_out, cv_bathy_out,  nb_att_z, vatt_info_z)
         ELSEIF (trim(ctype_z_out) == 'z' ) THEN
            !! depth vector copied on all grid-points
            CALL rd_vgrid(nk_out, cf_z_out, cv_z_out, depth_out(1,1,:))
            !WRITE(6,*) ''; WRITE(6,*) 'Output Depths ='; PRINT *, depth_out(1,1,:) ; WRITE(6,*) ''
            CALL GETVAR_ATTRIBUTES(cf_z_out, cv_z_out,  nb_att_z, vatt_info_z)
            DO ji=1,ni_out
               DO jj=1,nj_out
                  depth_out(ji,jj,:) = depth_out(1,1,:)
               ENDDO
            ENDDO
         ELSE
            PRINT*,''; WRITE(6,*) 'Not a valid output vertical coordinate' ; STOP
            !!
         ENDIF

         !RD fix this
         !         IF (trim(ctype_z_out) == 'z' ) THEN
         !            PRINT*,''; WRITE(6,*) 'Output Depths ='; PRINT *, depth_out(1,1,:) ; PRINT*,''
         !         ELSEIF ( trim(ctype_z_out) == 'sigma' ) THEN
         !            PRINT*,''; WRITE(6,*) 'Output on sigma coordinates' ; PRINT*,''
         !         ENDIF

      END IF

      !RD fix this
      !         CALL rd_vgrid(nk_out, cf_z_out, cv_z_out, depth_out)
      !         WRITE(6,*) ''; WRITE(6,*) 'Output Depths ='; PRINT *, depth_out ; WRITE(6,*) ''
      !         CALL GETVAR_ATTRIBUTES(cf_z_out, cv_z_out,  nb_att_z, vatt_info_z)

      !!  Getting target mask (mandatory doing 3D interpolation!)
      IF ( lmout .OR. l_int_3d ) THEN

         IF ( (l3d).AND.(jplev > 1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) '****************************************************************'
            WRITE(6,*) 'We do not know output mask yet, since it is at a given depth!'
            WRITE(6,*) '--> next version of SOSIE!'
            WRITE(6,*) 'So we do not mask output file'
            WRITE(6,*) '****************************************************************'
            WRITE(6,*) ''

            mask_out = 1

         ELSE

            IF ( TRIM(cf_lsm_out) == 'missing_value' ) THEN
               
               WRITE(6,*) 'Opening target land-sea mask from missing_value!'
               CALL CHECK_4_MISS(cf_x_out, cv_lsm_out, lmval, rmv, ca_missval)
               IF ( .NOT. lmval ) THEN
                  PRINT *, 'ERROR lili (get_trg_conf of mod_grids.f90) : '//TRIM(cv_lsm_out)//' has no missing value attribute!'
                  PRINT *, '      (in '//TRIM(cf_x_out)//')'
                  STOP
               END IF
               !!
               ALLOCATE ( z3d_tmp(ni_out,nj_out,nk_out) )
               if0 = 0 ; iv0 = 0   ! Read data field (at time 1 if time exists)
               IF ( l_int_3d ) THEN
                  CALL GETVAR_3D(if0, iv0, cf_x_out, cv_lsm_out, 1, 1,    z3d_tmp)
               ELSE
                  CALL GETVAR_2D(if0, iv0, cf_x_out, cv_lsm_out, 1, 1, 1, z3d_tmp(:,:,1))
               END IF
               mask_out = 1
               WHERE ( z3d_tmp == rmv ) mask_out = 0
               DEALLOCATE ( z3d_tmp )

            ELSE

               IF ( l_int_3d ) THEN
                  IF ( ((trim(cf_lsm_out) == '').OR.(trim(cv_lsm_out) == '')) ) THEN
                     WRITE(6,*) 'WARNING: no target 3D land-sea mask provided (cf_lsm_out)!'
                     mask_out = 1
                  ELSE
                     !! select coord type
                     IF ( trim(ctype_z_out) == 'sigma' ) THEN
                        WRITE(6,*) 'Opening 2D land-sea mask file on target grid: ', trim(cf_lsm_out)
                        !! read 2D mask for output and make it 3D
                        CALL GETMASK_2D(cf_lsm_out, cv_lsm_out, mask_out(:,:,1))
                        DO jz0=2,nk_out
                           mask_out(:,:,jz0) = mask_out(:,:,1)
                        ENDDO
                     ELSEIF ( trim(ctype_z_out) == 'z' ) THEN
                        WRITE(6,*) 'Opening 3D land-sea mask file on target grid: ',trim(cf_lsm_out)
                        WRITE(6,*) '             => name mask : ',trim(cv_lsm_out)
                        CALL GETMASK_3D(cf_lsm_out, cv_lsm_out, mask_out(:,:,:))
                        WRITE(6,*) ''
                     ELSE
                        STOP
                     ENDIF
                  END IF
               ELSE
                  WRITE(6,*) 'Opening 2D land-sea mask file on target grid: ', trim(cf_lsm_out)
                  CALL GETMASK_2D(cf_lsm_out, cv_lsm_out, mask_out(:,:,1))
               END IF
               WRITE(6,*) ''
            END IF

         END IF

      END IF

   END SUBROUTINE get_trg_conf







   SUBROUTINE rd_grid(iv, lreg, cfgrd, cvx, cvy, rlon, rlat)

      USE io_ezcdf

      !!   AUTHOR:
      !!   -------
      !!   Laurent Brodeau, 2015
      !!
      !!   may 2007: modified by Pierre Mathiot to detect and handle
      !!               regular 2D input longitude and regular 2D input
      !!               latitude when the input grid is declared as irregular.
      !!               Allowing akima interpolation.

      !! Arguments:
      INTEGER,                   INTENT(in)  :: iv !: iv = -1 means we're handling input grid
      !!                                                /= -1 means target grid
      LOGICAL,                   INTENT(in)  :: lreg
      CHARACTER(len=400),        INTENT(in)  :: cfgrd
      CHARACTER(len=80),         INTENT(in)  :: cvx, cvy
      REAL(8),   DIMENSION(:,:), INTENT(out) :: rlon, rlat

      !! Local
      LOGICAL :: lreg2d, l2dyreg_x, l2dyreg_y
      CHARACTER(len=8) :: cdomain, clreg
      INTEGER :: &
         &     idx, lx, ly, &
         &     ii, ij, if1, iv1, &
         &     ilx1, ily1, ilz1,   &
         &     ilx2, ily2, ilz2
      REAL(8),   DIMENSION(:,:), ALLOCATABLE :: zrlon, zrlat

      IF ( iv == -1 ) THEN
         cdomain = 'source'
         clreg   = 'lregin'
         idx = 1
      ELSE
         cdomain = 'target'
         clreg   = 'lregout'
         idx = 2
      END IF

      IF ( lreg ) THEN
         IF ( (size(rlon,2) /= size(rlat,2)).OR.(size(rlon,2) /= 1) ) THEN
            WRITE(6,*) 'ERROR 1: rd_grid (prepare.F90)!' ; STOP
         END IF
         lx = size(rlon,1) ; ly = size(rlat,1)
      ELSE
         IF ( (size(rlon,1) /= size(rlat,1)).OR.(size(rlon,1) /= size(rlat,1)) ) THEN
            WRITE(6,*) 'ERROR 2: rd_grid (prepare.F90)!' ; STOP
         END IF
         lx = size(rlon,1) ; ly = size(rlon,2)
      END IF
      rlon = 0. ; rlat = 0.

      !! Checking the dimension of longitude variable
      CALL DIMS(cfgrd, cvx, ilx1, ily1, ilz1, nrec)
      CALL DIMS(cfgrd, cvy, ilx2, ily2, ilz2, nrec)

      WRITE(6,*) ''


      IF (lreg) THEN
         !!         -------------------------------
         !!            R E G U L A R   G R I D :
         !!         -------------------------------
         WRITE(6,*) 'Opening regular grid ', trim(cvx), ' , ', trim(cvy), &
            &                 ' in file ', trim(cfgrd), ' .'

         l2dyreg_x = .FALSE.
         IF ( (ilz1 /= -1).or.(ily1 /= -1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) 'Warning! '//trim(cdomain)//' longitude ', trim(cvx), ' is not 1D!'
            WRITE(6,*) 'In the file ', trim(cfgrd)
            WRITE(6,*) ' => maybe you should specify '//trim(clreg)//' = F in the namelist?'
            WRITE(6,*) ' => anyway we assume you know what you are doing!'
            l2dyreg_x = .TRUE.
         END IF

         l2dyreg_y = .FALSE.
         IF ( (ilz2 /= -1).or.(ily2 /= -1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) 'Warning! '//trim(cdomain)//' latitude ', trim(cvy), ' is not 1D!'
            WRITE(6,*) 'In the file ', trim(cfgrd)
            WRITE(6,*) ' => maybe you should specify '//trim(clreg)//' = F in the namelist?'
            WRITE(6,*) ' => anyway we assume you know what you are doing!'
            l2dyreg_y = .TRUE.
         END IF

         IF ( l2dyreg_x .AND. (.NOT. l2dyreg_y) ) THEN
            WRITE(6,*) 'Error! '//trim(cdomain)//' longitude and latidude do not agree in shape (1D vs 2D)!' ; STOP
         END IF

         IF ( l2dyreg_x .AND. l2dyreg_y ) THEN
            WRITE(6,*) ' '
            WRITE(6,*) ' =================================================================================================='
            WRITE(6,*) ' *** Assuming that '//trim(cdomain)//' grid is regular even though longitude and latidude are 2D!'
            WRITE(6,*) '                     (because you set '//trim(clreg)//'=.TRUE. in the namelist)'
            WRITE(6,*) ' =================================================================================================='
            l_2d_grid_yet_regular(idx) = .TRUE.
            WRITE(6,*) ' '
         END IF


         IF ( l_2d_grid_yet_regular(idx) ) THEN
            !! Getting regular 2D grid :
            ALLOCATE ( zrlon(lx,ly), zrlat(lx,ly) )
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvx, 0, 0, 0, zrlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvy, 0, 0, 0, zrlat)
            !!
            !! Checking for Regular 2D longitude and latitude
            lreg2d = .TRUE.
            !! Checking if longitude array changes with latitude
            DO ij = 2, ly
               IF ( ABS(SUM(zrlon(:,ij)-zrlon(:,1))) > 1.e-6 ) THEN
                  lreg2d = .FALSE.
                  EXIT
               END IF
            END DO
            !! Checking if latitude array changes with longitude
            DO ii = 2, lx
               IF ( ABS(SUM(zrlat(ii,:)-zrlat(1,:))) > 1.e-6 ) THEN
                  lreg2d = .FALSE.
                  EXIT
               END IF
            END DO
            IF ( lreg2d ) THEN
               WRITE(6,*) ' *** OK! You were right, '//trim(cdomain)//' longitude and latitude are 2D but the grid is regular!'
               WRITE(6,*) ''
            ELSE
               WRITE(6,*) 'Error! We checked, and 2D '//trim(cdomain)//' longitude/latitude do not qualify for a regular grid!!!!'
               WRITE(6,*) '       => so set '//trim(clreg)//' to .FALSE. !'
               STOP
            END IF
            !!
            !! Giving values to 1D arrays:
            rlon(:,1) = zrlon(:,2)
            rlat(:,1) = zrlat(2,:)
            DEALLOCATE ( zrlon, zrlat)
            !!
         ELSE
            !!
            !! Normal case: Getting regular 1D grid :
            IF ( (ilx1 /= lx).or.(ilx2 /= ly) ) THEN
               WRITE(6,*) 'Error! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',  &
                  &   trim(cvy), &
                  &   ' do not agree in dimension with configuration dimension!'
               WRITE(6,*) 'In the file ', trim(cfgrd) ; STOP
            END IF
            CALL GETVAR_1D(cfgrd, cvx, rlon(:,1))
            CALL GETVAR_1D(cfgrd, cvy, rlat(:,1))
         END IF


      ELSEIF (.NOT. lreg) THEN
         !!         ----------------------------------
         !!            I R R E G U L A R   G R I D :
         !!         ----------------------------------
         WRITE(6,*) 'Opening irregular grid ', trim(cvx), ' , ', trim(cvy), &
            &                 ' in file ', trim(cfgrd), ' .'
         WRITE(6,*) ''

         IF (ilx1 /= ilx2) THEN

            IF ( (ily1 == ily2).and.(ilz1 == ilz2).and.(ily1 == -1).and.(ilz1 == -1) ) THEN
               WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',   &
                  &       trim(cvy), ' seem to be regular (1D), check namelist!'
               WRITE(6,*) 'In the file ', trim(cfgrd)
               WRITE(6,*) '       => so maybe try to set '//trim(clreg)//' to .TRUE. in the namelist!'
               STOP
            ELSE
               WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',  &
                  &   trim(cvy), ' do not agree in dimension!'
               WRITE(6,*) 'In the file ', trim(cfgrd) ; STOP
            END IF
         END IF

         IF ( (ily1 /= ily2).or.(ilz1 /= ilz2) ) THEN
            WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',    &
               &   trim(cvy), &
               &   ' do not agree in dimension!'
            WRITE(6,*) 'In the file ', trim(cfgrd); STOP
         END IF
         IF ( (ilx1 /= lx).or.(ily1 /= ly) ) THEN
            WRITE(6,*) 'Eror! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',   &
               &   trim(cvy),   &
               &   ' do not agree in dimension with configuration dimension!'
            WRITE(6,*) 'In the file ', trim(cfgrd); STOP
         END IF


         !! Getting source longitude array at level=1 (surface) and time =1 :
         !! -----------------------------------------------------------------
         !! If 3d dimension, we chose first level
         IF ( ilz1 /= -1 ) THEN
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvx, 0, 1, 0, rlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvy, 0, 1, 0, rlat)
         ELSE
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvx, 0, 0, 0, rlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D_R8(if1, iv1, cfgrd, cvy, 0, 0, 0, rlat)
         END IF

      END IF


   END SUBROUTINE rd_grid



   SUBROUTINE rd_vgrid(nk, cfgrd, cvz, vdepth)

      USE io_ezcdf

      INTEGER,                   INTENT(in) :: nk
      CHARACTER(len=400),        INTENT(in) :: cfgrd
      CHARACTER(len=80),         INTENT(in) :: cvz
      REAL(8), DIMENSION(nk) ,  INTENT(out) :: vdepth

      CALL GETVAR_1D(cfgrd, cvz, vdepth)

   END SUBROUTINE rd_vgrid





   SUBROUTINE know_dim_in()

      USE io_ezcdf

      INTEGER :: jk0, jrec

      nk_in = 1

      !! Determine input dimensions from input file :
      CALL DIMS(cf_in, cv_in, ni_in, nj_in, jk0, jrec)

      !PRINT *, 'mod_grids.f90, cf_in, cv_in =>', TRIM(cf_in), '', TRIM(cv_in)
      !PRINT *, 'ni_in, nj_in, jk0, jrec =>', ni_in, nj_in, jk0, jrec

      IF ( nj_in == -1 ) THEN
         WRITE(6,*) 'prepare.know_dim_in => ERROR! variable ',trim(cv_in),' should be at least 2D!!!'
         STOP
      END IF


      IF ( jplev == -1 ) THEN
         !! This is the overide case! Means 2D+T !
         !! Case when we want to read a 2D+T field but someone screwed up and the record
         !! dimension is not of 'UNLIMITED' type inside the netcdf file...
         !! So it just  overides good sence and force sosie to understand that
         !! your field to interpolate is 2D with a time record
         !! (usually the case if the time record dimension in your
         !! input file is not declared as UNLIMITED => bad! :(
         WRITE(6,*) ''
         Ntr = jk0
         WRITE(6,*) 'WARNING: know_dim_in of mod_grids.f90 !!!'
         WRITE(6,*) '   => we force input field "'//TRIM(cv_in)//'" to be 2D + time !!!'
         WRITE(6,*) '   => because you specified "jplev = -1" in the namelist!'
         WRITE(6,*) '   => the time-record dimension is therefore:', Ntr
         WRITE(6,*) ''

      ELSE

         IF ( jrec == -1) ltime=.FALSE.  ! no time records

         !! 3D variable
         IF ( jk0 > 0 ) THEN

            l3d = .TRUE.
            nk_in = jk0 ; WRITE(6,*) 'Number of level into input file =', nk_in

            IF (jplev /= 0) THEN
               IF ( (jplev <= 0).OR.(jplev > nk_in) ) THEN
                  WRITE(6,*) 'Level jplev is wrong! -> number of level =', nk_in
                  WRITE(6,*) 'In the namelist, jplev =', jplev ;  STOP
               END IF
               WRITE(6,*) 'We are going to interpolate level', jplev, ' of input file'
            ELSE
               WRITE(6,*) ''
               WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
               WRITE(6,*) 'WARNING: We are going to perform a 3D interpolation !!!'
               WRITE(6,*) '      => if this is not what you want, please specify '
               WRITE(6,*) '         which level (jplev) should be interpolated.'
               WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
               WRITE(6,*) ''
               l_int_3d = .TRUE.
            END IF

            IF ( jrec > 0 ) THEN
               !! 3D variable + time
               WRITE(6,*) 'Variable is 3D + time !'
               Ntr0 = jrec
            ELSE
               !! 3D variable without time
               WRITE(6,*) 'Variable is 3D without time axis'
               Ntr0 = 1
            END IF

         END IF


         !! 2D
         !! ~~~
         IF ( jk0 == -1 ) THEN

            l3d = .FALSE.

            IF ( jrec > 0 ) THEN
               !! 2D variable with time
               WRITE(6,*) 'Variable is 2D with time axis'
               Ntr0 = jrec
            ELSE
               !! 2D variable without time
               WRITE(6,*) 'Variable is 2D without time axis'
               Ntr0 = 1
            END IF

         END IF

         IF (.NOT. ltime) Ntr0 = 1

         Ntr = Ntr0;  j_start = 1 ; j_stop = Ntr0

      END IF

      j_start = 1 ; j_stop = Ntr


      IF ( (jt1 > 0).AND.(jt2 > 0) ) THEN
         Ntr = jt2 - jt1 + 1 ;  ;  j_start = jt1 ;  j_stop = jt2
         !! jrec is the time dimension of the input file, Ntr is the length requested by the user :
         IF ( ( Ntr > jrec ).OR.(jt1 < 1).OR.(jt1 > jt2).OR.(jt2 > jrec) ) THEN
            WRITE(6,*) ''; WRITE(6,*) 'Check jt1 and jt2 in the namelist:'
            WRITE(6, '("=> the time dimension of ",a," is ", i5)') TRIM(cv_in), jrec
            WRITE(6, '("   and you specified jt1, jt2 =",i5," ,",i5)') jt1, jt2
            STOP
         END IF
      END IF

      IF ( l_int_3d .AND. trim(ctype_z_in) == 'sigma' ) nk_in = ssig_in%Nlevels ! ugly but should work

   END SUBROUTINE know_dim_in






   SUBROUTINE know_dim_out
      !!
      USE io_ezcdf
      !!
      nk_out = 1
      !!

      IF ( TRIM(cmethod) == 'no_xy' ) THEN

         ni_out = ni_in
         nj_out = nj_in

      ELSE

         !! If 3D interpolation and building spherical grid, we use levels from source grid:
         IF ( l_int_3d .AND. (TRIM(cf_x_out)  == 'spheric') ) THEN
            cf_z_out = cf_z_in ;  cv_z_out = cv_z_in
         END IF
         !!
         IF (lregout) THEN
            IF ( TRIM(cf_x_out) == 'spheric') THEN
               WRITE(6,*) ''; WRITE(6,*) 'Building regular spherical output grid!'
               READ(cv_lon_out,*) dx ; READ(cv_lat_out,*) dy
               ni_out = INT(360./dx) ; nj_out = INT(180./dy)
               GOTO 100
            ELSE
               CALL DIMS(cf_x_out, cv_lon_out, ni_out, n1, n2, nrec)
               CALL DIMS(cf_x_out, cv_lat_out, nj_out, n1, n2, nrec)
            END IF
         ELSE
            CALL DIMS(cf_x_out, cv_lon_out, ni_out, nj_out, n1, nrec)
         END IF
         !!

      END IF

100   CONTINUE


      IF ( l_int_3d ) THEN
         !!
         WRITE(6,*) ''
         IF ( TRIM(cf_x_out)  == 'spheric' ) THEN
            WRITE(6,*) 'Since we are building our own spherical target grid,'
            WRITE(6,*) 'we are going to use source levels as target levels!'
         END IF
         !!
         WRITE(6,*) ''
         WRITE(6,*) ' => we read target levels in the following file:'
         PRINT *, TRIM(cf_z_out); WRITE(6,*) ''
         IF ( trim(ctype_z_out) == 'sigma' ) THEN
            nk_out = ssig_out%Nlevels
         ELSE
            CALL DIMS(cf_z_out, cv_z_out, nk_out, n1, n2, nrec)
         ENDIF
         WRITE(6,*) 'nk_out = ', nk_out ; WRITE(6,*) ''
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(6,*) 'Target grid dimension is', ni_out,'x',nj_out,'x',nk_out
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !!
      ELSE
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(6,*) 'Target grid dimension is', ni_out,'x', nj_out
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      END IF
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
      !!
   END SUBROUTINE know_dim_out
   !!
   !!
   !!
   !!
   !!
   SUBROUTINE FIX_LONG_1D(nx, vlon, n_inc, i_chg_x)
      !!
      !!
      !! Making 1D longitude understandable
      !! ==================================
      !!
      INTEGER,                  INTENT(in)    :: nx
      REAL(8), DIMENSION(nx),  INTENT(inout) :: vlon
      INTEGER,                  INTENT(out)   :: n_inc
      INTEGER,                  INTENT(out)   :: i_chg_x
      !!
      REAL(8), DIMENSION(:), ALLOCATABLE    :: x1d_b
      INTEGER :: jx



      IF ( (vlon(1) == -180.0).AND.(vlon(nx) == 180.0) ) THEN
         !! In that very special case we do the following trick to avoid having
         !! twice the point "lon. 180" in the final longitude array
         vlon(1)  = -179.99
         vlon(nx) =  179.99
         IF ( ewper /= 0 ) THEN
            ewper = 0
            WRITE(6,*) 'prepare.F90: FIX_LONG_1D => "ewper" forced to 0 !'; WRITE(6,*) ''
         END IF
      END IF


      !! We want positive and increasing longitude !
      WHERE ( vlon < 0. )  vlon = vlon + 360.
      !!
      i_chg_x = 0 ; n_inc = 1
      DO jx = 2, nx
         IF ( vlon(jx) < vlon(jx-1) ) THEN
            IF ( i_chg_x /= 0) THEN  ! another changing sign point has already been found!
               WRITE(6,*) 'Your longitude array is a mess!'
               WRITE(6,*) ' --> Fix it (positive and increasing)!'; STOP
            ELSE
               i_chg_x = jx
               n_inc = -1
            END IF
         END IF
      END DO
      !!
      !!
      IF ( n_inc == -1 ) THEN
         !!
         ALLOCATE( x1d_b(nx) )
         x1d_b = vlon
         !!
         DO jx = i_chg_x, nx
            vlon(jx - i_chg_x + 1) = x1d_b(jx)
         END DO
         !!
         DO jx = 1, i_chg_x - 1
            vlon(nx - i_chg_x + 1 + jx) = x1d_b(jx)
         END DO
         !!
         DEALLOCATE( x1d_b )
         !!
      END IF
      !!
   END SUBROUTINE FIX_LONG_1D
   !!
   !!
   SUBROUTINE FIX_LATD_1D(ny, vlat, n_inc)
      !!
      !! Making 1D latitude increasing with jy
      !! =====================================
      !!
      INTEGER,                  INTENT(in)    :: ny
      REAL(8), DIMENSION(ny), INTENT(inout) :: vlat
      INTEGER,                  INTENT(out)   :: n_inc
      !!
      REAL(8), DIMENSION(:), ALLOCATABLE :: y1d_b
      INTEGER :: jy
      !!
      n_inc = 1
      !!
      IF ( lregin .AND. ( vlat(1) > vlat(ny) ) ) THEN
         !!
         n_inc = -1 ; ALLOCATE( y1d_b(ny) )
         !!
         WRITE(6,*) ''; WRITE(6,*) 'Latitude does not seem to increase with j'
         WRITE(6,*) '--> We reverse 1D input latitude array!'
         y1d_b = vlat
         DO jy = 1, ny
            vlat(jy) = y1d_b(ny-jy+1)
         END DO
         !!
         DEALLOCATE( y1d_b )
         !!
      END IF
      !!
   END SUBROUTINE FIX_LATD_1D



   FUNCTION IS_ORCA_NORTH_FOLD( Xtest, cname_long )

      !!----------------------------------------------------------
      !! Tell if there is a a 2/point band overlaping folding att the north pole
      !! typical of the ORCA grid
      !!
      !!  0 => not an orca grid (or unknown one)
      !!  4 => North fold T-point pivot (ex: ORCA2)
      !!  6 => North fold F-point pivot (ex: ORCA1)
      !!----------------------------------------------------------

      IMPLICIT NONE
      ! Argument
      REAL(8), DIMENSION(:,:), INTENT(in)    :: Xtest
      CHARACTER(len=*), INTENT(in), OPTIONAL :: cname_long
      !INTEGER                               :: IS_ORCA_NORTH_FOLD
      TYPE(grid_type)                        :: IS_ORCA_NORTH_FOLD
      !!
      INTEGER :: nx, ny
      CHARACTER(len=128) :: cnlon

      !! We need all this 'cname_long' stuff because with our method, there is a
      !! confusion between "Grid_U with T-fold" and "Grid_V with F-fold"
      !! => so knowing the name of the longitude array (as in namelist, and hence as
      !!    in netcdf file) might help taking the righ decision !!! UGLY!!!
      !! => not implemented yet
      cnlon = 'nav_lon' !
      IF ( PRESENT(cname_long) ) cnlon = TRIM(cname_long)




      IS_ORCA_NORTH_FOLD%ifld_nord =  0
      IS_ORCA_NORTH_FOLD%cgrd_type = 'X'

      nx = SIZE(Xtest,1)
      ny = SIZE(Xtest,2)

      IF ( ny > 3 ) THEN ! (case if called with a 1D array, ignoring...)


         IF ( SUM( Xtest(2:nx/2,ny) - Xtest(nx:nx-nx/2+2:-1,ny-2) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 4 ! T-pivot, grid_T
            IS_ORCA_NORTH_FOLD%cgrd_type = 'T'
         END IF
         !---
         IF ( SUM( Xtest(2:nx/2,ny) - Xtest(nx-1:nx-nx/2+1:-1,ny-2) ) == 0. ) THEN
            IF (TRIM(cnlon)=='glamu') THEN
               IS_ORCA_NORTH_FOLD%ifld_nord = 4 ! T-pivot, grid_T
               IS_ORCA_NORTH_FOLD%cgrd_type = 'U'
            END IF
            !! LOLO: PROBLEM == 6, V !!!
         END IF
         !---
         IF ( SUM( Xtest(2:nx/2,ny) - Xtest(nx:nx-nx/2+2:-1,ny-3) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 4 ! T-pivot, grid_V
            IS_ORCA_NORTH_FOLD%cgrd_type = 'V'
         END IF


         IF ( SUM( Xtest(2:nx/2,ny) - Xtest(nx-1:nx-nx/2+1:-1,ny-1) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 6 ! F-pivot, grid_T
            IS_ORCA_NORTH_FOLD%cgrd_type = 'T'
         END IF
         IF ( SUM( Xtest(2:nx/2,ny) - Xtest(nx-2:nx-nx/2:-1,ny-1) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 6 ! F-pivot, grid_U
            IS_ORCA_NORTH_FOLD%cgrd_type = 'U'
         END IF
         !---
         IF ( SUM( Xtest(2:nx/2,ny) - Xtest(nx-1:nx-nx/2+1:-1,ny-2) ) == 0. ) THEN
            IF (TRIM(cnlon)=='glamv') THEN
               IS_ORCA_NORTH_FOLD%ifld_nord = 6 ! F-pivot, grid_V
               IS_ORCA_NORTH_FOLD%cgrd_type = 'V'
            END IF
            !! LOLO: PROBLEM == 4, U !!!
         END IF
         !---

      END IF

   END FUNCTION IS_ORCA_NORTH_FOLD





END MODULE MOD_GRIDS
