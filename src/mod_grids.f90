MODULE  MOD_GRIDS

   USE mod_conf
   USE mod_scoord
   USE io_ezcdf
   USE mod_manip, ONLY : degE_to_degWE

   IMPLICIT none

   PRIVATE

   INTERFACE CREATE_LSM
      MODULE PROCEDURE CREATE_LSM_2D, CREATE_LSM_3D
   END INTERFACE CREATE_LSM


   PUBLIC :: SRC_DOMAIN, &
      &      TRG_DOMAIN, &
      &      IS_ORCA_NORTH_FOLD, &
      &      TERMINATE, &
      &      CREATE_LSM

   INTEGER :: &
      &     ji, jj, jk, jt0, jz0, &
      &     n1, n2, n3, nrec,   &
      &     if0, iv0

   REAL(wpl) :: rmv, dx, dy
   LOGICAL :: lmval


CONTAINS


   SUBROUTINE SRC_DOMAIN()

      INTEGER :: jt

      !! Determine source dimensions from input file :
      CALL know_dim_src()

      !! Allocate source arrays with source dimensions :
      ALLOCATE ( data_src(ni_src,nj_src), mask_src(ni_src,nj_src,nk_src), data_src_b(ni_src,nj_src),    &
         &       mask_src_b(ni_src,nj_src,nk_src), vt0(Ntr0), vt(Ntr) )
      
      IF( l_save_drwn .OR. (ixtrpl_bot>0) ) ALLOCATE ( data_src_drowned(ni_src,nj_src,nk_src) )
      
      vt(:) = 0.

      IF( l_reg_src ) THEN
         ALLOCATE ( lon_src(ni_src,1),     lat_src(nj_src,1)     )
      ELSE
         ALLOCATE ( lon_src(ni_src,nj_src), lat_src(ni_src,nj_src) )
      END IF

      IF( l_itrp_3d ) ALLOCATE ( data3d_src(ni_src,nj_src,nk_src), depth_src(ni_src,nj_src,nk_src) )
      !! if source grid is in terrain-following, need array for source bathy
      IF( l_itrp_3d .AND. trim(ctype_z_src) == 'sigma' ) ALLOCATE ( bathy_src(ni_src,nj_src) )
      !! Filling time array :
      IF( ltime  )  THEN
         IF( lct ) THEN       ! time is being controlled
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
         ELSE    ! we get time vector and its attributes into source file
            CALL GETVAR_1D(cf_src, cv_t_src, vt0)
            vt(:) = vt0(j_start:j_stop)
            CALL GETVAR_ATTRIBUTES(cf_src, cv_t_src,  nb_att_t, vatt_info_t)
         END IF
      END IF

      CALL get_src_conf()

      max_lon_src = MAXVAL(lon_src);   max_lat_src = MAXVAL(lat_src)
      min_lon_src = MINVAL(lon_src);   min_lat_src = MINVAL(lat_src)

      PRINT *, ''
      gt_orca_src = IS_ORCA_NORTH_FOLD( lat_src , cname_long=trim(cv_lon_src) )
      i_orca_src  = gt_orca_src%ifld_nord
      IF( i_orca_src == 4 ) PRINT *, ' Source grid is an ORCA grid with north-pole T-point folding!'
      IF( i_orca_src == 6 ) PRINT *, ' Source grid is an ORCA grid with north-pole F-point folding!'
      PRINT *, ''

   END SUBROUTINE SRC_DOMAIN




   SUBROUTINE TRG_DOMAIN()

      REAL(8), DIMENSION(:,:), ALLOCATABLE :: xdum
      REAL(8) :: zz


      IF( (TRIM(cmethod) == 'no_xy').OR.(l_save_drwn ) ) THEN
         CALL GETVAR_ATTRIBUTES(cf_x_src, cv_lon_src, nb_att_lon_src, vatt_info_lon_src)
         CALL GETVAR_ATTRIBUTES(cf_x_src, cv_lat_src, nb_att_lat_src, vatt_info_lat_src)
      END IF


      IF( TRIM(cmethod) == 'no_xy' ) THEN
         l_reg_trg = l_reg_src
         !! Geting them from source file:
         cv_lon_trg = cv_lon_src
         cv_lat_trg = cv_lat_src
         ni_trg = ni_src
         nj_trg = nj_src
         nb_att_lon_trg = nb_att_lon_src
         nb_att_lat_trg = nb_att_lat_src
         vatt_info_lon_trg = vatt_info_lon_src
         vatt_info_lat_trg = vatt_info_lat_src
      END IF




      CALL know_dim_trg()

      !! Allocate target arrays with target dimensions :
      ALLOCATE ( mask_trg(ni_trg,nj_trg,nk_trg), data_trg(ni_trg,nj_trg), IGNORE(ni_trg,nj_trg) )

      IF( .NOT. lmout ) mask_trg(:,:,:) = 1

      IF( l_reg_trg ) THEN
         ALLOCATE ( lon_trg(ni_trg, 1) , lat_trg(nj_trg, 1), lon_trg_b(ni_trg, 1) )
      ELSE
         ALLOCATE ( lon_trg(ni_trg, nj_trg) , lat_trg(ni_trg, nj_trg), lon_trg_b(ni_trg, nj_trg) )
      END IF

      IF(l_itrp_3d) THEN
         ALLOCATE ( depth_src_trgt2d(ni_trg, nj_trg, nk_src), depth_trg(ni_trg,nj_trg,nk_trg) )
         ALLOCATE ( data3d_tmp(ni_trg, nj_trg, nk_src), data3d_trg(ni_trg,nj_trg,nk_trg) )
      END IF

      !! If target grid is terrain-following, then allocate for bathy_trg
      IF( l_itrp_3d .AND. trim(ctype_z_trg) == 'sigma' ) ALLOCATE ( bathy_trg(ni_trg,nj_trg) )
      IF( l_itrp_3d .AND. trim(ctype_z_trg) == 'sigma' ) ALLOCATE ( Cs_rho(nk_trg), Sc_rho(nk_trg) )

      IF( l_itrp_3d .AND. trim(ctype_z_trg) == 'sigma' ) CALL compute_scoord_2_4(ssig_trg,nk_trg,Cs_rho,Sc_rho)

      jj_ex_top = 0 ; jj_ex_btm = 0

      IF( TRIM(cmethod) == 'no_xy' ) THEN

         CALL get_trg_conf()

         !! The very few things we have to do if no 2D horizontal interpolation needed:
         max_lat_trg  = max_lat_src
         nlat_icr_trg = nlat_icr_src

         CALL rd_vgrid(nk_trg, cf_z_trg, cv_z_trg, depth_trg(1,1,:))
         WRITE(6,*) ''; WRITE(6,*) 'Target Depths ='; PRINT *, depth_trg(1,1,:) ; WRITE(6,*) ''
         CALL GETVAR_ATTRIBUTES(cf_z_trg, cv_z_trg,  nb_att_z_trg, vatt_info_z_trg)
         DO ji=1,ni_trg
            DO jj=1,nj_trg
               depth_trg(ji,jj,:) = depth_trg(1,1,:)
            ENDDO
         ENDDO

         lon_trg   = lon_src
         lon_trg_b = lon_src
         lat_trg   = lat_src

      ELSE

         !! Normal case => 2D horizontal interpolation will be performed !

         CALL get_trg_conf()

         IF( (nk_src > 1) .AND. (TRIM(ctype_z_src) == 'z') .AND. (nk_src == nk_trg) ) THEN
            WRITE(6,'("Mhhh interesting, target and vertical both have ", i4.4 ," vertical levels!")') nk_src
            zz = SUM( ( 1000.*(depth_trg(1,1,:) - depth_src(1,1,:)) )**2 )
            IF( zz < 5.0 ) THEN
               PRINT *, ' => well, they are actually identical!', zz
               l_identical_levels = .TRUE.
            ELSE
               PRINT *, ' => well, they are not identical...'
               PRINT *, '      => zz =', zz
            END IF
            PRINT *, ''
         END IF

         max_lat_trg = maxval(lat_trg) ;   min_lat_trg = minval(lat_trg) ;
         PRINT*,''
         WRITE(6,*) 'Latitude max on source grid =', max_lat_src
         WRITE(6,*) 'Latitude max on target grid =', max_lat_trg
         PRINT*,''
         WRITE(6,*) 'Latitude min on source grid =', min_lat_src
         WRITE(6,*) 'Latitude min on target grid =', min_lat_trg
         PRINT*,''

         !! Building IGNORE mask:
         IGNORE = 1
         IF( (.NOT.l_glob_lon_wize).OR.(.NOT.l_glob_lat_wize) ) ALLOCATE ( xdum(ni_trg,nj_trg) )
         IF( .NOT.l_glob_lon_wize ) THEN
            WRITE(*,'("  => going to disregard points of target domain with lon < ",f7.2," and lon > ",f7.2)') lon_min_1,lon_max_1
            IF( l_reg_trg ) THEN
               DO jj=1,nj_trg
                  xdum(:,jj) = lon_trg(:,1)
               END DO
            ELSE
               xdum = lon_trg
            END IF
            xdum = degE_to_degWE(xdum)
            IF( (lon_min_1 > 180.).OR.(lon_max_1 > 180.) ) THEN
               lon_min_1 = degE_to_degWE(lon_min_1)
               lon_max_1 = degE_to_degWE(lon_max_1)
            END IF
            WHERE ( xdum < lon_min_1 ) IGNORE=0
            WHERE ( xdum > lon_max_1 ) IGNORE=0
         END IF
         IF( .NOT.l_glob_lat_wize ) THEN
            WRITE(*,'("  => going to disregard points of target domain with lat < ",f7.2," and lat > ",f7.2)') min_lat_src, max_lat_src
            IF( l_reg_trg ) THEN
               DO ji=1,ni_trg
                  xdum(ji,:) = lat_trg(:,1)
               END DO
            ELSE
               xdum = lat_trg
            END IF
            WHERE ( xdum < min_lat_src ) IGNORE=0
            WHERE ( xdum > max_lat_src ) IGNORE=0
            PRINT *, ''
         END IF

         !LOLO:
         !CALL DUMP_FIELD(REAL(IGNORE(:,:),4), 'IGNORE.nc', 'lsm')
         !STOP

         IF( (.NOT.l_glob_lon_wize).OR.(.NOT.l_glob_lat_wize) ) DEALLOCATE ( xdum )


         !! Is target latitude increasing with j : 1 = yes | -1 = no
         nlat_icr_trg = 1

         IF( l_reg_trg ) THEN
            IF(lat_trg(1,1) > lat_trg(nj_trg,1)) nlat_icr_trg = -1
         ELSE
            IF(lat_trg(1,1) > lat_trg(1,nj_trg)) nlat_icr_trg = -1
         END IF

         !! Find first jj that exits max_lat_src --> a mettre ailleurs!
         IF( l_reg_trg ) THEN
            IF  ( max_lat_trg >  max_lat_src ) THEN
               jj_ex_top = (nlat_icr_trg + 1)/2 + (1 - nlat_icr_trg)/2*nj_trg
               DO WHILE ( lat_trg(jj_ex_top,1) < max_lat_src )
                  jj_ex_top = jj_ex_top + nlat_icr_trg
               END DO
            END IF
            !! Find first ji that exits max_lat_src --> a mettre ailleurs!
            IF  ( min_lat_trg <  min_lat_src ) THEN
               jj_ex_btm = (nlat_icr_trg + 1)/2*nj_trg + (1 - nlat_icr_trg)/2
               DO WHILE ( lat_trg(jj_ex_btm,1) > min_lat_src )
                  jj_ex_btm = jj_ex_btm - nlat_icr_trg
               END DO
            END IF
         END IF

         IF(jj_ex_top > 0) THEN
            jj_ex_top = jj_ex_top - nlat_icr_trg
            WRITE(6,*) 'jj_ex_top =', jj_ex_top, lat_trg(jj_ex_top,1)
         END IF

         IF(jj_ex_btm > 0) THEN
            jj_ex_btm = jj_ex_btm + nlat_icr_trg
            WRITE(6,*) 'jj_ex_btm =', jj_ex_btm, lat_trg(jj_ex_btm,1)
            !IF(jj_ex_btm > 0) jj_ex_btm = jj_ex_btm + 3*nlat_icr_trg !lolo Works but 3 points is too much!!!
         END IF

      END IF

      !! LOLO: masking target mask
      IF(jj_ex_btm > 0) THEN
         DO jj=(nlat_icr_trg + 1)/2+(1 - nlat_icr_trg)/2*nj_trg,jj_ex_btm,nlat_icr_trg
            mask_trg(:,jj,:) = 0
         END DO
         !debug: CALL DUMP_2D_FIELD(REAL(mask_trg(:,:,1),4), 'mask_trg.nc', 'lsm')
      END IF

      ! Type of target grid (only matters if ORCA grid...)
      gt_orca_trg = IS_ORCA_NORTH_FOLD( lon_trg , cname_long=trim(cv_lon_trg) )
      i_orca_trg = gt_orca_trg%ifld_nord
      c_orca_trg = gt_orca_trg%cgrd_type
      IF( i_orca_trg == 4 ) PRINT *, ' Target grid is an ORCA '//c_orca_trg//' grid with north-pole T-point folding!'
      IF( i_orca_trg == 6 ) PRINT *, ' Target grid is an ORCA '//c_orca_trg//' grid with north-pole F-point folding!'
      PRINT *, ''

   END SUBROUTINE TRG_DOMAIN





   SUBROUTINE TERMINATE()

      DEALLOCATE ( data_src, mask_src, data_src_b, lon_src, lat_src, mask_src_b, vt0, vt )

      DEALLOCATE ( mask_trg, data_trg, lat_trg, lon_trg, lon_trg_b )

      IF(l_itrp_3d) THEN
         DEALLOCATE ( data3d_src, depth_src )
         DEALLOCATE ( data3d_tmp, depth_trg, data3d_trg )
         DEALLOCATE ( depth_src_trgt2d )
         IF(TRIM(ctype_z_src) == 'sigma' ) DEALLOCATE ( bathy_src )
         IF(TRIM(ctype_z_trg) == 'sigma' ) DEALLOCATE ( bathy_trg )
      END IF

   END SUBROUTINE TERMINATE






   !! LOCAL SUBROUTINES
   !! ~~~~~~~~~~~~~~~~~

   SUBROUTINE get_src_conf()

      USE mod_manip
      !! Local :
      REAL    :: lon_min_2, lon_max_2
      LOGICAL :: l_loc1, l_loc2

      !! Getting grid on source domain:
      CALL rd_grid(-1, l_reg_src, cf_x_src, cv_lon_src, cv_lat_src, lon_src, lat_src)

      IF( l_itrp_3d ) THEN
         IF( trim(ctype_z_src) == 'sigma' ) THEN
            !! read source bathymetry
            CALL GETVAR_2D(if0,iv0,cf_bathy_src, cv_bathy_src, 0, 0, 0, bathy_src(:,:))
            !! compute 3D depth_src for source variable from bathy and sigma parameters
            CALL depth_from_scoord(ssig_src, bathy_src, ni_src, nj_src, nk_src, depth_src)
         ELSEIF( TRIM(ctype_z_src) == 'z' ) THEN
            !! in z case, the depth vector is copied at each grid-point
            CALL rd_vgrid(nk_src, cf_z_src, cv_z_src, depth_src(1,1,:))
            !WRITE(6,*) ''; WRITE(6,*) 'Source Depths ='; PRINT *, depth_src(1,1,:) ; WRITE(6,*) ''
            DO ji=1,ni_src
               DO jj=1,nj_src
                  depth_src(ji,jj,:) = depth_src(1,1,:)
               ENDDO
            ENDDO
            IF( l_save_drwn) CALL GETVAR_ATTRIBUTES(cf_z_src, cv_z_src,  nb_att_z_src, vatt_info_z_src)
         ELSE
            PRINT*,''; PRINT *, 'Not a valid source vertical coordinate' ; PRINT*,''
         ENDIF

         IF( TRIM(ctype_z_src) == 'z' ) THEN
            PRINT*,''; WRITE(6,*) 'Source has z coordinates and depth vector is:'; PRINT *, depth_src(1,1,:); PRINT*,''
         ELSEIF( TRIM(ctype_z_src) == 'sigma' ) THEN
            PRINT*,''; WRITE(6,*) 'Source has sigma coordinates and depth range is ', MINVAL(depth_src), &
               &                           ' to ', MAXVAL(depth_src) ; PRINT*,''
         ELSE
            PRINT*,''; WRITE(6,*) 'You should not see this' ; STOP
         ENDIF
      END IF

      lon_min_1 = MINVAL(lon_src)
      lon_max_1 = MAXVAL(lon_src)
      PRINT *, ' *** Minimum longitude on source domain before reorg. : ', REAL(lon_min_1,4)
      PRINT *, ' *** Maximum longitude on source domain before reorg. : ', REAL(lon_max_1,4)

      IF( TRIM(cmethod) /= 'no_xy' ) THEN
         IF( l_reg_src ) THEN
            !! Fixing source 1D longitude:
            CALL FIX_LONG_1D(ni_src, lon_src(:,1), nlon_icr_src, i_chg_lon)
            !! Fixing source 1D latitude:
            CALL FIX_LATD_1D(nj_src, lat_src(:,1), nlat_icr_src)
         ELSE
            WHERE ( lon_src < 0. )  lon_src = lon_src + 360.
         END IF
      END IF

      lon_min_2 = MINVAL(lon_src)
      lon_max_2 = MAXVAL(lon_src)
      PRINT *, ' *** Minimum longitude on source domain now: ', lon_min_2
      PRINT *, ' *** Maximum longitude on source domain now: ', lon_max_2

      ! lolo: IMPROVE! This is disgusting:
      l_loc1 = (lon_min_1 <  0.).AND.(lon_min_1 > -175.).AND.(lon_max_1 >  0. ).AND.(lon_max_1 <  175.)
      l_loc2 = (lon_min_2 >= 0.).AND.(lon_min_2 <   2.5).AND.(lon_max_2 >357.5).AND.(lon_max_2 <= 360.)
      IF( (.NOT. l_loc1).AND.(l_loc2) ) THEN
         l_glob_lon_wize = .TRUE.
         PRINT *, 'Looks like global setup (longitude-wise at least...)'
      ELSE
         PRINT *, 'Looks like regional setup (longitude-wise at least...)'
         l_glob_lon_wize = .FALSE.
         WRITE(*,'("  => going to disregard points of target domain with lon < ",f7.2," and lon > ",f7.2)') lon_min_1,lon_max_1
      END IF
      PRINT *, ''

      l_glob_lat_wize = .TRUE.
      IF( MAXVAL(lat_src) < 88. ) l_glob_lat_wize =.FALSE.

      !! Getting land-sea mask on source domain
      mask_src(:,:,:) = 1 ! by default everything is considered sea (helps for smoothing when no LSM)

      IF( l_drown_src ) THEN

         IF( TRIM(cf_lsm_src) == '' ) THEN
            WRITE(6,*) 'ERROR! if you want to "drown" input field (idrown[1]>0) then "cf_lsm_src"'
            WRITE(6,*) '       cannot be an empty string!'
            STOP
         END IF

         IF( (TRIM(cf_lsm_src) == 'missing_value' ).OR. &
            & (TRIM(cf_lsm_src) == 'val+')          .OR. &
            & (TRIM(cf_lsm_src) == 'val-')          .OR. &
            & (TRIM(cf_lsm_src) == 'value')         .OR. &
            & (TRIM(cf_lsm_src) == 'nan')           .OR. &
            & (TRIM(cf_lsm_src) =='NaN') )        THEN

            CALL CREATE_LSM( 'source', cf_lsm_src, cv_lsm_src, mask_src,   cf_fld=cf_src, cv_fld=cv_src )

         ELSE

            !! We are reading source land-sea mask from a file:
            !! -----------------------------------------------
            !! checking dimension:
            CALL DIMS(cf_lsm_src, cv_lsm_src, n1, n2, n3, nrec)

            IF( l3d .AND. ( jplev > 1 ) ) THEN
               IF( n3 == nk_src ) THEN
                  WRITE(6,*) 'Opening 3D land-sea mask on source grid for level', jplev
                  PRINT *, trim(cv_lsm_src)
                  !! if terrain-following, open the 2d mask, not sure interp one single level works
                  IF(TRIM(ctype_z_src) == 'sigma' ) THEN
                     CALL GETMASK_2D(cf_lsm_src, cv_lsm_src, mask_src(:,:,1))
                  ELSEIF(trim(ctype_z_src) == 'z' ) THEN
                     CALL GETMASK_2D(cf_lsm_src, cv_lsm_src, mask_src(:,:,1), jlev=jplev)
                  ELSE
                     STOP 'ERROR: should not be here! #1 (mod_grids.f90)'
                  ENDIF
               ELSE
                  WRITE(6,*) 'PROBLEM! You want to interpolate level', jplev
                  WRITE(6,*) 'but your source land-sea mask is not 3D!'
                  WRITE(6,*) '=> set "idrown" to "0,0" in the namelist'
                  WRITE(6,*) 'If you want to "DROWN" a level other than the surface,'
                  WRITE(6,*) 'please provide a 3D source land-sea mask'
                  STOP
               END IF
            ELSEIF( l_itrp_3d ) THEN
               !! RD: dims reads n3 = -1 on lsm, needs to force n3 to ROMS Nlevels
               !! RD: maybe there is a more elegant way to do that
               IF(trim(ctype_z_src) == 'sigma' ) n3 = ssig_src%Nlevels
               IF( n3 == nk_src ) THEN
                  !! if terrain-following, read 2D mask and copy it on all levels
                  IF(trim(ctype_z_src) == 'sigma' ) THEN
                     WRITE(6,*) 'Opening 2D land-sea mask file on source grid: ', trim(cf_lsm_src)
                     CALL GETMASK_2D(cf_lsm_src, cv_lsm_src, mask_src(:,:,1))
                     DO jz0=2,nk_src
                        mask_src(:,:,jz0) = mask_src(:,:,1)
                     ENDDO
                  ELSEIF(trim(ctype_z_src) == 'z' ) THEN
                     WRITE(6,*) 'Opening 3D land-sea mask file on source grid, ', trim(cv_lsm_src)
                     CALL GETMASK_3D(cf_lsm_src, cv_lsm_src, mask_src)
                  ELSE
                     STOP 'ERROR: should not be here! #2 (mod_grids.f90)'
                  ENDIF
               ELSE
                  WRITE(6,*) 'We need to open the 3D source land-sea mask,'
                  WRITE(6,*) 'but the vertical dimension of it does not match!'
                  STOP
               END IF
            ELSE
               WRITE(6,*) 'Opening land-sea mask file on source grid, ', trim(cv_lsm_src)
               CALL GETMASK_2D(cf_lsm_src, cv_lsm_src, mask_src(:,:,1))
            END IF
         END IF

      END IF  ! IF( l_drown_src ) THEN

      !! Need to modify the mask if lon or lat have been modified :
      IF( nlat_icr_src == -1 ) CALL FLIP_UD(mask_src)
      IF( nlon_icr_src == -1 ) CALL LONG_REORG_3D(i_chg_lon, mask_src)

   END SUBROUTINE get_src_conf



   SUBROUTINE get_trg_conf()

      IF( TRIM(cmethod) /= 'no_xy' ) THEN

         IF( (l_reg_trg).AND.(TRIM(cf_x_trg) == 'spheric') ) THEN

            !! Building target grid:
            READ(cv_lon_trg,*) dx ; READ(cv_lat_trg,*) dy
            cv_lon_trg = 'lon'           ; cv_lat_trg = 'lat'
            WRITE(6,*) '  * dx, dy =', dx, dy
            WRITE(6,*) '  * ni_trg, nj_trg =', ni_trg, nj_trg ;  PRINT*,''
            DO ji = 1, ni_trg
               lon_trg(ji,1) = dx/2.0 + dx*REAL(ji - 1 , 8)
            END DO
            DO jj = 1, nj_trg
               lat_trg(jj,1) = -90 + dy/2.0 + dy*REAL(jj - 1 , 8)
            END DO

            WRITE(6,*) ''; WRITE(6,*) 'Target Longitude array (deg.E):'; PRINT *, lon_trg; WRITE(6,*) ''
            WRITE(6,*) 'Target Latitude array (deg.N):';  PRINT *, lat_trg; WRITE(6,*) ''; WRITE(6,*) ''

         ELSE
            !! Getting target grid from netcdf file:
            CALL rd_grid(ivect, l_reg_trg, cf_x_trg, cv_lon_trg, cv_lat_trg, lon_trg, lat_trg)
            IF( l_reg_trg ) THEN
               WRITE(6,*) ''; WRITE(6,*) 'Target Longitude array (deg.E):'; PRINT *, lon_trg; WRITE(6,*) ''
               WRITE(6,*) 'Target Latitude array (deg.N):';  PRINT *, lat_trg; WRITE(6,*) ''; WRITE(6,*) ''
            END IF
         END IF

         !! Netcdf attributes for longitude and latitude:
         IF( TRIM(cf_x_trg) == 'spheric' ) THEN
            nb_att_lon_trg = 1
            vatt_info_lon_trg(:)%cname = 'null'
            vatt_info_lon_trg(1)%cname = 'units'
            vatt_info_lon_trg(1)%itype = 2 ! char
            vatt_info_lon_trg(1)%val_char = 'degrees_east'
            vatt_info_lon_trg(1)%ilength = LEN('degrees_east')
            nb_att_lat_trg = 1
            vatt_info_lat_trg(:)%cname = 'null'
            vatt_info_lat_trg(1)%cname = 'units'
            vatt_info_lat_trg(1)%itype = 2 ! char
            vatt_info_lat_trg(1)%val_char = 'degrees_west'
            vatt_info_lat_trg(1)%ilength = LEN('degrees_west')
            !!
         ELSE
            !! Geting them from target file:
            CALL GETVAR_ATTRIBUTES(cf_x_trg, cv_lon_trg, nb_att_lon_trg, vatt_info_lon_trg)
            CALL GETVAR_ATTRIBUTES(cf_x_trg, cv_lat_trg, nb_att_lat_trg, vatt_info_lat_trg)
         END IF

         lon_trg_b = lon_trg
         WHERE ( lon_trg < 0. ) lon_trg = lon_trg + 360.

      END IF ! IF( TRIM(cmethod) /= 'no_xy' )


      IF( l_itrp_3d ) THEN
         IF( trim(cf_x_trg)  == 'spheric') THEN
            cf_z_trg = cf_z_src
            cv_z_trg = cv_z_src         !Important
         END IF

         IF( TRIM(ctype_z_trg) == 'sigma' ) THEN

            !! read bathy for target grid
            CALL GETVAR_2D(if0,iv0,cf_bathy_trg, cv_bathy_trg, 0, 0, 0, bathy_trg(:,:))
            !! compute target depth on target grid from bathy_trg and ssig_trg params
            CALL DEPTH_FROM_SCOORD(ssig_trg, bathy_trg, ni_trg, nj_trg, ssig_trg%Nlevels, depth_trg)
            CALL GETVAR_ATTRIBUTES(cf_bathy_trg, cv_bathy_trg,  nb_att_z_trg, vatt_info_z_trg)

         ELSEIF(trim(ctype_z_trg) == 'z' ) THEN

            !! depth vector copied on all grid-points
            CALL rd_vgrid(nk_trg, cf_z_trg, cv_z_trg, depth_trg(1,1,:))
            CALL GETVAR_ATTRIBUTES(cf_z_trg, cv_z_trg, nb_att_z_trg, vatt_info_z_trg)
            DO ji=1,ni_trg
               DO jj=1,nj_trg
                  depth_trg(ji,jj,:) = depth_trg(1,1,:)
               ENDDO
            ENDDO
            WRITE(6,*) ''; WRITE(6,*) 'Target Depths ='; PRINT *, depth_trg(1,1,:) ; WRITE(6,*) ''

         ELSE
            PRINT*,''; WRITE(6,*) 'Not a valid target vertical coordinate' ; STOP
            !!
         ENDIF

         !RD fix this
         !         IF(trim(ctype_z_trg) == 'z' ) THEN
         !            PRINT*,''; WRITE(6,*) 'Target Depths ='; PRINT *, depth_trg(1,1,:) ; PRINT*,''
         !         ELSEIF( trim(ctype_z_trg) == 'sigma' ) THEN
         !            PRINT*,''; WRITE(6,*) 'Target on sigma coordinates' ; PRINT*,''
         !         ENDIF

      END IF

      !RD fix this
      !         CALL rd_vgrid(nk_trg, cf_z_trg, cv_z_trg, depth_trg)
      !         WRITE(6,*) ''; WRITE(6,*) 'Target Depths ='; PRINT *, depth_trg ; WRITE(6,*) ''
      !         CALL GETVAR_ATTRIBUTES(cf_z_trg, cv_z_trg,  nb_att_z_trg, vatt_info_z_trg)

      !!  Getting target mask (mandatory doing 3D interpolation!)
      IF( lmout .OR. l_itrp_3d ) THEN


         IF( (l3d).AND.(jplev > 1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) '****************************************************************'
            WRITE(6,*) 'We do not know target mask yet, since it is at a given depth!'
            WRITE(6,*) '--> next version of SOSIE!'
            WRITE(6,*) 'So we do not mask target file'
            WRITE(6,*) '****************************************************************'
            WRITE(6,*) ''

            mask_trg(:,:,:) = 1

         ELSE


            IF( TRIM(cf_lsm_trg) =='missing_value' ) THEN
               CALL CREATE_LSM( 'target', cf_lsm_trg, cv_lsm_trg, mask_trg,  cf_fld=cf_x_trg, cv_fld=cv_lsm_trg )

            ELSEIF( (TRIM(cf_lsm_trg) == 'val+')  .OR. &
               &     (TRIM(cf_lsm_trg) == 'val-')  .OR. &
               &     (TRIM(cf_lsm_trg) == 'value') .OR. &
               &     (TRIM(cf_lsm_trg) == 'nan')   .OR. &
               &     (TRIM(cf_lsm_trg) =='NaN') )     THEN
               WRITE(6,*) 'ERROR! ("get_trg_conf" of mod_grids.f90): CREATE_LSM does not support method "'//TRIM(cf_lsm_trg)//'" yet for target domain!'
               STOP
            ELSE

               IF( l_itrp_3d ) THEN
                  IF( ((TRIM(cf_lsm_trg) == '').OR.(TRIM(cv_lsm_trg) == '')) ) THEN
                     WRITE(6,*) 'WARNING: no target 3D land-sea mask provided (cf_lsm_trg)!'
                     mask_trg = 1
                  ELSE
                     !! select coord type
                     IF( TRIM(ctype_z_trg) == 'sigma' ) THEN
                        WRITE(6,*) 'Opening 2D land-sea mask file on target grid: ', trim(cf_lsm_trg)
                        !! read 2D mask for target and make it 3D
                        CALL GETMASK_2D(cf_lsm_trg, cv_lsm_trg, mask_trg(:,:,1))
                        DO jz0=2,nk_trg
                           mask_trg(:,:,jz0) = mask_trg(:,:,1)
                        ENDDO
                     ELSEIF( (TRIM(ctype_z_trg) == 'z').AND.(TRIM(cmethod) /= 'no_xy') ) THEN
                        WRITE(6,*) 'Opening 3D land-sea mask file on target grid: ',TRIM(cf_lsm_trg)
                        WRITE(6,*) '             => name mask : ',TRIM(cv_lsm_trg)
                        CALL GETMASK_3D(cf_lsm_trg, cv_lsm_trg, mask_trg(:,:,:))
                        WRITE(6,*) ''
                     ELSE
                        WRITE(6,*) ' We skip reading the 3D target mask because cmethod=="no_xy" ... '
                     ENDIF
                  END IF
               ELSE
                  WRITE(6,*) 'Opening 2D land-sea mask file on target grid: ', TRIM(cf_lsm_trg)
                  CALL GETMASK_2D(cf_lsm_trg, cv_lsm_trg, mask_trg(:,:,1))
               END IF
               WRITE(6,*) ''
            END IF

         END IF

      END IF

   END SUBROUTINE get_trg_conf







   SUBROUTINE rd_grid(iv, lreg, cfgrd, cvx, cvy, rlon, rlat)
      !!
      !! Arguments:
      INTEGER,                   INTENT(in)  :: iv !: iv = -1 means we're handling source grid
      !!                                                /= -1 means target grid
      LOGICAL,                   INTENT(in)  :: lreg
      CHARACTER(len=400),        INTENT(in)  :: cfgrd
      CHARACTER(len=80),         INTENT(in)  :: cvx, cvy
      REAL(8),   DIMENSION(:,:), INTENT(out) :: rlon, rlat

      !! Local
      LOGICAL :: lreg2d, l2dyreg_x, l2dyreg_y
      CHARACTER(len=8) :: cdomain, clreg
      INTEGER :: &
         &     idx, Nx, Ny, &
         &     ii, ij, if1, iv1, &
         &     ilx1, ily1, ilz1,   &
         &     ilx2, ily2, ilz2
      REAL(8),   DIMENSION(:,:), ALLOCATABLE :: zrlon, zrlat

      IF( iv == -1 ) THEN
         cdomain = 'source'
         clreg   = 'l_reg_src'
         idx = 1
      ELSE
         cdomain = 'target'
         clreg   = 'l_reg_trg'
         idx = 2
      END IF

      IF( lreg ) THEN
         IF( (size(rlon,2) /= size(rlat,2)).OR.(size(rlon,2) /= 1) ) THEN
            WRITE(6,*) 'ERROR 1: rd_grid (mod_grids.f90)!' ; STOP
         END IF
         Nx = size(rlon,1) ; Ny = size(rlat,1)
      ELSE
         IF( (size(rlon,1) /= size(rlat,1)).OR.(size(rlon,1) /= size(rlat,1)) ) THEN
            WRITE(6,*) 'ERROR 2: rd_grid (mod_grids.f90)!' ; STOP
         END IF
         Nx = size(rlon,1) ; Ny = size(rlon,2)
      END IF
      rlon = 0. ; rlat = 0.

      !! Checking the dimension of longitude and latitude:
      CALL DIMS(cfgrd, cvx, ilx1, ily1, ilz1, nrec)
      CALL DIMS(cfgrd, cvy, ilx2, ily2, ilz2, nrec)

      WRITE(6,*) ''


      IF(lreg) THEN
         !!         -------------------------------
         !!            R E G U L A R   G R I D :
         !!         -------------------------------
         WRITE(6,*) 'Opening regular grid ', trim(cvx), ' , ', trim(cvy), &
            &                 ' in file ', trim(cfgrd), ' .'

         l2dyreg_x = .FALSE.
         IF( (ilz1 /= -1).or.(ily1 /= -1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) 'Warning! '//trim(cdomain)//' longitude ', trim(cvx), ' is not 1D!'
            WRITE(6,*) 'In the file ', trim(cfgrd)
            WRITE(6,*) ' => maybe you should specify '//trim(clreg)//' = F in the namelist?'
            WRITE(6,*) ' => anyway we assume you know what you are doing!'
            l2dyreg_x = .TRUE.
         END IF

         l2dyreg_y = .FALSE.
         IF( (ilz2 /= -1).or.(ily2 /= -1) ) THEN
            WRITE(6,*) ''
            WRITE(6,*) 'Warning! '//trim(cdomain)//' latitude ', trim(cvy), ' is not 1D!'
            WRITE(6,*) 'In the file ', trim(cfgrd)
            WRITE(6,*) ' => maybe you should specify '//trim(clreg)//' = F in the namelist?'
            WRITE(6,*) ' => anyway we assume you know what you are doing!'
            l2dyreg_y = .TRUE.
         END IF

         IF( l2dyreg_x .AND. (.NOT. l2dyreg_y) ) THEN
            WRITE(6,*) 'ERROR! '//trim(cdomain)//' longitude and latidude do not agree in shape (1D vs 2D)!' ; STOP
         END IF

         IF( l2dyreg_x .AND. l2dyreg_y ) THEN
            WRITE(6,*) ' '
            WRITE(6,*) ' =================================================================================================='
            WRITE(6,*) ' *** Assuming that '//trim(cdomain)//' grid is regular even though longitude and latidude are 2D!'
            WRITE(6,*) '                     (because you set '//trim(clreg)//'=.TRUE. in the namelist)'
            WRITE(6,*) ' =================================================================================================='
            l_2d_grid_yet_regular(idx) = .TRUE.
            WRITE(6,*) ' '
         END IF


         IF( l_2d_grid_yet_regular(idx) ) THEN
            !! Getting regular 2D grid :
            ALLOCATE ( zrlon(Nx,Ny), zrlat(Nx,Ny) )
            CALL GETVAR_2D(if1, iv1, cfgrd, cvx, 0, 0, 0, zrlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D(if1, iv1, cfgrd, cvy, 0, 0, 0, zrlat)
            !!
            !! Checking for Regular 2D longitude and latitude
            lreg2d = .TRUE.
            !! Checking if longitude array changes with latitude
            DO ij = 2, Ny
               IF( ABS(SUM(zrlon(:,ij)-zrlon(:,1))) > 1.e-6 ) THEN
                  lreg2d = .FALSE.
                  EXIT
               END IF
            END DO
            !! Checking if latitude array changes with longitude
            DO ii = 2, Nx
               IF( ABS(SUM(zrlat(ii,:)-zrlat(1,:))) > 1.e-6 ) THEN
                  lreg2d = .FALSE.
                  EXIT
               END IF
            END DO
            IF( lreg2d ) THEN
               WRITE(6,*) ' *** OK! You were right, '//trim(cdomain)//' longitude and latitude are 2D but the grid is regular!'
               WRITE(6,*) ''
            ELSE
               WRITE(6,*) 'ERROR! We checked, and 2D '//trim(cdomain)//' longitude/latitude do not qualify for a regular grid!!!!'
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
            IF( (ilx1 /= Nx).or.(ilx2 /= Ny) ) THEN
               WRITE(6,*) 'ERROR! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',  &
                  &   trim(cvy), &
                  &   ' do not agree in dimension with configuration dimension!'
               WRITE(6,*) 'In the file ', trim(cfgrd) ; STOP
            END IF
            CALL GETVAR_1D(cfgrd, cvx, rlon(:,1))
            CALL GETVAR_1D(cfgrd, cvy, rlat(:,1))
         END IF


      ELSEIF(.NOT. lreg) THEN
         !!         ----------------------------------
         !!            I R R E G U L A R   G R I D :
         !!         ----------------------------------
         WRITE(6,*) 'Opening irregular grid ', TRIM(cvx),', ',TRIM(cvy),' in file:'
         WRITE(6,*) '   ', TRIM(cfgrd)
         !!WRITE(6,*) ' => ilx1, ilx2, ily1, ily2, ilz1, ilz2, Nx, Ny =', ilx1, ilx2, ily1, ily2, ilz1, ilz2, Nx, Ny
         WRITE(6,*) ''

         IF(ilx1 /= ilx2) THEN

            IF( (ily1 == ily2).and.(ilz1 == ilz2).and.(ily1 == -1).and.(ilz1 == -1) ) THEN
               WRITE(6,*) 'ERROR! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',   &
                  &       trim(cvy), ' seem to be regular (1D), check namelist!'
               WRITE(6,*) 'In the file ', trim(cfgrd)
               WRITE(6,*) '       => so maybe try to set '//trim(clreg)//' to .TRUE. in the namelist!'
               STOP
            ELSE
               WRITE(6,*) 'ERROR! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',  &
                  &   trim(cvy), ' do not agree in dimension!'
               WRITE(6,*) 'In the file ', trim(cfgrd) ; STOP
            END IF
         END IF

         IF( (ily1 /= ily2).or.(ilz1 /= ilz2) ) THEN
            WRITE(6,*) 'ERROR! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',    &
               &   trim(cvy), &
               &   ' do not agree in dimension!'
            WRITE(6,*) 'In the file ', trim(cfgrd); STOP
         END IF

         IF( (ilx1 /= Nx).OR.(ily1 /= Ny) ) THEN
            WRITE(6,*) 'ERROR! '//trim(cdomain)//' longitude ', trim(cvx), ' and latitude ',   &
               &   trim(cvy),   &
               &   ' do not agree in dimension with configuration dimension!'
            WRITE(6,*) 'In the file ', trim(cfgrd); STOP
         END IF


         !! Getting source longitude array at level=1 (surface) and time =1 :
         !! -----------------------------------------------------------------
         !! If 3d dimension, we chose first level
         IF( ilz1 /= -1 ) THEN
            CALL GETVAR_2D(if1, iv1, cfgrd, cvx, 0, 1, 0, rlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D(if1, iv1, cfgrd, cvy, 0, 1, 0, rlat)
         ELSE
            CALL GETVAR_2D(if1, iv1, cfgrd, cvx, 0, 0, 0, rlon) ; if1 = 0 ; iv1 = 0
            CALL GETVAR_2D(if1, iv1, cfgrd, cvy, 0, 0, 0, rlat)
         END IF

      END IF


   END SUBROUTINE rd_grid



   SUBROUTINE rd_vgrid(nk, cfgrd, cvz, vdepth)

      INTEGER,                   INTENT(in) :: nk
      CHARACTER(len=400),        INTENT(in) :: cfgrd
      CHARACTER(len=80),         INTENT(in) :: cvz
      REAL(4), DIMENSION(nk) ,  INTENT(out) :: vdepth

      CALL GETVAR_1D(cfgrd, cvz, vdepth)

   END SUBROUTINE rd_vgrid





   SUBROUTINE know_dim_src()

      INTEGER :: jk0, jrec

      nk_src = 1

      CALL DIMS(cf_src, cv_src, ni_src, nj_src, jk0, jrec)

      IF( nj_src == -1 ) THEN
         WRITE(6,*) 'ERROR (know_dim_src)! variable ',TRIM(cv_src),' should be at least 2D!!!'
         STOP
      END IF


      IF( jplev == -1 ) THEN
         !! This is the overide case! Means 2D+T !
         !! Case when we want to read a 2D+T field but someone screwed up and the record
         !! dimension is not of 'UNLIMITED' type inside the netcdf file...
         !! So it just  overides good sence and force sosie to understand that
         !! your field to interpolate is 2D with a time record
         !! (usually the case if the time record dimension in your
         !! source file is not declared as UNLIMITED => bad! :(
         WRITE(6,*) ''
         Ntr = jk0
         WRITE(6,*) 'WARNING: know_dim_src of mod_grids.f90 !!!'
         WRITE(6,*) '   => we force source field "'//TRIM(cv_src)//'" to be 2D + time !!!'
         WRITE(6,*) '   => because you specified "jplev = -1" in the namelist!'
         WRITE(6,*) '   => the time-record dimension is therefore:', Ntr
         WRITE(6,*) ''
         Ntr0 = Ntr

      ELSE

         IF( jrec == -1) ltime=.FALSE.  ! no time records

         !! 3D variable
         IF( jk0 > 0 ) THEN

            l3d = .TRUE.
            nk_src = jk0 ; WRITE(6,*) 'Number of level into source file =', nk_src

            IF(jplev /= 0) THEN
               IF( (jplev <= 0).OR.(jplev > nk_src) ) THEN
                  WRITE(6,*) 'Level jplev is wrong! -> number of level =', nk_src
                  WRITE(6,*) 'In the namelist, jplev =', jplev ;  STOP
               END IF
               WRITE(6,*) 'We are going to interpolate level', jplev, ' of source file'
            ELSE
               WRITE(6,*) ''
               WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
               WRITE(6,*) 'WARNING: We are going to perform a 3D interpolation !!!'
               WRITE(6,*) '      => if this is not what you want, please specify '
               WRITE(6,*) '         which level (jplev) should be interpolated.'
               WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
               WRITE(6,*) ''
               l_itrp_3d = .TRUE.
            END IF

            IF( jrec > 0 ) THEN
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
         IF( jk0 == -1 ) THEN

            l3d = .FALSE.

            IF( jrec > 0 ) THEN
               !! 2D variable with time
               WRITE(6,*) 'Variable is 2D with time axis'
               Ntr0 = jrec
            ELSE
               !! 2D variable without time
               WRITE(6,*) 'Variable is 2D without time axis'
               Ntr0 = 1
            END IF

         END IF

         IF(.NOT. ltime) Ntr0 = 1

         Ntr = Ntr0
         j_start = 1
         j_stop = Ntr0

      END IF

      j_start = 1 ; j_stop = Ntr


      IF( (jt1 > 0).AND.(jt2 > 0) ) THEN
         Ntr = jt2 - jt1 + 1 ;  ;  j_start = jt1 ;  j_stop = jt2
         !! jrec is the time dimension of the source file, Ntr is the length requested by the user :
         IF( (jplev>=0).AND.( (Ntr>jrec).OR.(jt1 < 1).OR.(jt1 > jt2).OR.(jt2 > jrec) ) ) THEN
            WRITE(6,*) ''; WRITE(6,*) 'Check jt1 and jt2 in the namelist:'
            WRITE(6, '("=> the time dimension of ",a," is ", i5)') TRIM(cv_src), jrec
            WRITE(6, '("   and you specified jt1, jt2 =",i5," ,",i5)') jt1, jt2
            STOP
         END IF
      END IF

      IF( l_itrp_3d .AND. trim(ctype_z_src) == 'sigma' ) nk_src = ssig_src%Nlevels ! ugly but should work

   END SUBROUTINE know_dim_src






   SUBROUTINE know_dim_trg
      !!
      nk_trg = 1
      !!
      IF( TRIM(cmethod) /= 'no_xy' ) THEN
         !! If 3D interpolation and building spherical grid, we use levels from source grid:
         IF( l_itrp_3d .AND. (TRIM(cf_x_trg)  == 'spheric') ) THEN
            cf_z_trg = cf_z_src ;  cv_z_trg = cv_z_src
         END IF
         !!
         IF(l_reg_trg) THEN
            IF( TRIM(cf_x_trg) == 'spheric') THEN
               WRITE(6,*) ''; WRITE(6,*) 'Building regular spherical target grid!'
               READ(cv_lon_trg,*) dx
               READ(cv_lat_trg,*) dy
               ni_trg = INT(360./dx)
               nj_trg = INT(180./dy)
               GOTO 100
            ELSE
               CALL DIMS(cf_x_trg, cv_lon_trg, ni_trg, n1, n2, nrec)
               CALL DIMS(cf_x_trg, cv_lat_trg, nj_trg, n1, n2, nrec)
            END IF
         ELSE
            CALL DIMS(cf_x_trg, cv_lon_trg, ni_trg, nj_trg, n1, nrec)
         END IF
         !!
      END IF ! IF( TRIM(cmethod) /= 'no_xy' )

100   CONTINUE

      IF( l_itrp_3d ) THEN
         !!
         WRITE(6,*) ''
         IF( TRIM(cf_x_trg)  == 'spheric' ) THEN
            WRITE(6,*) 'Since we are building our own spherical target grid,'
            WRITE(6,*) 'we are going to use source levels as target levels!'
         END IF
         !!
         WRITE(6,*) ''
         WRITE(6,*) ' => we read target levels in the following file:'
         PRINT *, TRIM(cf_z_trg); WRITE(6,*) ''
         IF( trim(ctype_z_trg) == 'sigma' ) THEN
            nk_trg = ssig_trg%Nlevels
         ELSE
            CALL DIMS(cf_z_trg, cv_z_trg, nk_trg, n1, n2, nrec)
         ENDIF
         WRITE(6,*) 'nk_trg = ', nk_trg ; WRITE(6,*) ''
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(6,*) 'Target grid dimension is', ni_trg,'x',nj_trg,'x',nk_trg
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         !!
      ELSE
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         WRITE(6,*) 'Target grid dimension is', ni_trg,'x', nj_trg
         WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      END IF
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
      !!
   END SUBROUTINE know_dim_trg
   !!
   !!
   !!
   !!
   !!
   SUBROUTINE FIX_LONG_1D(nx, vlon, n_icr, i_chg_x)
      !!
      !! Makes 1D longitude understandable
      !!
      INTEGER,                  INTENT(in)    :: nx
      REAL(8), DIMENSION(nx),  INTENT(inout) :: vlon
      INTEGER,                  INTENT(out)   :: n_icr
      INTEGER,                  INTENT(out)   :: i_chg_x
      !!
      REAL(8), DIMENSION(:), ALLOCATABLE    :: x1d_b
      INTEGER :: jx
      !!---------------------------------------------------------------
      IF( (vlon(1) == -180.0).AND.(vlon(nx) == 180.0) ) THEN
         !! to avoid having twice the point "lon. 180" in the final longitude array
         vlon(1)  = -179.99
         vlon(nx) =  179.99
         IF( ewper_src /= 0 ) THEN
            ewper_src = 0
            WRITE(6,*) 'mod_grids.f90: FIX_LONG_1D => "ewper_src" forced to 0 !'; WRITE(6,*) ''
         END IF
      END IF

      !! We want positive and increasing longitude !
      WHERE ( vlon < 0. )  vlon = vlon + 360.
      !!
      i_chg_x = 0 ; n_icr = 1
      DO jx = 2, nx
         IF( vlon(jx) < vlon(jx-1) ) THEN
            IF( i_chg_x /= 0) THEN  ! another changing sign point has already been found!
               WRITE(6,*) 'Your longitude array is a mess!'
               WRITE(6,*) ' --> Fix it (positive and increasing)!'; STOP
            ELSE
               i_chg_x = jx
               n_icr = -1
            END IF
         END IF
      END DO
      !!
      IF( n_icr == -1 ) THEN
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


   SUBROUTINE FIX_LATD_1D(ny, vlat, n_icr)
      !!
      !! Makes 1D latitude increasing with jy
      !!
      INTEGER,                  INTENT(in)    :: ny
      REAL(8), DIMENSION(ny), INTENT(inout) :: vlat
      INTEGER,                  INTENT(out)   :: n_icr
      !!
      REAL(8), DIMENSION(:), ALLOCATABLE :: y1d_b
      INTEGER :: jy
      !!
      n_icr = 1
      !!
      IF( l_reg_src .AND. ( vlat(1) > vlat(ny) ) ) THEN
         !!
         n_icr = -1 ; ALLOCATE( y1d_b(ny) )
         !!
         WRITE(6,*) ''; WRITE(6,*) 'Latitude does not seem to increase with j'
         WRITE(6,*) '--> We reverse 1D source latitude array!'
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
      IF( PRESENT(cname_long) ) cnlon = TRIM(cname_long)




      IS_ORCA_NORTH_FOLD%ifld_nord =  0
      IS_ORCA_NORTH_FOLD%cgrd_type = 'X'

      nx = SIZE(Xtest,1)
      ny = SIZE(Xtest,2)

      IF( ny > 3 ) THEN ! (case if called with a 1D array, ignoring...)


         IF( SUM( Xtest(2:nx/2,ny) - Xtest(nx:nx-nx/2+2:-1,ny-2) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 4 ! T-pivot, grid_T
            IS_ORCA_NORTH_FOLD%cgrd_type = 'T'
         END IF
         !---
         IF( SUM( Xtest(2:nx/2,ny) - Xtest(nx-1:nx-nx/2+1:-1,ny-2) ) == 0. ) THEN
            IF(TRIM(cnlon)=='glamu') THEN
               IS_ORCA_NORTH_FOLD%ifld_nord = 4 ! T-pivot, grid_T
               IS_ORCA_NORTH_FOLD%cgrd_type = 'U'
            END IF
            !! LOLO: PROBLEM == 6, V !!!
         END IF
         !---
         IF( SUM( Xtest(2:nx/2,ny) - Xtest(nx:nx-nx/2+2:-1,ny-3) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 4 ! T-pivot, grid_V
            IS_ORCA_NORTH_FOLD%cgrd_type = 'V'
         END IF


         IF( SUM( Xtest(2:nx/2,ny) - Xtest(nx-1:nx-nx/2+1:-1,ny-1) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 6 ! F-pivot, grid_T
            IS_ORCA_NORTH_FOLD%cgrd_type = 'T'
         END IF
         IF( SUM( Xtest(2:nx/2,ny) - Xtest(nx-2:nx-nx/2:-1,ny-1) ) == 0. ) THEN
            IS_ORCA_NORTH_FOLD%ifld_nord = 6 ! F-pivot, grid_U
            IS_ORCA_NORTH_FOLD%cgrd_type = 'U'
         END IF
         !---
         IF( SUM( Xtest(2:nx/2,ny) - Xtest(nx-1:nx-nx/2+1:-1,ny-2) ) == 0. ) THEN
            IF(TRIM(cnlon)=='glamv') THEN
               IS_ORCA_NORTH_FOLD%ifld_nord = 6 ! F-pivot, grid_V
               IS_ORCA_NORTH_FOLD%cgrd_type = 'V'
            END IF
            !! LOLO: PROBLEM == 4, U !!!
         END IF
         !---

      END IF

   END FUNCTION IS_ORCA_NORTH_FOLD


   SUBROUTINE CREATE_LSM_3D( cinfo, cmthd, cnumv, mask,  cf_fld, cv_fld, xfield )
      CHARACTER(len=6),  INTENT(in) :: cinfo  ! 'target' or 'source'
      CHARACTER(len=*),  INTENT(in) :: cmthd, cnumv  ! 'missing_value','value','val+','val-','nan','Nan'
      INTEGER(1), DIMENSION(:,:,:), INTENT(inout) :: mask
      !!       !! if xfield is present then don't read cv out of cf_fld
      !!         => so it's either cf_fld=xxx, cv_fld=xxx OR xfield=XXXX, not both !!!
      CHARACTER(len=*), OPTIONAL,            INTENT(in) :: cf_fld, cv_fld
      REAL(wpl), OPTIONAL, DIMENSION(:,:,:), INTENT(in) :: xfield

      LOGICAL :: l_use_field_array
      INTEGER :: ni, nj, nk
      REAL(wpl), DIMENSION(:,:,:), ALLOCATABLE :: z3d_tmp
      REAL    :: rval_thrshld

      IF( .NOT. ((cinfo == 'source').OR.(cinfo == 'target')) ) THEN
         PRINT *, 'ERROR (CREATE_LSM_3D of mod_grids.f90) : unknown "cinfo" => ', cinfo ; STOP
      END IF
      IF(PRESENT(cf_fld)) THEN
         IF(.NOT. PRESENT(cv_fld)) THEN
            PRINT *, 'ERROR (CREATE_LSM_3D of mod_grids.f90) : if "cf_fld" is specified so must be "cv_fld" !!!'; STOP
         END IF
      END IF
      IF( (PRESENT(cf_fld)).AND.(PRESENT(xfield)) ) THEN
         PRINT *, 'ERROR (CREATE_LSM_3D of mod_grids.f90) : if "cf_fld" and "cv_fld" are specified then "xfield" should not !!!'; STOP
      END IF

      ni = SIZE(mask,1)
      nj = SIZE(mask,2)
      nk = SIZE(mask,3)

      l_use_field_array = .FALSE.
      IF(PRESENT(xfield)) THEN
         l_use_field_array = .TRUE.
         IF( (SIZE(xfield,1)/=ni).OR.(SIZE(xfield,2)/=nj).OR.(SIZE(xfield,1)/=nk) ) THEN
            PRINT *, 'ERROR (CREATE_LSM_3D of mod_grids.f90) : if "mask" and "xfield" do not agree in shape!!!'; STOP
         END IF
      END IF

      ALLOCATE ( z3d_tmp(ni,nj,nk) )

      IF( l_use_field_array ) THEN
         z3d_tmp(:,:,:) = xfield(:,:,:)
      ELSE
         !! Will read data field (at time record # jt0 if time exists) :
         IF( ltime  ) jt0 = 1
         CALL GETVAR_3D(if0, iv0, cf_fld, cv_fld, Ntr,      jt0, z3d_tmp)
         DO jk =  1, nk
            DO jj =  1, nj
               DO ji =  1, ni
                  IF( ISNAN(z3d_tmp(ji,jj,jk)) ) z3d_tmp(ji,jj,jk) = rmissval
               END DO
            END DO
         END DO

      END IF

      mask(:,:,:) = 1

      IF( TRIM(cmthd) == 'missing_value' ) THEN
         WRITE(6,*) 'Opening land-sea mask "'//TRIM(cinfo)//'" from missing_value of source field "'//TRIM(cv_fld)//'"!'
         CALL CHECK_4_MISS(cf_fld, cv_fld, lmval, rmv, ca_missval)
         !PRINT *, 'LOLO: rmv =', rmv, ISNAN(rmv)

         IF( .NOT. lmval ) THEN
            PRINT *, 'ERROR (CREATE_LSM_3D of mod_grids.f90) : '//TRIM(cv_fld)//' has no missing value attribute!'
            PRINT *, '      (in '//TRIM(cf_fld)//')'
            STOP
         END IF

         IF( ISNAN(rmv) ) THEN
            WHERE ( z3d_tmp < -9990. ) mask = 0  ! NaN have been replaced with rmissval earlier !
         ELSE
            WHERE ( z3d_tmp == rmv   ) mask = 0
         END IF
         !!CALL DUMP_FIELD(REAL(mask,4), 'mask_trg.tmp', 'lsm')

      ELSEIF( (TRIM(cmthd) == 'val+').OR.(TRIM(cmthd) == 'val-').OR.(TRIM(cmthd) == 'value') ) THEN
         READ(cnumv,*) rval_thrshld
         IF(TRIM(cmthd) == 'val+')  WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from values >=', rval_thrshld
         IF(TRIM(cmthd) == 'val-')  WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from values <=', rval_thrshld
         IF(TRIM(cmthd) == 'value') WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from values ==', rval_thrshld
         IF(TRIM(cmthd) == 'val+') THEN
            WHERE ( z3d_tmp >= rval_thrshld ) mask = 0
         END IF
         IF(TRIM(cmthd) == 'val-') THEN
            WHERE ( z3d_tmp <= rval_thrshld ) mask = 0
         END IF
         IF(TRIM(cmthd) == 'value') THEN
            WHERE ( z3d_tmp == rval_thrshld ) mask = 0
         END IF

      ELSEIF((TRIM(cmthd)=='nan').OR.(TRIM(cmthd)=='NaN')) THEN
         WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from NaN values in '//TRIM(cv_src)//' at jt = ', jt0
         ! =>  NaN values have been flagged earlier and replaced by rmissval
         WHERE ( z3d_tmp < -9998. ) mask = 0
      ELSE
         WRITE(6,*) 'ERROR! (CREATE_LSM_3D of mod_grids.f90): Unknown value for "cmthd": '//TRIM(cmthd)//' !' ; STOP
      END IF
      DEALLOCATE ( z3d_tmp )
   END SUBROUTINE CREATE_LSM_3D





   SUBROUTINE CREATE_LSM_2D( cinfo, cmthd, cnumv, mask,  cf_fld, cv_fld, xfield )
      CHARACTER(len=6),  INTENT(in) :: cinfo  ! 'target' or 'source'
      CHARACTER(len=*),  INTENT(in) :: cmthd, cnumv  ! 'missing_value','value','val+','val-','nan','Nan'
      INTEGER(1), DIMENSION(:,:), INTENT(inout) :: mask
      !!       !! if xfield is present then don't read cv out of cf_fld
      !!         => so it's either cf_fld=xxx, cv_fld=xxx OR xfield=XXXX, not both !!!
      CHARACTER(len=*), OPTIONAL,            INTENT(in) :: cf_fld, cv_fld
      REAL(wpl), OPTIONAL, DIMENSION(:,:), INTENT(in) :: xfield
      !!
      LOGICAL :: l_use_field_array
      INTEGER :: ni, nj
      REAL(wpl), DIMENSION(:,:), ALLOCATABLE :: z2d_tmp
      REAL    :: rval_thrshld

      IF( .NOT. ((cinfo == 'source').OR.(cinfo == 'target')) ) THEN
         PRINT *, 'ERROR (CREATE_LSM_2D of mod_grids.f90) : unknown "cinfo" => ', cinfo ; STOP
      END IF
      IF(PRESENT(cf_fld)) THEN
         IF(.NOT. PRESENT(cv_fld)) THEN
            PRINT *, 'ERROR (CREATE_LSM_2D of mod_grids.f90) : if "cf_fld" is specified so must be "cv_fld" !!!'; STOP
         END IF
      END IF
      IF( (PRESENT(cf_fld)).AND.(PRESENT(xfield)) ) THEN
         PRINT *, 'ERROR (CREATE_LSM_2D of mod_grids.f90) : if "cf_fld" and "cv_fld" are specified then "xfield" should not !!!'; STOP
      END IF

      ni = SIZE(mask,1)
      nj = SIZE(mask,2)

      l_use_field_array = .FALSE.
      IF(PRESENT(xfield)) THEN
         l_use_field_array = .TRUE.
         IF( (SIZE(xfield,1)/=ni).OR.(SIZE(xfield,2)/=nj) ) THEN
            PRINT *, 'ERROR (CREATE_LSM_2D of mod_grids.f90) : if "mask" and "xfield" do not agree in shape!!!'; STOP
         END IF
      END IF

      ALLOCATE ( z2d_tmp(ni,nj) )

      IF( l_use_field_array ) THEN
         z2d_tmp(:,:) = xfield(:,:)
      ELSE
         !! Will read data field (at time record # jt0 if time exists) :
         IF( ltime  ) jt0 = 1
         jz0 = 0
         IF( l3d ) jz0 = jplev
         CALL GETVAR_2D(if0, iv0, cf_fld, cv_fld, Ntr, jz0, jt0, z2d_tmp)
         DO jj =  1, nj
            DO ji =  1, ni
               IF( ISNAN(z2d_tmp(ji,jj)) ) z2d_tmp(ji,jj) = rmissval
            END DO
         END DO
      END IF


      mask(:,:) = 1

      IF( TRIM(cmthd) == 'missing_value' ) THEN
         WRITE(6,*) 'Opening land-sea mask "'//TRIM(cinfo)//'" from missing_value of source field "'//TRIM(cv_fld)//'"!'
         CALL CHECK_4_MISS(cf_fld, cv_fld, lmval, rmv, ca_missval)
         IF( .NOT. lmval ) THEN
            PRINT *, 'ERROR (CREATE_LSM_2D of mod_grids.f90) : '//TRIM(cv_fld)//' has no missing value attribute!'
            PRINT *, '      (in '//TRIM(cf_fld)//')'
            STOP
         END IF

         IF( ISNAN(rmv) ) THEN
            WHERE ( z2d_tmp < -9990. ) mask = 0  ! NaN have been replaced with rmissval earlier !
         ELSE
            WHERE ( z2d_tmp == rmv )   mask = 0
         END IF

      ELSEIF( (TRIM(cmthd) == 'val+').OR.(TRIM(cmthd) == 'val-').OR.(TRIM(cmthd) == 'value') ) THEN
         READ(cnumv,*) rval_thrshld
         IF(TRIM(cmthd) == 'val+')  WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from values >=', rval_thrshld
         IF(TRIM(cmthd) == 'val-')  WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from values <=', rval_thrshld
         IF(TRIM(cmthd) == 'value') WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from values ==', rval_thrshld
         IF(TRIM(cmthd) == 'val+') THEN
            WHERE ( z2d_tmp >= rval_thrshld ) mask = 0
         END IF
         IF(TRIM(cmthd) == 'val-') THEN
            WHERE ( z2d_tmp <= rval_thrshld ) mask = 0
         END IF
         IF(TRIM(cmthd) == 'value') THEN
            WHERE ( z2d_tmp == rval_thrshld ) mask = 0
         END IF

      ELSEIF((TRIM(cmthd)=='nan').OR.(TRIM(cmthd)=='NaN')) THEN
         WRITE(6,*) ' Land-sea mask "'//TRIM(cinfo)//'" is defined from NaN values in '//TRIM(cv_src)//' at jt = ', jt0
         ! =>  NaN values have been flagged earlier and replaced by rmissval
         WHERE ( z2d_tmp < -9998. ) mask = 0

      ELSE
         WRITE(6,*) 'ERROR! (CREATE_LSM_2D of mod_grids.f90): Unknown value for "cmthd": '//TRIM(cmthd)//' !' ; STOP
      END IF
      DEALLOCATE ( z2d_tmp )

   END SUBROUTINE CREATE_LSM_2D



END MODULE MOD_GRIDS
