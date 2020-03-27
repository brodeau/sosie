PROGRAM SOSIE

   !!----------------------------------------------------------------------
   !!
   !!     SOSIE is Only a Surface Interpolation Environment
   !!     =================================================
   !!
   !!     Version 3.0, Janvier 2019
   !!
   !!
   !!    PURPOSE:
   !!    --------
   !!    * Reads a geophysical field on a given source grid in a netcdf file
   !!
   !!    * Remaps it on a given target grid using either the Akima method or the bilinear
   !!      methods in case the source grid is distorded (horizontally)
   !!
   !!    * Writes the remaped field in a netcf file
   !!
   !!
   !!    METHODS:
   !!    --------
   !!      * Horizontal interpolation
   !!          - Akima (1970) method: "A New Method of Interpolation and Smooth Surface Fitting Based
   !!                               on local Procedures, J.of Applied Comput. Math., 17, 589-602."
   !!          - Standard bilinear
   !!      * Verical interpolation (for 3D fields)
   !!
   !!
   !!
   !!     AUTHORS:
   !!     -------
   !!     Laurent Brodeau, 2019
   !!     Jean-Michel Brankart, Jean-Marc Molines, Pierre Mathiot
   !!
   !!            Contact: https://github.com/brodeau/sosie
   !!
   !!--------------------------------------------------------------------------

   USE io_ezcdf     !* routines for netcdf input/output
   USE mod_conf     !* important parameters, namelist, arrays, etc...
   USE mod_init     !* important parameters, namelist, arrays, etc...
   USE mod_grids    !* grided-domain module
   USE mod_interp


   IMPLICIT NONE

   CHARACTER(len=800) :: cextinf

   INTEGER    :: &
      &   ji, jj, &
      &   jt, jte, jct, jcz,  &
      &   idf_o, idf_i, idv_o, idv_i, &  !: netcdf ID for source/target file
      &   idf_id, idv_id  !: drowned source field...

   REAL(4) :: rfct_miss=1.


   !OPEN(UNIT=6, FORM='FORMATTED', RECL=512)  ! problem with Gfortan 4.8...

   WRITE(6,*)''
   WRITE(6,*)'=========================================================='
   WRITE(6,*)'            S  O  S  I  E    version ', trim(csosie_version)
   WRITE(6,*)'=========================================================='
   WRITE(6,*)''


   !! Fecthing command line arguments if any:
   CALL GET_ARGUMENTS()    ! MODULE mod_init

   !! Reading namelist:
   CALL READ_NMLST(1)      ! MODULE mod_init

   !! Setting defaults for some logical flags :
   l3d = .FALSE. ; l_itrp_3d = .FALSE. ; ltime = .TRUE.

   !! Getting source domain and metrics:
   CALL SRC_DOMAIN()       ! MODULE mod_grids
   mask_src_b = mask_src     ! making a backup of the land-sea mask on the source grid

   !! Getting target domain and metrics:
   CALL TRG_DOMAIN()       ! MODULE mod_grids

   !! Just to display some info of the current setting:
   CALL REMINDER()         ! MODULE mod_init


   IF ( Ntr == 0 ) THEN
      PRINT *, 'PROBLEM: something is wrong => Ntr = 0 !!!'; STOP
   END IF

   IF ( .NOT. lmout) rfct_miss = 0.

   cextinf =                'Horizontal grid read in '//TRIM(cf_x_trg)
   IF (l_itrp_3d) cextinf = TRIM(cextinf)//' / Vertical grid read in '//TRIM(cf_z_trg)
   cextinf = TRIM(cextinf)//' / Source field read in '//TRIM(cf_src)
   cextinf = TRIM(cextinf)//' / Interpolation method: '//TRIM(cmethod)


   !! Netcdf Atributes of the interpolated field as in input file:
   CALL GETVAR_ATTRIBUTES(cf_src, cv_src,  nb_att_F, vatt_info_F) ; !getting all attributes for treated field !lolo
   !! Overwritting some attributes given in the namelist if /= '':
   IF ( TRIM(cu_out)  /= '' ) CALL FORCE_ATTR('units',      cu_out, vatt_info_F)
   IF ( TRIM(cln_out) /= '' ) CALL FORCE_ATTR('long_name', cln_out, vatt_info_F)







   !!                -------------------------------
   !!                  M A I N   T I M E   L O O P
   !!                -------------------------------

   l_first_call_interp_routine = .TRUE.

   jct = 0 ;  jcz = 0
   IF ( ltime ) jct = 1
   IF (  l3d  ) jcz = 1

   DO jt = 1, Ntr

      jte = j_start + jt - 1     ! actual time record to extract

      WRITE(*,'(" *** treating record #",i5)') jte

      IF ( .not. l_itrp_3d ) THEN

         !! ================
         !! 2D INTERPOLATION
         !! ================

         data_src = 0.

         !! Read data 2D field at time jte :
         CALL GETVAR_2D(idf_i, idv_i, cf_src, cv_src, Ntr, jplev*jcz, jte*jct, data_src, jt1=j_start, jt2=j_stop)

         !IF ((TRIM(cf_lsm_src)=='nan').OR.(TRIM(cf_lsm_src)=='NaN')) THEN
         !! Replacing NaN with 0. to avoid some fuck-up later...
         DO jj =  1, nj_src
            DO ji =  1, ni_src
               IF ( ISNAN(data_src(ji,jj)) ) data_src(ji,jj) = -9999.
            END DO
         END DO
         !END IF


         CALL INTERP_2D(jt)

         !! => data_trg for current time step is ready to be written in netcdf file

         !! Print current record into netcdf file
         CALL P2D_T(idf_o, idv_o, Ntr, jt,    &
            &      lon_trg_b, lat_trg, vt, data_trg,    &
            &      cf_out, cv_lon_trg, cv_lat_trg, cv_t_out,    &
            &      cv_out, rfct_miss*REAL(rmiss_val,4), &
            &      attr_lon=vatt_info_lon_trg, attr_lat=vatt_info_lat_trg, attr_time=vatt_info_t, attr_F=vatt_info_F, &
            &      cextrainfo=cextinf)

         IF ( l_save_drwn ) THEN
            CALL P2D_T(idf_id, idv_id, Ntr, jt,    &
               &       lon_src, lat_src, vt, data_src_drowned(:,:,1),    &
               &       TRIM(cf_src)//'.drwn', cv_lon_src, cv_lat_src, cv_t_src,    &
               &       cv_src, 0., &
               &       attr_lon=vatt_info_lon_src, attr_lat=vatt_info_lat_src, attr_time=vatt_info_t, attr_F=vatt_info_F )
         END IF


      ELSE



         !! ================
         !! 3D INTERPOLATION
         !! ================

         data3d_src = 0.

         !! Read data 3D field at time jt :
         CALL GETVAR_3D(idf_i, idv_i, cf_src, cv_src, Ntr, jte*jct, data3d_src, jt1=j_start, jt2=j_stop)

         IF ((TRIM(cf_lsm_src)=='nan').OR.(TRIM(cf_lsm_src)=='NaN')) THEN
            !! Replacing NaN with 0. to avoid some fuck-up later...
            WHERE(mask_src==0) data3d_src = 0.
         END IF

         CALL INTERP_3D(jt)

         !! Print current record into netcdf file
         !!  => data3d_trg for current time step is ready to be written in netcdf file

         IF (trim(ctype_z_trg) == 'z') THEN
            CALL P3D_T(idf_o, idv_o, Ntr, jt, &
               &       lon_trg_b, lat_trg, REAL(depth_trg(1,1,:),8), vt, data3d_trg,                &
               &       cf_out, cv_lon_trg, cv_lat_trg, cv_z_out, cv_t_out, &
               &       cv_out, rfct_miss*REAL(rmiss_val,4), &
               &       attr_lon=vatt_info_lon_trg, attr_lat=vatt_info_lat_trg, attr_z=vatt_info_z_trg, &
               &       attr_time=vatt_info_t, attr_F=vatt_info_F, &
               &       cextrainfo=cextinf)

            IF ( l_save_drwn ) THEN
               CALL P3D_T(idf_id, idv_id, Ntr, jt,    &
                  &       lon_src, lat_src, REAL(depth_src(1,1,:),8), vt, data_src_drowned(:,:,:), &
                  &       TRIM(cf_src)//'.drwn', cv_lon_src, cv_lat_src, cv_z_src, cv_t_src,    &
                  &       cv_src, 0., &
                  &       attr_lon=vatt_info_lon_src, attr_lat=vatt_info_lat_src, attr_z=vatt_info_z_src, &
                  &       attr_time=vatt_info_t, attr_F=vatt_info_F )
            END IF




         ELSEIF (trim(ctype_z_trg) == 'sigma' ) THEN
            CALL P3D_T(idf_o, idv_o, Ntr, jt, &
               &       lon_trg_b, lat_trg, Sc_rho(:), vt, data3d_trg,                &
               &       cf_out, cv_lon_trg, cv_lat_trg, cv_z_out, cv_t_out, &
               &       cv_out, rfct_miss*REAL(rmiss_val,4), &
               &       attr_lon=vatt_info_lon_trg, attr_lat=vatt_info_lat_trg, attr_z=vatt_info_z_trg, &
               &       attr_time=vatt_info_t, attr_F=vatt_info_F, &
               &       cextrainfo=cextinf)
         ELSE
            PRINT *, 'Unknown vertical coordinate' ; STOP
         ENDIF

      END IF ! .not. l_itrp_3d

   END DO  ! end of time loop

   CALL TERMINATE() ! deleting and de-allocating arrays...

   PRINT *, ''; PRINT *, TRIM(cf_out), ' is created'; PRINT *, ''

   IF ( l_save_drwn ) PRINT *, 'Also saved drowned input field => ', TRIM(cf_src)//'.drwn'

   PRINT *, ''

   CLOSE(6)

END PROGRAM SOSIE
