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

   !$ USE omp_lib
   USE io_ezcdf     !* routines for netcdf input/output
   USE mod_conf     !* important parameters, namelist, arrays, etc...
   USE mod_init     !* important parameters, namelist, arrays, etc...
   USE mod_grids    !* grided-domain module
   USE mod_interp


   IMPLICIT NONE

   CHARACTER(len=800) :: cextinf

   INTEGER    :: &
      &   jt, jte, jct, jcz,  &
      &   idf_o, idf_i, idv_o, idv_i, &  !: netcdf ID for source/target file
      &   idf_id, idv_id  !: drowned source field...

   INTEGER :: noinc, ileft, ipl, jo

   REAL(4) :: rfct_miss=1.

   !! for timing:
   CHARACTER(LEN=8)      :: cldate
   CHARACTER(LEN=10)     :: cltime
   CHARACTER(LEN=5)      :: clzone
   INTEGER, DIMENSION(8), SAVE :: ivalue1, ivalue2

   INTEGER :: NTHREADS, TID


   !OPEN(UNIT=6, FORM='FORMATTED', RECL=512)  ! problem with Gfortan 4.8...

   WRITE(6,*)''
   WRITE(6,*)'=========================================================='
   WRITE(6,*)'            S  O  S  I  E    version ', trim(csosie_version)
   WRITE(6,*)'=========================================================='
   WRITE(6,*)''

   Nthrd = 1


   CALL DATE_AND_TIME(cldate,cltime,clzone,ivalue1)
   !PRINT '(i8,1x,8i4)', istp, ivalue(:)   ! DRAKKAR code : print time step for run monitoring
   WRITE(6,'(" ### staring at: ",8i4)') ivalue1(:)





   !! Fecthing command line arguments if any:
   CALL GET_ARGUMENTS()    ! MODULE mod_init


   !! ========== OpenMP only ==========================================
   !$ CALL OMP_SET_NUM_THREADS(Nthrd_fix)
   !$OMP PARALLEL PRIVATE(NTHREADS, TID)
   !$ TID = OMP_GET_THREAD_NUM()
   !$ IF (TID == 0) Nthrd = OMP_GET_NUM_THREADS()
   !$OMP END PARALLEL
   !$ WRITE(6,'(" ### Going to use ",i2," OpenMP threads!")') Nthrd
   !!===================================================================


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


   !! OMP decomposition
   !
   !$ PRINT *, ''
   !$ PRINT *, ' OMP thread decomposition along target grid X-axis:'
   !$ PRINT *, '  Nthrd =', Nthrd
   !$ PRINT *, '  ni_trg=', ni_trg
   !PRINT *, '  ni_trg/Nthrd , ni_trg%Nthrd =>', ni_trg/Nthrd , MOD(ni_trg,Nthrd)
   !PRINT *, '   noinc =', noinc

   !! ALONG I:
   ALLOCATE ( io1(Nthrd), io2(Nthrd), i_seg_s(Nthrd) )
   IF ( Nthrd > 1) THEN
      noinc = ni_trg/Nthrd
      ipl = 0
      io1(1) = 1
      io2(1) = io1(1) + noinc -1
      DO jo = 2, Nthrd
         ileft = Nthrd - jo + 1
         io1(jo) = io1(jo-1) + noinc + ipl
         IF ( ileft == MOD(ni_trg,Nthrd) ) ipl=1
         io2(jo) = io1(jo)   + noinc -1 + ipl
      END DO
   ELSE
      !! No OMP (Nthrd==1):
      io1(1) = 1
      io2(1) = ni_trg
      i_seg_s(1) = ni_trg
   END IF

   PRINT *, '  io1 =', io1
   PRINT *, '  io2 =', io2
   DO jo = 1, Nthrd
      i_seg_s(jo) = SIZE(lon_trg(io1(jo):io2(jo),0)) ! lazy!!! but trustworthy I guess...
      PRINT *, '  SIZE segment #',INT(jo,1),' => ', i_seg_s(jo)
   END DO
   PRINT *, ''




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

   ALLOCATE ( l_first_call_interp_routine(Nthrd) )
   l_first_call_interp_routine(:) = .true.

   ALLOCATE( l_always_first_call(Nthrd) )
   
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

         IF ((TRIM(cf_lsm_src)=='nan').OR.(TRIM(cf_lsm_src)=='NaN')) THEN
            !! Replacing NaN with 0. to avoid some fuck-up later...
            WHERE(mask_src(:,:,1)==0) data_src = 0.
         END IF


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
            PRINT *, 'Not! sosie.f90 l_save_drwn'
            !CALL P2D_T(idf_id, idv_id, Ntr, jt,    &
            !   !&       lon_src, lat_src, vt, data_src_drowned(:,:,1),    &
            !   &       lon_src, lat_src, vt, data_src(:,:,1),    &
            !   &       TRIM(cf_src)//'.drwn', cv_lon_src, cv_lat_src, cv_t_src,    &
            !   &       cv_src, 0., &
            !   &       attr_lon=vatt_info_lon_src, attr_lat=vatt_info_lat_src, attr_time=vatt_info_t, attr_F=vatt_info_F )
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

         CALL INTERP_3D()

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
               PRINT *, 'Not! sosie.f90 l_save_drwn'
               !CALL P3D_T(idf_id, idv_id, Ntr, jt,    &
               !   &       lon_src, lat_src, REAL(depth_src(1,1,:),8), vt, data_src_drowned(:,:,:), &
               !   &       TRIM(cf_src)//'.drwn', cv_lon_src, cv_lat_src, cv_z_src, cv_t_src,    &
               !   &       cv_src, 0., &
               !   &       attr_lon=vatt_info_lon_src, attr_lat=vatt_info_lat_src, attr_z=vatt_info_z_src, &
               !   &       attr_time=vatt_info_t, attr_F=vatt_info_F )
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


   CALL DATE_AND_TIME(cldate,cltime,clzone,ivalue2)
   WRITE(6,'(" ###  started at: ",8i4)') ivalue1(:)
   WRITE(6,'(" ### stopping at: ",8i4)') ivalue2(:)

   PRINT *, ''
   rfct_miss = ivalue1(8)/1000. + ivalue1(7) + ivalue1(6)*60. + ivalue1(5)*3600.
   rfct_miss = ivalue2(8)/1000. + ivalue2(7) + ivalue2(6)*60. + ivalue2(5)*3600. - rfct_miss
   WRITE(6,'(" ### Total time in seconds: ",f11.4)') rfct_miss




   CLOSE(6)

END PROGRAM SOSIE
