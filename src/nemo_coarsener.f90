PROGRAM NEMO_COARSENER

   USE netcdf
   USE io_ezcdf
   USE mod_manip
   USE mod_nemo
   USE crs
   USE crsdom

   IMPLICIT NONE

   !! ************************ Configurable part ****************************
   LOGICAL, PARAMETER :: &
      &   l_debug    = .FALSE., &
      &   l_drown_in = .FALSE. ! Not needed since we ignore points that are less than 1 point away from land... drown the field to avoid spurious values right at the coast!
   !!
   INTEGER,PARAMETER :: sp=SELECTED_REAL_KIND(6,37)      ! NF90_FLOAT

   REAL(wp), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER,PARAMETER :: numerr = 0
   INTEGER :: chunksize = 3200000, &
      &       deflate_level = 5
   !!
   INTEGER :: ji, jj, jv, nb_normal
   !!


   CHARACTER(LEN=nf90_max_name) :: filebase, suffix, attname, dimname, varname, time, date, zone, timestamp
   INTEGER :: ncid, ndims, nvars, natts, outid, idim, attid, ntype, varid,  dimlen, rbdims
   INTEGER :: nmax_unlimited
   INTEGER :: dimid, unlimitedDimId

   INTEGER :: nbdim, itype

   CHARACTER(LEN=nf90_max_name), DIMENSION(:),   ALLOCATABLE :: c_list_var_names
   INTEGER,                      DIMENSION(:),   ALLOCATABLE :: i_list_var_types, i_list_var_ndims, i_list_var_ids
   INTEGER,                      DIMENSION(:,:), ALLOCATABLE :: i_list_var_dim_ids
   LOGICAL, DIMENSION(:),   ALLOCATABLE :: l_var_is_done

   INTEGER, DIMENSION(2) :: local_sizes
   INTEGER, DIMENSION(:), ALLOCATABLE  :: global_sizes, rebuild_dims
   REAL(sp) :: ValMin, ValMax, InMin, InMax, rmdi

   INTEGER, DIMENSION(4) :: indimids
   INTEGER, DIMENSION(:), ALLOCATABLE  :: outdimids, outdimlens, indimlens, inncids
   REAL(wp), DIMENSION(:), ALLOCATABLE :: mdiVals

   CHARACTER(LEN=nf90_max_name), DIMENSION(:), ALLOCATABLE :: indimnames
   CHARACTER(LEN=nf90_max_name), DIMENSION(2) :: cdims

   LOGICAL :: l_3d, l_findDims = .FALSE.


   !! Coupe stuff:
   REAL(wp), DIMENSION(:), ALLOCATABLE :: Vt
   REAL(4), DIMENSION(:), ALLOCATABLE :: vdepth
   REAL(4), DIMENSION(:,:), ALLOCATABLE :: vdepth_b


   CHARACTER(len=8), DIMENSION(3), PARAMETER :: csurf_var1 = (/ 'sosstsst', 'sosaline', 'sossheig' /)


   !! Grid, default name :
   CHARACTER(len=80) :: &
      &    cv_in, cv_mm, &
      &    cv_t   = 'time_counter', &
      &    cv_lon = 'nav_lon',      & ! input grid longitude name, T-points
      &    cv_lat = 'nav_lat',      & ! input grid latitude name,  T-points
      &    cv_z   = 'nav_lev'         ! input grid latitude name,  T-points


   CHARACTER(len=128), DIMENSION(4)  :: vlist_coor

   CHARACTER(len=256)  :: cr, cmissval_in
   !CHARACTER(len=512)  :: cdir_home, cdir_out, cdir_tmpdir, cdum, cconf
   !!
   !!
   !!******************** End of conf for user ********************************
   !!
   !!               ** don't change anything below **
   !!
   LOGICAL :: l_exist, lmv_in, & ! input field has a missing value attribute
      &       l_coor_info=.false.

   !!
   !!
   CHARACTER(len=400)  :: &
      &    cf_in='', cf_mm='', cf_get_lat_lon='', cf_out=''
   !!
   INTEGER      :: &
      &    jarg, nb_coor, &
      &    id_v, &
      &    i0=0, j0=0, &
      &    ifi=0, ivi=0, &
      &    ifo=0, ivo=0, &
      &    Nt=0, nk=0, &
      &    ni1, nj1, nk1, &
      &    iargc
   !!
   REAL(sp), DIMENSION(1) :: r4
   REAL(wp), DIMENSION(1) :: r8

   INTEGER :: id_x, id_y, id_z, id_t, id_dim_extra, nlextra
   !!
   !INTEGER :: ji_min, ji_max, jj_min, jj_max, nib, njb


   REAL(4),    DIMENSION(:,:), ALLOCATABLE :: xlon, xlat, xdum_r4
   REAL(wp),   DIMENSION(:,:), ALLOCATABLE :: glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif
   REAL(wp),   DIMENSION(:,:), ALLOCATABLE :: e1t, e2t, e1u, e2u, e1v, e2v, e1f, e2f

   REAL(wp),   DIMENSION(:,:,:), ALLOCATABLE :: e3t, e3u, e3v

   INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: imaskt_crs, imasku_crs, imaskv_crs, imaskf_crs
   REAL(4),    DIMENSION(:,:), ALLOCATABLE :: xdum_r4_crs
   REAL(wp),    DIMENSION(:,:), ALLOCATABLE :: xdum_r8_crs
   REAL(wp),    DIMENSION(:),   ALLOCATABLE :: xr8

   INTEGER :: jt
   !!
   !REAL(wp) :: rt, rt0, rdt, &
   !   &       t_min_e, t_max_e, t_min_m, t_max_m, &
   !   &       alpha, beta, t_min, t_max
   !!
   CHARACTER(LEN=2), DIMENSION(6), PARAMETER :: &
      &            clist_opt = (/ '-h','-m','-i','-o','-x','-y' /)

   !REAL(wp) :: lon_min_1, lon_max_1, lon_min_2, lon_max_2, lat_min, lat_max, r_obs

   !REAL(wp) :: lon_min_trg, lon_max_trg, lat_min_trg, lat_max_trg

   INTEGER :: Nb_att_lon, Nb_att_lat, Nb_att_time, Nb_att_vin
   TYPE(var_attr), DIMENSION(nbatt_max) :: &
      &   v_att_list_lon, v_att_list_lat, v_att_list_time, v_att_list_vin

   REAL(4) :: rmissv_in
   !CALL GET_ENVIRONMENT_VARIABLE("HOME", cdir_home)
   !CALL GET_ENVIRONMENT_VARIABLE("TMPDIR", cdir_tmpdir)


   !cdir_out = TRIM(cdir_tmpdir)//'/EXTRACTED_BOXES' ! where to write data!
   !cdir_out = '.'



   !! Getting string arguments :
   !! --------------------------

   jarg = 0

   DO WHILE ( jarg < iargc() )

      jarg = jarg + 1
      CALL getarg(jarg,cr)

      SELECT CASE (trim(cr))

      CASE('-h')
         call usage()

      CASE('-m')
         CALL GET_MY_ARG('mesh_mask', cf_mm)

      CASE('-i')
         CALL GET_MY_ARG('input file', cf_in)

         !CASE('-v')
         !   CALL GET_MY_ARG('input file', cv_in)
         !
      CASE('-o')
         CALL GET_MY_ARG('input file', cf_out)

      CASE('-x')
         CALL GET_MY_ARG('longitude', cv_lon)

      CASE('-y')
         CALL GET_MY_ARG('latitude', cv_lat)

      CASE DEFAULT
         PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
         CALL usage()

      END SELECT

   END DO

   !IF ( (TRIM(cv_in) == '').OR.(TRIM(cf_in) == '') ) THEN
   IF (TRIM(cf_in) == '') THEN
      PRINT *, ''
      PRINT *, 'You must at least specify input file (-i) !!!'
      CALL usage()
   END IF

   IF ( TRIM(cf_out) == '' ) THEN
      PRINT *, ''
      PRINT *, 'You must at least specify output file (-o) !!!'
      CALL usage()
   END IF

   PRINT *, ''
   PRINT *, ''; PRINT *, 'Use "-h" for help'; PRINT *, ''
   PRINT *, ''

   PRINT *, ' * Input file = ', trim(cf_in)
   !PRINT *, '   => associated variable names = ', TRIM(cv_in)
   !PRINT *, '   => associated longitude/latitude/time = ', trim(cv_lon), ', ', trim(cv_lat)


   PRINT *, ''

   !! Name of config: lulu
   !idot = SCAN(cf_in, '/', back=.TRUE.)
   !cdum = cf_in(idot+1:)
   !idot = SCAN(cdum, '.', back=.TRUE.)
   !cconf = cdum(:idot-1)
   !PRINT *, ' *** CONFIG: cconf ='//TRIM(cconf) ; PRINT *, ''


   !! testing longitude and latitude
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INQUIRE(FILE=TRIM(cf_in), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: input file not found! ', TRIM(cf_in)
      call usage()
   END IF
   INQUIRE(FILE=TRIM(cf_mm), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      PRINT *, 'ERROR: mesh_mask file not found! ', TRIM(cf_mm)
      call usage()
   END IF


   !CALL coordinates_from_var_attr(cf_in, cv_in, nb_coor, vlist_coor)
   !IF ( nb_coor > 0 ) THEN
   !   l_coor_info = .TRUE.
   !   IF ( nb_coor /= 3 ) STOP 'ERROR: since this is the 2D version we expect nb_coor == 3'
   !   PRINT *, ''
   !   PRINT *, ' *** We update names of coordinates as follows:'
   !   cv_t   = TRIM(vlist_coor(1))
   !   cv_lat = TRIM(vlist_coor(2))
   !   cv_lon = TRIM(vlist_coor(3))

   !   PRINT *, '    cv_t   = ', TRIM(cv_t)
   !   PRINT *, '    cv_lat = ', TRIM(cv_lat)
   !   PRINT *, '    cv_lon = ', TRIM(cv_lon)
   !   PRINT *, ''
   !END IF



   ! there's always a 'nav_lon' in a NEMO file right?
   CALL DIMS(cf_in, 'nav_lon', ni1, nj1, nk1, Nt)



   !CALL DIMS(cf_in, cv_lat, ni2, nj2, nk, Nt)
   !IF ( (nj1==-1).AND.(nj2==-1) ) THEN
   !   ni = ni1 ; nj = ni2
   !   PRINT *, 'Grid is 1D: ni, nj =', ni, nj
   !   l_reg_src = .TRUE.
   !ELSE
   !   IF ( (ni1==ni2).AND.(nj1==nj2) ) THEN
   !      ni = ni1 ; nj = nj1
   !      PRINT *, 'Grid is 2D: ni, nj =', ni, nj
   !      l_reg_src = .FALSE.
   !   ELSE
   !      PRINT *, 'ERROR: problem with grid!' ; STOP
   !   END IF
   !END IF

   CALL DIMS(cf_mm, 'tmask', jpi, jpj, jpk, Nt)
   PRINT *, ' *** 3D domain size from meshmask is:'
   PRINT *, '     jpi, jpj, jpk =>', jpi, jpj, jpk
   PRINT *, ''

   IF ( (jpi/=ni1).OR.(jpj/=nj1) ) STOP 'Problem of shape between input field and mesh_mask!'


   jpiglo = jpi
   jpjglo = jpj



   !!LILI
   !! Need to know if there is at least one file with the vertical dimension in the input file!

   CALL check_nf90( NF90_OPEN( TRIM(cf_in), nf90_share, ncid ) )
   CALL check_nf90( NF90_INQUIRE( ncid, ndims, nvars, natts ) )

   PRINT *, ' ncid, ndims, nvars, natts =', ncid, ndims, nvars, natts

   CALL check_nf90( NF90_INQUIRE( ncid, unlimitedDimId = id_t ) )

   WRITE(numout,*) ''
   WRITE(numout,'(" *** input file contains ",i1," dimensions of which time has id ",i1," !")') ndims, id_t


   ALLOCATE(indimlens(ndims), indimnames(ndims), outdimlens(ndims))

   l_3d = .FALSE.

   DO idim = 1, ndims
      CALL check_nf90( nf90_inquire_dimension( ncid, idim, dimname, dimlen ) )
      indimlens(idim) = dimlen
      indimnames(idim) = dimname

      outdimlens(idim) = dimlen
      !! We're going to write a coarsened file so x & y dimension must be set accordingly!
      IF ( TRIM(dimname) == 'x' ) THEN
         id_x = idim
         IF ( .NOT. (dimlen == jpi) ) STOP 'PROBLEM #1'
      ELSEIF ( TRIM(dimname) == 'y' ) THEN
         id_y = idim
         IF ( .NOT. (dimlen == jpj) ) STOP 'PROBLEM #2'
      ELSEIF ( TRIM(dimname) == 'deptht' ) THEN
         id_z = idim
         IF ( .NOT. (dimlen == jpk) ) STOP 'PROBLEM #3'
         IF ( dimlen > 1 ) l_3d = .TRUE.
      END IF

      IF ( idim == id_t ) Nt = dimlen

      WRITE(numout,'(" *** dimension #",i1," => ",a," (len = ",i4.4,")")') idim, TRIM(indimnames(idim)), indimlens(idim)

   END DO

   nb_normal = 4 ; ! x,y,z,t

   IF ( .NOT. l_3d) THEN
      jpk = 1
      nb_normal = 3 ; ! x,y,t
      id_z = -1
   END IF

   IF ( ndims > nb_normal ) THEN
      !! There is this annoying extra dimension:
      IF ( ndims == nb_normal+1 ) THEN
         DO idim = 1, ndims
            IF ( .NOT. ANY( (/id_x,id_y,id_z,id_t/) == idim ) ) THEN
               id_dim_extra = idim
               EXIT
            END IF
         END DO
         nlextra = outdimlens(id_dim_extra)
      ELSE
         STOP 'PROBLEM: there is more than 1 extra weird dimension!!!'
      END IF
   END IF


   WRITE(numout,*) ''
   WRITE(numout,*) '  *** id_x =', id_x, jpi
   WRITE(numout,*) '  *** id_y =', id_y, jpj
   IF ( l_3d ) WRITE(numout,*) '  *** id_z =', id_z, jpk
   WRITE(numout,*) '  *** id_t =', id_t, Nt, ' time records!'
   WRITE(numout,*) '  *** id_dim_extra =', id_dim_extra, nlextra
   WRITE(numout,*) '' ; WRITE(numout,*) ''


   IF ( l_3d )  ALLOCATE ( vdepth(jpk) , vdepth_b(jpk,nlextra) )


   jpkm1 = MAX(jpk-1,1)


   !ni = ni1 ; jpj = ni1
   !! Source:
   ALLOCATE ( xlon(jpi,jpj), xlat(jpi,jpj), xdum_r4(jpi,jpj), &
      &       glamt(jpi,jpj), gphit(jpi,jpj), glamu(jpi,jpj), gphiu(jpi,jpj), glamv(jpi,jpj), gphiv(jpi,jpj), glamf(jpi,jpj), gphif(jpi,jpj),  &
      &       e1t(jpi,jpj), e2t(jpi,jpj), e1u(jpi,jpj), e2u(jpi,jpj), e1v(jpi,jpj), e2v(jpi,jpj), e1f(jpi,jpj), e2f(jpi,jpj) )

   ALLOCATE ( e3t(jpi,jpj,jpk) )




   !! Getting source land-sea mask:
   PRINT *, ''
   ALLOCATE ( tmask(jpi,jpj,jpk), umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk), fmask(jpi,jpj,jpk) )
   PRINT *, ' *** Reading land-sea mask'
   CALL GETMASK_3D(cf_mm, 'tmask', tmask, jz1=1, jz2=jpk)
   CALL GETMASK_3D(cf_mm, 'umask', umask, jz1=1, jz2=jpk)
   CALL GETMASK_3D(cf_mm, 'vmask', vmask, jz1=1, jz2=jpk)
   CALL GETMASK_3D(cf_mm, 'fmask', fmask, jz1=1, jz2=jpk)
   PRINT *, ' Done!'; PRINT *, ''










   !! Getting model longitude & latitude:
   ! Longitude array:
   !PRINT *, ''
   !PRINT *, ' *** Going to fetch longitude array:'
   !CALL GETVAR_ATTRIBUTES(cf_get_lat_lon, cv_lon,  Nb_att_lon, v_att_list_lon)
   !!PRINT *, '  => attributes are:', v_att_list_lon(:Nb_att_lon)
   !CALL GETVAR_2D(i0, j0, cf_get_lat_lon, cv_lon, 0, 0, 0, xlon) ; i0=0 ; j0=0
   !PRINT *, '  '//TRIM(cv_lon)//' sucessfully fetched!'; PRINT *, ''

   !! Latitude array:
   !PRINT *, ''
   !PRINT *, ' *** Going to fetch latitude array:'
   !CALL GETVAR_ATTRIBUTES(cf_get_lat_lon, cv_lat,  Nb_att_lat, v_att_list_lat)
   !!PRINT *, '  => attributes are:', v_att_list_lat(:Nb_att_lat)
   !CALL GETVAR_2D   (i0, j0, cf_get_lat_lon, cv_lat, 0, 0, 0, xlat)
   !i0=0 ; j0=0
   !PRINT *, '  '//TRIM(cv_lat)//' sucessfully fetched!'; PRINT *, ''

   !CALL CHECK_4_MISS(cf_in, cv_in, lmv_in, rmissv_in, cmissval_in)
   !IF ( .not. lmv_in ) rmissv_in = 0.

   !CALL GETVAR_ATTRIBUTES(cf_in, cv_t,  Nb_att_time, v_att_list_time)
   !PRINT *, '  => attributes of '//TRIM(cv_t)//' are:', v_att_list_time(:Nb_att_time)
   !CALL GETVAR_ATTRIBUTES(cf_in, cv_in,  Nb_att_vin, v_att_list_vin)
   !PRINT *, '  => attributes of '//TRIM(cv_in)//' are:', v_att_list_vin(:Nb_att_vin)


   ALLOCATE ( Vt(Nt) )
   CALL GETVAR_1D(cf_in, cv_t, Vt)
   PRINT *, 'Vt = ', Vt(:)






   !!LOLO:
   jperio = 0

   noso = -1 !lolo


   !! We need to fake a MPP with only 1 proc...
   jpni  = 1
   jpnj  = 1
   jpnij = 1
   narea = 1 ! (1 proc)
   ALLOCATE ( nfipproc(jpni,jpnj), mig(jpi), mjg(jpj), nimppt(jpnij), njmppt(jpnij) )
   nfipproc(:,:) = 0 !???
   jpizoom = 1
   jpjzoom = 1
   nimpp = 1
   njmpp =1
   DO ji = 1, jpi                 ! local domain indices ==> data domain indices
      mig(ji) = ji + jpizoom - 1 + nimpp - 1
   END DO
   DO jj = 1, jpj
      mjg(jj) = jj + jpjzoom - 1 + njmpp - 1
   END DO



   npolj  = 0
   nperio = 0


   nlci  = jpi
   nlcj  = jpj

   nimpp = 1
   njmpp = 1

   nldi   =   1
   nldj   =   1
   nlei   = jpi
   nlej   = jpj


   rfactx_r = 1. / nn_factx
   rfacty_r = 1. / nn_facty


   !---------------------------------------------------------
   ! 2. Define Global Dimensions of the coarsened grid
   !---------------------------------------------------------
   CALL crs_dom_def

   PRINT *, ' *** After crs_dom_def: jpi_crs, jpj_crs =', jpi_crs, jpj_crs

   PRINT *, 'TARGET coarsened horizontal domain, jpi_crs, jpj_crs =', jpi_crs, jpj_crs

   PRINT *, ' *** nn_factx, nn_facty', nn_factx, nn_facty

   PRINT *, ''



   CALL GETVAR_2D(i0, j0, cf_mm, 'glamt', 0, 0, 0, glamt) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'gphit', 0, 0, 0, gphit) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'glamu', 0, 0, 0, glamu) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'gphiu', 0, 0, 0, gphiu) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'glamv', 0, 0, 0, glamv) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'gphiv', 0, 0, 0, gphiv) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'glamf', 0, 0, 0, glamf) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'gphif', 0, 0, 0, gphif) ; i0=0 ; j0=0


   CALL GETVAR_2D(i0, j0, cf_mm, 'e1t', 0, 0, 0, e1t) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'e2t', 0, 0, 0, e2t) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'e1u', 0, 0, 0, e1u) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'e2u', 0, 0, 0, e2u) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'e1v', 0, 0, 0, e1v) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'e2v', 0, 0, 0, e2v) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'e1f', 0, 0, 0, e1f) ; i0=0 ; j0=0
   CALL GETVAR_2D(i0, j0, cf_mm, 'e2f', 0, 0, 0, e2f) ; i0=0 ; j0=0

   !! Should get from VVL! Not meshmask:
   CALL GETVAR_2D(i0, j0, cf_mm, 'e3t_0', 0, 1, 0, e3t(:,:,1)) ; i0=0 ; j0=0





   CALL DUMP_2D_FIELD(REAL(tmask(:,:,1),4), 'tmask.tmp', 'tmask' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(glamt(:,:)  ,4), 'glamt.tmp', 'glamt' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(gphit(:,:)  ,4), 'gphit.tmp', 'gphit' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e1t(:,:)  ,4), 'e1t.tmp', 'e1t' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e2t(:,:)  ,4), 'e2t.tmp', 'e2t' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)



   !---------------------------------------------------------
   ! 3. Mask and Mesh
   !---------------------------------------------------------
   !     Set up the masks and meshes
   !  3.a. Get the masks
   CALL crs_dom_msk




   CALL DUMP_2D_FIELD(REAL(tmask_crs(:,:,1),4), 'tmask_crs.tmp', 'tmask_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(umask_crs(:,:,1),4), 'umask_crs.tmp', 'umask_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(vmask_crs(:,:,1),4), 'vmask_crs.tmp', 'vmask_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(fmask_crs(:,:,1),4), 'fmask_crs.tmp', 'fmask_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)


   PRINT *, ''
   PRINT *, ' nrestx, nresty = ', nrestx, nresty
   PRINT *, ''

   gphit_crs = 0.0
   glamt_crs = 0.0
   gphiu_crs = 0.0
   glamu_crs = 0.0
   gphiv_crs = 0.0
   glamv_crs = 0.0
   gphif_crs = 0.0
   glamf_crs = 0.0



   IF ( nresty /= 0 .AND. nrestx /= 0 ) THEN
      CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs )
      CALL crs_dom_coordinates( gphiu, glamu, 'U', gphiu_crs, glamu_crs )
      CALL crs_dom_coordinates( gphiv, glamv, 'V', gphiv_crs, glamv_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
   ELSEIF ( nresty /= 0 .AND. nrestx == 0 ) THEN
      CALL crs_dom_coordinates( gphiu, glamu, 'T', gphit_crs, glamt_crs )
      CALL crs_dom_coordinates( gphiu, glamu, 'U', gphiu_crs, glamu_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'V', gphiv_crs, glamv_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
   ELSEIF ( nresty == 0 .AND. nrestx /= 0 ) THEN
      CALL crs_dom_coordinates( gphiv, glamv, 'T', gphit_crs, glamt_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'U', gphiu_crs, glamu_crs )
      CALL crs_dom_coordinates( gphiv, glamv, 'V', gphiv_crs, glamv_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
   ELSE
      CALL crs_dom_coordinates( gphif, glamf, 'T', gphit_crs, glamt_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'U', gphiu_crs, glamu_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'V', gphiv_crs, glamv_crs )
      CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )
   ENDIF


   CALL DUMP_2D_FIELD(REAL(glamt_crs(:,:),4), 'glamt_crs.tmp', 'glamt_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(gphit_crs(:,:),4), 'gphit_crs.tmp', 'gphit_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(glamu_crs(:,:),4), 'glamu_crs.tmp', 'glamu_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(gphiu_crs(:,:),4), 'gphiu_crs.tmp', 'gphiu_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(glamv_crs(:,:),4), 'glamv_crs.tmp', 'glamv_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(gphiv_crs(:,:),4), 'gphiv_crs.tmp', 'gphiv_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(glamf_crs(:,:),4), 'glamf_crs.tmp', 'glamf_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(gphif_crs(:,:),4), 'gphif_crs.tmp', 'gphif_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)




   e1t_crs = 0.0
   e2t_crs = 0.0
   e1u_crs = 0.0
   e2u_crs = 0.0
   e1v_crs = 0.0
   e2v_crs = 0.0
   e1f_crs = 0.0
   e2f_crs = 0.0

   CALL crs_dom_hgr( e1t, e2t, 'T', e1t_crs, e2t_crs )
   CALL crs_dom_hgr( e1u, e2u, 'U', e1u_crs, e2u_crs )
   CALL crs_dom_hgr( e1v, e2v, 'V', e1v_crs, e2v_crs )
   CALL crs_dom_hgr( e1f, e2f, 'F', e1f_crs, e2f_crs )

   CALL DUMP_2D_FIELD(REAL(e1t_crs(:,:),4), 'e1t_crs.tmp', 'e1t_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e2t_crs(:,:),4), 'e2t_crs.tmp', 'e2t_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e1u_crs(:,:),4), 'e1u_crs.tmp', 'e1u_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e2u_crs(:,:),4), 'e2u_crs.tmp', 'e2u_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e1v_crs(:,:),4), 'e1v_crs.tmp', 'e1v_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e2v_crs(:,:),4), 'e2v_crs.tmp', 'e2v_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e1f_crs(:,:),4), 'e1f_crs.tmp', 'e1f_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   CALL DUMP_2D_FIELD(REAL(e2f_crs(:,:),4), 'e2f_crs.tmp', 'e2f_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)




   PRINT *, ''
   PRINT *, ''
   PRINT *, ''



   !2. Read in the global dimensions from the input file and set up the output file






   !2.1 Set up the output file
   CALL check_nf90( nf90_create( TRIM(cf_out), nf90_netcdf4, outid, chunksize=chunksize ) )


   !2.2.0 Find out how many dimensions are required to be rebuilt and which ones they are
   !CALL check_nf90( nf90_inquire_attribute( ncid, nf90_global, 'DOMAIN_dimensions_ids', itype, rbdims, attid ) )
   !ALLOCATE(rebuild_dims(rbdims))
   !CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_dimensions_ids', rebuild_dims ) )
   !ALLOCATE(global_sizes(rbdims))
   !CALL check_nf90( nf90_get_att( ncid, nf90_global, 'DOMAIN_size_global', global_sizes ) )
   !IF (l_verbose) WRITE(numout,*) 'Size of global arrays: ', global_sizes


   !2.2.1 Copy the dimensions into the output file apart from rebuild_dims() which are dimensioned globally



   !! Defining dimention in output file:

   outdimlens(id_x) = jpi_crs
   outdimlens(id_y) = jpj_crs

   DO idim = 1, ndims
      CALL check_nf90( nf90_inquire_dimension( ncid, idim, dimname, dimlen ) )
      IF( idim == id_t ) THEN
         CALL check_nf90( nf90_def_dim( outid, dimname, nf90_unlimited, dimid) )
         nmax_unlimited = dimlen
      ELSE
         PRINT *, ' nf90_def_dim( outid, dimname, outdimlens(idim), dimid)'
         CALL check_nf90( nf90_def_dim( outid, dimname, outdimlens(idim), dimid) )
      ENDIF
   END DO






   ! nmax_unlimited is only used for time-chunking so we set it to be at least 1 to
   ! account for files with no record dimension or zero length record dimension(!)
   nmax_unlimited = max(nmax_unlimited,1)
   WRITE(numout,*) ''


   IF ( ndims > id_t ) THEN
      IF ( ndims == id_t+1 ) THEN
         id_dim_extra = ndims
         nlextra = outdimlens(ndims)
      ELSE
         STOP 'PROBLEM: there is more than 1 extra weird dimension!!!'
      END IF
   END IF






   !2.2.3 Copy the variable definitions and attributes into the output file.
   ALLOCATE( mdiVals(nvars), c_list_var_names(nvars), i_list_var_types(nvars), &
      &      i_list_var_ndims(nvars), i_list_var_ids(nvars), i_list_var_dim_ids(4,nvars), &
      &      l_var_is_done(nvars) )

   ALLOCATE ( xr8(nlextra) )





   i_list_var_dim_ids(:,:) = 0
   mdiVals(:)=0

   DO jv = 1, nvars

      id_v = jv ! lolo!
      indimids(:) = 0

      CALL check_nf90( nf90_inquire_variable( ncid, id_v, varname, ntype, ndims, indimids, natts ) )
      c_list_var_names(jv)     = TRIM(varname)
      i_list_var_types(jv)     = ntype
      i_list_var_ndims(jv)     = ndims
      i_list_var_dim_ids(:,jv) = indimids(:)
      !CALL check_nf90( nf90_inquire_variable( ncid, id_v) !, ndims=nd)
      i_list_var_ids(jv)   = id_v

      ALLOCATE(outdimids(ndims))
      DO idim = 1, ndims
         outdimids(idim) = indimids(idim)
      END DO



      CALL check_nf90( nf90_def_var( outid, varname, ntype, outdimids, varid, &
         deflate_level=deflate_level ) )
      DEALLOCATE(outdimids)
      WRITE(numout,*) 'Defining variable '//TRIM(varname)//'...'
      IF( natts > 0 ) THEN
         DO attid = 1, natts
            CALL check_nf90( nf90_inq_attname( ncid, varid, attid, attname ) )
            IF ( attname == "_FillValue" ) THEN
               CALL check_nf90( nf90_get_att( ncid, varid, attname, rmdi ) )
               mdiVals(jv)=rmdi
            ENDIF
            CALL check_nf90( nf90_copy_att( ncid, varid, attname, outid, varid ) )
         END DO
      ENDIF
   END DO



   !2.3 End definitions in output file and copy 1st file ncid to the inncids array

   CALL check_nf90( nf90_enddef( outid ) )
   !inncids(1) = ncid
   !WRITE(numout,*) 'Finished defining output file.'



   ALLOCATE ( xdum_r4_crs(jpi_crs,jpj_crs), xdum_r8_crs(jpi_crs,jpj_crs) )





   PRINT *, ''

   l_var_is_done(:) = .FALSE.

   ! A. Writing all variables that do not have a time record:
   !PRINT *, ' *** ID of unlimited (record) dimension is:', id_t

   DO jv = 1, nvars

      cv_in = TRIM(c_list_var_names(jv))
      id_v        = i_list_var_ids(jv)
      itype       = i_list_var_types(jv)
      nbdim       = i_list_var_ndims(jv)
      indimids(:) = i_list_var_dim_ids(:,jv)

      IF ( .NOT. ANY(indimids == id_t) ) THEN
         PRINT *, '  *** Variable '//TRIM(cv_in)//' does not have a time record, writing it before time loop!'

         ! LOLO add 1D like depth here !!!! IF( nbdim == 1 ) THEN

         IF( nbdim == 1 ) THEN
            IF ( (TRIM(cv_in)=='deptht').OR.(TRIM(cv_in)=='depthw') ) THEN
               CALL check_nf90( nf90_get_var( ncid,  id_v, vdepth ) )!lili
               CALL check_nf90( nf90_put_var( outid, id_v, vdepth ) )
               l_var_is_done(jv) = .TRUE.
            END IF


         ELSEIF( nbdim == 2 ) THEN
            IF ( (TRIM(cv_in)=='nav_lon').OR.(TRIM(cv_in)=='glamt') ) THEN
               ! It's beed coarsened earlier! with crs_dom_coordinates => !CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs )
               WRITE(numout,*) '   *** writing coarsened '//TRIM(cv_in)//' in '//TRIM(cf_out)
               CALL check_nf90( nf90_put_var( outid, jv, glamt_crs ) )
               l_var_is_done(jv) = .TRUE.
            ELSEIF ( (TRIM(cv_in)=='nav_lat').OR.(TRIM(cv_in)=='gphit') ) THEN
               ! It's beed coarsened earlier! with crs_dom_coordinates => !CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs )
               WRITE(numout,*) '   *** writing coarsened '//TRIM(cv_in)//' in '//TRIM(cf_out)
               CALL check_nf90( nf90_put_var( outid, jv, gphit_crs ) )
               l_var_is_done(jv) = .TRUE.
            ELSEIF ( (TRIM(cv_in)=='deptht_bounds').OR.(TRIM(cv_in)=='depthw_bounds') ) THEN
               !CALL check_nf90( nf90_get_var( ncid,  id_v, vdepth_b ) )
               !CALL check_nf90( nf90_put_var( outid, id_v, vdepth_b ) )
               l_var_is_done(jv) = .TRUE.

            ELSE
               WRITE(numout,*) 'ERROR: unknown variable without time)!!!'; STOP
            END IF !IF ( (TRIM(cv_in)=='nav_lon').OR.(TRIM(cv_in)=='glamt') )

         ELSE

            WRITE(numout,*) 'ERROR: unknown number of dimensions for variable without time)!!!'; STOP

         END IF ! IF( nbdim == 1 )

      END IF ! IF ( .NOT. ANY(indimids == id_t) )

   END DO  ! DO jv = 1, nvars


   STOP 'lala'


   !! Everything that depends on time record:

   e3t(:,:,:) = 1._wp

   DO jt=1, Nt

      PRINT *, ''

      DO jv = 1, nvars
         IF ( .NOT. l_var_is_done(jv) ) THEN

            PRINT *, ''
            cv_in = TRIM(c_list_var_names(jv))
            id_v        = i_list_var_ids(jv)
            itype       = i_list_var_types(jv)
            nbdim       = i_list_var_ndims(jv)
            indimids(:) = i_list_var_dim_ids(:,jv)

            PRINT *, ' ### Reading field '//TRIM(cv_in)//' at record #',jt
            PRINT *, ' * var ID         =>', id_v
            PRINT *, ' * var type       =>', itype
            PRINT *, ' * var nb of dims =>', nbdim
            PRINT *, ' * dim ids        =>', indimids(:)
            PRINT *, ''

            ! T
            IF( nbdim == 1 ) THEN
               IF ( .NOT. (indimids(1)==id_t) ) STOP 'ERROR: we should only have record-dependant vectors here!!!'
               SELECT CASE( itype )
               CASE( NF90_DOUBLE )
                  CALL check_nf90( nf90_get_var( ncid,  id_v, r8, start=(/jt/), count=(/1/)) )
                  CALL check_nf90( nf90_put_var( outid, id_v, r8, start=(/jt/), count=(/1/)) )
               CASE( NF90_FLOAT )
                  CALL check_nf90( nf90_get_var( ncid,  id_v, r4, start=(/jt/), count=(/1/)) )
                  CALL check_nf90( nf90_put_var( outid, id_v, r4, start=(/jt/), count=(/1/)) )
               CASE DEFAULT
                  WRITE(numout,*) 'ERROR: unknown netcdf type!!!'; STOP
               END SELECT
               WRITE(numout,*) '   *** '//TRIM(cv_in)//' written at record #',jt
               WRITE(numout,*) ''
            END IF


            ! 1D+T
            IF( nbdim == 2 ) THEN
               IF ( .NOT. ((indimids(1)==id_dim_extra).AND.(indimids(2)==id_t)) ) STOP 'ERROR: we do not know wthat this 2D field is! '
               SELECT CASE( itype )
               CASE( NF90_DOUBLE )
                  CALL check_nf90( NF90_GET_VAR(ncid,  id_v, xr8, start=(/1,jt/), count=(/nlextra,1/)) )
                  CALL check_nf90( NF90_PUT_VAR(outid, id_v, xr8, start=(/1,jt/), count=(/nlextra,1/)) )

               CASE DEFAULT
                  WRITE(numout,*) 'ERROR: unknown netcdf type!!!'; STOP
               END SELECT


               WRITE(numout,*) '   *** just read '//TRIM(cv_in)//' at record #',jt




               !WRITE(numout,*) ' Mhhh... We should not have 2D variables in here... =>', TRIM(cv_in); STOP
            END IF

            !! Reading field at record jt:

            ! 2D+T
            IF( nbdim == 3 ) THEN
               IF ( .NOT. ((indimids(1)==id_x).AND.(indimids(2)==id_y).AND.(indimids(3)==id_t)) ) STOP 'ERROR: we do not know wthat this 3D field is! '
               SELECT CASE( itype )
               CASE( NF90_FLOAT )
                  CALL check_nf90( NF90_GET_VAR(ncid, id_v, xdum_r4, start=(/1,1,jt/), count=(/jpi,jpj,1/)) )
               CASE DEFAULT
                  WRITE(numout,*) 'ERROR: unknown netcdf type!!! We should not have 2D+T arrays as double precision!'; STOP
               END SELECT
               WRITE(numout,*) '   *** just read '//TRIM(cv_in)//' at record #',jt

               ! Time for coarsening! Use different routines depending on the type of field!
               IF ( (TRIM(cv_in)=='sossheig').OR.(TRIM(cv_in)=='sosstsst').OR.(TRIM(cv_in)=='sosaline') ) THEN
                  WRITE(numout,*) '   *** Coarsening '//TRIM(cv_in)//' with crs_dom_ope / VOL / T !'
                  CALL crs_dom_ope( REAL(xdum_r4,8) , 'VOL', 'T', REAL(tmask,8), xdum_r8_crs , p_e12=e1t*e2t, p_e3=e3t, psgn=1.0_wp )
                  WRITE(numout,*) '   *** writing '//TRIM(cv_in)//' in '//TRIM(cf_out)
                  xdum_r4_crs = xdum_r8_crs
1                 WHERE ( tmask_crs(:,:,1) == 0 ) xdum_r4_crs = 1.e+20
                  CALL check_nf90( nf90_put_var( outid, id_v,  xdum_r4_crs, start=(/1,1,jt/), count=(/jpi_crs,jpj_crs,1/)) )
                  WRITE(numout,*) '   ***  written!'
               END IF

            END IF !IF( nbdim == 3 )

            !WRITE(numout,*) 'ERROR: unknown number of dimmensions!!!'





            !CALL check_nf90( NF90_GET_VAR(ncid, id_v, xdum_r4_crs, start=(/1,1/), count=(/lx,ly/))


         END IF ! IF ( .NOT. l_var_is_done(jv) )


      END DO !DO jv = 1, nvars





      !IF ( lmv_in ) THEN
      !WHERE ( tmask_crs(:,:,1) == 0 ) xdum_r4_crs = rmissv_in
      !END IF



   END DO !DO jt=1, Nt


   CALL check_nf90( nf90_close( outid ) )
   !STOP'lili'


   STOP 'LOLO'





CONTAINS






   SUBROUTINE GET_MY_ARG(cname, cvalue)
      CHARACTER(len=*), INTENT(in)    :: cname
      CHARACTER(len=*), INTENT(inout) :: cvalue
      !!
      IF ( jarg + 1 > iargc() ) THEN
         PRINT *, 'ERROR: Missing ',trim(cname),' name!' ; call usage()
      ELSE
         jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         IF ( ANY(clist_opt == trim(cr)) ) THEN
            PRINT *, 'ERROR: Missing',trim(cname),' name!'; call usage()
         ELSE
            cvalue = trim(cr)
         END IF
      END IF
   END SUBROUTINE GET_MY_ARG





   SUBROUTINE usage()
      !!
      !OPEN(UNIT=6, FORM='FORMATTED', RECL=512)
      !!
      WRITE(6,*) ''
      WRITE(6,*) '   List of command line options:'
      WRITE(6,*) '   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(6,*) ''
      WRITE(6,*) ' -m <mesh_mask.nc>    => file containing grid metrics of model'
      WRITE(6,*) ''
      WRITE(6,*) ' -i <input_file.nc>   => input file to coarsen'
      WRITE(6,*) ''
      !WRITE(6,*) ' -v <name_field>      => variable to coarsen'
      !WRITE(6,*) ''
      WRITE(6,*) ' -o <output_file.nc>  => file to be created'
      WRITE(6,*) ''
      WRITE(6,*) '    Optional:'
      WRITE(6,*) ' -h                   => Show this message'
      WRITE(6,*) ''
      WRITE(6,*) ' -x  <name>           => Specify longitude name in input file (default: '//TRIM(cv_lon)//')'
      WRITE(6,*) ''
      WRITE(6,*) ' -y  <name>           => Specify latitude  name in input file  (default: '//TRIM(cv_lon)//')'
      WRITE(6,*) ''
      WRITE(6,*) ''
      !!
      STOP
      !!
   END SUBROUTINE usage


   SUBROUTINE check_nf90(status, errorFlag)
      !---------------------------------------------------------------------
      !  Checks return code from nf90 library calls and warns if needed
      !  If errorFlag is present then it just increments this flag (OMP use)
      !
      !---------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: status
      INTEGER, INTENT(INOUT), OPTIONAL :: errorFlag
      !---------------------------------------------------------------------
      IF( status /= nf90_noerr ) THEN
         WRITE(numerr,*) 'ERROR! : '//TRIM(nf90_strerror(status))
         IF( PRESENT( errorFlag ) ) THEN
            errorFlag = errorFlag + status
         ELSE
            WRITE(numerr,*) "*** NEMO offline coarsening failed ***"
            WRITE(numerr,*)
            STOP 4
         ENDIF
      ENDIF
   END SUBROUTINE check_nf90






END PROGRAM NEMO_COARSENER




!!
