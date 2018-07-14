PROGRAM NEMO_COARSENER

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

   REAL(wp), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER :: ji, jj
   !!

   !! Coupe stuff:
   REAL(wp), DIMENSION(:), ALLOCATABLE :: Vt


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
      &    i0=0, j0=0, &
      &    ifi=0, ivi=0, &
      &    ifo=0, ivo=0, &
      &    Nt=0, nk=0, &
      &    ni1, nj1, &
      &    iargc
   !!

   !!
   !INTEGER :: ji_min, ji_max, jj_min, jj_max, nib, njb


   INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: imaskt, imasku, imaskv, imaskf
   REAL(4),    DIMENSION(:,:), ALLOCATABLE :: xlon, xlat, xdum_r4
   REAL(wp),   DIMENSION(:,:), ALLOCATABLE :: glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif
   REAL(wp),   DIMENSION(:,:), ALLOCATABLE :: e1t, e2t, e1u, e2u, e1v, e2v, e1f, e2f
   
   REAL(wp),   DIMENSION(:,:,:), ALLOCATABLE :: e3t, e3u, e3v

   INTEGER(1), DIMENSION(:,:), ALLOCATABLE :: imaskt_crs, imasku_crs, imaskv_crs, imaskf_crs
   REAL(4),    DIMENSION(:,:), ALLOCATABLE :: xdum_r4_crs
   REAL(wp),    DIMENSION(:,:), ALLOCATABLE :: xdum_r8_crs
   !REAL(wp),    DIMENSION(:,:), ALLOCATABLE :: glamt_crs, gphit_crs


   !!
   INTEGER :: jt
   !!
   !REAL(wp) :: rt, rt0, rdt, &
   !   &       t_min_e, t_max_e, t_min_m, t_max_m, &
   !   &       alpha, beta, t_min, t_max
   !!
   CHARACTER(LEN=2), DIMENSION(7), PARAMETER :: &
      &            clist_opt = (/ '-h','-m','-i','-v','-o','-x','-y' /)

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

      CASE('-v')
         CALL GET_MY_ARG('input file', cv_in)

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

   IF ( (trim(cv_in) == '').OR.(trim(cf_in) == '') ) THEN
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




   CALL coordinates_from_var_attr(cf_in, cv_in, nb_coor, vlist_coor)
   IF ( nb_coor > 0 ) THEN
      l_coor_info = .TRUE.
      IF ( nb_coor /= 3 ) STOP 'ERROR: since this is the 2D version we expect nb_coor == 3'
      PRINT *, ''
      PRINT *, ' *** We update names of coordinates as follows:'
      cv_t   = TRIM(vlist_coor(1))
      cv_lat = TRIM(vlist_coor(2))
      cv_lon = TRIM(vlist_coor(3))

      PRINT *, '    cv_t   = ', TRIM(cv_t)
      PRINT *, '    cv_lat = ', TRIM(cv_lat)
      PRINT *, '    cv_lon = ', TRIM(cv_lon)
      PRINT *, ''
   END IF



   !cf_get_lat_lon = cf_mm
   cf_get_lat_lon = cf_in

   CALL DIMS(cf_get_lat_lon, cv_lon, ni1, nj1, nk, Nt)
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

   jpiglo = ni1
   jpjglo = nj1

   jpi = jpiglo
   jpj = jpjglo
   jpk = 1 ! 2D
   jpkm1 = MAX(jpk-1,1)

   CALL DIMS(cf_in, cv_in, jpi, jpj, nk, Nt)
   PRINT *, ' *** input field: jpi, jpj, nk, Nt =>', jpi, jpj, nk, Nt
   PRINT *, ''

   IF ( (jpi/=ni1).OR.(jpj/=nj1) ) STOP 'Problem of shape between input field and mesh_mask!'

   !ni = ni1 ; jpj = ni1
   !! Source:
   ALLOCATE ( xlon(jpi,jpj), xlat(jpi,jpj), xdum_r4(jpi,jpj), &
      &       glamt(jpi,jpj), gphit(jpi,jpj), glamu(jpi,jpj), gphiu(jpi,jpj), glamv(jpi,jpj), gphiv(jpi,jpj), glamf(jpi,jpj), gphif(jpi,jpj),  &
      &       e1t(jpi,jpj), e2t(jpi,jpj), e1u(jpi,jpj), e2u(jpi,jpj), e1v(jpi,jpj), e2v(jpi,jpj), e1f(jpi,jpj), e2f(jpi,jpj),  &
      &       imaskt(jpi,jpj),imasku(jpi,jpj),imaskv(jpi,jpj),imaskf(jpi,jpj) )

   ALLOCATE ( e3t(jpi,jpj,jpk) )


   

   !! Getting source land-sea mask:
   PRINT *, '';
   PRINT *, ' *** Reading land-sea mask'
   CALL GETMASK_2D(cf_mm, 'tmask', imaskt)
   CALL GETMASK_2D(cf_mm, 'umask', imasku)
   CALL GETMASK_2D(cf_mm, 'vmask', imaskv)
   CALL GETMASK_2D(cf_mm, 'fmask', imaskf)
   PRINT *, ' Done!'; PRINT *, ''



   !! Target:
   !! Coarsening stuff:
   !jpi_crs = INT( (jpi - 2) / nn_factx ) + 2
   !jpj_crs = INT( (jpj - MOD(jpj, nn_facty)) / nn_facty ) + 3

   !jpiglo_crs = jpi_crs
   !jpjglo_crs = jpj_crs






   !! Getting model longitude & latitude:
   ! Longitude array:
   PRINT *, ''
   PRINT *, ' *** Going to fetch longitude array:'
   CALL GETVAR_ATTRIBUTES(cf_get_lat_lon, cv_lon,  Nb_att_lon, v_att_list_lon)
   !PRINT *, '  => attributes are:', v_att_list_lon(:Nb_att_lon)
   CALL GETVAR_2D(i0, j0, cf_get_lat_lon, cv_lon, 0, 0, 0, xlon) ; i0=0 ; j0=0
   PRINT *, '  '//TRIM(cv_lon)//' sucessfully fetched!'; PRINT *, ''

   ! Latitude array:
   PRINT *, ''
   PRINT *, ' *** Going to fetch latitude array:'
   CALL GETVAR_ATTRIBUTES(cf_get_lat_lon, cv_lat,  Nb_att_lat, v_att_list_lat)
   !PRINT *, '  => attributes are:', v_att_list_lat(:Nb_att_lat)
   CALL GETVAR_2D   (i0, j0, cf_get_lat_lon, cv_lat, 0, 0, 0, xlat)
   i0=0 ; j0=0
   PRINT *, '  '//TRIM(cv_lat)//' sucessfully fetched!'; PRINT *, ''

   CALL CHECK_4_MISS(cf_in, cv_in, lmv_in, rmissv_in, cmissval_in)
   IF ( .not. lmv_in ) rmissv_in = 0.

   CALL GETVAR_ATTRIBUTES(cf_in, cv_t,  Nb_att_time, v_att_list_time)
   !PRINT *, '  => attributes of '//TRIM(cv_t)//' are:', v_att_list_time(:Nb_att_time)
   CALL GETVAR_ATTRIBUTES(cf_in, cv_in,  Nb_att_vin, v_att_list_vin)
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

   !STOP 'LULU'

   ALLOCATE ( tmask(jpi,jpj,jpk), umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk), fmask(jpi,jpj,jpk) )

   tmask(:,:,1) = imaskt(:,:)
   umask(:,:,1) = imasku(:,:)
   vmask(:,:,1) = imaskv(:,:)
   fmask(:,:,1) = imaskf(:,:)


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


   ALLOCATE ( xdum_r4_crs(jpi_crs,jpj_crs), xdum_r8_crs(jpi_crs,jpj_crs) )

   DO jt=1, Nt

      PRINT *, ''
      PRINT *, ' Reading field '//TRIM(cv_in)//' at record #',jt

      CALL GETVAR_2D   (ifi, ivi, cf_in, cv_in, Nt, 0, jt, xdum_r4)

      xdum_r4 = xdum_r4*REAL(imaskt,4)

      IF ( TRIM(cv_in) == 'sossheig' ) THEN
         e3t(:,:,:) = 1._wp
         CALL crs_dom_ope( REAL(xdum_r4,8) , 'VOL', 'T', REAL(tmask,8), xdum_r8_crs , p_e12=e1t*e2t, p_e3=e3t, psgn=1.0_wp )
      END IF
      
      !IF ( ANY( csurf_var1 == TRIM(cv_in) ) ) THEN
      !   CALL crs_dom_ope( REAL(xdum_r4,8), 'VOL', 'T', REAL(tmask,8), xdum_r8_crs, p_e12=e1t*e2t, p_e3=e3t, psgn=1.0_wp )         
      !ELSE
      !   PRINT *, 'Unknown variable: ', TRIM(cv_in) ; PRINT *, ''
      !END IF
      
      xdum_r4_crs = REAL(xdum_r8_crs,4)
      
      !IF ( jt == 1 ) THEN
      !   imaskt(:,:) = 1
      !   WHERE ( xdum_r4 > 10000. ) imaskt = 0
      !   imaskt_crs(:,:) = imaskt(1:jpi,1:jpj)
      !END IF

      !xdum_r4 = xdum_r4*REAL(imaskt,4)
      !xdum_r4 = xdum_r4*xdum_r4


      !xdum_r4_crs(:,:) = xdum_r4(1:jpi,1:jpj)



      IF ( lmv_in ) THEN
         WHERE ( tmask_crs(:,:,1) == 0 ) xdum_r4_crs = rmissv_in
      END IF

      CALL P2D_T( ifo, ivo, Nt, jt, glamt_crs, gphit_crs, Vt, xdum_r4_crs, cf_out, &
         &        cv_lon, cv_lat, cv_t, cv_in, rmissv_in,     &
         &        attr_lon=v_att_list_lon, attr_lat=v_att_list_lat, attr_time=v_att_list_time, &
         &        attr_F=v_att_list_vin, l_add_valid_min_max=.false. )


   END DO




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
      WRITE(6,*) ' -v <name_field>      => variable to coarsen'
      WRITE(6,*) ''
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




END PROGRAM NEMO_COARSENER




!!
