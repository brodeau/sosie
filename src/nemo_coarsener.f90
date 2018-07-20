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
      &   l_debug    = .TRUE., &
      &   l_write_crs_mm = .TRUE. !: write the coarsened mesh_mask.nc ! need 3D to be on!!!
   !!
   INTEGER,PARAMETER :: sp=SELECTED_REAL_KIND(6,37)      ! NF90_FLOAT

   REAL(wp), PARAMETER :: res = 0.1  ! resolution in degree
   !!
   INTEGER,PARAMETER :: numerr = 0
   INTEGER :: chunksize = 3200000, &
      &       deflate_level = 5
   !!
   INTEGER :: ji, jj, jk, jt, jv, nb_normal
   !!


   CHARACTER(LEN=nf90_max_name) :: cnm_att, cnm_dim, cnm_var
   INTEGER :: idf_src, ndims, nvars, natts, idf_trg, idim, attid, ntype, varid,  dimlen
   INTEGER :: nmax_unlimited
   INTEGER :: dimid

   INTEGER :: idf_mm_crs, idd_x, idd_y, idd_z, idd_t, &
      &       idv_nlon, idv_nlat, idv_nlev, &
      &       idv_tmsk, idv_umsk, idv_vmsk, idv_fmsk 
   
   INTEGER :: nbdim, itype

   CHARACTER(LEN=nf90_max_name), DIMENSION(:),   ALLOCATABLE :: c_list_var_names
   INTEGER,                      DIMENSION(:),   ALLOCATABLE :: i_list_var_types, i_list_var_ndims, i_list_var_ids
   INTEGER,                      DIMENSION(:,:), ALLOCATABLE :: i_list_var_dim_ids
   LOGICAL, DIMENSION(:),   ALLOCATABLE :: l_var_is_done

   REAL(sp) :: rmdi

   INTEGER, DIMENSION(4) :: indimids
   INTEGER, DIMENSION(:), ALLOCATABLE  :: outdimids, outdimlens, indimlens
   REAL(wp), DIMENSION(:), ALLOCATABLE :: mdiVals

   CHARACTER(LEN=nf90_max_name), DIMENSION(:), ALLOCATABLE :: indimnames

   LOGICAL :: l_3d, l_exist


   !! Coupe stuff:
   REAL(wp), DIMENSION(:), ALLOCATABLE :: Vt
   REAL(4), DIMENSION(:), ALLOCATABLE :: vdepth
   REAL(4), DIMENSION(:,:), ALLOCATABLE :: vdepth_b


   CHARACTER(len=8), DIMENSION(3), PARAMETER :: csurf_var1 = (/ 'sosstsst', 'sosaline', 'sossheig' /)


   !! Grid, default name :
   CHARACTER(len=80) :: &
      &    cv_in, &
      &    cv_t   = 'time_counter', &
      &    cv_lon = 'nav_lon',      & ! input grid longitude name, T-points
      &    cv_lat = 'nav_lat'         ! input grid latitude name,  T-points




   CHARACTER(len=1)   :: cgp=''   !: grid point (T/U/V/W)
   CHARACTER(len=256) :: cr

   CHARACTER(len=400)  :: &
      &    cf_in='', cf_mm='', cf_out=''

   INTEGER      :: &
      &    jarg, iargc, &
      &    id_v, &
      &    i0=0, j0=0, &
      &    Nt=0, &
      &    ni1, nj1, nk1

   REAL(sp), DIMENSION(1) :: r4
   REAL(wp), DIMENSION(1) :: r8

   INTEGER :: id_x, id_y, id_z, id_t, id_dim_axbnd, nlaxbnd
   !!
   !INTEGER :: ji_min, ji_max, jj_min, jj_max, nib, njb


   REAL(4),  DIMENSION(:,:),   ALLOCATABLE :: x2d_r4, x2d_r4_crs
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: x2d_r8_crs
   REAL(4),  DIMENSION(:,:,:), ALLOCATABLE :: x3d_r4, x3d_r4_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: x3d_r8_crs, zfse3t, zfse3u, zfse3v, zfse3w, e3t_max_crs
   REAL(wp), DIMENSION(:),     ALLOCATABLE :: xr8

   CHARACTER(LEN=2), DIMENSION(5), PARAMETER :: clist_opt = (/ '-h','-m','-i','-P','-o' /)

   !CALL GET_ENVIRONMENT_VARIABLE("HOME", cdir_home)
   !CALL GET_ENVIRONMENT_VARIABLE("TMPDIR", cdir_tmpdir)
   !cdir_out = TRIM(cdir_tmpdir)//'/EXTRACTED_BOXES' ! where to write data!
   !cdir_out = '.'


   OPEN(UNIT=numout, FORM='FORMATTED', RECL=512)



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

      CASE('-P')
         CALL GET_MY_ARG('grid point', cgp)

      CASE('-o')
         CALL GET_MY_ARG('input file', cf_out)

      CASE('-x')
         CALL GET_MY_ARG('longitude', cv_lon)

      CASE('-y')
         CALL GET_MY_ARG('latitude', cv_lat)

      CASE DEFAULT
         WRITE(numout,*) 'Unknown option: ', trim(cr) ; WRITE(numout,*) ''
         CALL usage()

      END SELECT

   END DO

   IF (TRIM(cf_in) == '') THEN
      WRITE(numout,*) ''
      WRITE(numout,*) 'Specify the input file with the "-i" switch !!!'
      CALL usage()
   END IF

   IF (TRIM(cgp) == '') THEN
      WRITE(numout,*) ''
      WRITE(numout,*) 'Specify the C-grid point with the "-P" switch !!!'
      CALL usage()
   END IF

   IF ( TRIM(cf_out) == '' ) THEN
      WRITE(numout,*) ''
      WRITE(numout,*) 'You must at least specify output file (-o) !!!'
      CALL usage()
   END IF

   WRITE(numout,*) ''
   WRITE(numout,*) ''; WRITE(numout,*) 'Use "-h" for help'; WRITE(numout,*) ''
   WRITE(numout,*) ''

   WRITE(numout,*) ' * Input file = ', trim(cf_in)


   WRITE(numout,*) ''

   !! Name of config: lulu
   !idot = SCAN(cf_in, '/', back=.TRUE.)
   !cdum = cf_in(idot+1:)
   !idot = SCAN(cdum, '.', back=.TRUE.)
   !cconf = cdum(:idot-1)
   !WRITE(numout,*) ' *** CONFIG: cconf ='//TRIM(cconf) ; WRITE(numout,*) ''


   !! testing longitude and latitude
   !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   INQUIRE(FILE=TRIM(cf_in), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      WRITE(numout,*) 'ERROR: input file not found! ', TRIM(cf_in)
      call usage()
   END IF
   INQUIRE(FILE=TRIM(cf_mm), EXIST=l_exist )
   IF ( .NOT. l_exist ) THEN
      WRITE(numout,*) 'ERROR: mesh_mask file not found! ', TRIM(cf_mm)
      call usage()
   END IF

   ! there's always a 'nav_lon' in a NEMO file right?
   CALL DIMS(cf_in, 'nav_lon', ni1, nj1, nk1, Nt)

   CALL DIMS(cf_mm, 'tmask', jpi, jpj, jpk, Nt)
   WRITE(numout,*) ' *** 3D domain size from meshmask is:'
   WRITE(numout,*) '     jpi, jpj, jpk =>', jpi, jpj, jpk
   WRITE(numout,*) ''

   IF ( (jpi/=ni1).OR.(jpj/=nj1) ) STOP 'Problem of shape between input field and mesh_mask!'


   jpiglo = jpi
   jpjglo = jpj



   !!LILI
   !! Need to know if there is at least one file with the vertical dimension in the input file!

   CALL check_nf90( NF90_OPEN( TRIM(cf_in), nf90_share, idf_src ) )
   CALL check_nf90( NF90_INQUIRE( idf_src, ndims, nvars, natts ) )

   WRITE(numout,*) ' idf_src, ndims, nvars, natts =', idf_src, ndims, nvars, natts

   CALL check_nf90( NF90_INQUIRE( idf_src, unlimitedDimId = id_t ) )

   WRITE(numout,*) ''
   WRITE(numout,'(" *** input file contains ",i1," dimensions of which time has id ",i1," !")') ndims, id_t


   ALLOCATE(indimlens(ndims), indimnames(ndims), outdimlens(ndims))

   l_3d = .FALSE.

   DO idim = 1, ndims
      CALL check_nf90( nf90_inquire_dimension( idf_src, idim, cnm_dim, dimlen ) )
      indimlens(idim) = dimlen
      indimnames(idim) = cnm_dim

      outdimlens(idim) = dimlen
      !! We're going to write a coarsened file so x & y dimension must be set accordingly!
      IF ( TRIM(cnm_dim) == 'x' ) THEN
         id_x = idim
         IF ( .NOT. (dimlen == jpi) ) STOP 'PROBLEM #1'
      ELSEIF ( TRIM(cnm_dim) == 'y' ) THEN
         id_y = idim
         IF ( .NOT. (dimlen == jpj) ) STOP 'PROBLEM #2'
      ELSEIF ( TRIM(cnm_dim) == 'deptht' ) THEN
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
      !! There is this annoying axbnd dimension:
      IF ( ndims == nb_normal+1 ) THEN
         DO idim = 1, ndims
            IF ( .NOT. ANY( (/id_x,id_y,id_z,id_t/) == idim ) ) THEN
               id_dim_axbnd = idim
               EXIT
            END IF
         END DO
         nlaxbnd = outdimlens(id_dim_axbnd)
      ELSE
         STOP 'PROBLEM: there is more than 1 axbnd weird dimension!!!'
      END IF
   END IF


   WRITE(numout,*) ''
   WRITE(numout,*) '  *** id_x =', id_x, jpi
   WRITE(numout,*) '  *** id_y =', id_y, jpj
   IF ( l_3d ) WRITE(numout,*) '  *** id_z =', id_z, jpk
   WRITE(numout,*) '  *** id_t =', id_t, Nt, ' time records!'
   WRITE(numout,*) '  *** id_dim_axbnd =', id_dim_axbnd, nlaxbnd
   WRITE(numout,*) '' ; WRITE(numout,*) ''





   jpkm1 = MAX(jpk-1,1)


   !ni = ni1 ; jpj = ni1
   !! Source:
   ALLOCATE ( x2d_r4(jpi,jpj) )
   IF (( cgp == 'T' ).OR.(l_write_crs_mm)) ALLOCATE ( e3t_0(jpi,jpj,jpk), glamt(jpi,jpj), gphit(jpi,jpj), e1e2t(jpi,jpj), e1t(jpi,jpj), e2t(jpi,jpj) )
   IF (( cgp == 'U' ).OR.(l_write_crs_mm)) ALLOCATE ( e3u_0(jpi,jpj,jpk), glamu(jpi,jpj), gphiu(jpi,jpj), e1u(jpi,jpj), e2u(jpi,jpj) )
   IF (( cgp == 'V' ).OR.(l_write_crs_mm)) ALLOCATE ( e3v_0(jpi,jpj,jpk), glamv(jpi,jpj), gphiv(jpi,jpj), e1v(jpi,jpj), e2v(jpi,jpj) )

   IF ( l_write_crs_mm ) ALLOCATE ( e3w_0(jpi,jpj,jpk), glamf(jpi,jpj), gphif(jpi,jpj), e1f(jpi,jpj), e2f(jpi,jpj) )

   IF ( l_3d )  ALLOCATE ( vdepth(jpk) , vdepth_b(nlaxbnd,jpk) , x3d_r4(jpi,jpj,jpk) ) !, e3t(jpi,jpj,jpk) )

   
   IF ( l_3d .AND. l_write_crs_mm ) CALL GETVAR_1D(cf_mm, 'nav_lev', vdepth)


   


   !! Getting source land-sea mask:
   WRITE(numout,*) ''
   ALLOCATE ( tmask(jpi,jpj,jpk), umask(jpi,jpj,jpk), vmask(jpi,jpj,jpk), fmask(jpi,jpj,jpk) )
   WRITE(numout,*) ' *** Reading land-sea mask'
   IF (( cgp == 'T' ).OR.(l_write_crs_mm)) CALL GETMASK_3D(cf_mm, 'tmask', tmask, jz1=1, jz2=jpk)
   IF (( cgp == 'U' ).OR.(l_write_crs_mm)) CALL GETMASK_3D(cf_mm, 'umask', umask, jz1=1, jz2=jpk)
   IF (( cgp == 'V' ).OR.(l_write_crs_mm)) CALL GETMASK_3D(cf_mm, 'vmask', vmask, jz1=1, jz2=jpk)
   IF ( l_write_crs_mm )                   CALL GETMASK_3D(cf_mm, 'fmask', fmask, jz1=1, jz2=jpk)
   WRITE(numout,*) ' Done!'; WRITE(numout,*) ''



   ALLOCATE ( Vt(Nt) )
   CALL GETVAR_1D(cf_in, cv_t, Vt)
   WRITE(numout,*) 'Vt = ', Vt(:)





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

   jphgr_msh = 0

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

   WRITE(numout,*) ' *** After crs_dom_def: jpi_crs, jpj_crs =', jpi_crs, jpj_crs
   WRITE(numout,*) 'TARGET coarsened horizontal domain, jpi_crs, jpj_crs =', jpi_crs, jpj_crs
   WRITE(numout,*) ' *** nn_factx, nn_facty', nn_factx, nn_facty
   WRITE(numout,*) ''


   IF (( cgp == 'T' ).OR.(l_write_crs_mm)) THEN
      CALL GETVAR_2D(i0, j0, cf_mm, 'glamt', 0, 0, 0, glamt) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'gphit', 0, 0, 0, gphit) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e1t', 0, 0, 0, e1t) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e2t', 0, 0, 0, e2t) ; i0=0 ; j0=0
      e1e2t(:,:) = e1t(:,:)*e2t(:,:)
      CALL GETVAR_3D(i0, j0, cf_mm, 'e3t_0', 0, 1,    e3t_0(:,:,:), jz1=1, jz2=jpk)
   END IF

   IF (( cgp == 'U' ).OR.(l_write_crs_mm)) THEN
      CALL GETVAR_2D(i0, j0, cf_mm, 'glamu', 0, 0, 0, glamu) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'gphiu', 0, 0, 0, gphiu) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e1u', 0, 0, 0, e1u) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e2u', 0, 0, 0, e2u) ; i0=0 ; j0=0
      CALL GETVAR_3D(i0, j0, cf_mm, 'e3u_0', 0, 1,    e3u_0(:,:,:), jz1=1, jz2=jpk)
   END IF

   IF (( cgp == 'V' ).OR.(l_write_crs_mm)) THEN
      CALL GETVAR_2D(i0, j0, cf_mm, 'glamv', 0, 0, 0, glamv) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'gphiv', 0, 0, 0, gphiv) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e1v', 0, 0, 0, e1v) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e2v', 0, 0, 0, e2v) ; i0=0 ; j0=0
      CALL GETVAR_3D(i0, j0, cf_mm, 'e3v_0', 0, 1,    e3v_0(:,:,:), jz1=1, jz2=jpk)
   END IF

   IF ( l_write_crs_mm ) THEN
      CALL GETVAR_2D(i0, j0, cf_mm, 'glamf', 0, 0, 0, glamf) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'gphif', 0, 0, 0, gphif) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e1f', 0, 0, 0, e1f) ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, 'e2f', 0, 0, 0, e2f) ; i0=0 ; j0=0
      CALL GETVAR_3D(i0, j0, cf_mm, 'e3w_0', 0, 1,    e3w_0(:,:,:), jz1=1, jz2=jpk)
   END IF


   !CALL DUMP_FIELD(REAL(tmask(:,:,1),4), 'tmask.tmp', 'tmask' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   !CALL DUMP_FIELD(REAL(glamt(:,:)  ,4), 'glamt.tmp', 'glamt' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   !CALL DUMP_FIELD(REAL(gphit(:,:)  ,4), 'gphit.tmp', 'gphit' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   !CALL DUMP_FIELD(REAL(e1t(:,:)  ,4), 'e1t.tmp', 'e1t' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
   !CALL DUMP_FIELD(REAL(e2t(:,:)  ,4), 'e2t.tmp', 'e2t' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)





   !---------------------------------------------------------
   ! 3. Mask and Mesh
   !---------------------------------------------------------

   !  3.a. Get the masks
   tmask_crs(:,:,:) = 0 ; umask_crs(:,:,:) = 0 ; vmask_crs(:,:,:) = 0 ; fmask_crs(:,:,:) = 0

   CALL crs_dom_msk

   IF ( l_debug ) THEN 
      IF (( cgp == 'T' ).OR.(l_write_crs_mm)) CALL DUMP_FIELD(REAL(tmask_crs,4), 'tmask_crs.tmp', 'tmask_crs' )
      IF (( cgp == 'U' ).OR.(l_write_crs_mm)) CALL DUMP_FIELD(REAL(umask_crs,4), 'umask_crs.tmp', 'umask_crs' )
      IF (( cgp == 'V' ).OR.(l_write_crs_mm)) CALL DUMP_FIELD(REAL(vmask_crs,4), 'vmask_crs.tmp', 'vmask_crs' )
      IF ( l_write_crs_mm )                   CALL DUMP_FIELD(REAL(fmask_crs,4), 'fmask_crs.tmp', 'fmask_crs' )
   END IF


   !  3.b. Get the coordinates
   !      Odd-numbered reduction factor, center coordinate on T-cell
   !      Even-numbered reduction factor, center coordinate on U-,V- faces or f-corner.
   !
   gphit_crs = 0.0 ; glamt_crs = 0.0 ; gphiu_crs = 0.0 ; glamu_crs = 0.0
   gphiv_crs = 0.0 ; glamv_crs = 0.0 ; gphif_crs = 0.0 ; glamf_crs = 0.0

   IF (( cgp == 'T' ).OR.(l_write_crs_mm)) CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs )
   IF (( cgp == 'U' ).OR.(l_write_crs_mm)) CALL crs_dom_coordinates( gphiu, glamu, 'U', gphiu_crs, glamu_crs )
   IF (( cgp == 'V' ).OR.(l_write_crs_mm)) CALL crs_dom_coordinates( gphiv, glamv, 'V', gphiv_crs, glamv_crs )
   IF ( l_write_crs_mm )                   CALL crs_dom_coordinates( gphif, glamf, 'F', gphif_crs, glamf_crs )



   !  3.c. Get the horizontal mesh

   !      3.c.1 Horizontal scale factors

   e1t_crs(:,:) = 0.0 ; e2t_crs(:,:) = 0.0 ; e1u_crs(:,:) = 0.0 ; e2u_crs(:,:) = 0.0
   e1v_crs(:,:) = 0.0 ; e2v_crs(:,:) = 0.0 ; e1f_crs(:,:) = 0.0 ; e2f_crs(:,:) = 0.0

   IF (( cgp == 'T' ).OR.(l_write_crs_mm)) CALL crs_dom_hgr( e1t, e2t, 'T', e1t_crs, e2t_crs )
   IF (( cgp == 'U' ).OR.(l_write_crs_mm)) CALL crs_dom_hgr( e1u, e2u, 'U', e1u_crs, e2u_crs )
   IF (( cgp == 'V' ).OR.(l_write_crs_mm)) CALL crs_dom_hgr( e1v, e2v, 'V', e1v_crs, e2v_crs )
   IF ( l_write_crs_mm )                   CALL crs_dom_hgr( e1f, e2f, 'F', e1f_crs, e2f_crs )

   WHERE(e1t_crs == 0._wp) e1t_crs=r_inf
   WHERE(e1u_crs == 0._wp) e1u_crs=r_inf
   WHERE(e1v_crs == 0._wp) e1v_crs=r_inf
   WHERE(e1f_crs == 0._wp) e1f_crs=r_inf
   WHERE(e2t_crs == 0._wp) e2t_crs=r_inf
   WHERE(e2u_crs == 0._wp) e2u_crs=r_inf
   WHERE(e2v_crs == 0._wp) e2v_crs=r_inf
   WHERE(e2f_crs == 0._wp) e2f_crs=r_inf

   e1e2t_crs(:,:) = e1t_crs(:,:) * e2t_crs(:,:)







   IF ( l_debug ) THEN
      IF (( cgp == 'T' ).OR.(l_write_crs_mm)) THEN
         CALL DUMP_FIELD(REAL(glamt_crs(:,:),4), 'glamt_crs.tmp', 'glamt_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(gphit_crs(:,:),4), 'gphit_crs.tmp', 'gphit_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(e1t_crs(:,:),4), 'e1t_crs.tmp', 'e1t_crs' )
         CALL DUMP_FIELD(REAL(e2t_crs(:,:),4), 'e2t_crs.tmp', 'e2t_crs' )
         CALL DUMP_FIELD(REAL(e1e2t_crs(:,:),4), 'e1e2t_crs.tmp', 'e1e2t_crs' )
      END IF
      IF (( cgp == 'U' ).OR.(l_write_crs_mm)) THEN
         CALL DUMP_FIELD(REAL(glamu_crs(:,:),4), 'glamu_crs.tmp', 'glamu_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(gphiu_crs(:,:),4), 'gphiu_crs.tmp', 'gphiu_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(e1u_crs(:,:),4), 'e1u_crs.tmp', 'e1u_crs' )
         CALL DUMP_FIELD(REAL(e2u_crs(:,:),4), 'e2u_crs.tmp', 'e2u_crs' )
      END IF
      IF (( cgp == 'V' ).OR.(l_write_crs_mm)) THEN
         CALL DUMP_FIELD(REAL(glamv_crs(:,:),4), 'glamv_crs.tmp', 'glamv_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(gphiv_crs(:,:),4), 'gphiv_crs.tmp', 'gphiv_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(e1v_crs(:,:),4), 'e1v_crs.tmp', 'e1v_crs' )
         CALL DUMP_FIELD(REAL(e2v_crs(:,:),4), 'e2v_crs.tmp', 'e2v_crs' )
      END IF
      IF ( l_write_crs_mm ) THEN
         CALL DUMP_FIELD(REAL(glamf_crs(:,:),4), 'glamf_crs.tmp', 'glamf_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(gphif_crs(:,:),4), 'gphif_crs.tmp', 'gphif_crs' ) !,  xlon, xlat, cv_lo, cv_la,  rfill)
         CALL DUMP_FIELD(REAL(e1f_crs(:,:),4), 'e1f_crs.tmp', 'e1f_crs' )
         CALL DUMP_FIELD(REAL(e2f_crs(:,:),4), 'e2f_crs.tmp', 'e2f_crs' )
      END IF
   END IF




   !      3.c.2 Coriolis factor
   IF ( l_write_crs_mm ) THEN
      SELECT CASE( jphgr_msh )   ! type of horizontal mesh

      CASE ( 0, 1, 4 )           ! mesh on the sphere
         ff_crs(:,:) = 2. * omega * SIN( rad * gphif_crs(:,:) )
      CASE DEFAULT
         IF(lwp)    WRITE(numout,*) 'nemo_coarsener.f90. crs_init. Only jphgr_msh = 0, 1 or 4 supported'
      END SELECT

      IF ( l_debug ) CALL DUMP_FIELD(REAL(ff_crs(:,:),4), 'ff_crs.tmp', 'ff_crs' )
   END IF


   !    3.d.1 mbathy ( vertical k-levels of bathymetry )
   IF ( l_write_crs_mm ) THEN
      IF ( l_3d ) THEN
         CALL crs_dom_bat
         IF ( l_debug ) CALL DUMP_FIELD(REAL(mbathy_crs(:,:),4), 'mbathy_crs.tmp', 'mbathy_crs' )
      END IF
   END IF


   IF ( cgp == 'T' ) THEN
      ALLOCATE (zfse3t(jpi,jpj,jpk))
      zfse3t(:,:,:) = e3t_0(:,:,:) !fse3t(:,:,:)
   END IF
   IF ( cgp == 'U' ) THEN
      ALLOCATE (zfse3u(jpi,jpj,jpk))
      zfse3u(:,:,:) = e3u_0(:,:,:)
   END IF
   IF ( cgp == 'V' ) THEN
      ALLOCATE (zfse3v(jpi,jpj,jpk))
      zfse3v(:,:,:) = e3v_0(:,:,:)
   END IF
   IF ( cgp == 'W' ) THEN
      ALLOCATE (zfse3w(jpi,jpj,jpk))
      zfse3w(:,:,:) = e3w_0(:,:,:)
   END IF

   

   !    3.d.2   Surfaces
   e2e3u_crs(:,:,:)=0._wp
   e2e3u_msk(:,:,:)=0._wp
   e1e3v_crs(:,:,:)=0._wp
   e1e3v_msk(:,:,:)=0._wp
   CALL crs_dom_sfc( REAL(tmask,8), 'W', e1e2w_crs, e1e2w_msk, p_e1=e1t, p_e2=e2t    )
   IF ( cgp == 'U' ) CALL crs_dom_sfc( REAL(umask,8), 'U', e2e3u_crs, e2e3u_msk, p_e2=e2u, p_e3=zfse3u )
   IF ( cgp == 'V' ) CALL crs_dom_sfc( REAL(vmask,8), 'V', e1e3v_crs, e1e3v_msk, p_e1=e1v, p_e3=zfse3v )

   IF ( l_debug ) THEN
      CALL DUMP_FIELD(REAL(e1e2w_crs,4), 'e1e2w_crs.tmp', 'e1e2w_crs' )
   END IF







   !    3.d.3   Vertical scale factors
   !
   IF ( cgp == 'T' ) CALL crs_dom_e3( e1t, e2t, zfse3t, p_sfc_3d_crs=e1e2w_crs, cd_type='T', p_mask=REAL(tmask,wp), p_e3_crs=e3t_0_crs, p_e3_max_crs=e3t_max_0_crs)
   IF ( cgp == 'W' ) CALL crs_dom_e3( e1t, e2t, zfse3w, p_sfc_3d_crs=e1e2w_crs, cd_type='W', p_mask=REAL(tmask,wp), p_e3_crs=e3w_0_crs, p_e3_max_crs=e3w_max_0_crs)
   IF ( cgp == 'U' ) CALL crs_dom_e3( e1u, e2u, zfse3u, p_sfc_2d_crs=e2u_crs  , cd_type='U', p_mask=REAL(umask,wp), p_e3_crs=e3u_0_crs, p_e3_max_crs=e3u_max_0_crs)
   IF ( cgp == 'V' ) CALL crs_dom_e3( e1v, e2v, zfse3v, p_sfc_2d_crs=e1v_crs  , cd_type='V', p_mask=REAL(vmask,wp), p_e3_crs=e3v_0_crs, p_e3_max_crs=e3v_max_0_crs)

   WHERE(e3t_max_0_crs == 0._wp) e3t_max_0_crs=r_inf
   WHERE(e3u_max_0_crs == 0._wp) e3u_max_0_crs=r_inf
   WHERE(e3v_max_0_crs == 0._wp) e3v_max_0_crs=r_inf
   WHERE(e3w_max_0_crs == 0._wp) e3w_max_0_crs=r_inf

   ht_0_crs(:,:)=0._wp
   DO jk = 1, jpk
      ht_0_crs(:,:)=ht_0_crs(:,:)+e3t_0_crs(:,:,jk)*tmask_crs(:,:,jk)
   ENDDO

   IF ( l_debug ) THEN
      IF ( cgp == 'T' ) CALL DUMP_FIELD(REAL(e3t_0_crs,4), 'e3t_0_crs.tmp', 'e3t_0_crs' )
      IF ( cgp == 'U' ) CALL DUMP_FIELD(REAL(e3u_0_crs,4), 'e3u_0_crs.tmp', 'e3u_0_crs' )
      IF ( cgp == 'V' ) CALL DUMP_FIELD(REAL(e3v_0_crs,4), 'e3v_0_crs.tmp', 'e3v_0_crs' )
      IF ( cgp == 'W' ) CALL DUMP_FIELD(REAL(e3w_0_crs,4), 'e3w_0_crs.tmp', 'e3w_0_crs' )
      IF ( cgp == 'T' ) CALL DUMP_FIELD(REAL(ht_0_crs,4),  'ht_0_crs.tmp' , 'ht_0_crs' )
   END IF



   !    3.d.3   Vertical depth (meters)
   !cbr: il semblerait que p_e3=... ne soit pas utile ici !!!!!!!!!
   !CALL crs_dom_ope( gdept_0, 'MAX', 'T', tmask, gdept_0_crs, p_e3=zfse3t, psgn=1.0 )
   !CALL crs_dom_ope( gdepw_0, 'MAX', 'W', tmask, gdepw_0_crs, p_e3=zfse3w, psgn=1.0 )
   !DO jk = 1, jpk
   !   DO ji = 1, jpi_crs
   !      DO jj = 1, jpj_crs
   !         IF( gdept_0_crs(ji,jj,jk) .LE. 0._wp ) gdept_0_crs(ji,jj,jk) = gdept_1d(jk)
   !         IF( gdepw_0_crs(ji,jj,jk) .LE. 0._wp ) gdepw_0_crs(ji,jj,jk) = gdepw_1d(jk)
   !      ENDDO
   !   ENDDO
   !ENDDO


   !!.... 3 to 6 skipped !!! => check crsini in NEMO !!!


   !---------------------------------------------------------
   ! 7. Finish and clean-up
   !---------------------------------------------------------

   !needed later... DEALLOCATE ( zfse3t, zfse3u, zfse3v, zfse3w )

   IF ( cgp == 'U' ) DEALLOCATE ( zfse3u )
   IF ( cgp == 'V' ) DEALLOCATE ( zfse3v )
   IF ( cgp == 'W' ) DEALLOCATE ( zfse3w )

   ! -------------------------- crs init done ! ----------------------------


   WRITE(numout,*) ''
   WRITE(numout,*) ''
   WRITE(numout,*) ''

   
   !! Writing the mesh_mask if required:

   IF ( l_write_crs_mm .AND. l_3d ) THEN !lili
      CALL check_nf90( NF90_CREATE( 'mesh_mask_crs.nc', nf90_netcdf4, idf_mm_crs, chunksize=chunksize ) )
      CALL check_nf90( NF90_DEF_DIM(idf_mm_crs, 'x', jpi_crs,        idd_x) )
      CALL check_nf90( NF90_DEF_DIM(idf_mm_crs, 'y', jpj_crs,        idd_y) )
      CALL check_nf90( NF90_DEF_DIM(idf_mm_crs, 'z', jpk,            idd_z) )
      CALL check_nf90( NF90_DEF_DIM(idf_mm_crs, 't', nf90_unlimited, idd_t) )
      !
      CALL check_nf90( NF90_DEF_VAR(idf_mm_crs, 'nav_lon', NF90_FLOAT, (/idd_x,idd_y/), idv_nlon, deflate_level=9) )
      CALL check_nf90( NF90_DEF_VAR(idf_mm_crs, 'nav_lat', NF90_FLOAT, (/idd_x,idd_y/), idv_nlat, deflate_level=9) )
      CALL check_nf90( NF90_DEF_VAR(idf_mm_crs, 'nav_lev', NF90_FLOAT, (/idd_z/),       idv_nlev, deflate_level=9) )
      !
      CALL check_nf90( NF90_DEF_VAR(idf_mm_crs, 'tmask', NF90_BYTE, (/idd_x,idd_y,idd_z,idd_t/), idv_tmsk, deflate_level=9) )
      CALL check_nf90( NF90_DEF_VAR(idf_mm_crs, 'umask', NF90_BYTE, (/idd_x,idd_y,idd_z,idd_t/), idv_umsk, deflate_level=9) )
      CALL check_nf90( NF90_DEF_VAR(idf_mm_crs, 'vmask', NF90_BYTE, (/idd_x,idd_y,idd_z,idd_t/), idv_vmsk, deflate_level=9) )
      CALL check_nf90( NF90_DEF_VAR(idf_mm_crs, 'fmask', NF90_BYTE, (/idd_x,idd_y,idd_z,idd_t/), idv_fmsk, deflate_level=9) )

      CALL check_nf90( NF90_ENDDEF( idf_mm_crs) )

      CALL check_nf90( NF90_PUT_VAR(idf_mm_crs, idv_nlon, REAL(glamt_crs,4)) )
      CALL check_nf90( NF90_PUT_VAR(idf_mm_crs, idv_nlat, REAL(gphit_crs,4)) )
      CALL check_nf90( NF90_PUT_VAR(idf_mm_crs, idv_nlev, vdepth           ) )
      CALL check_nf90( NF90_PUT_VAR(idf_mm_crs, idv_tmsk, tmask_crs, start=(/1,1,1,1/), count=(/jpi_crs,jpj_crs,jpk,1/)) )
      CALL check_nf90( NF90_PUT_VAR(idf_mm_crs, idv_umsk, umask_crs, start=(/1,1,1,1/), count=(/jpi_crs,jpj_crs,jpk,1/)) )
      CALL check_nf90( NF90_PUT_VAR(idf_mm_crs, idv_vmsk, vmask_crs, start=(/1,1,1,1/), count=(/jpi_crs,jpj_crs,jpk,1/)) )
      CALL check_nf90( NF90_PUT_VAR(idf_mm_crs, idv_fmsk, fmask_crs, start=(/1,1,1,1/), count=(/jpi_crs,jpj_crs,jpk,1/)) )
      CALL check_nf90( NF90_CLOSE( idf_mm_crs ) )

      STOP 'mesh_mask_crs.nc'

   END IF





   







   !2.1 Set up the output file
   CALL check_nf90( nf90_create( TRIM(cf_out), nf90_netcdf4, idf_trg, chunksize=chunksize ) )


   !! Defining dimention in output file:
   outdimlens(id_x) = jpi_crs ! => overwriting values taken from source domain !
   outdimlens(id_y) = jpj_crs ! =>         "

   DO idim = 1, ndims
      CALL check_nf90( nf90_inquire_dimension( idf_src, idim, cnm_dim, dimlen ) )
      IF( idim == id_t ) THEN
         CALL check_nf90( nf90_def_dim( idf_trg, cnm_dim, nf90_unlimited, dimid) )
         nmax_unlimited = dimlen
      ELSE
         WRITE(numout,*) ' nf90_def_dim( idf_trg, cnm_dim, outdimlens(idim), dimid)'
         CALL check_nf90( nf90_def_dim( idf_trg, cnm_dim, outdimlens(idim), dimid) )
      ENDIF
   END DO

   ! nmax_unlimited is only used for time-chunking so we set it to be at least 1 to
   ! account for files with no record dimension or zero length record dimension(!)
   nmax_unlimited = max(nmax_unlimited,1)
   WRITE(numout,*) ''


   IF ( ndims > id_t ) THEN
      IF ( ndims == id_t+1 ) THEN
         id_dim_axbnd = ndims
         nlaxbnd = outdimlens(ndims)
      ELSE
         STOP 'PROBLEM: there is more than 1 axbnd weird dimension!!!'
      END IF
   END IF






   !2.2.3 Copy the variable definitions and attributes into the output file.
   ALLOCATE( mdiVals(nvars), c_list_var_names(nvars), i_list_var_types(nvars), &
      &      i_list_var_ndims(nvars), i_list_var_ids(nvars), i_list_var_dim_ids(4,nvars), &
      &      l_var_is_done(nvars) )

   ALLOCATE ( xr8(nlaxbnd) )





   i_list_var_dim_ids(:,:) = 0
   mdiVals(:)=0

   DO jv = 1, nvars

      id_v = jv ! lolo!
      indimids(:) = 0

      CALL check_nf90( NF90_INQUIRE_VARIABLE( idf_src, id_v, cnm_var, ntype, ndims, indimids, natts ) )
      c_list_var_names(jv)     = TRIM(cnm_var)
      i_list_var_types(jv)     = ntype
      i_list_var_ndims(jv)     = ndims
      i_list_var_dim_ids(:,jv) = indimids(:)
      !CALL check_nf90( nf90_inquire_variable( idf_src, id_v) !, ndims=nd)
      i_list_var_ids(jv)   = id_v

      ALLOCATE(outdimids(ndims))
      DO idim = 1, ndims
         outdimids(idim) = indimids(idim)
      END DO



      CALL check_nf90( nf90_def_var( idf_trg, cnm_var, ntype, outdimids, varid, &
         deflate_level=deflate_level ) )
      DEALLOCATE(outdimids)
      WRITE(numout,*) 'Defining variable '//TRIM(cnm_var)//'...'
      IF( natts > 0 ) THEN
         DO attid = 1, natts
            CALL check_nf90( NF90_INQ_ATTNAME( idf_src, varid, attid, cnm_att ) )
            IF ( cnm_att == "_FillValue" ) THEN
               CALL check_nf90( NF90_GET_ATT( idf_src, varid, cnm_att, rmdi ) )
               mdiVals(jv)=rmdi
            ENDIF
            CALL check_nf90( NF90_COPY_ATT( idf_src, varid, cnm_att, idf_trg, varid ) )
         END DO
      ENDIF
   END DO



   !2.3 End definitions in output file and copy 1st file idf_src to the inidf_srcs array

   CALL check_nf90( NF90_ENDDEF( idf_trg ) )
   !inidf_srcs(1) = idf_src
   !WRITE(numout,*) 'Finished defining output file.'

   ALLOCATE ( x2d_r4_crs(jpi_crs,jpj_crs), x2d_r8_crs(jpi_crs,jpj_crs) )
   IF ( l_3d ) ALLOCATE ( x3d_r4_crs(jpi_crs,jpj_crs,jpk), x3d_r8_crs(jpi_crs,jpj_crs,jpk), e3t_max_crs(jpi_crs,jpj_crs,jpk) )





   WRITE(numout,*) ''

   l_var_is_done(:) = .FALSE.

   ! A. Writing all variables that do not have a time record:

   DO jv = 1, nvars

      cv_in = TRIM(c_list_var_names(jv))
      id_v        = i_list_var_ids(jv)
      itype       = i_list_var_types(jv)
      nbdim       = i_list_var_ndims(jv)
      indimids(:) = i_list_var_dim_ids(:,jv)

      IF ( .NOT. ANY(indimids == id_t) ) THEN
         WRITE(numout,*) '  *** Variable '//TRIM(cv_in)//': no time record, writing before time loop!'

         ! LOLO add 1D like depth here !!!! IF( nbdim == 1 ) THEN

         IF( nbdim == 1 ) THEN
            IF ( (TRIM(cv_in)=='deptht').OR.(TRIM(cv_in)=='depthw') ) THEN
               CALL check_nf90( nf90_get_var( idf_src,  id_v, vdepth ) )
               CALL check_nf90( nf90_put_var( idf_trg, id_v, vdepth ) )
               l_var_is_done(jv) = .TRUE.
            ELSE
               WRITE(numout,*) 'ERROR: unknown 1D variable without time!!!'; STOP
            END IF

         ELSEIF( nbdim == 2 ) THEN
            IF ( (TRIM(cv_in)=='nav_lon').OR.(TRIM(cv_in)=='glamt') ) THEN
               ! It's beed coarsened earlier! with crs_dom_coordinates => !CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs )
               WRITE(numout,*) '      ==> writing coarsened '//TRIM(cv_in)//' in '//TRIM(cf_out)
               CALL check_nf90( nf90_put_var( idf_trg, id_v, glamt_crs ) )
               l_var_is_done(jv) = .TRUE.
            ELSEIF ( (TRIM(cv_in)=='nav_lat').OR.(TRIM(cv_in)=='gphit') ) THEN
               ! It's beed coarsened earlier! with crs_dom_coordinates => !CALL crs_dom_coordinates( gphit, glamt, 'T', gphit_crs, glamt_crs )
               WRITE(numout,*) '      ==> writing coarsened '//TRIM(cv_in)//' in '//TRIM(cf_out)
               CALL check_nf90( nf90_put_var( idf_trg, id_v, gphit_crs ) )
               l_var_is_done(jv) = .TRUE.

            ELSEIF ( (TRIM(cv_in)=='deptht_bounds').OR.(TRIM(cv_in)=='depthw_bounds') ) THEN
               CALL check_nf90( nf90_get_var( idf_src,  id_v, vdepth_b ) )
               CALL check_nf90( nf90_put_var( idf_trg, id_v, vdepth_b ) )
               l_var_is_done(jv) = .TRUE.

            ELSE
               WRITE(numout,*) 'ERROR: unknown 2D variable without time!!!'; STOP
            END IF !IF ( (TRIM(cv_in)=='nav_lon').OR.(TRIM(cv_in)=='glamt') )

         ELSE

            WRITE(numout,*) 'ERROR: unknown number of dimensions for variable without time)!!!'; STOP

         END IF ! IF( nbdim == 1 )

         WRITE(numout,*) ''
      END IF ! IF ( .NOT. ANY(indimids == id_t) )

   END DO  ! DO jv = 1, nvars





   !! B. Will read, coarsen (if / x,y), and write, everything that depends on time record:

   !e3t_0(:,:,:) = 1._wp

   DO jt=1, Nt

      WRITE(numout,*) ''

      DO jv = 1, nvars
         IF ( .NOT. l_var_is_done(jv) ) THEN

            WRITE(numout,*) ''
            cv_in = TRIM(c_list_var_names(jv))
            id_v        = i_list_var_ids(jv)
            itype       = i_list_var_types(jv)
            nbdim       = i_list_var_ndims(jv)
            indimids(:) = i_list_var_dim_ids(:,jv)

            WRITE(numout,*) ' ### Reading field '//TRIM(cv_in)//' at record #',jt
            WRITE(numout,*) ' * var ID         =>', id_v
            WRITE(numout,*) ' * var type       =>', itype
            WRITE(numout,*) ' * var nb of dims =>', nbdim
            WRITE(numout,*) ' * dim ids        =>', indimids(:)
            WRITE(numout,*) ''

            ! T
            IF( nbdim == 1 ) THEN
               IF ( .NOT. (indimids(1)==id_t) ) STOP 'ERROR: we should only have record-dependant vectors here!!!'
               SELECT CASE( itype )
               CASE( NF90_DOUBLE )
                  CALL check_nf90( nf90_get_var( idf_src,  id_v, r8, start=(/jt/), count=(/1/)) )
                  CALL check_nf90( nf90_put_var( idf_trg, id_v, r8, start=(/jt/), count=(/1/)) )
               CASE( NF90_FLOAT )
                  CALL check_nf90( nf90_get_var( idf_src,  id_v, r4, start=(/jt/), count=(/1/)) )
                  CALL check_nf90( nf90_put_var( idf_trg, id_v, r4, start=(/jt/), count=(/1/)) )
               CASE DEFAULT
                  WRITE(numout,*) 'ERROR: unknown netcdf type!!!'; STOP
               END SELECT
               WRITE(numout,*) '   *** '//TRIM(cv_in)//' written at record #',jt
               WRITE(numout,*) ''
            END IF


            ! 1D+T
            IF( nbdim == 2 ) THEN
               IF ( .NOT. ((indimids(1)==id_dim_axbnd).AND.(indimids(2)==id_t)) ) STOP 'ERROR: we do not know wthat this 2D field is! '
               SELECT CASE( itype )
               CASE( NF90_DOUBLE )
                  CALL check_nf90( NF90_GET_VAR(idf_src,  id_v, xr8, start=(/1,jt/), count=(/nlaxbnd,1/)) )
                  CALL check_nf90( NF90_PUT_VAR(idf_trg, id_v, xr8, start=(/1,jt/), count=(/nlaxbnd,1/)) )

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
                  CALL check_nf90( NF90_GET_VAR(idf_src, id_v, x2d_r4, start=(/1,1,jt/), count=(/jpi,jpj,1/)) )
               CASE DEFAULT
                  WRITE(numout,*) 'ERROR: unknown netcdf type!!! We should not have 2D+T arrays as double precision!'; STOP
               END SELECT
               WRITE(numout,*) '   *** just read '//TRIM(cv_in)//' at record #',jt

               ! Time for coarsening! Use different routines depending on the type of field!
               IF ( (TRIM(cv_in)=='sossheig').OR.(TRIM(cv_in)=='sosstsst').OR.(TRIM(cv_in)=='sosaline') ) THEN
                  WRITE(numout,*) '   *** Coarsening '//TRIM(cv_in)//' with crs_dom_ope / VOL / T !'
                  CALL crs_dom_ope( REAL(x2d_r4,8) , 'VOL', 'T', REAL(tmask,8), x2d_r8_crs , p_e12=e1e2t, p_e3=zfse3t, psgn=1.0_wp )
                  WRITE(numout,*) '   *** writing '//TRIM(cv_in)//' in '//TRIM(cf_out)
                  x2d_r4_crs = x2d_r8_crs
                  WHERE ( tmask_crs(:,:,1) == 0 ) x2d_r4_crs = 1.e+20
                  CALL check_nf90( nf90_put_var( idf_trg, id_v,  x2d_r4_crs, start=(/1,1,jt/), count=(/jpi_crs,jpj_crs,1/)) )
                  WRITE(numout,*) '   ***  written!'
               END IF

            END IF !IF( nbdim == 3 )



            IF( nbdim == 4 ) THEN

               PRINT *, ' indimids =', indimids(:)


               IF ( .NOT. ((indimids(1)==id_x).AND.(indimids(2)==id_y).AND.(indimids(3)==id_z).AND.(indimids(4)==id_t)) ) STOP 'ERROR: we do not know wthat this 4D field is! '
               SELECT CASE( itype )
               CASE( NF90_FLOAT )
                  CALL check_nf90( NF90_GET_VAR(idf_src, id_v, x3d_r4, start=(/1,1,1,jt/), count=(/jpi,jpj,jpk,1/)) )
               CASE DEFAULT
                  WRITE(numout,*) 'ERROR: unknown netcdf type!!! We should not have 3D+T arrays as double precision!'; STOP
               END SELECT
               WRITE(numout,*) '   *** just read '//TRIM(cv_in)//' at record #',jt

               ! Time for coarsening! Use different routines depending on the type of field!
               IF ( (TRIM(cv_in)=='e3t') ) THEN
                  !! Should come prior to any other 3D field because we need actual e3t!? (VVL)
                  WRITE(numout,*) '   *** zfse3t is updated with field '//TRIM(cv_in)//' !!!'
                  zfse3t(:,:,:) = x3d_r4(:,:,:)
                  WRITE(numout,*) '   *** Coarsening '//TRIM(cv_in)//' with crs_dom_e3 / T !'
                  CALL crs_dom_e3( e1t, e2t, zfse3t, p_sfc_3d_crs=e1e2w_crs, cd_type='T', p_mask=REAL(tmask,8), p_e3_crs=x3d_r8_crs, p_e3_max_crs=e3t_max_crs )
                  x3d_r4_crs = REAL(x3d_r8_crs,4)
                  WHERE ( tmask_crs == 0 ) x3d_r4_crs = 1.e+20
                  CALL check_nf90( nf90_put_var( idf_trg, id_v,  x3d_r4_crs, start=(/1,1,1,jt/), count=(/jpi_crs,jpj_crs,jpk,1/)) )
                  WRITE(numout,*) '   ***  written!'
               END IF

               IF ( (TRIM(cv_in)=='votemper').OR.(TRIM(cv_in)=='vosaline') ) THEN
                  WRITE(numout,*) '   *** Coarsening '//TRIM(cv_in)//' with crs_dom_ope / VOL / T !'
                  CALL crs_dom_ope( REAL(x3d_r4,8) , 'VOL', 'T', REAL(tmask,8), x3d_r8_crs , p_e12=e1e2t, p_e3=zfse3t, psgn=1.0_wp )
                  WRITE(numout,*) '   *** writing '//TRIM(cv_in)//' in '//TRIM(cf_out)
                  x3d_r4_crs = REAL(x3d_r8_crs,4)
                  WHERE ( tmask_crs == 0 ) x3d_r4_crs = 1.e+20
                  CALL check_nf90( nf90_put_var( idf_trg, id_v,  x3d_r4_crs, start=(/1,1,1,jt/), count=(/jpi_crs,jpj_crs,jpk,1/)) )
                  WRITE(numout,*) '   ***  written!'
               END IF

            END IF !IF( nbdim == 4 )




            !WRITE(numout,*) 'ERROR: unknown number of dimmensions!!!'





            !CALL check_nf90( NF90_GET_VAR(idf_src, id_v, x2d_r4_crs, start=(/1,1/), count=(/lx,ly/))


         END IF ! IF ( .NOT. l_var_is_done(jv) )


      END DO !DO jv = 1, nvars





   END DO !DO jt=1, Nt


   CALL check_nf90( nf90_close( idf_trg ) )
   CALL check_nf90( nf90_close( idf_src ) )




   WRITE(numout,*) ''
   WRITE(numout,*) 'Over!'
   WRITE(numout,*) ''




CONTAINS






   SUBROUTINE GET_MY_ARG(cname, cvalue)
      CHARACTER(len=*), INTENT(in)    :: cname
      CHARACTER(len=*), INTENT(inout) :: cvalue
      !!
      IF ( jarg + 1 > iargc() ) THEN
         WRITE(numout,*) 'ERROR: Missing ',trim(cname),' name!' ; call usage()
      ELSE
         jarg = jarg + 1 ;  CALL getarg(jarg,cr)
         IF ( ANY(clist_opt == trim(cr)) ) THEN
            WRITE(numout,*) 'ERROR: Missing',trim(cname),' name!'; call usage()
         ELSE
            cvalue = trim(cr)
         END IF
      END IF
   END SUBROUTINE GET_MY_ARG





   SUBROUTINE usage()
      !!
      WRITE(numout,*) ''
      WRITE(numout,*) '   List of command line options:'
      WRITE(numout,*) '   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      WRITE(numout,*) ''
      WRITE(numout,*) ' -m <mesh_mask.nc>    => file containing grid metrics of model'
      WRITE(numout,*) ''
      WRITE(numout,*) ' -i <input_file.nc>   => input file to coarsen'
      WRITE(numout,*) ''
      WRITE(numout,*) ' -P <T/U/V/W>         => C-grid point on which fields in input file are given'
      WRITE(numout,*) ''
      WRITE(numout,*) ' -o <output_file.nc>  => file to be created'
      WRITE(numout,*) ''
      WRITE(numout,*) '    Optional:'
      WRITE(numout,*) ' -h                   => Show this message'
      WRITE(numout,*) ''
      WRITE(numout,*) ''
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
