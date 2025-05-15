PROGRAM CORR_VECT

   USE mod_conf
   USE mod_init
   USE io_ezcdf
   USE mod_grids,  ONLY: IS_ORCA_NORTH_FOLD
   USE mod_nemotools
   USE mod_drown
   
   !!========================================================================
   !! Purpose :  correct vector components 'uraw' and 'vraw' directly
   !!            interpolated on an irregular grid
   !! ---------  Only for ORCA family of grids!!!!
   !!
   !! Author :   Laurent Brodeau
   !! --------
   !!
   !!========================================================================

   IMPLICIT NONE

   LOGICAL, PARAMETER :: ldebug = .false.

   !! Grid :
   CHARACTER(len=80), PARAMETER   :: &
      &    cv_glamt     = 'glamt',   &   ! input grid longitude name, T-points
      &    cv_gphit     = 'gphit',   &   ! input grid latitude name,  T-points
      &    cv_glamu     = 'glamu',   &   ! input grid longitude name, U-points
      &    cv_gphiu     = 'gphiu',   &   ! input grid latitude name,  U-points
      &    cv_glamv     = 'glamv',   &   ! input grid longitude name, V-points
      &    cv_gphiv     = 'gphiv',   &   ! input grid latitude name,  V-points
      &    cv_glamf     = 'glamf',   &   ! input grid longitude name, F-points
      &    cv_gphif     = 'gphif',   &   ! input grid latitude name,  F-points
      &    cv_depth     = 'deptht'       !  depth at T-points (U-points and V-points too)


   CHARACTER(LEN=400)  :: cn_xtr_x, cn_xtr_y, cextra_x, cextra_y

   CHARACTER(len=3)    :: cdum
   CHARACTER(len=3)    :: cgrid_trg='T'
   CHARACTER(len=80)   :: cv_time_0 = 'none'
   CHARACTER(len=800)  :: cr, cf_mm, cf_ang, cnmlst_x, cnmlst_y
   CHARACTER(len=64)   :: cNP_lat   = '89.7'

   CHARACTER(len=80) :: &
      &    cv_u_in='0', cv_v_in='0', &
      &    cv_out_U    = 'uraw',     &   ! raw U name
      &    cv_out_V    = 'vraw'          ! raw V name

   CHARACTER(len=800)  :: &
      &    cf_raw_U,    &  ! file containing u_raw
      &    cf_raw_V,    &  ! file containing v_raw
      &    cf_out_U, cf_out_V,    &
      &    cufilin = 'none',  cvfilin = 'none'

   CHARACTER(len=80)  ::  &
      &    cv_rot_U ,  &  ! output name for U corrected
      &    cv_rot_V       ! output name for V corrected


   TYPE(grid_type) :: gt_orca

   INTEGER      :: &
      &    jarg, i3d, nbc, &
      &    iorca=0,      &
      &    nlext=3, &
      &    i0, j0, &
      &    ni, nj, nk, nk1, nk2, &
      &    ni1, nj1, Ntr1,       &
      &    ni2, nj2, Ntr2,       &
      &    ni_g, nj_g, nk_g,    &
      &    iargc,         &
      &    idf_u, idv_u, idf_v, idv_v, &
      &    id_f1, id_v1, &
      &    id_f2, id_v2

   INTEGER :: ipb, ipc, jiNP0, jjNP1, jjNP2

   INTEGER(1), DIMENSION(:,:),   ALLOCATABLE :: mask_np
   INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mask_t, mask_u, mask_v

   REAL(4), DIMENSION(:,:), ALLOCATABLE :: ztmp4


   REAL(4), DIMENSION(:,:,:), ALLOCATABLE :: U_c, V_c

   REAL(8), DIMENSION(:,:), ALLOCATABLE ::      &
      &    XCOS_t, XSIN_t, XCOS_u, XSIN_u, XCOS_v, XSIN_v, XCOS_f, XSIN_f, &
      &    zU, zV,  ztmp, Xdum, Ydum, &
      &    xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f

   REAL(8), DIMENSION(:), ALLOCATABLE ::   vtime, vdepth
   REAL(8) :: rNp_lat=90.
   REAL    :: rmv_eff

   INTEGER :: jt, jk

   LOGICAL :: &
      &    l_read_angles = .FALSE., &
      &    l_inv    = .FALSE., &
      &    l_anlt   = .FALSE., & !: analytic input field (U=1, V=0) DEBUG...
      &    l_3d_inv = .FALSE., & !: will treat 3d files in inverse mode...
      &    lmout_x, lmout_y,   &
      &    lexist,             & 
      &    lNPfix = .FALSE.

   CHARACTER(LEN=2), DIMENSION(12), PARAMETER :: &
      &            clist_opt = (/ '-I','-h','-m','-A','-G','-p','-f','-i','-v','-t','-1','-P' /)

   WRITE(6,*)''
   WRITE(6,*)'=========================================================='
   WRITE(6,*)'            S  O  S  I  E    version ', trim(csosie_version)
   WRITE(6,*)''
   WRITE(6,*)'            Vector rotation for distorted mapping       '
   WRITE(6,*)'=========================================================='
   WRITE(6,*)''


   !! Getting string arguments :
   !! --------------------------

   jarg = 0

   DO WHILE ( jarg < iargc() )
      !!
      jarg = jarg + 1
      CALL getarg(jarg,cr)
      !!
      !!
      SELECT CASE (trim(cr))
         !!
         !!
      CASE('-h')
         call usage_corr_vect()
         !!
         !!
      CASE('-I')
         l_inv = .TRUE.
         !!
      CASE('-f')
         IF( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
            WRITE(6,*)'ERROR: Missing namelist name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: ', trim(cr), ' is definitively not the name of the namelist!'
               call usage_corr_vect()
            ELSE
               cf_nml_sosie = trim(cr)
            END IF
         END IF
         !!
      CASE('-m')
         IF( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
            WRITE(6,*)'ERROR: Missing mesh_mask file name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: ', trim(cr), ' is definitively not the name of the mesh_mask file!'
               call usage_corr_vect()
            ELSE
               cf_mm = trim(cr)
            END IF
         END IF
         !!
      CASE('-A')
         IF( jarg + 1 > iargc() ) THEN ! checking that there is at least an other argument following
            WRITE(6,*)'ERROR: Missing angle file name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: ', trim(cr), ' is definitively not the name of the angle file!'
               call usage_corr_vect()
            ELSE
               cf_ang = TRIM(cr)
               l_read_angles = .TRUE.
            END IF
         END IF
         !!
      CASE('-G')
         IF( jarg + 1 > iargc() ) THEN
            WRITE(6,*)'ERROR: Missing grid type to write to ("T" or "UV"?) !!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == TRIM(cr)) ) THEN
               WRITE(6,*)'ERROR: Missing grid type to write to ("T" or "UV"?) !!'; CALL usage_corr_vect()
            ELSE
               cgrid_trg = TRIM(cr)
            END IF
         END IF
         !!
      CASE('-i')
         IF( jarg + 2 > iargc() ) THEN
            WRITE(6,*)'ERROR: Missing input file names!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: Missing input file name 1 !' ; call usage_corr_vect()
            ELSE
               cufilin = trim(cr)
            END IF
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: Missing input file name 2 !' ; call usage_corr_vect()
            ELSE
               cvfilin = trim(cr)
            END IF
         END IF
         !!
         !!
      CASE('-v')
         IF( jarg + 2 > iargc() ) THEN
            WRITE(6,*)'ERROR: Missing input file names!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: Missing input variable name 1 !' ; call usage_corr_vect()
            ELSE
               cv_u_in = trim(cr)
            END IF
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: Missing input variable name 2 !' ; call usage_corr_vect()
            ELSE
               cv_v_in = trim(cr)
            END IF
         END IF
         !!
      CASE('-t')
         IF( jarg + 1 > iargc() ) THEN
            WRITE(6,*)'ERROR: Missing time name!' ; call usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == trim(cr)) ) THEN
               WRITE(6,*)'ERROR: Missing time name!'; call usage_corr_vect()
            ELSE
               cv_time_0 = trim(cr)
            END IF
         END IF
         !!
      CASE('-1')
         l_anlt = .TRUE.
         !!
      CASE('-P')
         IF( jarg + 1 > iargc() ) THEN
            WRITE(6,*)'ERROR: Missing threshold latitude for NorthPole fix => -P <lat_lim_deg> !!' ; CALL usage_corr_vect()
         ELSE
            jarg = jarg + 1 ;  CALL getarg(jarg,cr)
            IF( ANY(clist_opt == TRIM(cr)) ) THEN
               WRITE(6,*)'ERROR: Missing threshold latitude for NorthPole fix => -P <lat_lim_deg> !!' ; CALL usage_corr_vect()
            ELSE
               lNPfix = .TRUE.
               cNP_lat  = TRIM(cr)
            END IF
         END IF
         !!
      CASE DEFAULT
         WRITE(6,*)'Unknown option: ', TRIM(cr) ; WRITE(6,*)''
         CALL usage_corr_vect()
         !!
      END SELECT
      !!
   END DO


   IF( (.NOT. l_inv).AND.( TRIM(cufilin) /= 'none' ) ) THEN
      WRITE(6,*)'ERROR: specify the "-I" switch if you want to perform inverse correction!'
      STOP
   END IF

   IF( (cgrid_trg /= 'U,V').AND.(TRIM(cgrid_trg) /= 'T') ) THEN
      WRITE(6,*)'ERROR: unknown target grid point type: '//trim(cgrid_trg)//'!'
      CALL usage_corr_vect()
   END IF

   WRITE(6,*)''; WRITE(6,*)'Use "-h" for help'; WRITE(6,*)''
   WRITE(6,*)''


   IF( l_inv ) THEN

      IF( (TRIM(cv_u_in)=='0').OR.(TRIM(cv_u_in)=='0') ) THEN
         WRITE(6,*)'ERROR: you must specify the names of vector components with the "-v" switch!'
         STOP
      END IF
      WRITE(6,*)' * Vector files to unrotate = ', trim(cufilin), ' , ', trim(cvfilin)
      WRITE(6,*)'   => associated variable names = ', TRIM(cv_u_in), ' , ', TRIM(cv_v_in)
      IF( trim(cv_time_0) == 'none' ) THEN
         WRITE(6,*)'ERROR: you must specify the name of time variable with the "-t" switch!'; STOP
      END IF
      WRITE(6,*)'   => time variable name = ', trim(cv_time_0)
   END IF

   WRITE(6,*)' * mesh_mask file to use = ', TRIM(cf_mm)





   IF(.NOT. l_inv) THEN

      !!     N O R M A L   C O R R E C T I O N

      cnmlst_x = TRIM(cf_nml_sosie)//'_x'
      cnmlst_y = TRIM(cf_nml_sosie)//'_y'

      WRITE(6,*)' * namelists we expect => ', TRIM(cnmlst_x)//' and '//TRIM(cnmlst_y)
      WRITE(6,*)''

      INQUIRE(FILE=trim(cf_mm), EXIST=lexist )
      IF( .NOT. lexist ) THEN
         WRITE(*,'("The mesh_mask file ",a," file was not found!")') trim(cf_mm)
         CALL usage_corr_vect()
      END IF


      IF( l_read_angles ) THEN
         INQUIRE(FILE=TRIM(cf_ang), EXIST=lexist )
         IF( .NOT. lexist ) THEN
            WRITE(*,'("The mesh_mask file ",a," file was not found!")') TRIM(cf_ang)
            CALL usage_corr_vect()
         END IF
      END IF


      !! Namelist of X component:
      INQUIRE(FILE=TRIM(cnmlst_x), EXIST=lexist )
      IF( .NOT. lexist ) THEN
         WRITE(*,'("The namelist file ",a," file was not found!")') TRIM(cnmlst_x)
         CALL usage_corr_vect()
      END IF
      WRITE(6,*)''
      cf_nml_sosie = TRIM(cnmlst_x)
      CALL READ_NMLST(2)
      lmout_x  = lmout
      cv_rot_U = cv_out
      cn_xtr_x = cextra


      !! Namelist of Y component:
      INQUIRE(FILE=TRIM(cnmlst_y), EXIST=lexist )
      IF( .NOT. lexist ) THEN
         WRITE(*,'("The namelist file ",a," file was not found!")') TRIM(cnmlst_y)
         CALL usage_corr_vect()
      END IF
      WRITE(6,*)''
      cf_nml_sosie = TRIM(cnmlst_y)
      CALL READ_NMLST(2)
      lmout_y  = lmout
      cv_rot_V = cv_out
      cn_xtr_y = cextra


      IF( trim(cgrid_trg) == 'T' ) THEN
         WRITE(6,*)' *** Will save on T-grid points.'
         cextra_x = 'gridT_'//TRIM(cn_xtr_x)
         cextra_y = 'gridT_'//TRIM(cn_xtr_y)
      ELSEIF( cgrid_trg == 'U,V' ) THEN
         WRITE(6,*)' *** Will save on U,V-grid points.'
         cextra_x = 'gridU_'//TRIM(cn_xtr_x)
         cextra_y = 'gridV_'//TRIM(cn_xtr_y)
      ELSE
         WRITE(6,*)'ERROR: "cgrid_trg" value unknown: ', TRIM(cgrid_trg) ; STOP
      END IF
      WRITE(6,*)''

      IF( l_reg_trg ) THEN
         WRITE(6,*)'Vector correction only makes sense if your target grid is distorded!'
         WRITE(6,*)'  => check "l_reg_trg" into the namelist...' ; WRITE(6,*)''; STOP
      END IF


      cf_raw_U = trim(cd_out)//'/'//trim(cv_out_U)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cn_xtr_x)//'.nc'
      cf_raw_V = trim(cd_out)//'/'//trim(cv_out_V)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cn_xtr_y)//'.nc'

      WRITE(6,*)'The two input pre-interpolated needed files are :'
      PRINT *, trim(cf_raw_U) ;   PRINT *, trim(cf_raw_V) ;

      cf_out_U = TRIM(cd_out)//'/'//TRIM(cv_rot_U)//'_'//TRIM(csource)//'-' &
         &   //TRIM(ctarget)//'_'//TRIM(cextra_x)//'.nc'
      cf_out_V = TRIM(cd_out)//'/'//TRIM(cv_rot_V)//'_'//TRIM(csource)//'-' &
         &   //TRIM(ctarget)//'_'//TRIM(cextra_y)//'.nc'

      WRITE(6,*)'' ;   WRITE(6,*)'output files :'
      PRINT *, trim(cf_out_U) ;   PRINT *, trim(cf_out_V) ; WRITE(6,*)'' ; WRITE(6,*)''

      WRITE(6,*)'File containing x and y raw components of vector to be treated :'
      PRINT *, trim(cf_raw_U)
      PRINT *, trim(cf_raw_V) ; WRITE(6,*)''
      WRITE(6,*)'File containing grid :'
      PRINT *, TRIM(cf_mm) ; WRITE(6,*)''
      IF( l_read_angles ) THEN
         WRITE(6,*)'File containing angles :'
         PRINT *, TRIM(cf_ang) ; WRITE(6,*)''
      END IF

      !! North Pole fix?
      IF(lNPfix) THEN
         READ(cNP_lat,'(f)') rNP_lat
         PRINT *, ' *** We are going to correct the North-Pole region from:',REAL(rNP_lat,4),'deg.N!'
      END IF
      !   !
      !   ipc = INDEX(TRIM(cNP_fckp),',') !,BACK=.TRUE.)
      !   cr = cNP_fckp(1:ipc-1)
         !   READ(cr,'(i)') jiNP0
      !   ipb = ipc
      !   !
      !   ipc = INDEX(TRIM(cNP_fckp(ipb+1:)),',') + ipb
      !   cr = cNP_fckp(ipb+1:ipc-1)
      !   READ(cr,'(i)') jjNP1
      !   !
      !   cr = TRIM(cNP_fckp(ipc+1:))
      !   READ(cr,'(i)') jjNP2
      !   !
      !   PRINT *, ' ==> ji0, jj1, jj2 =', jiNP0, jjNP1, jjNP2!      !
      !IF( jjNP2 < jjNP1 ) THEN
      !      WRITE(6,*)'ERROR: (North Pole fix) => jjNP2 must be >= jjNP1 !!!'
      !      STOP
      !   END IF
      !END IF
      
      !! Geting array dimension and testing...
      !! -------------------------------------

      CALL DIMS(cf_raw_U, cv_out_U, ni1, nj1, nk1, Ntr1)
      CALL DIMS(cf_raw_V, cv_out_V, ni2, nj2, nk2, Ntr2)

      CALL DIMS(cf_mm, cv_glamt, ni_g, nj_g, nk_g, Ntr)

      !! testing ni agreement :
      IF( (ni1 /= ni2).or.(ni1 /= ni_g).or.(ni2 /= ni_g) ) THEN
         WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : the 3 files dont agree for x length.' ; STOP
      END IF

      !! testing nj agreement :
      IF( (nj1 /= nj2).or.(nj1 /= nj_g).or.(nj2 /= nj_g) ) THEN
         WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : the 3 files dont agree for y length.'; STOP
      END IF

      ni = ni1 ; nj = nj1

      IF( l_read_angles ) THEN
         CALL DIMS(cf_ang, 'cost', ni_g, nj_g, nk_g, Ntr)
         PRINT *, 'ni, nj, ni_g, nj_g', ni, nj, ni_g, nj_g
         IF( (nj /= nj_g).OR.(ni /= ni_g) ) THEN
            WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : `Angle file` does not agree in shape with setup!'; STOP
         END IF
      END IF



      !! testing nk agreement :
      IF( nk1 /= nk2 ) THEN
         WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : u and v files dont agree for z length.'
         STOP
      END IF

      nk = nk1


      !! testing nt agreement :
      IF( Ntr1 /= Ntr2 ) THEN
         WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : u and v files dont agree for time length.'
         STOP
      END IF

      Ntr = Ntr1


      IF( nk > 1 ) THEN
         i3d = 1
      ELSE
         nk  = 1
         i3d = 0
      END IF

      WRITE(*,'("Space dimension is : ",i4," x",i4," x",i4)') ni, nj, nk
      WRITE(*,'(i4," time records for u and v")') Ntr
      WRITE(6,*)''


      !! Allocations :
      ALLOCATE ( ztmp4(ni,nj) ,  &
         &    U_c(ni,nj,nk) , V_c(ni,nj,nk), mask_t(ni,nj,nk), &
         &    xlon_t(ni,nj)  , xlat_t(ni,nj), Xdum(ni,nj), Ydum(ni,nj), &
         &    xlon_u(ni,nj)  , xlat_u(ni,nj) , xlon_v(ni,nj)  , xlat_v(ni,nj) , &
         &    xlon_f(ni,nj) , xlat_f(ni,nj), ztmp(ni,nj), &
         &    vtime(Ntr)   )
      IF( cgrid_trg == 'U,V' ) ALLOCATE ( mask_u(ni,nj,nk) , mask_v(ni,nj,nk) )

      IF( lNPfix ) ALLOCATE ( mask_np(ni,nj)  )

      ALLOCATE (XCOS_t(ni,nj) , XSIN_t(ni,nj) , XCOS_u(ni,nj) , XSIN_u(ni,nj) , XCOS_v(ni,nj) , XSIN_v(ni,nj) , &
         &      XCOS_f(ni,nj) , XSIN_f(ni,nj) , zU(ni,nj) , zV(ni,nj) )


      IF( i3d == 1 ) THEN
         ALLOCATE ( vdepth(nk) )
         CALL GETVAR_1D(cf_mm, cv_depth, vdepth)
         IF( cgrid_trg == 'U,V' ) THEN
            CALL GETMASK_3D(cf_mm, 'umask', mask_u(:,:,:))
            CALL GETMASK_3D(cf_mm, 'vmask', mask_v(:,:,:))
         ELSE
            CALL GETMASK_3D(cf_mm, 'tmask', mask_t(:,:,:))
         END IF
      ELSE
         IF( cgrid_trg == 'U,V' ) THEN
            CALL GETMASK_2D(cf_mm, 'umask', mask_u(:,:,1), jlev=1)
            CALL GETMASK_2D(cf_mm, 'vmask', mask_v(:,:,1), jlev=1)
         ELSE
            CALL GETMASK_2D(cf_mm, 'tmask', mask_t(:,:,1), jlev=1)
         END IF
      END IF

      !!  Getting longitude and latitude form grid file :
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamt, 1, 1, 1, xlon_t)
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphit, 1, 1, 1, xlat_t)
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamu, 1, 1, 1, xlon_u)
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiu, 1, 1, 1, xlat_u)
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamv, 1, 1, 1, xlon_v)
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiv, 1, 1, 1, xlat_v)
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamf, 1, 1, 1, xlon_f)
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphif, 1, 1, 1, xlat_f)
      
      WRITE(6,*)''

      !! Is this a known ORCA grid (just for info now, not used!):
      gt_orca = IS_ORCA_NORTH_FOLD( xlat_t )
      iorca = gt_orca%ifld_nord
      IF( iorca == 4 ) WRITE(6,*)' Grid is an ORCA grid with north-pole T-point folding!'
      IF( iorca == 6 ) WRITE(6,*)' Grid is an ORCA grid with north-pole F-point folding!'
      WRITE(6,*)' Grid is not an ORCA grid!'
      WRITE(6,*)''

      WRITE(6,*)''
      IF( l_read_angles ) THEN
         WRITE(6,*)' *** Reading COS and SIN of rotation angles at t,u,v,f points'
         WRITE(6,*)'     => in file "',TRIM(cf_ang),'"'
         CALL GETVAR_2D(i0, j0, cf_ang, 'cost', 1, 1, 1, XCOS_t)
         CALL GETVAR_2D(i0, j0, cf_ang, 'sint', 1, 1, 1, XSIN_t)
         CALL GETVAR_2D(i0, j0, cf_ang, 'cosu', 1, 1, 1, XCOS_u)
         CALL GETVAR_2D(i0, j0, cf_ang, 'sinu', 1, 1, 1, XSIN_u)
         CALL GETVAR_2D(i0, j0, cf_ang, 'cosv', 1, 1, 1, XCOS_v)
         CALL GETVAR_2D(i0, j0, cf_ang, 'sinv', 1, 1, 1, XSIN_v)
         CALL GETVAR_2D(i0, j0, cf_ang, 'cosf', 1, 1, 1, XCOS_f)
         CALL GETVAR_2D(i0, j0, cf_ang, 'sinf', 1, 1, 1, XSIN_f)
         WRITE(6,*)''; WRITE(6,*)' *** Done reading!'
      ELSE
         !!  Getting cosine and sine corresponding to the angle of the local distorsion of the grid:
         WRITE(6,*)' *** Computing COS and SIN of rotation angles at t,u,v,f points with `ANGLE2()`'
         CALL ANGLE2( iorca, xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f, &
            &        XCOS_t, XSIN_t, XCOS_u, XSIN_u, XCOS_v, XSIN_v, XCOS_f, XSIN_f )
      END IF
      WRITE(6,*)''; WRITE(6,*)''

      XCOS_t(:,1) = 0.
      XSIN_t(:,1) = 0.
      
      IF( ldebug ) THEN
         CALL DUMP_FIELD(REAL(XCOS_t,4), 'cost_angle.nc', 'cost')
         CALL DUMP_FIELD(REAL(XSIN_t,4), 'sint_angle.nc', 'sint')
         CALL DUMP_FIELD(REAL(XCOS_u,4), 'cosu_angle.nc', 'cosu')
         CALL DUMP_FIELD(REAL(XSIN_u,4), 'sinu_angle.nc', 'sinu')
         CALL DUMP_FIELD(REAL(XCOS_v,4), 'cosv_angle.nc', 'cosv')
         CALL DUMP_FIELD(REAL(XSIN_v,4), 'sinv_angle.nc', 'sinv')
         CALL DUMP_FIELD(REAL(XCOS_f,4), 'cosf_angle.nc', 'cosf')
         CALL DUMP_FIELD(REAL(XSIN_f,4), 'sinf_angle.nc', 'sinf')
         CALL DUMP_FIELD( ATAN(XSIN_t/XCOS_t)/rad, 'ANGLE_t.nc', 'angle_t_deg' )
      END IF

      IF( lNPfix ) THEN
         mask_np(:,:) = 1.
         WHERE( xlat_t > rNP_lat ) mask_np = 0.
         IF(ldebug) CALL DUMP_FIELD( REAL(mask_np,4) , 'MASK_NP.nc', 'mask' )
      END IF
      
      !!  Getting time from the u_raw file or the namelist :
      IF( lct ) THEN       ! time is being controlled
         DO jt = 1, Ntr
            vtime(jt) = t0 + t_stp*REAL(jt)
         END DO
         nb_att_t = 1
         vatt_info_t(:)%cname = 'null'
         vatt_info_t(1)%cname = 'units'
         vatt_info_t(1)%itype = 2 ! char
         vatt_info_t(1)%val_char = 'unknown'
         vatt_info_t(1)%ilength = LEN('unknown')
      ELSE                  ! we use time from input file
         CALL GETVAR_1D(cf_raw_U, cv_t_out, vtime)
         CALL GETVAR_ATTRIBUTES(cf_raw_U, cv_t_out, nb_att_t, vatt_info_t)
      END IF



      DO jt = 1, Ntr

         WRITE(6,*)''; WRITE(6,*)'Time step =', jt ; WRITE(6,*)''

         DO jk = 1, nk

            IF( l_anlt ) THEN
               ! Analytical field for debugging purposes...
               zU = 1.
               zV = 0.
            ELSE
               !! Getting uncorrected U on grid T:
               CALL  GETVAR_2D(idf_u, idv_u, cf_raw_U, cv_out_U, Ntr, jk*i3d, jt, ztmp4, Nk=nk)
               zU = ztmp4

               !! Getting uncorrected V on grid T:
               CALL  GETVAR_2D(idf_v, idv_v, cf_raw_V, cv_out_V, Ntr, jk*i3d, jt, ztmp4, Nk=nk)
               zV = ztmp4
            END IF


            
            !! Getting U and V on T-grid:
            !! Correcting U (East-North to i-component) :
            Xdum = XCOS_t*zU + XSIN_t*zV
            !! Correcting V (East-North to j-component) :
            Ydum = XCOS_t*zV - XSIN_t*zU


            IF( TRIM(cgrid_trg) == 'T' ) THEN
               !! to T-grid, nothing to do...
               IF( iorca > 0 ) CALL lbc_lnk( iorca, Xdum, 'T', -1.0_8 )
               U_c(:,:,jk) = REAL(Xdum , 4)
               IF( iorca > 0 ) CALL lbc_lnk( iorca, Ydum, 'T', -1.0_8 )
               V_c(:,:,jk) = REAL(Ydum , 4)
               !!
            ELSEIF( cgrid_trg == 'U,V' ) THEN
               !! to U-V grid, need to interpolate
               ztmp(1:ni-1,:) = 0.5*(Xdum(1:ni-1,:)+Xdum(2:ni,:))
               IF( iorca > 0 ) CALL lbc_lnk( iorca, ztmp, 'U', -1.0_8 )
               U_c(:,:,jk) = REAL(ztmp , 4)
               !!
               ztmp(:,1:nj-1) = 0.5*(Ydum(:,1:nj-1)+Ydum(:,2:nj))
               IF( iorca > 0 ) CALL lbc_lnk( iorca, ztmp, 'V', -1.0_8 )
               V_c(:,:,jk) = REAL(ztmp , 4)
               !!
            END IF

            IF( lNPfix ) THEN
               CALL DROWN( ewper_trg, U_c(:,:,jk), mask_np(:,:), nb_inc=100 )
               CALL DROWN( ewper_trg, V_c(:,:,jk), mask_np(:,:), nb_inc=100 )
            END IF
            
         END DO ! jk


         IF( lmout_x .AND. lmout_y ) THEN
            IF( cgrid_trg == 'U,V' ) THEN
               WHERE ( mask_u == 0 ) U_c = rmiss_val
               WHERE ( mask_v == 0 ) V_c = rmiss_val
            ELSE
               WHERE ( mask_t == 0 ) U_c = rmiss_val
               WHERE ( mask_t == 0 ) V_c = rmiss_val
            END IF
         ELSE
            rmiss_val = 0.
         END IF

         !! Boundaries: #fixme
         IF( ewper_trg >= 1 ) THEN
            U_c(1,:,:) = U_c(ni-ewper_trg+1,:,:)
            V_c(1,:,:) = V_c(ni-ewper_trg+1,:,:)
         ELSE
            U_c(1,:,:) = rmiss_val
            V_c(1,:,:) = rmiss_val
         END IF
         U_c(:, 1,:) = rmiss_val
         V_c(:, 1,:) = rmiss_val
         U_c(:,nj,:) = rmiss_val
         V_c(:,nj,:) = rmiss_val



         IF( ABS(rmiss_val) < 1.e-6 ) THEN
            !! ==> rmis_val==0 !
            rmv_eff = -9999.
         ELSE
            rmv_eff = rmiss_val
         ENDIF


         IF( i3d == 1 ) THEN

            !! 3D:
            IF( cgrid_trg == 'U,V' ) THEN
               CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vdepth, vtime, U_c(:,:,:),  &
                  &    cf_out_U, 'nav_lon_u', 'nav_lat_u', cv_depth, cv_t_out, cv_rot_U,      &
                  &    rmv_eff, attr_t=vatt_info_t)
               CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vdepth, vtime, V_c(:,:,:),  &
                  &    cf_out_V, 'nav_lon_v', 'nav_lat_v', cv_depth, cv_t_out, cv_rot_V,      &
                  &    rmv_eff, attr_t=vatt_info_t)
            ELSE
               CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_t, xlat_t, vdepth, vtime, U_c(:,:,:),  &
                  &    cf_out_U, 'nav_lon', 'nav_lat', cv_depth, cv_t_out, cv_rot_U,      &
                  &    rmv_eff, attr_t=vatt_info_t)
               CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_t, xlat_t, vdepth, vtime, V_c(:,:,:),  &
                  &    cf_out_V, 'nav_lon', 'nav_lat', cv_depth, cv_t_out, cv_rot_V,      &
                  &    rmv_eff, attr_t=vatt_info_t)
            END IF

         ELSE

            !! 2D:
            IF( cgrid_trg == 'U,V' ) THEN
               CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u,         vtime, U_c(:,:,1), &
                  &    cf_out_U, 'nav_lon_u', 'nav_lat_u', cv_t_out, cv_rot_U,       &
                  &    rmv_eff, attr_t=vatt_info_t)
               CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v,         vtime, V_c(:,:,1), &
                  &    cf_out_V, 'nav_lon_v', 'nav_lat_v', cv_t_out, cv_rot_V,   &
                  &    rmv_eff, attr_t=vatt_info_t)
            ELSE
               CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_t, xlat_t,         vtime, U_c(:,:,1), &
                  &    cf_out_U, 'nav_lon', 'nav_lat', cv_t_out, cv_rot_U,       &
                  &    rmv_eff, attr_t=vatt_info_t)
               CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_t, xlat_t,         vtime, V_c(:,:,1), &
                  &    cf_out_V, 'nav_lon', 'nav_lat', cv_t_out, cv_rot_V,   &
                  &    rmv_eff, attr_t=vatt_info_t)
            END IF

         END IF

      END DO ! jt











































   ELSE
      !! !!     I N V E R S E   C O R R E C T I O N

      WRITE(6,*)'Will unrotate vector fields given on an irregular grid!'
      !!
      !IF( l_reg_src ) THEN
      !   WRITE(6,*)'Reverse vector correction only makes sense if your source grid is distorded!'
      !   WRITE(6,*)'  => check "l_reg_src" into the namelist...' ; WRITE(6,*)''; STOP
      !END IF
      !!
      !!
      !!
      cf_raw_U = trim(cd_out)//'/'//trim(cv_out_U)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.nc'
      cf_raw_V = trim(cd_out)//'/'//trim(cv_out_V)//'_'//trim(csource)//'-' &
         &   //trim(ctarget)//'_'//trim(cextra)//'.nc'
      !!
      WRITE(6,*)'The two unrotated raw files will be produced :'
      PRINT *, trim(cf_raw_U) ;   PRINT *, trim(cf_raw_V) ;
      !!
      WRITE(6,*)'' ;   WRITE(6,*)'Input files :'
      PRINT *, trim(cufilin) ;   PRINT *, trim(cvfilin) ; WRITE(6,*)'' ; WRITE(6,*)''
      !!
      WRITE(6,*)'File containing input grid :'
      PRINT *, trim(cf_mm) ; WRITE(6,*)''


      !! Creating name for unrotated output file:
      nbc = LEN_TRIM(cufilin)
      cdum = cufilin(nbc-2:nbc)
      IF( cdum /= '.nc' ) THEN
         !IF( cdum == 'nc4' ) THEN
         !   cfext = 'nc4'
         !   nlext = 4
         !ELSE
         WRITE(6,*)'Unknow file extension for ',TRIM(cufilin) ; STOP
         !END IF
      END IF

      cf_raw_U = cufilin(1:nbc-nlext)
      cf_raw_U = trim(cf_raw_U)//'_unrotated.nc'

      nbc = LEN_TRIM(cvfilin)
      cf_raw_V = cvfilin(1:nbc-nlext)
      cf_raw_V = TRIM(cf_raw_V)//'_unrotated.nc'

      cv_out_U = trim(cv_u_in)//'_unrotated'
      cv_out_V = trim(cv_v_in)//'_unrotated'


      !! Geting array dimension and testing...
      !! -------------------------------------

      CALL DIMS(cufilin, cv_u_in, ni1, nj1, nk1, Ntr1)
      CALL DIMS(cvfilin, cv_v_in, ni2, nj2, nk2, Ntr2)

      CALL DIMS(cf_mm,   cv_glamt, ni_g, nj_g, nk_g, Ntr)

      !! testing ni agreement :
      IF( (ni1 /= ni2).or.(ni1 /= ni_g).or.(ni2 /= ni_g) ) THEN
         WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : the 3 files dont agree for x length.'
         STOP
      END IF

      !! testing nj agreement :
      IF( (nj1 /= nj2).or.(nj1 /= nj_g).or.(nj2 /= nj_g) ) THEN
         WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : the 3 files dont agree for y length.'
         STOP
      END IF


      nk = 1 ; i3d = 0

      !! testing 3D and nk agreement :

      IF( (nk1 /= nk2) ) THEN
         WRITE(6,*)'SHAPE CONSISTENCY ERROR:! : the 2 files dont agree for number of levels.'; STOP
      END IF


      IF( nk1 /= -1 ) THEN
         l_3d_inv = .TRUE.
         nk = nk1
         i3d = 1
         WRITE(6,*)''; WRITE(6,*)'Will perform 3D un-rotating!!!' ; WRITE(6,*)''
      END IF


      !! testing nt agreement :
      IF( Ntr1 /= Ntr2 ) THEN
         WRITE(6,*)'CONSISTENCY ERROR: u and v files dont agree for time length.'
         STOP
      END IF

      ni = ni1 ; nj = nj1 ; Ntr = Ntr1

      WRITE(*,'("Dimension is : ",i4," x",i4)') ni, nj
      WRITE(6,*)'Number of levels to treat =>', nk
      WRITE(*,'(i4," time records for u and v")') Ntr
      WRITE(6,*)''



      !! Allocations :
      !! -------------
      ALLOCATE (ztmp4(ni,nj) , &
         &     U_c(ni,nj,nk) , V_c(ni,nj,nk),                              &
         &     xlon_t(ni,nj) , xlat_t(ni,nj) , xlon_u(ni,nj) , xlat_u(ni,nj), xlon_v(ni,nj) , xlat_v(ni,nj) , &
         &     xlon_f(ni,nj) , xlat_f(ni,nj) , vtime(Ntr) , mask_u(ni,nj,nk) , mask_v(ni,nj,nk)  )
      !!
      ALLOCATE (XCOS_t(ni,nj) , XSIN_t(ni,nj) , XCOS_u(ni,nj) , XSIN_u(ni,nj) , XCOS_v(ni,nj) , XSIN_v(ni,nj) , &
         &      XCOS_f(ni,nj) , XSIN_f(ni,nj) , zU(ni,nj) , zV(ni,nj) )


      !!  Getting longitude and latitude form grid file :
      !! ------------------------------------------------
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamt, 1, 1, 1, xlon_t)  ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphit, 1, 1, 1, xlat_t)  ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamu, 1, 1, 1, xlon_u)  ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiu, 1, 1, 1, xlat_u)  ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamv, 1, 1, 1, xlon_v)  ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphiv, 1, 1, 1, xlat_v)  ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, cv_glamf, 1, 1, 1, xlon_f)  ; i0=0 ; j0=0
      CALL GETVAR_2D(i0, j0, cf_mm, cv_gphif, 1, 1, 1, xlat_f)  ; i0=0 ; j0=0

      IF( l_3d_inv ) THEN
         ALLOCATE ( vdepth(nk) )
         CALL GETVAR_1D(cf_mm, cv_depth, vdepth)
         CALL GETMASK_3D(cf_mm, 'umask', mask_u(:,:,:))
         CALL GETMASK_3D(cf_mm, 'vmask', mask_v(:,:,:))
      ELSE
         CALL GETMASK_2D(cf_mm, 'umask', mask_u(:,:,1))
         CALL GETMASK_2D(cf_mm, 'vmask', mask_v(:,:,1))
      END IF

      !!  Getting cos and sin of the grid distorsion angle:
      !! --------------------------------------------------
      CALL ANGLE2( iorca, xlon_t, xlat_t, xlon_u, xlat_u, xlon_v, xlat_v, xlon_f, xlat_f, &
         &         XCOS_t, XSIN_t, XCOS_u, XSIN_u, XCOS_v, XSIN_v, XCOS_f, XSIN_f )



      !!  Getting time from the u_raw file or the namelist :
      !! ---------------------------------------------------
      IF( lct ) THEN       ! time is being controlled
         DO jt = 1, Ntr
            vtime(jt) = t0 + t_stp*REAL(jt)
         END DO
         nb_att_t = 1
         vatt_info_t(:)%cname = 'null'
         vatt_info_t(1)%cname = 'units'
         vatt_info_t(1)%itype = 2 ! char
         vatt_info_t(1)%val_char = 'unknown'
         vatt_info_t(1)%ilength = LEN('unknown')

      ELSE                  ! we use time from input file
         CALL GETVAR_1D(cufilin, cv_time_0, vtime)
         CALL GETVAR_ATTRIBUTES(cufilin, cv_time_0, nb_att_t, vatt_info_t)
      END IF




      DO jt = 1, Ntr

         WRITE(6,*)''
         WRITE(6,*)''
         WRITE(6,*)''; WRITE(6,*)'Time step =', jt ; WRITE(6,*)''

         DO jk = 1, nk



            WRITE(6,*)'  *** Un-rotating level =', jk

            !! Getting U :
            !! -----------
            CALL  GETVAR_2D(idf_u, idv_u, cufilin, cv_u_in, Ntr, jk*i3d, jt, ztmp4, Nk=nk)
            zU = ztmp4



            !! Getting V :
            !! -----------
            CALL  GETVAR_2D(idf_v, idv_v, cvfilin, cv_v_in, Ntr, jk*i3d, jt, ztmp4, Nk=nk)
            zV = ztmp4

            !! Unrotating U :
            !! --------------
            !#LB: zU is on U-points!
            U_c(:,:,jk) = REAL(XCOS_t*zU - XSIN_t*zV , 4) ! note the '-' sign --> reverse correction
            !#LB: on what grid points is this supposed to be??? U or T? Presently the rest of the code assume it's U !!!


            !! Unrotating V :
            !! --------------
            !#LB: zV is on V-points!
            V_c(:,:,jk) = REAL(XCOS_t*zV + XSIN_t*zU , 4) ! note the + sign for reverse correction
            !#LB: on what grid points is this supposed to be??? V or T? Presently the rest of the code assume it's V !!!


         END DO  ! jk

         !#LB: need to fix the following if I'm wrong about U,V or T
         !#LB:  + also INCLUDE the possibility to chose 'T' or 'U,V' for the grid on which to save, just as for the rotation case (lines 522->539)

         WHERE ( mask_u == 0 ) U_c(:,:,1:nk) = rmissval
         WHERE ( mask_v == 0 ) V_c(:,:,1:nk) = rmissval

         IF( l_3d_inv ) THEN

            CALL P3D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vdepth, vtime, U_c(:,:,:),  &
               &    cf_raw_U, 'nav_lon_u', 'nav_lat_u', cv_depth, cv_time_0, cv_out_U,       &
               &    rmissval, attr_t=vatt_info_t)

            CALL P3D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vdepth, vtime, V_c(:,:,:),  &
               &    cf_raw_V, 'nav_lon_v', 'nav_lat_v', cv_depth, cv_time_0, cv_out_V,       &
               &    rmissval, attr_t=vatt_info_t)

         ELSE

            CALL P2D_T(id_f1, id_v1, Ntr, jt, xlon_u, xlat_u, vtime, U_c(:,:,1),     &
               &    cf_raw_U, 'nav_lon_u', 'nav_lat_u', cv_time_0, cv_out_U,       &
               &    rmissval, attr_t=vatt_info_t)

            CALL P2D_T(id_f2, id_v2, Ntr, jt, xlon_v, xlat_v, vtime, V_c(:,:,1), &
               &    cf_raw_V, 'nav_lon_v', 'nav_lat_v', cv_time_0, cv_out_V,   &
               &    rmissval, attr_t=vatt_info_t)

         END IF


      END DO

      WRITE(6,*)''
      WRITE(6,*)'Files created:'; WRITE(6,*)'   ', trim(cf_raw_U); WRITE(6,*)'   ', trim(cf_raw_V)
      WRITE(6,*)''

   END IF ! on l_inv


   WRITE(6,*)'Done!'


CONTAINS


   SUBROUTINE usage_corr_vect()
      !!
      WRITE(6,*)''
      WRITE(6,*)'   List of command line options:'
      WRITE(6,*)''
      WRITE(6,*)' -I   => will perform inverse correction, ie un-rotate a ORCA grid'
      WRITE(6,*)''
      WRITE(6,*)''
      WRITE(6,*)'  *** MANDATORY for both normal and inverse mode:'
      WRITE(6,*)''
      WRITE(6,*)' -m  <mesh_mask_file> => Specify which mesh_mask file to use'
      WRITE(6,*)''
      WRITE(6,*)''
      WRITE(6,*)'  ***  MANDATORY for normal mode (no -I switch) :'
      WRITE(6,*)''
      WRITE(6,*)' Note: In "normal" aka "rotation" mode, corr_vect.x will read U and V'
      WRITE(6,*)'       bluntly interpolated from a regular grid onto the T-grid points'
      WRITE(6,*)'       of the TARGET NEMO grid, and rotate them according to the local'
      WRITE(6,*)'       distortion of the TARGET NEMO grid.'
      WRITE(6,*)'       Result is the rotated vector on T-grid points, unless the "-G U,V"'
      WRITE(6,*)'       switch is used!'
      WRITE(6,*)'       If "-G U,V" is used rotated vector is interpolated on U and V points.'
      WRITE(6,*)''
      WRITE(6,*)' -f  <namelist_prefix> => common file name prefix for the 2 namelists to be read'
      WRITE(6,*)'                       => expects to find <namelist_prefix>_x & <namelist_prefix>_y !'
      WRITE(6,*)'                         [no namelist needed when inverse correction]'
      WRITE(6,*)''
      WRITE(6,*)''
      WRITE(6,*)'  *** MANDATORY for INVERSE MODE (-I switch) :'
      WRITE(6,*)''
      WRITE(6,*)' Note: In "inverse" aka "de-rotation" mode, corr_vect.x will read U and V'
      WRITE(6,*)'      of a vector native of NEMO (like current field output from NEMO)'
      WRITE(6,*)'      and will "un-rotate" it so each component can then be inter-'
      WRITE(6,*)'      polated on a regular grid.'
      WRITE(6,*)'      Input vector components, U,V, are provided on their respective'
      WRITE(6,*)'      U- & V-grid points.'
      WRITE(6,*)'      Result is the "de-rotated" vector on T-grid points of NEMO grid.'
      WRITE(6,*)''
      WRITE(6,*)' -i <x.nc> <y.nc>     => unrotate vector fields given in these 2 files'
      WRITE(6,*)'                         to the same grid'
      WRITE(6,*)''
      WRITE(6,*)' -v  <nameU> <nameV>  => names for x and y components in intput files'
      WRITE(6,*)''
      WRITE(6,*)' -t <time_name>       => name of time variable in <x.nc> and <y.nc>'
      WRITE(6,*)''
      WRITE(6,*)''
      WRITE(6,*)'  *** MISC options :'
      WRITE(6,*)''
      WRITE(6,*)' -A  <angle_file>     => if file name provided, will read sin,cos of angle rather than compute them!'
      WRITE(6,*)''
      WRITE(6,*)' -G  <T/U,V>          => Specify grid points where to save the (un)rotated vector'
      WRITE(6,*)'                       ==>    T    points "-G T" DEFAULT !'
      WRITE(6,*)'                       ==> U and V points "-G U,V"'
      WRITE(6,*)''
      WRITE(6,*)'  *** Misc "ad-hoc" fix options:'
      WRITE(6,*)''
      WRITE(6,*)' -P <lat_lim> => fix fields beyond this latitude (ex: -P 89.7)'
      WRITE(6,*)''
      WRITE(6,*)''
      WRITE(6,*)' -h                   => Show this message'
      WRITE(6,*)''
      STOP
      !!
   END SUBROUTINE usage_corr_vect
   !!
END PROGRAM CORR_VECT
