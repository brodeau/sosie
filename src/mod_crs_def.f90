MODULE mod_crs_def
   !!======================================================================
   !!                         ***  MODULE crs_dom  ***
   !!        Declare the coarse grid domain and other public variables
   !!        then allocate them if needed.
   !!======================================================================
   !!  History     2012-06  Editing  (J. Simeon, G. Madec, C. Ethe, C. Calone) Original code
   !!----------------------------------------------------------------------
   !USE par_oce
   !USE dom_oce
   !USE in_out_manager
   USE mod_nemo
   
   IMPLICIT NONE
   PUBLIC

   !PUBLIC crs_dom_alloc  ! Called from crsini.F90
   !PUBLIC crs_dom_alloc2  ! Called from crsini.F90
   !PUBLIC dom_grid_glo
   PUBLIC dom_grid_crs


   !!Known from dom_oce.F90:
   !!    nlcit , nlcjt    !: dimensions of every subdomain


   

   
   ! Domain variables
   INTEGER  ::  jpiglo_crs ,   &             !: 1st dimension of global coarse grid domain
      jpjglo_crs                   !: 2nd dimension of global coarse grid domain
   INTEGER  ::  jpi_crs ,   &                !: 1st dimension of local coarse grid domain
      jpj_crs                      !: 2nd dimension of local coarse grid domain
   INTEGER  ::  jpi_full ,  &                !: 1st dimension of local parent grid domain
      jpj_full                     !: 2nd dimension of local parent grid domain

   INTEGER  ::  nistr , njstr
   INTEGER  ::  niend , njend

   INTEGER  ::  jpi_crsm1, jpj_crsm1         !: loop indices
   INTEGER  ::  jpiglo_crsm1, jpjglo_crsm1   !: loop indices
   INTEGER  ::  nperio_full, nperio_crs      !: jperio of parent and coarse grids
   INTEGER  ::  npolj_full, npolj_crs        !: north fold mark
   INTEGER  ::  jpiglo_full, jpjglo_full     !: jpiglo / jpjglo
   INTEGER  ::  npiglo, npjglo               !: jpjglo
   INTEGER  ::  nlci_full, nlcj_full         !: i-, j-dimension of local or sub domain on parent grid
   INTEGER  ::  nldi_full, nldj_full         !: starting indices of internal sub-domain on parent grid
   INTEGER  ::  nlei_full, nlej_full         !: ending indices of internal sub-domain on parent grid
   INTEGER  ::  nlci_crs, nlcj_crs           !: i-, j-dimension of local or sub domain on coarse grid
   INTEGER  ::  nldi_crs, nldj_crs           !: starting indices of internal sub-domain on coarse grid
   INTEGER  ::  nlei_crs, nlej_crs           !: ending indices of internal sub-domain on coarse grid

   INTEGER  ::  narea_full, narea_crs        !: node
   INTEGER  ::  jpnij_full, jpnij_crs        !: =jpni*jpnj, the pe decomposition
   INTEGER  ::  jpim1_full, jpjm1_full       !:
   INTEGER  ::  nimpp_full, njmpp_full       !: global position of point (1,1) of subdomain on parent grid
   INTEGER  ::  nimpp_crs, njmpp_crs         !: set to 1,1 for now .  Valid only for monoproc
   INTEGER  ::  nreci_full, nrecj_full
   INTEGER  ::  nreci_crs, nrecj_crs
   !cc
   INTEGER ::   noea_full, nowe_full        !: index of the local neighboring processors in
   INTEGER ::   noso_full, nono_full        !: east, west, south and north directions
   INTEGER ::   npne_full, npnw_full        !: index of north east and north west processor
   INTEGER ::   npse_full, npsw_full        !: index of south east and south west processor
   INTEGER ::   nbne_full, nbnw_full        !: logical of north east & north west processor
   INTEGER ::   nbse_full, nbsw_full        !: logical of south east & south west processor
   INTEGER ::   nidom_full                  !: ???
   INTEGER ::   nproc_full                  !:number for local processor
   INTEGER ::   nbondi_full, nbondj_full    !: mark of i- and j-direction local boundaries
   INTEGER ::   noea_crs, nowe_crs          !: index of the local neighboring processors in
   INTEGER ::   noso_crs, nono_crs          !: east, west, south and north directions
   INTEGER ::   npne_crs, npnw_crs          !: index of north east and north west processor
   INTEGER ::   npse_crs, npsw_crs          !: index of south east and south west processor
   INTEGER ::   nbne_crs, nbnw_crs          !: logical of north east & north west processor
   INTEGER ::   nbse_crs, nbsw_crs          !: logical of south east & south west processor
   INTEGER ::   nidom_crs                   !: ???
   INTEGER ::   nproc_crs                   !:number for local processor
   INTEGER ::   nbondi_crs, nbondj_crs      !: mark of i- and j-direction local boundaries


   INTEGER, DIMENSION(:), ALLOCATABLE :: mis_crs, mie_crs, mis2_crs, mie2_crs  ! starting and ending i-indices of parent subset
   INTEGER, DIMENSION(:), ALLOCATABLE :: mjs_crs, mje_crs, mjs2_crs, mje2_crs ! starting and ending  j-indices of parent subset
   INTEGER, DIMENSION(:), ALLOCATABLE :: mjg_crs, mig_crs
   INTEGER, DIMENSION(:), ALLOCATABLE :: mi0_crs, mi1_crs, mj0_crs, mj1_crs
   INTEGER  :: mxbinctr, mybinctr            ! central point in grid box
   INTEGER, DIMENSION(:), ALLOCATABLE ::   nlcit_crs, nlcit_full  !: dimensions of every subdomain
   INTEGER, DIMENSION(:), ALLOCATABLE ::   nldit_crs, nldit_full     !: first, last indoor index for each i-domain
   INTEGER, DIMENSION(:), ALLOCATABLE ::   nleit_crs, nleit_full    !: first, last indoor index for each j-domain
   INTEGER, DIMENSION(:), ALLOCATABLE ::   nimppt_crs, nimppt_full    !: first, last indoor index for each j-domain
   !INTEGER, DIMENSION(:), ALLOCATABLE ::   nlcjt_crs, nlcjt_full  !: dimensions of every subdomain
   INTEGER, DIMENSION(:), ALLOCATABLE ::   nldjt_crs, nldjt_full     !: first, last indoor index for each i-domain
   INTEGER, DIMENSION(:), ALLOCATABLE ::   nlejt_crs, nlejt_full    !: first, last indoor index for each j-domain
   INTEGER, DIMENSION(:), ALLOCATABLE ::   njmppt_crs, njmppt_full    !: first, last indoor index for each j-domain


   ! Masks
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: tmask_crs, umask_crs, vmask_crs, fmask_crs
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE :: tmask_i_crs, rnfmsk_crs, tpol_crs, fpol_crs

   !    REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: tmask_i_crs, tpol, fpol

   ! Scale factors
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1t_crs, e2t_crs, e1e2t_crs ! horizontal scale factors grid type T
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1u_crs, e2u_crs ! horizontal scale factors grid type U
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1v_crs, e2v_crs ! horizontal scale factors grid type V
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1f_crs, e2f_crs ! horizontal scale factors grid type F
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_crs, e3u_crs, e3v_crs, e3f_crs, e3w_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_max_crs, e3u_max_crs, e3v_max_crs, e3f_max_crs, e3w_max_crs

   ! Surface
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e1e2w_crs, e2e3u_crs, e1e3v_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e1e2w_msk, e2e3u_msk, e1e3v_msk
   ! vertical scale factors
   ! Coordinates
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: gphit_crs, glamt_crs, gphif_crs, glamf_crs
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: gphiu_crs, glamu_crs, gphiv_crs, glamv_crs
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: ff_crs
   INTEGER,  DIMENSION(:,:),   ALLOCATABLE :: mbathy_crs, mbkt_crs, mbku_crs, mbkv_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: gdept_crs, gdepu_crs, gdepv_crs, gdepw_crs

   ! Weights
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: facsurfv, facsurfu, facvol_t, facvol_w
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ocean_volume_crs_t, ocean_volume_crs_w, bt_crs, r1_bt_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: crs_surfu_wgt, crs_surfv_wgt, crs_surfw_wgt, crs_volt_wgt

   ! CRS Namelist
   INTEGER           :: nn_factx   = 3       !: reduction factor of x-dimension of the parent grid
   INTEGER           :: nn_facty   = 3       !: reduction factor of y-dimension of the parent grid
   INTEGER           :: nn_binref  = 0       !: 0 = binning starts north fold (equator could be asymmetric)
   !: 1 = binning centers at equator (north fold my have artifacts)
   !:    for even reduction factors, equator placed in bin biased south
   INTEGER           :: nn_msh_crs = 1       !: Organization of mesh mask output
   !: 0 = no mesh mask output
   !: 1 = unified mesh mask output
   !: 2 = 2 separate mesh mask output
   !: 3 = 3 separate mesh mask output
   INTEGER           :: nn_crs_kz    =    0       !: type of Kz coarsening ( =0->VOL ; =1->MAX ; =2->MIN)
   LOGICAL           :: ln_crs_wn    = .FALSE.    !: coarsening wn or computation using horizontal divergence
   !
   INTEGER           :: nrestx, nresty       !: for determining odd or even reduction factor


   ! Grid reduction factors
   REAL(wp)     ::  rfactx_r                !: inverse of x-dim reduction factor
   REAL(wp)     ::  rfacty_r                !: inverse of y-dim reduction factor
   REAL(wp)     ::  rfactxy

   ! Physical and dynamical ocean fields for output or passing to TOP, time-mean fields

   ! Direction of lateral diffusion


   !! $Id: crs.F90 5215 2015-04-15 16:11:56Z nicolasmartin $
CONTAINS

   INTEGER FUNCTION crs_dom_alloc()
      !!-------------------------------------------------------------------
      !!                     *** FUNCTION crs_dom_alloc ***
      !!  ** Purpose :   Allocate public crs arrays
      !!-------------------------------------------------------------------
      !! Local variables
      INTEGER, DIMENSION(17) :: ierr

      ierr(:) = 0

      ! Set up bins for coarse grid, horizontal only.
      ALLOCATE( mis2_crs(jpiglo_crs), mie2_crs(jpiglo_crs),  &
         &       mjs2_crs(jpjglo_crs), mje2_crs(jpjglo_crs),  &
         &       mi0_crs (jpiglo_crs), mi1_crs (jpiglo_crs),  &
         &       mj0_crs (jpjglo_crs), mj1_crs (jpjglo_crs),  &
         &       mig_crs (jpi_crs)   , mjg_crs (jpj_crs)   ,  STAT=ierr(1) )


      ! Set up Mask and Mesh
      ALLOCATE( tmask_crs(jpi_crs,jpj_crs,jpk) , fmask_crs(jpi_crs,jpj_crs,jpk) ,  &
         &      umask_crs(jpi_crs,jpj_crs,jpk) , vmask_crs(jpi_crs,jpj_crs,jpk) , STAT=ierr(2))

      ALLOCATE( tmask_i_crs(jpi_crs,jpj_crs)   , rnfmsk_crs(jpi_crs,jpj_crs), &
         &         tpol_crs(jpiglo_crs,jpjglo_crs), fpol_crs(jpiglo_crs,jpjglo_crs), STAT=ierr(3) )

      ALLOCATE( gphit_crs(jpi_crs,jpj_crs) , glamt_crs(jpi_crs,jpj_crs) , &
         &      gphiu_crs(jpi_crs,jpj_crs) , glamu_crs(jpi_crs,jpj_crs) , &
         &      gphiv_crs(jpi_crs,jpj_crs) , glamv_crs(jpi_crs,jpj_crs) , &
         &      gphif_crs(jpi_crs,jpj_crs) , glamf_crs(jpi_crs,jpj_crs) , &
         &      ff_crs(jpi_crs,jpj_crs)    , STAT=ierr(4))

      ALLOCATE( e1t_crs(jpi_crs,jpj_crs) , e2t_crs(jpi_crs,jpj_crs) , &
         &      e1u_crs(jpi_crs,jpj_crs) , e2u_crs(jpi_crs,jpj_crs) , &
         &      e1v_crs(jpi_crs,jpj_crs) , e2v_crs(jpi_crs,jpj_crs) , &
         &      e1f_crs(jpi_crs,jpj_crs) , e2f_crs(jpi_crs,jpj_crs) , &
         &      e1e2t_crs(jpi_crs,jpj_crs), STAT=ierr(5))

      ALLOCATE( e3t_crs(jpi_crs,jpj_crs,jpk)    , e3w_crs(jpi_crs,jpj_crs,jpk)    , &
         &      e3u_crs(jpi_crs,jpj_crs,jpk)    , e3v_crs(jpi_crs,jpj_crs,jpk)    , &
         &      e3f_crs(jpi_crs,jpj_crs,jpk)    , e1e2w_msk(jpi_crs,jpj_crs,jpk)  , &
         &      e2e3u_msk(jpi_crs,jpj_crs,jpk)  , e1e3v_msk(jpi_crs,jpj_crs,jpk)  , &
         &      e1e2w_crs(jpi_crs,jpj_crs,jpk)  , e2e3u_crs(jpi_crs,jpj_crs,jpk)  , &
         &      e1e3v_crs(jpi_crs,jpj_crs,jpk)  , e3t_max_crs(jpi_crs,jpj_crs,jpk), &
         &      e3w_max_crs(jpi_crs,jpj_crs,jpk), e3u_max_crs(jpi_crs,jpj_crs,jpk), &
         &      e3v_max_crs(jpi_crs,jpj_crs,jpk), STAT=ierr(6))


      !ALLOCATE( facsurfv(jpi_crs,jpj_crs,jpk), facsurfu(jpi_crs,jpj_crs,jpk) , &
      !   &      facvol_t(jpi_crs,jpj_crs,jpk), facvol_w(jpi_crs,jpj_crs,jpk) , &
      !   &      ocean_volume_crs_t(jpi_crs,jpj_crs,jpk), ocean_volume_crs_w(jpi_crs,jpj_crs,jpk), &
      !   &      bt_crs(jpi_crs,jpj_crs,jpk)  , r1_bt_crs(jpi_crs,jpj_crs,jpk) , STAT=ierr(7))


      ALLOCATE( crs_surfu_wgt(jpi_crs,jpj_crs,jpk), crs_surfv_wgt(jpi_crs,jpj_crs,jpk) , &
         &      crs_surfw_wgt(jpi_crs,jpj_crs,jpk), crs_volt_wgt(jpi_crs,jpj_crs,jpk) , STAT=ierr(8))


      ALLOCATE( mbathy_crs(jpi_crs,jpj_crs), mbkt_crs(jpi_crs,jpj_crs) , &
         &      mbku_crs(jpi_crs,jpj_crs)  , mbkv_crs(jpi_crs,jpj_crs) , STAT=ierr(9))

      ALLOCATE( gdept_crs(jpi_crs,jpj_crs,jpk), gdepu_crs(jpi_crs,jpj_crs,jpk) , &
         &      gdepv_crs(jpi_crs,jpj_crs,jpk), gdepw_crs(jpi_crs,jpj_crs,jpk) , STAT=ierr(10) )



      !ALLOCATE( nmln_crs(jpi_crs,jpj_crs) , hmld_crs(jpi_crs,jpj_crs) , &
      !   &      hmlp_crs(jpi_crs,jpj_crs) , hmlpt_crs(jpi_crs,jpj_crs) , STAT=ierr(14) )

      !! LOLO: only mpp stuff: not needed!
      ! ALLOCATE( nimppt_crs(jpnij) , nlcit_crs(jpnij) , nldit_crs(jpnij) , nleit_crs(jpnij), &
      !    &  nimppt_full(jpnij) , nlcit_full(jpnij) , nldit_full(jpnij) , nleit_full(jpnij),   &
      !    njmppt_crs(jpnij) , nlcjt_crs(jpnij) , nldjt_crs(jpnij) , nlejt_crs(jpnij), &
      !    &  njmppt_full(jpnij) , nlcjt_full(jpnij) , nldjt_full(jpnij) , nlejt_full(jpnij)  , STAT=ierr(15) )
      

      crs_dom_alloc = MAXVAL(ierr)

   END FUNCTION crs_dom_alloc

   INTEGER FUNCTION crs_dom_alloc2()
      !!-------------------------------------------------------------------
      !!                     *** FUNCTION crs_dom_alloc ***
      !!  ** Purpose :   Allocate public crs arrays
      !!-------------------------------------------------------------------
      !! Local variables
      INTEGER, DIMENSION(1) :: ierr

      ierr(:) = 0

      ALLOCATE( mjs_crs(nlej_crs) , mje_crs(nlej_crs), mis_crs(nlei_crs) , mie_crs(nlei_crs), STAT=ierr(1) )
      crs_dom_alloc2 = MAXVAL(ierr)

   END FUNCTION crs_dom_alloc2

   SUBROUTINE dom_grid_glo
      !!--------------------------------------------------------------------
      !!                       ***  MODULE dom_grid_glo  ***
      !!
      !! ** Purpose : +Return back to parent grid domain
      !!---------------------------------------------------------------------

      !                         Return to parent grid domain

   END SUBROUTINE dom_grid_glo

   SUBROUTINE dom_grid_crs
      !!--------------------------------------------------------------------
      !!                       ***  MODULE dom_grid_crs  ***
      !!
      !! ** Purpose :  Save the parent grid information & Switch to coarse grid domain
      !!---------------------------------------------------------------------

      !
      !                        Switch to coarse grid domain
      jpi    = jpi_crs
      jpj    = jpj_crs
      jpim1  = jpi_crsm1
      jpjm1  = jpj_crsm1
      nperio = nperio_crs

      npolj  = npolj_crs
      jpiglo = jpiglo_crs
      jpjglo = jpjglo_crs


      !! LOLO: mpp stuff:
      !nlci   = nlci_crs
      !nlcj   = nlcj_crs
      !nldi   = nldi_crs
      !nlei   = nlei_crs
      !nlej   = nlej_crs
      !nldj   = nldj_crs
      !nimpp  = nimpp_crs
      !njmpp  = njmpp_crs

      !nlcit(:)  = nlcit_crs(:)
      !nldit(:)  = nldit_crs(:)
      !nleit(:)  = nleit_crs(:)
      !nimppt(:) = nimppt_crs(:)
      !nlcjt(:)  = nlcjt_crs(:)
      !nldjt(:)  = nldjt_crs(:)
      !nlejt(:)  = nlejt_crs(:)
      !njmppt(:) = njmppt_crs(:)
      !LOLO.
      
      !
   END SUBROUTINE dom_grid_crs


   !!======================================================================

END MODULE mod_crs_def
