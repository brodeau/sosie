MODULE crs
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
   !PUBLIC dom_grid_glo
   PUBLIC dom_grid_crs

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

   INTEGER ::  nfsloop_full,nfeloop_full
   INTEGER ::  nfsloop_crs ,nfeloop_crs

   INTEGER, DIMENSION(:), ALLOCATABLE :: mis_crs, mie_crs, mis2_crs, mie2_crs  ! starting and ending i-indices of parent subset
   INTEGER, DIMENSION(:), ALLOCATABLE :: mjs_crs, mje_crs, mjs2_crs, mje2_crs  ! starting and ending  j-indices of parent subset
   INTEGER, DIMENSION(:), ALLOCATABLE :: mjg_crs, mig_crs
   INTEGER, DIMENSION(:), ALLOCATABLE :: mi0_crs, mi1_crs, mj0_crs, mj1_crs
   INTEGER                            :: mxbinctr, mybinctr                    ! central point in grid box
   INTEGER, DIMENSION(:), ALLOCATABLE :: nlcit_crs, nlcit_full                 ! dimensions of every subdomain
   INTEGER, DIMENSION(:), ALLOCATABLE :: nldit_crs, nldit_full                 ! first, last indoor index for each i-domain
   INTEGER, DIMENSION(:), ALLOCATABLE :: nleit_crs, nleit_full                 ! first, last indoor index for each j-domain
   INTEGER, DIMENSION(:), ALLOCATABLE :: nimppt_crs, nimppt_full               ! first, last indoor index for each j-domain
   INTEGER, DIMENSION(:), ALLOCATABLE :: nlcjt_crs, nlcjt_full                 ! dimensions of every subdomain
   INTEGER, DIMENSION(:), ALLOCATABLE :: nldjt_crs, nldjt_full                 ! first, last indoor index for each i-domain
   INTEGER, DIMENSION(:), ALLOCATABLE :: nlejt_crs, nlejt_full                 ! first, last indoor index for each j-domain
   INTEGER, DIMENSION(:), ALLOCATABLE :: njmppt_crs, njmppt_full               ! first, last indoor index for each j-domain

   INTEGER, DIMENSION(:,:), ALLOCATABLE ::   nfiimpp_full
   INTEGER, DIMENSION(:,:), ALLOCATABLE ::   nfiimpp_crs

   ! Masks
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE,SAVE :: tmask_crs, umask_crs, vmask_crs, fmask_crs
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE,SAVE :: tmask_i_crs, rnfmsk_crs, tpol_crs, fpol_crs

   ! Scale factors
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1t_crs, e2t_crs, e1e2t_crs ! horizontal scale factors grid type T
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1u_crs, e2u_crs            ! horizontal scale factors grid type U
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1v_crs, e2v_crs            ! horizontal scale factors grid type V
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE :: e1f_crs, e2f_crs            ! horizontal scale factors grid type F

   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE :: ht_0_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_0_crs, e3u_0_crs, e3v_0_crs, e3f_0_crs, e3w_0_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_max_0_crs, e3u_max_0_crs, e3v_max_0_crs, e3f_max_0_crs, e3w_max_0_crs

   !#if defined key_vvl
   !   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_b_crs, e3u_b_crs, e3v_b_crs, e3f_b_crs, e3w_b_crs
   !   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_n_crs, e3u_n_crs, e3v_n_crs, e3f_n_crs, e3w_n_crs
   !   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_a_crs, e3u_a_crs, e3v_a_crs, e3f_a_crs, e3w_a_crs
   !   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e3t_max_n_crs, e3u_max_n_crs, e3v_max_n_crs, e3f_max_n_crs, e3w_max_n_crs
   !#endif

   ! Surface
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e1e2w_crs, e2e3u_crs, e1e3v_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: e1e2w_msk, e2e3u_msk, e1e3v_msk

   ! Coordinates
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE,SAVE :: gphit_crs, glamt_crs, gphif_crs, glamf_crs
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE,SAVE :: gphiu_crs, glamu_crs, gphiv_crs, glamv_crs
   REAL(wp), DIMENSION(:,:),   ALLOCATABLE,SAVE :: ff_crs
   INTEGER,  DIMENSION(:,:),   ALLOCATABLE,SAVE :: mbathy_crs, mbkt_crs, mbku_crs, mbkv_crs

   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE,SAVE :: gdept_0_crs, gdepu_0_crs, gdepv_0_crs, gdepw_0_crs
   !#if defined key_vvl
   !   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE,SAVE :: gdept_n_crs, gdepu_n_crs, gdepv_n_crs, gdepw_n_crs
   !#endif

   ! Weights
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: facsurfv, facsurfu, facvol_t, facvol_w
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: ocean_volume_crs_t, ocean_volume_crs_w, bt_crs, r1_bt_crs
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: crs_surfu_wgt, crs_surfv_wgt, crs_surfw_wgt, crs_volt_wgt

   ! Namelist
   INTEGER           :: nn_factx   = 3       !: reduction factor of x-dimension of the parent grid
   INTEGER           :: nn_facty   = 3       !: reduction factor of y-dimension of the parent grid
   INTEGER           :: nn_msh_crs = 1       !: Organization of mesh mask output
   !: 0 = no mesh mask output
   !: 1 = unified mesh mask output
   !: 2 = 2 separate mesh mask output
   !: 3 = 3 separate mesh mask output
   INTEGER           :: nn_crs_kz    =    0       !: type of Kz coarsening ( =0->VOL ; =1->MAX ; =2->MIN)
   LOGICAL           :: ln_crs_wn    = .FALSE.    !: coarsening wn or computation using horizontal divergence
   LOGICAL, PUBLIC   :: ln_crs_top   = .FALSE.    !:coarsening online for the bio
   !

   ! Grid reduction factors
   REAL(wp)     ::  rfactx_r                !: inverse of x-dim reduction factor
   REAL(wp)     ::  rfacty_r                !: inverse of y-dim reduction factor
   REAL(wp)     ::  rfactxy
   INTEGER      :: nrestx, nresty           !: for determining odd or even reduction factor
   INTEGER, DIMENSION(:), ALLOCATABLE      :: nfactx,nfacty


   ! Physical and dynamical ocean fields for output or passing to TOP, time-mean fields
   REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE      :: tsb_crs,tsn_crs,tsa_crs,rab_crs_n
   REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE      :: un_crs, vn_crs, wn_crs
   REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE      :: ub_crs, vb_crs
   REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE      :: hdivb_crs , hdivn_crs
   REAL(wp), DIMENSION(:,:)    , ALLOCATABLE      :: sshb_crs, sshn_crs , ssha_crs
   REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE      :: rhop_crs,rhd_crs,rn2_crs,rb2_crs
   REAL(wp), DIMENSION(:,:)    , ALLOCATABLE      :: gru_crs, grv_crs
   REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE      :: gtsu_crs, gtsv_crs
   !
   ! Surface fluxes to pass to TOP
   REAL(wp), PUBLIC, DIMENSION(:,:)  , ALLOCATABLE :: qsr_crs, fr_i_crs, wndm_crs
   REAL(wp), PUBLIC, DIMENSION(:,:)  , ALLOCATABLE :: emp_crs, emp_b_crs, sfx_crs
   REAL(wp), PUBLIC, DIMENSION(:,:)  , ALLOCATABLE :: fmmflx_crs
   REAL(wp), PUBLIC, DIMENSION(:,:)  , ALLOCATABLE :: utau_crs, vtau_crs, taum_crs
   REAL(wp), PUBLIC, DIMENSION(:,:)  , ALLOCATABLE :: rnf_crs,rnf_b_crs,h_rnf_crs
   INTEGER , PUBLIC, DIMENSION(:,:)  , ALLOCATABLE :: nk_rnf_crs
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: trc_i_crs,trc_o_crs
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: sbc_trc_crs, sbc_trc_b_crs

   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: uslp_crs, wslpi_crs          !: i_slope at U- and W-points
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: vslp_crs, wslpj_crs          !: j-slope at V- and W-points

   ! Horizontal diffusion
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   r_fact_lap_crs

   ! Vertical diffusion
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)  ::  avt_crs           !: vert. diffusivity coef. [m2/s] at w-point for temp
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:)  ::  en_crs            !: vert. diffusivity coef. [m2/s] at w-point for temp
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)    ::  avtb_2d_crs       !: vert. diffusivity coef. [m2/s] at w-point for temp

   ! Mixing and Mixed Layer Depth
   INTEGER,  PUBLIC, DIMENSION(:,:) , ALLOCATABLE ::  nmln_crs
   REAL(wp), PUBLIC, DIMENSION(:,:) , ALLOCATABLE :: hmlp_crs , hmlpt_crs , hmld_crs

   ! Direction of lateral diffusion


   !! $Id: crs.F90 7795 2017-03-15 08:04:30Z cbricaud $
CONTAINS

   INTEGER FUNCTION crs_dom_alloc()
      !!-------------------------------------------------------------------
      !!                     *** FUNCTION crs_dom_alloc ***
      !!  ** Purpose :   Allocate public crs arrays
      !!-------------------------------------------------------------------
      !! Local variables
      INTEGER, DIMENSION(15) :: ierr

      ierr(:) = 0

      ! Set up bins for coarse grid, horizontal only.
      ALLOCATE( mis2_crs(jpiglo_crs), mie2_crs(jpiglo_crs),  &
         &       mjs2_crs(jpjglo_crs), mje2_crs(jpjglo_crs),  &
         &       mi0_crs (jpiglo_crs), mi1_crs (jpiglo_crs),  &
         &       mj0_crs (jpjglo_crs), mj1_crs (jpjglo_crs),  &
         &       mig_crs (jpi_crs)   , mjg_crs (jpj_crs)   ,  &
         &       mis_crs (jpi_crs)   , mie_crs (jpi_crs)   ,  &
         &       mjs_crs (jpj_crs)   , mje_crs (jpj_crs)   ,  &
         &       nfactx  (jpi_crs)   , nfacty  (jpj_crs)   ,  &
         &       nimppt_crs(jpnij) , nlcit_crs(jpnij) , nldit_crs(jpnij) , nleit_crs(jpnij) , &
         &       nimppt_full(jpnij), nlcit_full(jpnij), nldit_full(jpnij), nleit_full(jpnij), &
         &       njmppt_crs(jpnij) , nlcjt_crs(jpnij) , nldjt_crs(jpnij) , nlejt_crs(jpnij) , &
         &       njmppt_full(jpnij), nlcjt_full(jpnij), nldjt_full(jpnij), nlejt_full(jpnij), &
         &       nfiimpp_full(jpni,jpnj) , nfiimpp_crs(jpni,jpnj) , STAT=ierr(1) )

      ! Set up Mask and Mesh
      PRINT *, ' ### crs_dom_alloc@crs => allocating tmask_crs:',jpi_crs,jpj_crs,jpk
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

      ALLOCATE( e3t_0_crs(jpi_crs,jpj_crs,jpk)    , e3w_0_crs(jpi_crs,jpj_crs,jpk)    , &
         &      e3u_0_crs(jpi_crs,jpj_crs,jpk)    , e3v_0_crs(jpi_crs,jpj_crs,jpk)    , &
         &           ht_0_crs(jpi_crs,jpj_crs),                                     &
         !#if defined key_vvl
         !         &      e3t_b_crs(jpi_crs,jpj_crs,jpk)    , e3w_b_crs(jpi_crs,jpj_crs,jpk)    , &
         !         &      e3u_b_crs(jpi_crs,jpj_crs,jpk)    , e3v_b_crs(jpi_crs,jpj_crs,jpk)    , &
         !         &      e3t_n_crs(jpi_crs,jpj_crs,jpk)    , e3w_n_crs(jpi_crs,jpj_crs,jpk)    , &
         !         &      e3u_n_crs(jpi_crs,jpj_crs,jpk)    , e3v_n_crs(jpi_crs,jpj_crs,jpk)    , &
         !         &      e3t_a_crs(jpi_crs,jpj_crs,jpk)    , e3w_a_crs(jpi_crs,jpj_crs,jpk)    , &
         !         &      e3u_a_crs(jpi_crs,jpj_crs,jpk)    , e3v_a_crs(jpi_crs,jpj_crs,jpk)    , &
         !#endif
         &      e1e2w_msk(jpi_crs,jpj_crs,jpk)  , &
         &      e2e3u_msk(jpi_crs,jpj_crs,jpk)  , e1e3v_msk(jpi_crs,jpj_crs,jpk)  , &
         &      e1e2w_crs(jpi_crs,jpj_crs,jpk)  , e2e3u_crs(jpi_crs,jpj_crs,jpk)  , &
         &      e1e3v_crs(jpi_crs,jpj_crs,jpk)  , &
         &      e3t_max_0_crs(jpi_crs,jpj_crs,jpk), e3w_max_0_crs(jpi_crs,jpj_crs,jpk) , &
         &      e3u_max_0_crs(jpi_crs,jpj_crs,jpk), e3v_max_0_crs(jpi_crs,jpj_crs,jpk) , &
         !#if defined key_vvl
         !         &      e3t_max_n_crs(jpi_crs,jpj_crs,jpk), e3w_max_n_crs(jpi_crs,jpj_crs,jpk) , &
         !         &      e3u_max_n_crs(jpi_crs,jpj_crs,jpk), e3v_max_n_crs(jpi_crs,jpj_crs,jpk) , &
         !#endif
         &      STAT=ierr(6))


      !ALLOCATE( facsurfv(jpi_crs,jpj_crs,jpk), facsurfu(jpi_crs,jpj_crs,jpk) , &
      !   &      facvol_t(jpi_crs,jpj_crs,jpk), facvol_w(jpi_crs,jpj_crs,jpk) , &
      !   &      ocean_volume_crs_t(jpi_crs,jpj_crs,jpk), ocean_volume_crs_w(jpi_crs,jpj_crs,jpk), &
      !   &      bt_crs(jpi_crs,jpj_crs,jpk)  , r1_bt_crs(jpi_crs,jpj_crs,jpk) , STAT=ierr(7))


      !ALLOCATE( crs_surfu_wgt(jpi_crs,jpj_crs,jpk), crs_surfv_wgt(jpi_crs,jpj_crs,jpk) , &
      !   &      crs_surfw_wgt(jpi_crs,jpj_crs,jpk), crs_volt_wgt(jpi_crs,jpj_crs,jpk) , STAT=ierr(8))


      ALLOCATE( mbathy_crs(jpi_crs,jpj_crs), mbkt_crs(jpi_crs,jpj_crs) , &
         &      mbku_crs(jpi_crs,jpj_crs)  , mbkv_crs(jpi_crs,jpj_crs) , STAT=ierr(9))

      ALLOCATE( gdept_0_crs(jpi_crs,jpj_crs,jpk), gdepu_0_crs(jpi_crs,jpj_crs,jpk) , &
         &      gdepv_0_crs(jpi_crs,jpj_crs,jpk), gdepw_0_crs(jpi_crs,jpj_crs,jpk) , &
         !#if defined key_vvl
         !         &      gdept_n_crs(jpi_crs,jpj_crs,jpk), gdepu_n_crs(jpi_crs,jpj_crs,jpk) , &
         !         &      gdepv_n_crs(jpi_crs,jpj_crs,jpk), gdepw_n_crs(jpi_crs,jpj_crs,jpk) , &
         !#endif
         & STAT=ierr(10))


      !      ALLOCATE( ub_crs(jpi_crs,jpj_crs,jpk) , vb_crs(jpi_crs,jpj_crs,jpk) , &
      !         &      un_crs(jpi_crs,jpj_crs,jpk) , vn_crs(jpi_crs,jpj_crs,jpk)    ,  wn_crs(jpi_crs,jpj_crs,jpk) ,&
      !         &      hdivb_crs(jpi_crs,jpj_crs,jpk) , hdivn_crs(jpi_crs,jpj_crs,jpk) , &
      !         &      rhop_crs(jpi_crs,jpj_crs,jpk)  , &
      !         &      rb2_crs(jpi_crs,jpj_crs,jpk) ,rn2_crs(jpi_crs,jpj_crs,jpk) , &
      !         &      rhd_crs(jpi_crs,jpj_crs,jpk)   , rab_crs_n(jpi_crs,jpj_crs,jpk,jpts) , &
      !         &      avtb_2d_crs(jpi_crs,jpj_crs), &
      !         &      gtsu_crs(jpi_crs,jpj_crs,jpts) ,gtsv_crs(jpi_crs,jpj_crs,jpts) , &
      !         gru_crs(jpi_crs,jpj_crs) ,grv_crs(jpi_crs,jpj_crs) , STAT=ierr(11))

      !      ALLOCATE( sshb_crs(jpi_crs,jpj_crs), sshn_crs(jpi_crs,jpj_crs),  ssha_crs(jpi_crs,jpj_crs), &
      !         &     qsr_crs(jpi_crs ,jpj_crs), wndm_crs(jpi_crs,jpj_crs), utau_crs(jpi_crs,jpj_crs) , &
      !         &     vtau_crs(jpi_crs,jpj_crs), taum_crs(jpi_crs,jpj_crs),  &
      !         &     rnf_crs (jpi_crs,jpj_crs), rnf_b_crs(jpi_crs ,jpj_crs), nk_rnf_crs(jpi_crs ,jpj_crs), h_rnf_crs(jpi_crs ,jpj_crs), &
      !         &     emp_crs (jpi_crs,jpj_crs), emp_b_crs(jpi_crs,jpj_crs), &
      !         &     sbc_trc_crs (jpi_crs,jpj_crs,jpts), sbc_trc_b_crs(jpi_crs,jpj_crs,jpts), &
      !         &     trc_i_crs (jpi_crs,jpj_crs,jpts), trc_o_crs(jpi_crs,jpj_crs,jpts), &
      !         &     fr_i_crs(jpi_crs,jpj_crs), sfx_crs(jpi_crs ,jpj_crs), fmmflx_crs(jpi_crs ,jpj_crs),  STAT=ierr(12)  )
      !
      !#if defined key_traldf_c3d
      !      ALLOCATE( ahtt_crs(jpi_crs,jpj_crs,jpk) , ahtu_crs(jpi_crs,jpj_crs,jpk) , &
      !         & ahtv_crs(jpi_crs,jpj_crs,jpk) , ahtw_crs(jpi_crs,jpj_crs,jpk) , &
      !#elif defined key_traldf_c2d
      !         ALLOCATE( ahtt_crs(jpi_crs,jpj_crs    ) , ahtu_crs(jpi_crs,jpj_crs    ) , &
      !         & ahtv_crs(jpi_crs,jpj_crs    ) , ahtw_crs(jpi_crs,jpj_crs    ) , &
      !#elif defined key_traldf_c1d
      !         ALLOCATE( ahtt_crs(        jpk) , ahtu_crs(        jpk) , ahtv_crs(        jpk) , ahtw_crs(        jpk) , &
      !#endif
      !         & r_fact_lap_crs(jpi_crs,jpj_crs,jpk) , STAT=ierr(13) )

      !      ALLOCATE( tsb_crs(jpi_crs,jpj_crs,jpk,jpts), tsn_crs(jpi_crs,jpj_crs,jpk,jpts), tsa_crs(jpi_crs,jpj_crs,jpk,jpts),  &
      !         en_crs(jpi_crs,jpj_crs,jpk),   avt_crs(jpi_crs,jpj_crs,jpk),    &
      !# if defined key_zdfddm
      !         &      avs_crs(jpi_crs,jpj_crs,jpk),    &
      !# endif
      !         &      STAT=ierr(14) )!
      !
      !      ALLOCATE( nmln_crs(jpi_crs,jpj_crs) , hmld_crs(jpi_crs,jpj_crs) , &
      !         &      hmlp_crs(jpi_crs,jpj_crs) , hmlpt_crs(jpi_crs,jpj_crs) , STAT=ierr(15) )

      crs_dom_alloc = MAXVAL(ierr)

   END FUNCTION crs_dom_alloc

   SUBROUTINE dom_grid_glo
      !!--------------------------------------------------------------------
      !!                       ***  MODULE dom_grid_glo  ***
      !!
      !! ** Purpose : +Return back to parent grid domain
      !!---------------------------------------------------------------------

      !                         Return to parent grid domain
      jpi    = jpi_full
      jpj    = jpj_full
      jpim1  = jpim1_full
      jpjm1  = jpjm1_full

      npolj  = npolj_full
      jpiglo = jpiglo_full
      jpjglo = jpjglo_full

      nlci   = nlci_full
      nlcj   = nlcj_full
      nldi   = nldi_full
      nldj   = nldj_full
      nlei   = nlei_full
      nlej   = nlej_full
      nimpp  = nimpp_full
      njmpp  = njmpp_full

      !nlcit(:)  = nlcit_full(:)
      !nldit(:)  = nldit_full(:)
      !nleit(:)  = nleit_full(:)
      !nimppt(:) = nimppt_full(:)
      !nlcjt(:)  = nlcjt_full(:)
      !nldjt(:)  = nldjt_full(:)
      !nlejt(:)  = nlejt_full(:)
      !njmppt(:) = njmppt_full(:)

      !nfsloop = nfsloop_full
      !nfeloop = nfeloop_full

      nfiimpp(:,:) = nfiimpp_full(:,:)

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

      npolj  = npolj_crs
      jpiglo = jpiglo_crs
      jpjglo = jpjglo_crs


      nlci   = nlci_crs
      nlcj   = nlcj_crs
      nldi   = nldi_crs
      nlei   = nlei_crs
      nlej   = nlej_crs
      nldj   = nldj_crs
      nimpp  = nimpp_crs
      njmpp  = njmpp_crs

      !nlcit(:)  = nlcit_crs(:)
      !nldit(:)  = nldit_crs(:)
      !nleit(:)  = nleit_crs(:)
      !nimppt(:) = nimppt_crs(:)
      !nlcjt(:)  = nlcjt_crs(:)
      !nldjt(:)  = nldjt_crs(:)
      !nlejt(:)  = nlejt_crs(:)
      !njmppt(:) = njmppt_crs(:)

      !nfsloop = nfsloop_crs
      !nfeloop = nfeloop_crs

      nfiimpp(:,:) = nfiimpp_crs(:,:)

      !
   END SUBROUTINE dom_grid_crs


   !!======================================================================

END MODULE crs
