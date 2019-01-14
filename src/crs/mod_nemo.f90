MODULE mod_nemo

   IMPLICIT NONE
   PUBLIC

   !INTEGER, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   !INTEGER, PARAMETER ::   wp = dp                              !: working precision
   INTEGER, PARAMETER ::   wp = 4 !LOLO

   INTEGER ::   numout          =    6      !: logical unit for output print; Set to stdout to ensure any early

   LOGICAL, PARAMETER :: lk_mpp = .TRUE. ! NO MPP!!!

   REAL(wp), PARAMETER ::   rpi = 3.141592653589793_wp             !: pi
   REAL(wp), PARAMETER ::   rad = 3.141592653589793_wp / 180._wp   !: conversion from degre into radian
   REAL(wp), PARAMETER :: omega = 7.292115083046062E-005        !: earth rotation parameter           [s-1]

   INTEGER, PARAMETER ::   jpreci = 1   !: number of columns for overlap
   INTEGER, PARAMETER ::   jprecj = 1   !: number of rows    for overlap

   ! global or zoom domain size                      !!! * computational domain *
   INTEGER       ::    jpiglo           !: 1st dimension of global domain --> i
   INTEGER       ::    jpjglo           !: 2nd    -                  -    --> j

   INTEGER            ::   jpni         !: number of processors following i
   INTEGER            ::   jpnj         !: number of processors following j
   INTEGER            ::   jpnij        !: nb of local domain = nb of processors ( <= jpni x jpnj )


   INTEGER ::   noso, nono        !: east, west, south and north directions

   INTEGER ::   nlci, nldi, nlei  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER ::   nlcj, nldj, nlej  !: i-dimensions of the local subdomain and its first and last indoor indices
   INTEGER ::   noea, nowe        !: index of the local neighboring processors in

   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) ::   mig        !: local  ==> global domain i-index
   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) ::   mjg        !: local  ==> global domain j-index

   INTEGER ::   nimpp, njmpp      !: i- & j-indexes for mpp-subdomain left bottom

   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) ::   nlcit , nlcjt    !: dimensions of every subdomain
   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) ::   nldit , nldjt    !: first, last indoor index for each i-domain
   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) ::   nleit , nlejt    !: first, last indoor index for each j-domain
   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:,:) :: nfiimpp, nfipproc, nfilcit


   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) ::   ibonit, ibonjt   !: i-, j- processor neighbour existence
   INTEGER, ALLOCATABLE, SAVE, DIMENSION(:) ::   nimppt, njmppt   !: i-, j-indexes for each processor
   !! eNATL4:
   !&    jpiglo=559, &
   !   &    jpjglo=313


   INTEGER ::   nproc             !: number for local processor
   INTEGER ::   narea             !: number for local area

   LOGICAL, PARAMETER :: lwp = .true.

   INTEGER ::  jpi, jpj, jpim1, jpjm1, jpk, jpkm1

   INTEGER :: nperio, jperio
   INTEGER :: npolj

   !INTEGER :: & ! needed even if no mpp...
   !   &   nimpp_crs ,  &
   !   &   njmpp_crs  , &
   !   &   nlci_crs   , &
   !   &   nlcj_crs   , &
   !   &   nldi_crs   , &
   !   &   nldj_crs   , &
   !   &   nlei_crs   , &
   !   &   nlej_crs

   INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: tmask, umask, vmask, fmask

   REAL(4),    DIMENSION(:,:,:), ALLOCATABLE :: e3t_0, e3u_0, e3v_0, e3w_0 !LOLO
   REAL(wp),   DIMENSION(:,:,:), ALLOCATABLE :: e3t, e3u, e3v
   REAL(wp),   DIMENSION(:,:), ALLOCATABLE :: glamt, gphit, glamu, gphiu, glamv, gphiv, glamf, gphif
   REAL(wp),   DIMENSION(:,:), ALLOCATABLE :: e1t, e2t, e1u, e2u, e1v, e2v, e1f, e2f, e1e2t


   ! zoom starting position
   INTEGER       ::   jpizoom          !: left bottom (i,j) indices of the zoom
   INTEGER       ::   jpjzoom          !: in data domain indices

   INTEGER                                  ::   nsndto, nfsloop, nfeloop

   INTEGER :: jphgr_msh

END MODULE mod_nemo
