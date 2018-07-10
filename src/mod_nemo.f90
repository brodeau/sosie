MODULE mod_nemo

   IMPLICIT NONE
   PUBLIC

   INTEGER, PUBLIC, PARAMETER ::   dp = SELECTED_REAL_KIND(12,307)   !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER ::   wp = dp                              !: working precision

   INTEGER ::   numout          =    6      !: logical unit for output print; Set to stdout to ensure any early
   
   INTEGER, PUBLIC, PARAMETER ::   jpreci = 1   !: number of columns for overlap
   INTEGER, PUBLIC, PARAMETER ::   jprecj = 1   !: number of rows    for overlap

   ! global or zoom domain size                      !!! * computational domain *
   INTEGER       ::    jpiglo           !: 1st dimension of global domain --> i
   INTEGER       ::    jpjglo           !: 2nd    -                  -    --> j

   INTEGER, PUBLIC            ::   jpni         !: number of processors following i
   INTEGER, PUBLIC            ::   jpnj         !: number of processors following j
   INTEGER, PUBLIC            ::   jpnij        !: nb of local domain = nb of processors ( <= jpni x jpnj )


   INTEGER, PUBLIC ::   noso, nono        !: east, west, south and north directions

   INTEGER, PUBLIC ::   nlci, nldi, nlei  !: i-dimensions of the local subdomain and its first and last indoor indices    
   INTEGER, PUBLIC ::   nlcj, nldj, nlej  !: i-dimensions of the local subdomain and its first and last indoor indices 

   INTEGER, PUBLIC ::   nimpp, njmpp      !: i- & j-indexes for mpp-subdomain left bottom

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nlcit , nlcjt    !: dimensions of every subdomain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nldit , nldjt    !: first, last indoor index for each i-domain
   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nleit , nlejt    !: first, last indoor index for each j-domain

   INTEGER, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:) ::   nimppt, njmppt   !: i-, j-indexes for each processor
   !! eNATL4:
   !&    jpiglo=559, &
   !   &    jpjglo=313


   INTEGER, PUBLIC ::   nproc             !: number for local processor
   INTEGER, PUBLIC ::   narea             !: number for local area 
   
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




END MODULE mod_nemo
