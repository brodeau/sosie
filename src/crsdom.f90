MODULE crsdom
   !!===================================================================
   !!                  ***  crs.F90 ***
   !!  Purpose: Interface for calculating quantities from a
   !!           higher-resolution grid for the coarse grid.
   !!
   !!  Method:

   !!  References:  Aumont, O., J.C. Orr, D. Jamous, P. Monfray
   !!               O. Marti and G. Madec, 1998. A degradation
   !!               approach to accelerate simulations to steady-state
   !!               in a 3-D tracer transport model of the global ocean.
   !!               Climate Dynamics, 14:101-116.
   !!  History:
   !!       Original.   May 2012.  (J. Simeon, C. Calone, G. Madec, C. Ethe)
   !!===================================================================

   !USE dom_oce        ! ocean space and time domain and to get jperio
   USE mod_nemo
   USE crs             ! domain for coarse grid
   !USE in_out_manager
   !USE par_kind
   !USE crslbclnk
   !USE lib_mpp
   !USE ieee_arithmetic

   IMPLICIT NONE

   PRIVATE

   PUBLIC crs_dom_ope
   PUBLIC crs_dom_e3, crs_dom_sfc, crs_dom_msk, crs_dom_hgr, crs_dom_coordinates
   PUBLIC crs_dom_facvol, crs_dom_def, crs_dom_bat

   INTERFACE crs_dom_ope
      MODULE PROCEDURE crs_dom_ope_3d, crs_dom_ope_2d
   END INTERFACE crs_dom_ope

   REAL(wp),PUBLIC :: r_inf = 1e+7 !cbr 1e+36

   !! Substitutions
!#  include "domzgr_substitute.h90"

CONTAINS


   SUBROUTINE crs_dom_msk
      !!===================================================================
      !
      !
      !
      !!===================================================================
      INTEGER  ::  ji, jj, jk                   ! dummy loop indices
      INTEGER  ::  ijis,ijie,ijjs,ijje
      REAL(wp) ::  zmask
      !!-------------------------------------------------------------------

      ! Initialize
      tmask_crs(:,:,:) = 0.0
      vmask_crs(:,:,:) = 0.0
      umask_crs(:,:,:) = 0.0
      fmask_crs(:,:,:) = 0.0
      !
      DO jk = 1, jpkm1
         DO ji = nldi_crs, nlei_crs

            ijis = mis_crs(ji)
            ijie = mie_crs(ji)

            DO jj = nldj_crs, nlej_crs

               ijjs = mjs_crs(jj)
               ijje = mje_crs(jj)

               zmask = 0.0
               zmask = SUM( tmask(ijis:ijie,ijjs:ijje,jk) )
               IF ( zmask > 0.0 ) tmask_crs(ji,jj,jk) = 1.0

               zmask = 0.0
               zmask = SUM( vmask(ijis:ijie,ijje     ,jk) )
               IF ( zmask > 0.0 ) vmask_crs(ji,jj,jk) = 1.0

               zmask = 0.0
               zmask = SUM( umask(ijie     ,ijjs:ijje,jk) )
               IF ( zmask > 0.0 ) umask_crs(ji,jj,jk) = 1.0

               fmask_crs(ji,jj,jk) = fmask(ijie,ijje,jk)

            ENDDO
         ENDDO
      ENDDO

      !CALL crs_lbc_lnk( tmask_crs, 'T', 1.0 )
      !CALL crs_lbc_lnk( vmask_crs, 'V', 1.0 )
      !CALL crs_lbc_lnk( umask_crs, 'U', 1.0 )
      !CALL crs_lbc_lnk( fmask_crs, 'F', 1.0 )
      !
   END SUBROUTINE crs_dom_msk

   SUBROUTINE crs_dom_coordinates( p_gphi, p_glam, cd_type, p_gphi_crs, p_glam_crs )
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crs_coordinates ***
      !! ** Purpose :  Determine the coordinates for the coarse grid
      !!
      !! ** Method  :  From the parent grid subset, search for the central
      !!               point.  For an odd-numbered reduction factor,
      !!               the coordinate will be that of the central T-cell.
      !!               For an even-numbered reduction factor, of a non-square
      !!               coarse grid box, the coordinate will be that of
      !!               the east or north face or more likely.  For a square
      !!               coarse grid box, the coordinate will be that of
      !!               the central f-corner.
      !!
      !! ** Input   :  p_gphi = parent grid gphi[t|u|v|f]
      !!               p_glam = parent grid glam[t|u|v|f]
      !!               cd_type  = grid type (T,U,V,F)
      !! ** Output  :  p_gphi_crs = coarse grid gphi[t|u|v|f]
      !!               p_glam_crs = coarse grid glam[t|u|v|f]
      !!
      !! History. 1 Jun.
      !!----------------------------------------------------------------
      !! Arguments
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(in)  :: p_gphi  ! Parent grid latitude
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(in)  :: p_glam  ! Parent grid longitude
      CHARACTER(len=1),                     INTENT(in)  :: cd_type   ! grid type (T,U,V,F)
      REAL(wp), DIMENSION(jpi_crs,jpj_crs), INTENT(out) :: p_gphi_crs  ! Coarse grid latitude
      REAL(wp), DIMENSION(jpi_crs,jpj_crs), INTENT(out) :: p_glam_crs  ! Coarse grid longitude

      !! Local variables
      INTEGER :: ji, jj, jk                   ! dummy loop indices
      INTEGER :: iji, ijj
      INTEGER  :: ir,jr
      !!----------------------------------------------------------------
      p_gphi_crs(:,:)=0._wp
      p_glam_crs(:,:)=0._wp


      SELECT CASE ( cd_type )
      CASE ( 'T' )
         DO jj =  nldj_crs, nlej_crs
            ijj = mjs_crs(jj) + + INT(0.5*nfacty(jj))
            DO ji = nldi_crs, nlei_crs
               iji = mis_crs(ji) + INT(0.5*nfactx(ji))
               p_gphi_crs(ji,jj) = p_gphi(iji,ijj)
               p_glam_crs(ji,jj) = p_glam(iji,ijj)
            ENDDO
         ENDDO
      CASE ( 'U' )
         DO jj =  nldj_crs, nlej_crs
            ijj = mjs_crs(jj) + INT(0.5*nfacty(jj))
            DO ji = nldi_crs, nlei_crs
               iji = mie_crs(ji)
               p_gphi_crs(ji,jj) = p_gphi(iji,ijj)
               p_glam_crs(ji,jj) = p_glam(iji,ijj)

            ENDDO
         ENDDO
      CASE ( 'V' )
         DO jj =  nldj_crs, nlej_crs
            ijj = mje_crs(jj)
            DO ji = nldi_crs, nlei_crs
               iji = mis_crs(ji) + INT(0.5*nfactx(ji))
               p_gphi_crs(ji,jj) = p_gphi(iji,ijj)
               p_glam_crs(ji,jj) = p_glam(iji,ijj)
            ENDDO
         ENDDO
      CASE ( 'F' )
         DO jj =  nldj_crs, nlej_crs
            ijj = mje_crs(jj)
            DO ji = nldi_crs, nlei_crs
               iji = mie_crs(ji)
               p_gphi_crs(ji,jj) = p_gphi(iji,ijj)
               p_glam_crs(ji,jj) = p_glam(iji,ijj)
            ENDDO
         ENDDO
      END SELECT
      WRITE(narea+1000-1,*)"end glam_crs gphi_crs ",p_glam_crs(1,2),p_gphi_crs(1,2)

      ! Retroactively add back the boundary halo cells.
      !????      !CALL crs_lbc_lnk( p_gphi_crs, cd_type, 1.0 )
      !????      !CALL crs_lbc_lnk( p_glam_crs, cd_type, 1.0 )
      WRITE(narea+1000-1,*)"end1 glam_crs gphi_crs ",p_glam_crs(1,2),p_gphi_crs(1,2)
      !
   END SUBROUTINE crs_dom_coordinates

   SUBROUTINE crs_dom_hgr( p_e1, p_e2, cd_type, p_e1_crs, p_e2_crs )
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crs_dom_hgr ***
      !!
      !! ** Purpose :  Get coarse grid horizontal scale factors and unmasked fraction
      !!
      !! ** Method  :  For grid types T,U,V,Fthe 2D scale factors of
      !!               the coarse grid are the sum of the east or north faces of the
      !!               parent grid subset comprising the coarse grid box.
      !!               - e1,e2 Scale factors
      !!                 Valid arguments:
      !! ** Inputs  : p_e1, p_e2  = parent grid e1 or e2 (t,u,v,f)
      !!              cd_type     = grid type (T,U,V,F) for scale factors; for velocities (U or V)
      !! ** Outputs : p_e1_crs, p_e2_crs  = parent grid e1 or e2 (t,u,v,f)
      !!
      !! History.     4 Jun.  Write for WGT and scale factors only
      !!----------------------------------------------------------------
      !!
      !!  Arguments
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(in)  :: p_e1     ! Parent grid U,V scale factors (e1)
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(in)  :: p_e2     ! Parent grid U,V scale factors (e2)
      CHARACTER(len=1)                    , INTENT(in)  :: cd_type  ! grid type U,V

      REAL(wp), DIMENSION(jpi_crs,jpj_crs), INTENT(out) :: p_e1_crs ! Coarse grid box 2D quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs), INTENT(out) :: p_e2_crs ! Coarse grid box 2D quantity

      !! Local variables
      INTEGER :: ji, jj, jk     ! dummy loop indices
      INTEGER :: ijis,ijie,ijjs,ijje
      INTEGER :: i1, j1

      !!----------------------------------------------------------------
      ! Initialize

      DO ji = nldi_crs, nlei_crs

         ijis = mis_crs(ji)
         ijie = mie_crs(ji)

         DO jj = nldj_crs, nlej_crs

            ijjs = mjs_crs(jj)
            ijje = mje_crs(jj)

            i1=INT(0.5*nfactx(ji))
            j1=INT(0.5*nfacty(jj))

            ! Only for a factro 3 coarsening
            SELECT CASE ( cd_type )
            CASE ( 'T' )
               p_e1_crs(ji,jj) = REAL(nn_factx,wp)*p_e1(ijis+i1,ijjs+j1)
               p_e2_crs(ji,jj) = REAL(nn_facty,wp)*p_e2(ijis+i1,ijjs+j1)
            CASE ( 'U' )
               p_e1_crs(ji,jj) = REAL(nn_factx,wp)*p_e1(ijis+i1,ijjs+j1)
               p_e2_crs(ji,jj) = REAL(nn_facty,wp)*p_e2(ijie   ,ijjs+j1)

            CASE ( 'V' )
               p_e1_crs(ji,jj) = REAL(nn_factx,wp)*p_e1(ijis+i1,ijje   )
               p_e2_crs(ji,jj) = REAL(nn_facty,wp)*p_e2(ijis+i1,ijjs+j1)
            CASE ( 'F' )
               p_e1_crs(ji,jj) = REAL(nn_factx,wp)*p_e1(ijis+i1,ijje   )
               p_e2_crs(ji,jj) = REAL(nn_facty,wp)*p_e2(ijie   ,ijjs+j1)
            END SELECT
         ENDDO
      ENDDO


      !CALL crs_lbc_lnk( p_e1_crs, cd_type, 1.0 )
      !CALL crs_lbc_lnk( p_e2_crs, cd_type, 1.0 )

   END SUBROUTINE crs_dom_hgr


   SUBROUTINE crs_dom_facvol( p_mask, cd_type, p_e1, p_e2, p_e3, p_fld1_crs, p_fld2_crs )
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crsfun_wgt ***
      !! ** Purpose :  Three applications.
      !!               1) SUM. Get coarse grid horizontal scale factors and unmasked fraction
      !!               2) VOL. Get coarse grid box volumes
      !!               3) WGT. Weighting multiplier for volume-weighted and/or
      !!                       area-weighted averages.
      !!                       Weights (i.e. the denominator) calculated here
      !!                       to avoid IF-tests and division.
      !! ** Method  :  1) SUM.  For grid types T,U,V,F (and W) the 2D scale factors of
      !!               the coarse grid are the sum of the east or north faces of the
      !!               parent grid subset comprising the coarse grid box.
      !!               The fractions of masked:total surface (3D) on the east,
      !!               north and top faces is, optionally, also output.
      !!               - Top face area sum
      !!                 Valid arguments: cd_type, cd_op='W', p_pmask, p_e1, p_e2
      !!               - Top face ocean surface fraction
      !!                 Valid arguments: cd_type, cd_op='W', p_pmask, p_e1, p_e2
      !!               - e1,e2 Scale factors
      !!                 Valid arguments:
      !!               2) VOL.  For grid types W and T, the coarse grid box
      !!               volumes are output. Also optionally, the fraction of
      !!               masked:total volume of the parent grid subset is output (i.e. facvol).
      !!               3) WGT. Based on the grid type, the denominator is pre-determined here to
      !!               perform area- or volume- weighted averages,
      !!               to avoid IF-tests and divisions.
      !! ** Inputs  : p_e1, p_e2  = parent grid e1 or e2 (t,u,v,f)
      !!              p_pmask     = parent grid mask (T,U,V,F)
      !!              cd_type     = grid type (T,U,V,F) for scale factors; for velocities (U or V)
      !!              cd_op       = applied operation (SUM, VOL, WGT)
      !!              p_fse3      = (Optional) parent grid vertical level thickness (fse3u or fse3v)
      !! ** Outputs : p_cfield2d_1 = (Optional) 2D field on coarse grid
      !!              p_cfield2d_2 = (Optional) 2D field on coarse grid
      !!              p_cfield3d_1 = (Optional) 3D field on coarse grid
      !!              p_cfield3d_2 = (Optional) 3D field on coarse grid
      !!
      !! History.     4 Jun.  Write for WGT and scale factors only
      !!----------------------------------------------------------------
      !!
      !!  Arguments
      CHARACTER(len=1),                         INTENT(in)  :: cd_type  ! grid type U,V
      REAL(wp), DIMENSION(jpi,jpj,jpk)        , INTENT(in)  :: p_mask  ! Parent grid U,V mask
      REAL(wp), DIMENSION(jpi,jpj)            , INTENT(in)  :: p_e1     ! Parent grid U,V scale factors (e1)
      REAL(wp), DIMENSION(jpi,jpj)            , INTENT(in)  :: p_e2     ! Parent grid U,V scale factors (e2)
      REAL(wp), DIMENSION(jpi,jpj,jpk)        , INTENT(in)  :: p_e3     ! Parent grid vertical level thickness (fse3u, fse3v)

      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(out) :: p_fld1_crs ! Coarse grid box 3D quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(out) :: p_fld2_crs ! Coarse grid box 3D quantity

      !! Local variables
      REAL(wp)                                :: zdAm
      INTEGER                                 :: ji, jj, jk
      INTEGER :: ijis,ijie,ijjs,ijje

      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE     :: zvol, zmask
      !!----------------------------------------------------------------
      
      ALLOCATE ( zvol(jpi,jpj,jpk), zmask(jpi,jpj,jpk) )

      p_fld1_crs(:,:,:) = 0.0
      p_fld2_crs(:,:,:) = 0.0

      DO jk = 1, jpk
         zvol (:,:,jk) = p_e1(:,:) * p_e2(:,:) * p_e3(:,:,jk)
         zmask(:,:,jk) = p_mask(:,:,jk)
      ENDDO

      DO jk = 1, jpk
         DO ji = nldi_crs, nlei_crs

            ijis = mis_crs(ji)
            ijie = mie_crs(ji)

            DO jj = nldj_crs, nlej_crs

               ijjs = mjs_crs(jj)
               ijje = mje_crs(jj)

               p_fld1_crs(ji,jj,jk) =  SUM( zvol(ijis:ijie,ijjs:ijje,jk) )
               zdAm                 =  SUM( zvol(ijis:ijie,ijjs:ijje,jk) * zmask(ijis:ijie,ijjs:ijje,jk) )
               p_fld2_crs(ji,jj,jk) = zdAm / p_fld1_crs(ji,jj,jk)
            ENDDO
         ENDDO
      ENDDO
      !                                             !  Retroactively add back the boundary halo cells.
      !CALL crs_lbc_lnk( p_fld1_crs, cd_type, 1.0 )
      !CALL crs_lbc_lnk( p_fld2_crs, cd_type, 1.0 )
      !
      DEALLOCATE ( zvol, zmask )
      !
   END SUBROUTINE crs_dom_facvol


   SUBROUTINE crs_dom_ope_3d( p_fld, cd_op, cd_type, p_mask, p_fld_crs, p_e12, p_e3, p_surf_crs, p_mask_crs, psgn )
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crsfun_UV ***
      !! ** Purpose :  Average, area-weighted, of U or V on the east and north faces
      !!
      !! ** Method  :  The U and V velocities (3D) are determined as the area-weighted averages
      !!               on the east and north faces, respectively,
      !!               of the parent grid subset comprising the coarse grid box.
      !!               In the case of the V and F grid, the last jrow minus 1 is spurious.
      !! ** Inputs  : p_e1_e2     = parent grid e1 or e2 (t,u,v,f)
      !!              cd_type     = grid type (T,U,V,F) for scale factors; for velocities (U or V)
      !!              psgn        = sign change over north fold (See lbclnk.F90)
      !!              p_pmask     = parent grid mask (T,U,V,F) for scale factors;
      !!                                       for velocities (U or V)
      !!              p_fse3      = parent grid vertical level thickness (fse3u or fse3v)
      !!              p_pfield    = U or V on the parent grid
      !!              p_surf_crs  = (Optional) Coarse grid weight for averaging
      !! ** Outputs : p_cfield3d = 3D field on coarse grid
      !!
      !! History.  29 May.  completed draft.
      !!            4 Jun.  Revision for WGT
      !!            5 Jun.  Streamline for area-weighted average only ; separate scale factor and weights.
      !!----------------------------------------------------------------
      !!
      !!  Arguments
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in)        :: p_fld   ! T, U, V or W on parent grid
      CHARACTER(len=*),                         INTENT(in)           :: cd_op    ! Operation SUM, MAX or MIN
      CHARACTER(len=1),                         INTENT(in)           :: cd_type    ! grid type U,V
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in)           :: p_mask    ! Parent grid T,U,V mask
      REAL(wp), DIMENSION(jpi,jpj),             INTENT(in), OPTIONAL :: p_e12    ! Parent grid T,U,V scale factors (e1 or e2)
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in), OPTIONAL :: p_e3     ! Parent grid vertical level thickness (fse3u, fse3v)
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(in), OPTIONAL :: p_surf_crs ! Coarse grid area-weighting denominator
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(in), OPTIONAL :: p_mask_crs    ! Coarse grid T,U,V maska
      REAL(wp),                                 INTENT(in)           :: psgn    ! sign

      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(out)          :: p_fld_crs ! Coarse grid box 3D quantity

      !! Local variables
      INTEGER  :: ji, jj, jk
      INTEGER  :: ijis, ijie, ijjs, ijje
      INTEGER  :: ini, inj
      REAL(wp) :: zflcrs, zsfcrs
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zsurf, zsurfmsk, zmask,ztabtmp
      INTEGER  :: ir,jr
      REAL(wp), DIMENSION(nn_factx,nn_facty):: ztmp
      REAL(wp), DIMENSION(nn_factx*nn_facty):: ztmp1
      REAL(wp), DIMENSION(:), ALLOCATABLE   :: ztmp2
      INTEGER , DIMENSION(1)  :: zdim1
      REAL(wp) :: zmin,zmax
      !!----------------------------------------------------------------

      p_fld_crs(:,:,:) = 0.0

      SELECT CASE ( cd_op )

      CASE ( 'VOL' )

         ALLOCATE ( zsurf(jpi,jpj,jpk), zsurfmsk(jpi,jpj,jpk) )

         SELECT CASE ( cd_type )

         CASE( 'T', 'W','U','V' )
            DO jk = 1, jpk
               zsurf   (:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk) *  p_mask(:,:,jk)
               zsurfmsk(:,:,jk) =  zsurf(:,:,jk)
            ENDDO
            !
            DO jk = 1, jpk
               DO jj  = nldj_crs,nlej_crs
                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  DO ji = nldi_crs, nlei_crs

                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)

                     zflcrs = SUM( p_fld(ijis:ijie,ijjs:ijje,jk) * zsurfmsk(ijis:ijie,ijjs:ijje,jk) )
                     zsfcrs = SUM(                                 zsurfmsk(ijis:ijie,ijjs:ijje,jk) )

                     p_fld_crs(ji,jj,jk) = zflcrs
                     IF( zsfcrs /= 0.0 )  p_fld_crs(ji,jj,jk) = zflcrs / zsfcrs
                  ENDDO
               ENDDO
            ENDDO
            !
         CASE DEFAULT
            STOP
         END SELECT

         DEALLOCATE ( zsurf, zsurfmsk )

      CASE ( 'LOGVOL' )

         ALLOCATE ( zsurf(jpi,jpj,jpk), zsurfmsk(jpi,jpj,jpk), ztabtmp(jpi,jpj,jpk) )

         ztabtmp(:,:,:)=0._wp
         WHERE(p_fld* p_mask .NE. 0._wp ) ztabtmp =  LOG10(p_fld * p_mask)*p_mask
         ztabtmp = ztabtmp * p_mask

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            DO jk = 1, jpk
               zsurf   (:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk) *  p_mask(:,:,jk)
               zsurfmsk(:,:,jk) =  zsurf(:,:,jk)
            ENDDO
            !
            DO jk = 1, jpk
               DO jj  = nldj_crs,nlej_crs
                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  DO ji = nldi_crs, nlei_crs
                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)
                     zflcrs = SUM( ztabtmp(ijis:ijie,ijjs:ijje,jk) * zsurfmsk(ijis:ijie,ijjs:ijje,jk) )
                     zsfcrs = SUM(                                   zsurfmsk(ijis:ijie,ijjs:ijje,jk) )
                     p_fld_crs(ji,jj,jk) = zflcrs
                     IF( zsfcrs /= 0.0 )  p_fld_crs(ji,jj,jk) = zflcrs / zsfcrs
                     p_fld_crs(ji,jj,jk) = 10 ** ( p_fld_crs(ji,jj,jk) *  p_mask_crs(ji,jj,jk) ) * p_mask_crs(ji,jj,jk)
                  ENDDO
               ENDDO
            ENDDO
         CASE DEFAULT
            STOP
         END SELECT

         DEALLOCATE ( zsurf, zsurfmsk ,ztabtmp )

      CASE ( 'MED' )

         ALLOCATE ( zsurf(jpi,jpj,jpk), zsurfmsk(jpi,jpj,jpk) )

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )
            DO jk = 1, jpk
               zsurf   (:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk) *  p_mask(:,:,jk)
               zsurfmsk(:,:,jk) =  zsurf(:,:,jk)
            ENDDO
            !
            DO jk = 1, jpk

               DO jj  = nldj_crs,nlej_crs

                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  inj  = ijje-ijjs+1

                  DO ji = nldi_crs, nlei_crs
                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)
                     ini  = ijie-ijis+1

                     ztmp(1:ini,1:inj)= p_fld(ijis:ijie,ijjs:ijje,jk)
                     zdim1(1) = nn_factx*nn_facty
                     ztmp1(:) = RESHAPE( ztmp(:,:) , zdim1 )
                     CALL PIKSRT(nn_factx*nn_facty,ztmp1)

                     ir=0
                     jr=1
                     DO WHILE( jr .LE. nn_factx*nn_facty )
                        IF( ztmp1(jr) == 0. ) THEN
                           ir=jr
                           jr=jr+1
                        ELSE
                           EXIT
                        ENDIF
                     ENDDO
                     IF( ir .LE. nn_factx*nn_facty-1 )THEN
                        ALLOCATE( ztmp2(nn_factx*nn_facty-ir) )
                        ztmp2(1:nn_factx*nn_facty-ir) = ztmp1(1+ir:nn_factx*nn_facty)
                        jr=INT( 0.5 * REAL(nn_factx*nn_facty-ir,wp) )+1
                        p_fld_crs(ji,jj,jk) = ztmp2(jr)
                        DEALLOCATE( ztmp2 )
                     ELSE
                        p_fld_crs(ji,jj,jk) = 0._wp
                     ENDIF

                  ENDDO
               ENDDO
            ENDDO
         CASE DEFAULT
            STOP
         END SELECT

         DEALLOCATE ( zsurf, zsurfmsk )

      CASE ( 'SUM' )

         ALLOCATE ( zsurfmsk(jpi,jpj,jpk) )

         IF( PRESENT( p_e3 ) ) THEN
            DO jk = 1, jpk
               zsurfmsk(:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk) * p_mask(:,:,jk)
            ENDDO
         ELSE
            DO jk = 1, jpk
               zsurfmsk(:,:,jk) =  p_e12(:,:) * p_mask(:,:,jk)
            ENDDO
         ENDIF

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            DO jk = 1, jpk
               DO jj  = nldj_crs,nlej_crs
                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  DO ji = nldi_crs, nlei_crs
                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)

                     p_fld_crs(ji,jj,jk) = SUM( p_fld(ijis:ijie,ijjs:ijje,jk) * zsurfmsk(ijis:ijie,ijjs:ijje,jk) )
                  ENDDO
               ENDDO
            ENDDO

         CASE( 'V' )


            DO jk = 1, jpk
               DO jj  = nldj_crs,nlej_crs
                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  DO ji = nldi_crs, nlei_crs
                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)

                     p_fld_crs(ji,jj,jk) = SUM( p_fld(ijis:ijie,ijje,jk) * zsurfmsk(ijis:ijie,ijje,jk) )
                  ENDDO
               ENDDO
            ENDDO

         CASE( 'U' )

            DO jk = 1, jpk
               DO jj  = nldj_crs,nlej_crs
                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  DO ji = nldi_crs, nlei_crs
                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)

                     p_fld_crs(ji,jj,jk) = SUM( p_fld(ijie,ijjs:ijje,jk) * zsurfmsk(ijie,ijjs:ijje,jk) )
                  ENDDO
               ENDDO
            ENDDO

         END SELECT

         IF( PRESENT( p_surf_crs ) ) THEN
            WHERE ( p_surf_crs /= 0.0 ) p_fld_crs(:,:,:) = p_fld_crs(:,:,:) / p_surf_crs(:,:,:)
         ENDIF

         DEALLOCATE ( zsurfmsk )

      CASE ( 'MAX' )    !  search the max of unmasked grid cells

         ALLOCATE ( zmask(jpi,jpj,jpk) )

         DO jk = 1, jpk
            zmask(:,:,jk) = p_mask(:,:,jk)
         ENDDO

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            DO jk = 1, jpk
               DO jj  = nldj_crs,nlej_crs
                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  DO ji = nldi_crs, nlei_crs
                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)
                     p_fld_crs(ji,jj,jk) = MAXVAL( p_fld(ijis:ijie,ijjs:ijje,jk) * zmask(ijis:ijie,ijjs:ijje,jk) - &
                        & ( ( 1._wp - zmask(ijis:ijie,ijjs:ijje,jk))* r_inf )                )
                  ENDDO
               ENDDO
            ENDDO

         CASE( 'V' )
            PRINT *, 'ctl_stop ERROR: MAX OPERATOR and V CASE not available'; STOP

         CASE( 'U' )
            PRINT *, 'ctl_stop ERROR: MAX operator and U case not available'; STOP

         END SELECT

         DEALLOCATE ( zmask )

      CASE ( 'MIN' )      !   Search the min of unmasked grid cells

         ALLOCATE ( zmask(jpi,jpj,jpk) )
         DO jk = 1, jpk
            zmask(:,:,jk) = p_mask(:,:,jk)
         ENDDO

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            DO jk = 1, jpk
               DO jj  = nldj_crs,nlej_crs
                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)
                  DO ji = nldi_crs, nlei_crs
                     ijis = mis_crs(ji)
                     ijie = mie_crs(ji)

                     p_fld_crs(ji,jj,jk) = MINVAL( p_fld(ijis:ijie,ijjs:ijje,jk) * zmask(ijis:ijie,ijjs:ijje,jk) + &
                        & ( 1._wp - zmask(ijis:ijie,ijjs:ijje,jk)* r_inf )                )
                  ENDDO
               ENDDO
            ENDDO


         CASE( 'V' )
            PRINT *, 'ctl_stop ERROR: MIN operator and V case not available'; STOP

         CASE( 'U' )
            PRINT *, 'ctl_stop ERROR: MIN operator and U case not available'; STOP

         END SELECT
         !
         DEALLOCATE ( zmask )
         !
      END SELECT
      !
      !CALL crs_lbc_lnk( p_fld_crs, cd_type, psgn )
      !
   END SUBROUTINE crs_dom_ope_3d

   SUBROUTINE crs_dom_ope_2d( p_fld, cd_op, cd_type, p_mask, p_fld_crs, p_e12, p_e3, p_surf_crs, p_mask_crs, psgn )
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crsfun_UV ***
      !! ** Purpose :  Average, area-weighted, of U or V on the east and north faces
      !!
      !! ** Method  :  The U and V velocities (3D) are determined as the area-weighted averages
      !!               on the east and north faces, respectively,
      !!               of the parent grid subset comprising the coarse grid box.
      !!               In the case of the V and F grid, the last jrow minus 1 is spurious.
      !! ** Inputs  : p_e1_e2     = parent grid e1 or e2 (t,u,v,f)
      !!              cd_type     = grid type (T,U,V,F) for scale factors; for velocities (U or V)
      !!              psgn        = sign change over north fold (See lbclnk.F90)
      !!              p_pmask     = parent grid mask (T,U,V,F) for scale factors;
      !!                                       for velocities (U or V)
      !!              p_fse3      = parent grid vertical level thickness (fse3u or fse3v)
      !!              p_pfield    = U or V on the parent grid
      !!              p_surf_crs  = (Optional) Coarse grid weight for averaging
      !! ** Outputs : p_cfield3d = 3D field on coarse grid
      !!
      !! History.  29 May.  completed draft.
      !!            4 Jun.  Revision for WGT
      !!            5 Jun.  Streamline for area-weighted average only ; separate scale factor and weights.
      !!----------------------------------------------------------------
      !!
      !!  Arguments
      REAL(wp), DIMENSION(jpi,jpj),             INTENT(in)           :: p_fld   ! T, U, V or W on parent grid
      CHARACTER(len=*),                         INTENT(in)           :: cd_op    ! Operation SUM, MAX or MIN
      CHARACTER(len=1),                         INTENT(in)           :: cd_type    ! grid type U,V
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in)           :: p_mask    ! Parent grid T,U,V mask
      REAL(wp), DIMENSION(jpi,jpj),             INTENT(in), OPTIONAL :: p_e12    ! Parent grid T,U,V scale factors (e1 or e2)
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in), OPTIONAL :: p_e3     ! Parent grid vertical level thickness (fse3u, fse3v)
      REAL(wp), DIMENSION(jpi_crs,jpj_crs)    , INTENT(in), OPTIONAL :: p_surf_crs ! Coarse grid area-weighting denominator
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(in), OPTIONAL :: p_mask_crs    ! Coarse grid T,U,V mask
      REAL(wp),                                 INTENT(in)           :: psgn

      REAL(wp), DIMENSION(jpi_crs,jpj_crs)    , INTENT(out)          :: p_fld_crs ! Coarse grid box 3D quantity

      !! Local variables
      INTEGER  :: ji, jj, jk                 ! dummy loop indices
      INTEGER ::  ijis, ijie, ijjs, ijje
      REAL(wp) :: zflcrs, zsfcrs
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: zsurfmsk

      !!----------------------------------------------------------------

      p_fld_crs(:,:) = 0.0

      SELECT CASE ( TRIM(cd_op) )

      CASE ( 'VOL' )

         ALLOCATE ( zsurfmsk(jpi,jpj) )
         zsurfmsk(:,:) =  p_e12(:,:) * p_e3(:,:,1) * p_mask(:,:,1)

         DO jj  = nldj_crs,nlej_crs
            ijjs = mjs_crs(jj)
            ijje = mje_crs(jj)
            DO ji = nldi_crs, nlei_crs
               ijis = mis_crs(ji)
               ijie = mie_crs(ji)

               zflcrs = SUM( p_fld(ijis:ijie,ijjs:ijje) * zsurfmsk(ijis:ijie,ijjs:ijje) )
               zsfcrs = SUM(                              zsurfmsk(ijis:ijie,ijjs:ijje) )

               p_fld_crs(ji,jj) = zflcrs
               IF( zsfcrs /= 0.0 )  p_fld_crs(ji,jj) = zflcrs / zsfcrs
            ENDDO
         ENDDO
         DEALLOCATE ( zsurfmsk )
         !

      CASE ( 'SUM' )

         ALLOCATE ( zsurfmsk(jpi,jpj) )
         IF( PRESENT( p_e3 ) ) THEN
            zsurfmsk(:,:) =  p_e12(:,:) * p_e3(:,:,1) * p_mask(:,:,1)
         ELSE
            zsurfmsk(:,:) =  p_e12(:,:) * p_mask(:,:,1)
         ENDIF

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            DO jj  = nldj_crs,nlej_crs
               ijjs = mjs_crs(jj)
               ijje = mje_crs(jj)
               DO ji = nldi_crs, nlei_crs
                  ijis = mis_crs(ji)
                  ijie = mie_crs(ji)
                  p_fld_crs(ji,jj) = SUM( p_fld(ijis:ijie,ijjs:ijje) * zsurfmsk(ijis:ijie,ijjs:ijje) )
               ENDDO
            ENDDO

         CASE( 'V' )

            DO jj  = nldj_crs,nlej_crs
               ijjs = mjs_crs(jj)
               ijje = mje_crs(jj)
               DO ji = nldi_crs, nlei_crs
                  ijis = mis_crs(ji)
                  ijie = mie_crs(ji)
                  p_fld_crs(ji,jj) = SUM( p_fld(ijis:ijie,ijje) * zsurfmsk(ijis:ijie,ijje) )
               ENDDO
            ENDDO

         CASE( 'U' )

            DO jj  = nldj_crs,nlej_crs
               ijjs = mjs_crs(jj)
               ijje = mje_crs(jj)
               DO ji = nldi_crs, nlei_crs
                  ijis = mis_crs(ji)
                  ijie = mie_crs(ji)
                  p_fld_crs(ji,jj) = SUM( p_fld(ijie,ijjs:ijje) * zsurfmsk(ijie,ijjs:ijje) )
               ENDDO
            ENDDO

         END SELECT

         IF( PRESENT( p_surf_crs ) ) THEN
            WHERE ( p_surf_crs /= 0.0 ) p_fld_crs(:,:) = p_fld_crs(:,:) / p_surf_crs(:,:)
         ENDIF

         DEALLOCATE ( zsurfmsk )

      CASE ( 'MAX' )

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            DO jj  = nldj_crs,nlej_crs
               ijjs = mjs_crs(jj)
               ijje = mje_crs(jj)
               DO ji = nldi_crs, nlei_crs
                  ijis = mis_crs(ji)
                  ijie = mie_crs(ji)
                  p_fld_crs(ji,jj) = MAXVAL( p_fld(ijis:ijie,ijjs:ijje) * p_mask(ijis:ijie,ijjs:ijje,1) - &
                     & ( 1._wp - p_mask(ijis:ijie,ijjs:ijje,1) )                    )
               ENDDO
            ENDDO

         CASE( 'V' )
            PRINT *, 'ctl_stop ERROR: MAX operator and V case not available'; STOP

         CASE( 'U' )
            PRINT *, 'ctl_stop ERROR: MAX operator and U case not available'; STOP

         END SELECT

      CASE ( 'MIN' )      !   Search the min of unmasked grid cells

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            DO jj  = nldj_crs,nlej_crs
               ijjs = mjs_crs(jj)
               ijje = mje_crs(jj)
               DO ji = nldi_crs, nlei_crs
                  ijis = mis_crs(ji)
                  ijie = mie_crs(ji)
                  p_fld_crs(ji,jj) = MINVAL( p_fld(ijis:ijie,ijjs:ijje) * p_mask(ijis:ijie,ijjs:ijje,1) + &
                     & ( 1._wp - p_mask(ijis:ijie,ijjs:ijje,1) )                    )
               ENDDO
            ENDDO

         CASE( 'V' )
            PRINT *, 'ctl_stop ERROR: MIN operator and V case not available'; STOP

         CASE( 'U' )
            PRINT *, 'ctl_stop ERROR: MIN operator and U case not available'; STOP

         END SELECT
         !
      END SELECT
      !
      !CALL crs_lbc_lnk( p_fld_crs, cd_type, psgn )
      !
   END SUBROUTINE crs_dom_ope_2d

   SUBROUTINE crs_dom_e3( p_e1, p_e2, p_e3, p_sfc_2d_crs,  p_sfc_3d_crs, cd_type, p_mask, p_e3_crs, p_e3_max_crs)
      !!----------------------------------------------------------------
      !!
      !!
      !!
      !!
      !!----------------------------------------------------------------
      !!  Arguments
      CHARACTER(len=1),                         INTENT(in)          :: cd_type           ! grid type T, W ( U, V, F)
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in)          :: p_mask            ! Parent grid T mask
      REAL(wp), DIMENSION(jpi,jpj)    ,         INTENT(in)          :: p_e1, p_e2        ! 2D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in)          :: p_e3              ! 3D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi_crs,jpj_crs)    , INTENT(in),OPTIONAL :: p_sfc_2d_crs      ! Coarse grid box east or north face quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(in),OPTIONAL :: p_sfc_3d_crs      ! Coarse grid box east or north face quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout)       :: p_e3_crs          ! Coarse grid box east or north face quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout)       :: p_e3_max_crs      ! Coarse grid box east or north face quantity

      !! Local variables
      INTEGER ::  ji, jj, jk                   ! dummy loop indices
      INTEGER ::  ijis, ijie, ijjs, ijje
      REAL(wp) :: ze3crs

      !!----------------------------------------------------------------
      p_e3_crs    (:,:,:) = 0._wp
      p_e3_max_crs(:,:,:) = 0._wp


      SELECT CASE ( cd_type )

      CASE ('T')

         DO jk = 1, jpk
            DO ji = nldi_crs, nlei_crs

               ijis = mis_crs(ji)
               ijie = mie_crs(ji)

               DO jj = nldj_crs, nlej_crs

                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)

                  p_e3_max_crs(ji,jj,jk) = MAXVAL( p_e3(ijis:ijie,ijjs:ijje,jk) * p_mask(ijis:ijie,ijjs:ijje,jk) )

                  ze3crs = SUM( p_e1(ijis:ijie,ijjs:ijje) * p_e2(ijis:ijie,ijjs:ijje) * p_e3(ijis:ijie,ijjs:ijje,jk) * p_mask(ijis:ijie,ijjs:ijje,jk) )
                  IF( p_sfc_3d_crs(ji,jj,jk) .NE. 0._wp )p_e3_crs(ji,jj,jk) = ze3crs / p_sfc_3d_crs(ji,jj,jk)

               ENDDO
            ENDDO
         ENDDO

      CASE ('U')

         DO jk = 1, jpk
            DO ji = nldi_crs, nlei_crs

               ijis = mis_crs(ji)
               ijie = mie_crs(ji)

               DO jj = nldj_crs, nlej_crs

                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)

                  p_e3_max_crs(ji,jj,jk) = MAXVAL( p_e3(ijie,ijjs:ijje,jk) * p_mask(ijie,ijjs:ijje,jk) )

                  ze3crs = SUM( p_e2(ijie,ijjs:ijje) * p_e3(ijie,ijjs:ijje,jk) * p_mask(ijie,ijjs:ijje,jk) )
                  IF( p_sfc_2d_crs(ji,jj) .NE. 0._wp )p_e3_crs(ji,jj,jk) = ze3crs / p_sfc_2d_crs(ji,jj)
               ENDDO
            ENDDO
         ENDDO

      CASE ('V')

         DO jk = 1, jpk
            DO ji = nldi_crs, nlei_crs

               ijis = mis_crs(ji)
               ijie = mie_crs(ji)

               DO jj = nldj_crs, nlej_crs

                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)

                  p_e3_max_crs(ji,jj,jk) = MAXVAL( p_e3(ijis:ijie,ijje,jk) * p_mask(ijis:ijie,ijje,jk) )

                  ze3crs = SUM( p_e1(ijis:ijie,ijje) * p_e3(ijis:ijie,ijje,jk) * p_mask(ijis:ijie,ijje,jk) )
                  IF( p_sfc_2d_crs(ji,jj) .NE. 0._wp )p_e3_crs(ji,jj,jk) = ze3crs / p_sfc_2d_crs(ji,jj)

               ENDDO
            ENDDO
         ENDDO

      CASE ('W')

         DO jk = 1, jpk
            DO ji = nldi_crs, nlei_crs

               ijis = mis_crs(ji)
               ijie = mie_crs(ji)

               DO jj = nldj_crs, nlej_crs

                  ijjs = mjs_crs(jj)
                  ijje = mje_crs(jj)

                  p_e3_max_crs(ji,jj,jk) = MAXVAL( p_e3(ijis:ijie,ijjs:ijje,jk) * p_mask(ijis:ijie,ijjs:ijje,jk) )

                  ze3crs = SUM( p_e1(ijis:ijie,ijjs:ijje) * p_e2(ijis:ijie,ijjs:ijje) * p_e3(ijis:ijie,ijjs:ijje,jk) * p_mask(ijis:ijie,ijjs:ijje,jk) )
                  IF( p_sfc_3d_crs(ji,jj,jk) .NE. 0._wp )p_e3_crs(ji,jj,jk) = ze3crs / p_sfc_3d_crs(ji,jj,jk)

               ENDDO
            ENDDO
         ENDDO

      END SELECT

      !CALL crs_lbc_lnk( p_e3_max_crs, cd_type, 1.0 )
      !CALL crs_lbc_lnk( p_e3_crs    , cd_type, 1.0 )

   END SUBROUTINE crs_dom_e3

   SUBROUTINE crs_dom_sfc(p_mask, cd_type, p_surf_crs, p_surf_crs_msk,  p_e1, p_e2, p_e3 )
      !!=========================================================================================
      !!
      !!
      !!=========================================================================================
      !!  Arguments
      CHARACTER(len=1),                         INTENT(in)           :: cd_type      ! grid type T, W ( U, V, F)
      REAL(wp), DIMENSION(jpi,jpj,jpk)        , INTENT(in)           :: p_mask       ! Parent grid T mask
      REAL(wp), DIMENSION(jpi,jpj)            , INTENT(in), OPTIONAL :: p_e1, p_e2         ! 3D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi,jpj,jpk)        , INTENT(in), OPTIONAL :: p_e3         ! 3D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(out)          :: p_surf_crs ! Coarse grid box east or north face quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(out)          :: p_surf_crs_msk ! Coarse grid box east or north face quantity

      !! Local variables
      INTEGER  :: ji, jj, jk                   ! dummy loop indices
      INTEGER  :: ijis,ijie,ijjs,ijje
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zsurf, zsurfmsk
      !!----------------------------------------------------------------
      ! Initialize
      p_surf_crs(:,:,:)=0._wp
      p_surf_crs_msk(:,:,:)=0._wp

      ALLOCATE ( zsurf(jpi,jpj,jpk), zsurfmsk(jpi,jpj,jpk) )
      !
      SELECT CASE ( cd_type )

      CASE ('W')
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e1(:,:) * p_e2(:,:)
         ENDDO

      CASE ('V')
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e1(:,:) * p_e3(:,:,jk)
         ENDDO

      CASE ('U')
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e2(:,:) * p_e3(:,:,jk)
         ENDDO

      CASE DEFAULT
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e1(:,:) * p_e2(:,:)
         ENDDO
      END SELECT

      DO jk = 1, jpk
         zsurfmsk(:,:,jk) = zsurf(:,:,jk) * p_mask(:,:,jk)
      ENDDO

      SELECT CASE ( cd_type )

      CASE ('W')

         DO jk = 1, jpk
            DO jj = nldj_crs,nlej_crs
               ijjs=mjs_crs(jj)
               ijje=mje_crs(jj)
               DO ji = nldi_crs,nlei_crs
                  ijis=mis_crs(ji)
                  ijie=mie_crs(ji)
                  p_surf_crs    (ji,jj,jk) =  SUM(zsurf   (ijis:ijie,ijjs:ijje,jk) )
                  p_surf_crs_msk(ji,jj,jk) =  SUM(zsurfmsk(ijis:ijie,ijjs:ijje,jk) )
               ENDDO
            ENDDO
         ENDDO

      CASE ('U')

         DO jk = 1, jpk
            DO jj = nldj_crs,nlej_crs
               ijjs=mjs_crs(jj)
               ijje=mje_crs(jj)
               DO ji = nldi_crs,nlei_crs
                  ijis=mis_crs(ji)
                  ijie=mie_crs(ji)
                  p_surf_crs    (ji,jj,jk) =  SUM(zsurf   (ijie,ijjs:ijje,jk) )
                  p_surf_crs_msk(ji,jj,jk) =  SUM(zsurfmsk(ijie,ijjs:ijje,jk) )
               ENDDO
            ENDDO
         ENDDO

      CASE ('V')

         DO jk = 1, jpk
            DO jj = nldj_crs,nlej_crs
               ijjs=mjs_crs(jj)
               ijje=mje_crs(jj)
               DO ji = nldi_crs,nlei_crs
                  ijis=mis_crs(ji)
                  ijie=mie_crs(ji)
                  p_surf_crs    (ji,jj,jk) =  SUM(zsurf   (ijis:ijie,ijje,jk) )
                  p_surf_crs_msk(ji,jj,jk) =  SUM(zsurfmsk(ijis:ijie,ijje,jk) )
               ENDDO
            ENDDO
         ENDDO

      END SELECT

      !CALL crs_lbc_lnk( p_surf_crs    , cd_type , 1. )
      !CALL crs_lbc_lnk( p_surf_crs_msk, cd_type , 1. )

      DEALLOCATE ( zsurf, zsurfmsk )

   END SUBROUTINE crs_dom_sfc

   SUBROUTINE crs_dom_def
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crs_dom_def ***
      !! ** Purpose :  Three applications.
      !!               1) Define global domain indice of the croasening grid
      !!               2) Define local domain indice of the croasening grid
      !!               3) Define the processor domain indice for a croasening grid
      !!----------------------------------------------------------------
      INTEGER :: ji,jj,jk,ijjgloT,ijis,ijie,ijjs,ijje,jn      ! dummy indices
      INTEGER :: ierr                                ! allocation error status
      INTEGER :: iproci,iprocj,iproc,iprocno,iprocso,iimppt_crs
      INTEGER :: ii_start,ii_end,ij_start,ij_end
      !!----------------------------------------------------------------


      !==============================================================================================
      ! Define global and local domain sizes
      !==============================================================================================
      SELECT CASE ( jperio )

      CASE ( 0, 1 )
         jpiglo_crs   = INT( (jpiglo - 2) / nn_factx ) + 2
         jpjglo_crs   = INT( (jpjglo - 2) / nn_facty ) + 2

      CASE ( 3, 4 )    !
         jpiglo_crs   = INT( (jpiglo - 2) / nn_factx ) + 2
         jpjglo_crs   = INT( (jpjglo - MOD(jpjglo, nn_facty)) / nn_facty ) + 2

      END SELECT

      jpiglo_crsm1 = jpiglo_crs - 1
      jpjglo_crsm1 = jpjglo_crs - 1

      jpi_crs = ( jpiglo_crs   - 2 * jpreci + (jpni-1) ) / jpni + 2 * jpreci
      jpj_crs = ( jpjglo_crsm1 - 2 * jprecj + (jpnj-1) ) / jpnj + 2 * jprecj

      jpi_crsm1   = jpi_crs - 1
      jpj_crsm1   = jpj_crs - 1

      npolj_crs   = npolj

      ierr = crs_dom_alloc()          ! allocate most coarse grid arrays

      !==============================================================================================
      ! Define processor domain indices
      !==============================================================================================

         nimpp_crs  = 1
         njmpp_crs  = 1
         nlci_crs   = jpi_crs
         nlcj_crs   = jpj_crs
         nldi_crs   = 1
         nldj_crs   = 1
         nlei_crs   = jpi_crs
         nlej_crs   = jpj_crs


      nistr = mis_crs(2)  ;   niend = mis_crs(nlci_crs - 1)
      njstr = mjs_crs(3)  ;   njend = mjs_crs(nlcj_crs - 1)
      !
      !
   END SUBROUTINE crs_dom_def

   SUBROUTINE crs_dom_bat
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crs_dom_bat ***
      !! ** Purpose :  coarsenig bathy
      !!----------------------------------------------------------------
      !!
      !!  local variables
      INTEGER  :: ji,jj,jk      ! dummy indices
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE :: zmbk
      !!----------------------------------------------------------------

      ALLOCATE ( zmbk(jpi_crs,jpj_crs) )

      mbathy_crs(:,:) = jpkm1
      mbkt_crs(:,:) = 1
      mbku_crs(:,:) = 1
      mbkv_crs(:,:) = 1


      DO jj = 1, jpj_crs
         DO ji = 1, jpi_crs
            jk = 0
            DO WHILE( tmask_crs(ji,jj,jk+1) > 0.)
               jk = jk + 1
            ENDDO
            mbathy_crs(ji,jj) = float( jk )
         ENDDO
      ENDDO

      !#LOLO_BUG! ALLOCATE ( zmbk(jpi_crs,jpj_crs) )

      zmbk(:,:) = 0.0
      zmbk(:,:) = REAL( mbathy_crs(:,:), wp ) ;   !CALL crs_lbc_lnk(zmbk,'T',1.0)   ;   mbathy_crs(:,:) = INT( zmbk(:,:) )


      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    crsini : mbkt is ocean bottom k-index of T-, U-, V- and W-levels '
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~~~'
      !
      mbkt_crs(:,:) = MAX( mbathy_crs(:,:) , 1 )    ! bottom k-index of T-level (=1 over land)
      !                                     ! bottom k-index of W-level = mbkt+1

      DO jj = 1, jpj_crsm1                      ! bottom k-index of u- (v-) level
         DO ji = 1, jpi_crsm1
            mbku_crs(ji,jj) = MIN(  mbkt_crs(ji+1,jj  ) , mbkt_crs(ji,jj)  )
            mbkv_crs(ji,jj) = MIN(  mbkt_crs(ji  ,jj+1) , mbkt_crs(ji,jj)  )
         END DO
      END DO

      ! convert into REAL to use lbc_lnk ; impose a min value of 1 as a zero can be set in lbclnk
      zmbk(:,:) = 1.e0;
      zmbk(:,:) = REAL( mbku_crs(:,:), wp )   ;   !CALL crs_lbc_lnk(zmbk,'U',1.0) ;
      mbku_crs  (:,:) = MAX( INT( zmbk(:,:) ), 1 )
      zmbk(:,:) = REAL( mbkv_crs(:,:), wp )   ;   !CALL crs_lbc_lnk(zmbk,'V',1.0)
      mbkv_crs  (:,:) = MAX( INT( zmbk(:,:) ), 1 )
      !
      DEALLOCATE ( zmbk )
      !
   END SUBROUTINE crs_dom_bat

   SUBROUTINE PIKSRT(N,ARR)
      INTEGER                  ,INTENT(IN) :: N
      REAL(kind=8),DIMENSION(N),INTENT(INOUT) :: ARR

      INTEGER      :: i,j
      REAL(kind=8) :: a
      !!----------------------------------------------------------------

      DO j=2, N
         a=ARR(j)
         DO i=j-1,1,-1
            IF(ARR(i)<=a) goto 10
            ARR(i+1)=ARR(i)
         ENDDO
         i=0
10       ARR(i+1)=a
      ENDDO
      RETURN

   END SUBROUTINE PIKSRT


END MODULE crsdom
