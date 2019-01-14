MODULE crsdom
   !!===================================================================
   !!                  ***  crs.F90 ***
   !!  Purpose: Interface for calculating quantities from a
   !!           higher-resolution grid for the coarse grid.
   !!
   !!  Method:
   !!
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
      PRINT *, ' ### crs_dom_msk CTRL: nldi_crs, nlei_crs =', nldi_crs, nlei_crs
      PRINT *, ' ### crs_dom_msk CTRL: nldj_crs, nlej_crs =', nldj_crs, nlej_crs
      PRINT *, ' ### crs_dom_msk CTRL: mje_crs, mjs_crs =', mje_crs, mjs_crs
      PRINT *, ' ### crs_dom_msk CTRL: mie_crs, mis_crs =', mie_crs, mis_crs

      DO jk = 1, jpkm1
         DO ji = nldi_crs, nlei_crs

            ijis = mis_crs(ji)
            ijie = mie_crs(ji)

            IF ( (ijis > 0).AND.(ijie > 0) ) THEN
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
            END IF !IF ( (ijis > 0).AND.(ijie > 0) )
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

      !PRINT *, ' LOLO (crsdom.f90), jpi_full =', jpi_full
      SELECT CASE ( cd_type )
      CASE ( 'T' )
         DO jj =  nldj_crs, nlej_crs
            ijj = mjs_crs(jj) + INT(0.5*nfacty(jj))
            DO ji = nldi_crs, nlei_crs
               iji = mis_crs(ji) + INT(0.5*nfactx(ji))
               IF ( (iji > 0).AND.(iji<=jpi_full) ) THEN !LOLO
                  p_gphi_crs(ji,jj) = p_gphi(iji,ijj)
                  p_glam_crs(ji,jj) = p_glam(iji,ijj)
               END IF
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

            IF ( (ijis+i1 > 0).AND.(ijis+i1 <= jpi_full) ) THEN !LOLO
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
            END IF
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


      !PRINT *, ' LOLO: this is crs_dom_ope_3d !!!, jpk = ', jpk
      !PRINT *, ' LOLO: pe3 column:', p_e3(jpi/2,jpj/2,:)

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

                     IF ( ijis > 0 ) THEN !LOLO
                        zflcrs = SUM( p_fld(ijis:ijie,ijjs:ijje,jk) * zsurfmsk(ijis:ijie,ijjs:ijje,jk) )
                        zsfcrs = SUM(                                 zsurfmsk(ijis:ijie,ijjs:ijje,jk) )

                        p_fld_crs(ji,jj,jk) = zflcrs
                        IF( zsfcrs /= 0.0 )  p_fld_crs(ji,jj,jk) = zflcrs / zsfcrs
                     END IF
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

               IF ( ijis > 0 ) THEN !LOLO
                  zflcrs = SUM( p_fld(ijis:ijie,ijjs:ijje) * zsurfmsk(ijis:ijie,ijjs:ijje) )
                  zsfcrs = SUM(                              zsurfmsk(ijis:ijie,ijjs:ijje) )

                  p_fld_crs(ji,jj) = zflcrs
                  IF( zsfcrs /= 0.0 )  p_fld_crs(ji,jj) = zflcrs / zsfcrs
               END IF
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

                  IF ( ijis > 0 ) THEN !LOLO
                     p_e3_max_crs(ji,jj,jk) = MAXVAL( p_e3(ijis:ijie,ijjs:ijje,jk) * p_mask(ijis:ijie,ijjs:ijje,jk) )

                     ze3crs = SUM( p_e1(ijis:ijie,ijjs:ijje) * p_e2(ijis:ijie,ijjs:ijje) * p_e3(ijis:ijie,ijjs:ijje,jk) * p_mask(ijis:ijie,ijjs:ijje,jk) )
                     IF( p_sfc_3d_crs(ji,jj,jk) .NE. 0._wp )p_e3_crs(ji,jj,jk) = ze3crs / p_sfc_3d_crs(ji,jj,jk)
                  END IF
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
                  IF ( ijis > 0 ) THEN !LOLO
                     p_surf_crs    (ji,jj,jk) =  SUM(zsurf   (ijis:ijie,ijjs:ijje,jk) )
                     p_surf_crs_msk(ji,jj,jk) =  SUM(zsurfmsk(ijis:ijie,ijjs:ijje,jk) )
                  END IF
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

      !PRINT *, ''
      !PRINT *, ' ### mod_crs.f90 ### jpiglo, jpjglo = ', jpiglo, jpjglo
      !PRINT *, ' ### mod_crs.f90 ### jpi, jpj       = ', jpi, jpj


      !PRINT *, ' ### mod_crs.f90 ### jperio, nn_factx, nn_facty', jperio, nn_factx, nn_facty
      !PRINT *, ' ### mod_crs.f90 ### jpreci, jprecj', jpreci, jprecj
      !PRINT *, ''

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


      !PRINT *, ' ### mod_crs.f90 ### jpiglo_crs, jpjglo_crs = ', jpiglo_crs, jpjglo_crs
      !PRINT *, ' ### mod_crs.f90 ### jpi_crs, jpj_crs       = ', jpi_crs, jpj_crs
      !PRINT *, ''

      ierr = crs_dom_alloc()          ! allocate most coarse grid arrays

      !==============================================================================================
      ! Define processor domain indices
      !==============================================================================================
      IF( .NOT. lk_mpp ) THEN

         nimpp_crs  = 1
         njmpp_crs  = 1
         nlci_crs   = jpi_crs
         nlcj_crs   = jpj_crs
         nldi_crs   = 1
         nldj_crs   = 1
         nlei_crs   = jpi_crs
         nlej_crs   = jpj_crs

      ELSE

         nimpp_crs  = 1
         njmpp_crs  = 1
         nlci_crs   = jpi_crs
         nlcj_crs   = jpj_crs
         nldi_crs   = 1
         nldj_crs   = 1
         nlei_crs   = jpi_crs
         nlej_crs   = jpj_crs

         !==============================================================================================
         ! mpp_ini2
         !==============================================================================================

         !order of local domain in i and j directions
         DO ji = 1 , jpni
            DO jj = 1 ,jpnj
               IF( nfipproc(ji,jj)  == narea-1 )THEN
                  iproci=ji
                  iprocj=jj
               ENDIF
            ENDDO
         ENDDO

         !==========================================================================
         ! check
         !==========================================================================
         !CALL FLUSH(narea+1000-1)
         !WRITE(narea+1000-1,*)"nfipproc(ji,jj),narea :",nfipproc(iproci,iprocj),narea
         !WRITE(narea+1000-1,*)"proc i,j ",iproci,iprocj
         !WRITE(narea+1000-1,*)"jpni  jpnj jpnij ",jpni,jpnj,jpnij
         !WRITE(narea+1000-1,*)"nperio jperio ",nperio,jperio
         !WRITE(narea+1000-1,*)"nowe noea",nowe,noea
         !WRITE(narea+1000-1,*)"noso nono",noso,nono
         !WRITE(narea+1000-1,*)"nbondi nbondj ",nbondi,nbondj
         !WRITE(narea+1000-1,*)"jpiglo jpjglo ",jpiglo,jpjglo
         !WRITE(narea+1000-1,*)"jpi jpj ",jpi,jpj
         !WRITE(narea+1000-1,*)"nbondi nbondj",nbondi,nbondj
         !WRITE(narea+1000-1,*)"nimpp njmpp ",nimpp,njmpp
         !WRITE(narea+1000-1,*)"loc jpi nldi,nlei,nlci ",jpi, nldi        ,nlei         ,nlci
         !WRITE(narea+1000-1,*)"glo jpi nldi,nlei      ",jpi, nldi+nimpp-1,nlei+nimpp-1
         !WRITE(narea+1000-1,*)"loc jpj nldj,nlej,nlcj ",jpj, nldj        ,nlej         ,nlcj
         !WRITE(narea+1000-1,*)"glo jpj nldj,nlej      ",jpj, nldj+njmpp-1,nlej+njmpp-1
         !WRITE(narea+1000-1,*)"jpiglo_crs jpjglo_crs ",jpiglo_crs,jpjglo_crs
         !WRITE(narea+1000-1,*)"jpi_crs jpj_crs ",jpi_crs,jpj_crs
         !WRITE(narea+1000-1,*)"glamt gphit ",glamt(1,1),gphit(jpi,jpj),glamt(jpi,jpj),gphit(jpi,jpj)
         !WRITE(narea+1000-1,*)"min max tmask ",MINVAL(tmask),MAXVAL(tmask)
         !CALL FLUSH(narea+1000-1)
         !==========================================================================
         ! coarsened domain: dimensions along I
         !==========================================================================

         !------------------------------------------------------------------------------------
         !I-1 fill mis2_crs and mie2_crs: arrays to switch from physic grid to coarsened grid
         !------------------------------------------------------------------------------------

         ! !--------!--------!--------!
         ! !        !        !        !
         ! !        !        !        !
         ! !        !        !        !
         ! !--------!--------!--------!
         ! !        !        !        !
         ! !        ! ji,jj  !        !
         ! !        !        !        !
         ! !--------!--------!--------!
         ! !        !        !        !
         ! !        !        !        !
         ! !        !        !        !
         ! !--------!--------!--------!
         !  mis2_crs(ji)      mie2_crs(ji)


         SELECT CASE ( jperio )

         CASE ( 0, 1 )

            DO ji=1,jpiglo_crs
               ijis=nn_factx*(ji-1)+1
               ijie=nn_factx*(ji-1)+3
               mis2_crs(ji)=ijis
               mie2_crs(ji)=ijie
            ENDDO

         CASE ( 3, 4 )    !   3, 4 : T-Pivot at North Fold: make correspondinf the pivot points of the 2 grids

            DO ji=1,jpiglo_crs
               ijis=nn_factx*(ji-1)-2
               ijie=nn_factx*(ji-1)
               mis2_crs(ji)=ijis
               mie2_crs(ji)=ijie
               !WRITE(narea+1000-1,*)"glo crs",ji,ijis,ijie,ijis-nimpp+1,ijie-nimpp+1
            ENDDO

         CASE DEFAULT
            WRITE(numout,*) 'crs_init. Only jperio = 0, 1, 3, 4 supported; narea: ',narea
         END SELECT

         !-------------------------------------------------------------------------------
         ! I-2 find the first CRS cell which is inside the physic grid inner domain
         !-------------------------------------------------------------------------------
         ! ijis           : global indice of the first CRS cell which inside the physic grid inner domain
         ! mis2_crs(ijis) : global indice of the bottom-left physic cell corresponding to ijis cell
         ! ii_start       : local  ndice of the bottom-left physic cell corresponding to ijis cell

         ji=1
         DO WHILE( mis2_crs(ji) - nimpp + 1 .LT. 1 )
            ji=ji+1
            IF( ji==jpiglo_crs )EXIT
         END DO

         ijis=ji
         ii_start = mis2_crs(ijis)-nimpp+1
         !WRITE(narea+1000-1,*)"start ",ijis,mis2_crs(ijis),ii_start ; CALL FLUSH(narea+1000-1)

         !----------------------------------------------------------------------------------------------
         ! I-3 compute nldi_crs and compute mis2_crs and mie2_crs for the first cell of the local domain
         !---------------------------------------------------------------------------------------------
         nldi_crs = 2
         IF( nowe == -1 .AND. ( (jperio==3 .OR. jperio==4 ) .OR. ( (jperio==0 .OR. jperio==1 ) .AND. iproci .NE. 1 )) )THEN

            mie2_crs(ijis-1) = mis2_crs(ijis)-1

            SELECT CASE(ii_start)
            CASE(1)
               nldi_crs=2
               mie2_crs(ijis-1) = -1
               mis2_crs(ijis-1) = -1
            CASE(2)
               nldi_crs=2
               mis2_crs(ijis-1) = mie2_crs(ijis-1)
            CASE(3)
               nldi_crs=2
               mis2_crs(ijis-1) = mie2_crs(ijis-1) -1
            CASE DEFAULT
               WRITE(narea+8000-1,*)"WRONG VALUE FOR iistart ",ii_start
            END SELECT

         ENDIF

         !----------------------------------------------------------------------------------------------
         ! I-4 compute nimpp_crs
         !---------------------------------------------------------------------------------------------
         nimpp_crs = ijis-1
         IF( nimpp==1 )nimpp_crs=1

         !-------------------------------------------------------------------------------
         ! I-5 find the last CRS cell which is inside the physic grid inner domain
         !-------------------------------------------------------------------------------
         ! ijie           : global indice of the last CRS cell which inside the physic grid inner domain

         ji=jpiglo_crs
         DO WHILE( mie2_crs(ji) - nimpp + 1 .GT. nlci )
            ji=ji-1
            IF( ji==1 )EXIT
         END DO
         ijie=ji
         !WRITE(narea+1000-1,*)"end ",ijie ; CALL FLUSH(narea+1000-1)

         !-------------------------------------------------------------------------------
         ! I-6 compute nlei_crs and nlci_crs
         !-------------------------------------------------------------------------------
         nlei_crs=ijie-nimpp_crs+1
         !nlci_crs=nlei_crs ! cbr ???? +jpreci
         !IF( iproci == jpni ) THEN ; nlci_crs=nlei_crs ! cbr ???? +jpreci
         !ELSE                      ; nlci_crs=nlei_crs+1
         !ENDIF
         !cbr???? IF( iproci == jpni )nlei_crs=nlci_crs

         !-------------------------------------------------------------------------------
         ! I-7 local to global and global to local indices for CRS grid
         !-------------------------------------------------------------------------------
         DO ji = 1, jpi_crs
            mig_crs(ji) = ji + nimpp_crs - 1
            !WRITE(narea+1000-1,*)"crs loctoglo",ji,mig_crs(ji) ; CALL FLUSH(narea+1000-1)
         ENDDO
         DO ji = 1, jpiglo_crs
            mi0_crs(ji) = MAX( 1, MIN( ji - nimpp_crs + 1 , jpi_crs + 1 ) )
            mi1_crs(ji) = MAX( 0, MIN( ji - nimpp_crs + 1 , jpi_crs     ) )
         ENDDO

         !---------------------------------------------------------------------------------------
         ! I-8 CRS to physic grid: local indices mis_crs and mie_crs and local coarsening factor
         !---------------------------------------------------------------------------------------
         DO ji = 1, nlei_crs
            IF( mig_crs(ji) .GT. jpiglo_crs )WRITE(narea+1000-1,*)"BUG1 " ; CALL FLUSH(narea+1000-1)
            mis_crs(ji) = mis2_crs(mig_crs(ji)) - nimpp + 1
            mie_crs(ji) = mie2_crs(mig_crs(ji)) - nimpp + 1
            !IF( iproci == jpni  .AND. ji == nlei_crs )THEN
            !   mie_crs(ji) = nlei
            !   mie2_crs(mig_crs(ji)) = mig(nlei)
            !ENDIF
            nfactx(ji)  = mie_crs(ji)-mis_crs(ji)+1
         ENDDO

         !---------
         !cbr
         IF( iproci == 1 ) THEN !lolo
            nldi_crs=1
            mis_crs(1) = 1
            mie_crs(1) = 1
            mis2_crs(1) = 1
            mie2_crs(1) = 1
         ENDIF

         IF( iproci == jpni ) THEN
            nlei_crs=jpiglo_crs-nimpp_crs+1
            nlci_crs=nlei_crs
            mis_crs(nlei_crs) = 1
            mie_crs(nlei_crs) = 1
            mis2_crs(nlei_crs) = 1
            mie2_crs(nlei_crs) = 1
            nfactx(nlei_crs)=0
         ELSE
            nlci_crs=nlei_crs+1
         ENDIF

         !WRITE(narea+1000-1,*)"loc crs jpi nldi,nlei,nlci ",jpi_crs, nldi_crs            ,nlei_crs             ,nlci_crs
         !CALL FLUSH(narea+1000-1)
         !WRITE(narea+1000-1,*)"glo crs jpi nldi,nlei      ",jpi_crs, nldi_crs+nimpp_crs-1,nlei_crs+nimpp_crs-1
         !CALL FLUSH(narea+1000-1)
         !DO ji = 1, jpi_crs
         !   WRITE(narea+1000-1,'(A15,7(1X,I4.4))')"crs test i",ji,ji+nimpp_crs-1,mis_crs(ji),mie_crs(ji), &
         !        &                                        mis_crs(ji)+nimpp-1,mie_crs(ji)+nimpp-1,nfactx(ji)
         !ENDDO

         !==========================================================================
         ! coarsened domain: dimensions along J
         !==========================================================================

         !-----------------------------------------------------------------------------------
         !J-1 fill mjs2_crs and mje2_crs: arrays to switch from physic grid to coarsened grid
         !-----------------------------------------------------------------------------------

         ! !--------!--------!--------!
         ! !        !        !        !
         ! !        !        !        !
         ! !        !        !        ! mje2_crs(jj)
         ! !--------!--------!--------!
         ! !        !        !        !
         ! !        ! ji,jj  !        !
         ! !        !        !        !
         ! !--------!--------!--------!
         ! !        !        !        !
         ! !        !        !        ! mjs2_crs(jj)
         ! !        !        !        !
         ! !--------!--------!--------!

         SELECT CASE ( jperio )

         CASE ( 0, 1 )    !

            DO jj=1,jpjglo_crs
               ijjs=nn_facty*(jj-1)+1
               ijje=nn_facty*(jj-1)+3
               mjs2_crs(jj)=ijjs
               mje2_crs(jj)=ijje
            ENDDO

            PRINT *, 'LOLO/crsdom.f90: mjs2_crs, mje2_crs ='
            PRINT *, mje2_crs ; PRINT *, ''
            PRINT *, mjs2_crs ; PRINT *, ''
            PRINT *, ''


         CASE ( 3, 4 )    !   3, 4 : T-Pivot at North Fold

            DO jj=1,jpjglo_crs
               ijjs=nn_facty*(jj)-5
               ijje=nn_facty*(jj)-3
               mjs2_crs(jj)=ijjs
               mje2_crs(jj)=ijje
            ENDDO

         CASE DEFAULT
            WRITE(numout,*) 'crs_init. Only jperio = 0, 1, 3, 4 supported; narea: ',narea
         END SELECT

         !-------------------------------------------------------------------------------
         ! J-2 find the first CRS cell which is inside the physic grid inner domain
         !-------------------------------------------------------------------------------
         ! ijjs           : global indice of the first CRS cell which inside the physic grid inner domain
         ! mis2_crs(ijjs) : global indice of the bottom-left physic cell corresponding to ijis cell
         ! ij_start       : local  ndice of the bottom-left physic cell corresponding to ijis cell

         jj=1
         DO WHILE( mjs2_crs(jj) - njmpp + 1 .LT. 1 )
            jj=jj+1
            IF( jj==jpjglo_crs )EXIT
         END DO

         ijjs=jj
         ij_start = mjs2_crs(ijjs)-njmpp+1

         !----------------------------------------------------------------------------------------------
         ! J-3 compute nldj_crs and compute mjs2_crs and mje2_crs for the first cell of the local domain
         !---------------------------------------------------------------------------------------------
         nldj_crs = 2

         IF( jperio==3 .OR. jperio==4 )THEN

            IF( noso == -1 )THEN

               mje2_crs(ijjs-1) = mjs2_crs(ijjs)-1

               SELECT CASE(ij_start)
               CASE(1)
                  nldj_crs=2
                  mje2_crs(ijjs-1) = -1
                  mjs2_crs(ijjs-1) = -1
               CASE(2)
                  nldj_crs=2
                  mjs2_crs(ijjs-1) = mje2_crs(ijjs-1)
               CASE(3)
                  nldj_crs=2
                  mjs2_crs(ijjs-1) = mje2_crs(ijjs-1) -1
               CASE DEFAULT
                  WRITE(narea+8000-1,*)"WRONG VALUE FOR iistart ",ii_start
               END SELECT

            ENDIF
         ENDIF
         IF( jperio==0 .OR. jperio==1 )THEN
            IF( njmpp==1 )nldj_crs=1
         END IF

         !----------------------------------------------------------------------------------------------
         ! J-4 compute nimpp_crs
         !---------------------------------------------------------------------------------------------
         njmpp_crs = ijjs-1
         IF( njmpp==1 )njmpp_crs=1

         !-------------------------------------------------------------------------------
         ! J-5 find the last CRS cell which is inside the physic grid inner domain
         !-------------------------------------------------------------------------------
         ! ijje           : global indice of the last CRS cell which inside the physic grid inner domain

         jj=jpjglo_crs
         DO WHILE( mje2_crs(jj) - njmpp + 1 .GT. nlcj )
            jj=jj-1
            IF( jj==1 )EXIT
         END DO
         ijje=jj

         !-------------------------------------------------------------------------------
         ! J-6 compute nlej_crs and nlcj_crs
         !-------------------------------------------------------------------------------
         nlej_crs=ijje-njmpp_crs+1
         nlcj_crs=nlej_crs+jprecj

         IF( iprocj == jpnj )THEN
            IF( jperio==3 .OR. jperio==4 )THEN
               nlej_crs=jpj_crs
               nlcj_crs=nlej_crs
            ELSE
               nlej_crs= nlej_crs+1
               nlcj_crs= nlcj_crs+1
            ENDIF
         ENDIF


         ! LOLO: fix overide because we don't have MPP:
         nimpp_crs  = 1
         njmpp_crs  = 1
         nlci_crs   = jpi_crs
         nlcj_crs   = jpj_crs
         nldi_crs   = 1
         nldj_crs   = 1
         nlei_crs   = jpi_crs
         nlej_crs   = jpj_crs
         !LOLO.

         !PRINT *, ' ### mod_crs.f90 ### nlei_crs, nlej_crs', nlei_crs, nlej_crs
         !PRINT *, ''


         !WRITE(narea+1000-1,*)"loc crs jpj nldj,nlej,nlcj ",jpj_crs, nldj_crs            ,nlej_crs             ,nlcj_crs
         !CALL FLUSH(narea+1000-1)
         !WRITE(narea+1000-1,*)"glo crs jpj nldj,nlej      ",jpj_crs, nldj_crs+njmpp_crs-1,nlej_crs+njmpp_crs-1
         !CALL FLUSH(narea+1000-1)
         !-------------------------------------------------------------------------------
         ! J-7 local to global and global to local indices for CRS grid
         !-------------------------------------------------------------------------------
         DO jj = 1, jpj_crs
            mjg_crs(jj) = jj + njmpp_crs - 1
         ENDDO
         DO jj = 1, jpjglo_crs
            mj0_crs(jj) = MAX( 1, MIN( jj - njmpp_crs + 1 , jpj_crs + 1 ) )
            mj1_crs(jj) = MAX( 0, MIN( jj - njmpp_crs + 1 , jpj_crs     ) )
         ENDDO

         !---------------------------------------------------------------------------------------
         ! J-8 CRS to physic grid: local indices mis_crs and mie_crs and local coarsening factor
         !---------------------------------------------------------------------------------------

         !PRINT *, 'mjg = ', mjg(:)
         !PRINT *, ' SIZE(mjg,1) = ', SIZE(mjg,1)
         !STOP'LOLO! mjg must be allocated!!!'

         DO jj = 1, nlej_crs
            mjs_crs(jj) = mjs2_crs(mjg_crs(jj)) - njmpp + 1
            mje_crs(jj) = mje2_crs(mjg_crs(jj)) - njmpp + 1
            IF( iprocj == jpnj  .AND. jj == nlej_crs )THEN
               mjs_crs(jj) = nlej
               mjs2_crs(mjg_crs(jj)) = mjg(nlej)
               mje_crs(jj) = nlej
               mje2_crs(mjg_crs(jj)) = mjg(nlej)
            ENDIF
            nfacty(jj)  = mje_crs(jj)-mjs_crs(jj)+1
            !WRITE(narea+1000-1,*)"test J",jj,mjg_crs(jj),mjs_crs(jj),mje_crs(jj),mjs_crs(jj)+njmpp-1,mje_crs(jj)+njmpp-1,nfacty(jj)
         ENDDO

         !WRITE(narea+1000-1,*)"loc crs jpj nldj,nlej,nlcj ",jpj_crs, nldj_crs ,nlej_crs,nlcj_crs ; CALL FLUSH(narea+1000-1)
         !WRITE(narea+1000-1,*)"glo crs jpj nldj,nlej      ",jpj_crs, nldj_crs+njmpp_crs-1,nlej_crs+njmpp_crs-1 ;CALL FLUSH(narea+1000-1)
         !DO jj = 1, jpj_crs
         !   WRITE(narea+1000-1,'(A15,7(1X,I4.4))')"crs test j",jj,jj+njmpp_crs-1,mjs_crs(jj),mje_crs(jj), &
         !        & mjs_crs(jj)+njmpp-1,mje_crs(jj)+njmpp-1,nfacty(jj)
         !ENDDO

         !==========================================================================
         ! send local start and end indices to all procs
         !==========================================================================

         nldit_crs(:)=0 ; nleit_crs(:)=0 ; nlcit_crs(:)=0 ; nimppt_crs(:)=0
         nldjt_crs(:)=0 ; nlejt_crs(:)=0 ; nlcjt_crs(:)=0 ; njmppt_crs(:)=0

         !CALL mppgatheri((/nlci_crs/),0,nlcit_crs) ; CALL mppgatheri((/nlcj_crs/),0,nlcjt_crs)
         !CALL mppgatheri((/nldi_crs/),0,nldit_crs) ; CALL mppgatheri((/nldj_crs/),0,nldjt_crs)
         !CALL mppgatheri((/nlei_crs/),0,nleit_crs) ; CALL mppgatheri((/nlej_crs/),0,nlejt_crs)
         !CALL mppgatheri((/nimpp_crs/),0,nimppt_crs) ; CALL mppgatheri((/njmpp_crs/),0,njmppt_crs)

         DO jj = 1 ,jpnj
            DO ji = 1 , jpni
               jn=nfipproc(ji,jj)+1
               IF( jn .GE. 1 )THEN
                  nfiimpp_crs(ji,jj)=nimppt_crs(jn)
               ELSE
                  nfiimpp_crs(ji,jj) = ANINT( REAL( (nfiimpp(ji,jj) + 1 ) / nn_factx, wp ) ) + 1
               ENDIF
            ENDDO
         ENDDO

         !nogather=T
         nfsloop_crs = 1
         nfeloop_crs = nlci_crs
         DO jn = 2,jpni-1
            IF( nfipproc(jn,jpnj) .eq. (narea - 1) )THEN
               IF (nfipproc(jn - 1 ,jpnj) .eq. -1 )THEN
                  nfsloop_crs = nldi_crs
               ENDIF
               IF( nfipproc(jn + 1,jpnj) .eq. -1 )THEN
                  nfeloop_crs = nlei_crs
               ENDIF
            ENDIF
         END DO

         !==========================================================================
         ! check
         !==========================================================================
         !WRITE(narea+1000-1,*)"mpp_crs ",nimpp_crs,njmpp_crs
         !WRITE(narea+1000-1,*)"loc crs jpi nldi,nlei,nlci ",jpi_crs, nldi_crs            ,nlei_crs             ,nlci_crs
         !WRITE(narea+1000-1,*)"glo crs jpi nldi,nlei      ",jpi_crs, nldi_crs+nimpp_crs-1,nlei_crs+nimpp_crs-1
         !WRITE(narea+1000-1,*)"loc crs jpj nldj,nlej,nlcj ",jpj_crs, nldj_crs            ,nlej_crs             ,nlcj_crs
         !WRITE(narea+1000-1,*)"glo crs jpj nldj,nlej      ",jpj_crs, nldj_crs+njmpp_crs-1,nlej_crs+njmpp_crs-1

         !==========================================================================
         ! Save the parent grid information
         !==========================================================================
         IF( jpizoom /= 1 .OR. jpjzoom /= 1)    STOP  !cbr mettre un ctlstp et ailleurs ( crsini )
         jpi_full    = jpi
         jpj_full    = jpj
         jpim1_full  = jpim1
         jpjm1_full  = jpjm1
         npolj_full  = npolj
         jpiglo_full = jpiglo
         jpjglo_full = jpjglo

         nlcj_full   = nlcj
         nlci_full   = nlci
         nldi_full   = nldi
         nldj_full   = nldj
         nlei_full   = nlei
         nlej_full   = nlej
         nimpp_full  = nimpp
         njmpp_full  = njmpp

         !PRINT *, 'nimppt = ', nimppt(:)
         !PRINT *, ' SIZE(nimppt,1) = ', SIZE(nimppt,1)
         !STOP'LOLO! mjg must be allocated!!!'

         !nlcit_full(:)  = nlcit(:)
         !nldit_full(:)  = nldit(:)
         !nleit_full(:)  = nleit(:)
         !nimppt_full(:) = nimppt(:)
         !nlcjt_full(:)  = nlcjt(:)
         !nldjt_full(:)  = nldjt(:)
         !nlejt_full(:)  = nlejt(:)
         !njmppt_full(:) = njmppt(:)

         !nfsloop_full = nfsloop
         !nfeloop_full = nfeloop

         !nfiimpp_full(:,:) = nfiimpp(:,:)


         !==========================================================================
         ! control
         !==========================================================================
         CALL dom_grid_crs  !swich from mother grid to coarsened grid

         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'crs_init : coarse grid dimensions'
            WRITE(numout,*) '~~~~~~~   coarse domain global j-dimension           jpjglo = ', jpjglo
            WRITE(numout,*) '~~~~~~~   coarse domain global i-dimension           jpiglo = ', jpiglo
            WRITE(numout,*) '~~~~~~~   coarse domain local  i-dimension              jpi = ', jpi
            WRITE(numout,*) '~~~~~~~   coarse domain local  j-dimension              jpj = ', jpj
            WRITE(numout,*)
            WRITE(numout,*) ' nproc  = '     , nproc
            WRITE(numout,*) ' nlci   = '     , nlci
            WRITE(numout,*) ' nlcj   = '     , nlcj
            WRITE(numout,*) ' nldi   = '     , nldi
            WRITE(numout,*) ' nldj   = '     , nldj
            WRITE(numout,*) ' nlei   = '     , nlei
            WRITE(numout,*) ' nlej   = '     , nlej
            WRITE(numout,*) ' nlei_full='    , nlei_full
            WRITE(numout,*) ' nldi_full='    , nldi_full
            WRITE(numout,*) ' nimpp  = '     , nimpp
            WRITE(numout,*) ' njmpp  = '     , njmpp
            WRITE(numout,*) ' njmpp_full  = ', njmpp_full
            WRITE(numout,*)
         ENDIF

         CALL dom_grid_glo ! switch from coarsened grid to mother grid

         nrestx = MOD( nn_factx, 2 )   ! check if even- or odd- numbered reduction factor
         nresty = MOD( nn_facty, 2 )

         IF( nresty == 0 )THEN
            IF( npolj == 3 ) npolj_crs = 5
            IF( npolj == 5 ) npolj_crs = 3
         ENDIF

         rfactxy = nn_factx * nn_facty

      ENDIF ! lk_mpp
      !
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
      zmbk(:,:) = REAL( mbathy_crs(:,:), wp ) ;   !CALL crs_lbc_lnk(zmbk,'T',1.0)
      mbathy_crs(:,:) = INT( zmbk(:,:) )


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
      REAL(wp),DIMENSION(N),INTENT(INOUT) :: ARR

      INTEGER      :: i,j
      REAL(wp) :: a
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
