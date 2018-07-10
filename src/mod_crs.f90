MODULE mod_crs
   !!===================================================================
   !!                  ***  crs.F90 ***
   !!  Purpose: Interface for calculating quantities from a
   !!           higher-resolution grid for the coarse grid.
   !!
   !!  Method:  Given the user-defined reduction factor,
   !!           the averaging bins are set:
   !!           - nn_binref = 0, starting from the north
   !!           to the south in the model interior domain,
   !!           in this way the north fold and redundant halo cells
   !!           could be handled in a consistent manner and
   !!           the irregularities of bin size can be handled
   !!           more naturally by the presence of land
   !!           in the southern boundary.  Thus the southernmost bin
   !!           could be of an irregular bin size.
   !!           Information on the parent grid is retained, specifically,
   !!           each coarse grid cell's volume and ocean surface
   !!           at the faces, relative to the parent grid.
   !!           - nn_binref = 1 (not yet available), starting
   !!           at a centralized bin at the equator, being only
   !!           truly centered for odd-numbered j-direction reduction
   !!           factors.
   !!  References:  Aumont, O., J.C. Orr, D. Jamous, P. Monfray
   !!               O. Marti and G. Madec, 1998. A degradation
   !!               approach to accelerate simulations to steady-state
   !!               in a 3-D tracer transport model of the global ocean.
   !!               Climate Dynamics, 14:101-116.
   !!  History:
   !!       Original.   May 2012.  (J. Simeon, C. Calone, G. Madec, C. Ethe)
   !!===================================================================

   !USE dom_oce        ! ocean space and time domain and to get jperio
   !USE wrk_nemo       ! work arrays
   USE mod_nemo
   USE mod_crs_def            ! domain for coarse grid/ 'crs' in NEMO
   !USE in_out_manager
   !USE par_kind
   !USE crslbclnk
   !USE lib_mpp


   IMPLICIT NONE

   PRIVATE

   PUBLIC crs_dom_ope
   !PUBLIC crs_dom_e3, crs_dom_sfc, crs_dom_msk, crs_dom_hgr, crs_dom_coordinates
   !PUBLIC crs_dom_facvol, crs_dom_def, crs_dom_bat

   INTERFACE crs_dom_ope
      MODULE PROCEDURE crs_dom_ope_3d, crs_dom_ope_2d
   END INTERFACE crs_dom_ope

   REAL(wp) :: r_inf = 1e+36

   !! Substitutions
   !!#  include "domzgr_substitute.h90"

   !! $Id: mod_crs.F90 5302 2015-05-28 07:11:24Z smasson $
CONTAINS


   SUBROUTINE crs_dom_msk

      INTEGER  ::  ji, jj, jk                   ! dummy loop indices
      INTEGER  ::  ijie,ijis,ijje,ijjs,ij,je_2
      REAL(wp) ::  zmask

      ! Initialize

      tmask_crs(:,:,:) = 0.0
      vmask_crs(:,:,:) = 0.0
      umask_crs(:,:,:) = 0.0
      fmask_crs(:,:,:) = 0.0


      IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
         IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
            je_2 = mje_crs(2)   ;  ij = je_2
         ENDIF
      ELSE
         je_2 = mje_crs(2)      ;  ij = mjs_crs(2)
      ENDIF
      DO jk = 1, jpkm1
         DO ji = 2, nlei_crs
            ijis = mis_crs(ji)  ;  ijie = mie_crs(ji)
            !
            zmask = 0.0
            zmask = SUM( tmask(ijis:ijie,ij:je_2,jk) )
            IF ( zmask > 0.0 ) tmask_crs(ji,2,jk) = 1.0

            zmask = 0.0
            zmask = SUM( vmask(ijis:ijie,je_2     ,jk) )
            IF ( zmask > 0.0 ) vmask_crs(ji,2,jk) = 1.0

            zmask = 0.0
            zmask = SUM(umask(ijie,ij:je_2,jk))
            IF ( zmask > 0.0 ) umask_crs(ji,2,jk) = 1.0

            fmask_crs(ji,je_2,jk) = fmask(ijie,2,jk)
         ENDDO
      ENDDO
      !
      DO jk = 1, jpkm1
         DO ji = 2, nlei_crs
            ijis = mis_crs(ji)     ;   ijie = mie_crs(ji)
            DO jj = 3, nlej_crs
               ijjs = mjs_crs(jj)  ;   ijje = mje_crs(jj)

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

      !LOLO: ADD:
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
      INTEGER :: ijis, ijjs


      SELECT CASE ( cd_type )
      CASE ( 'T' )
         DO jj =  nldj_crs, nlej_crs
            ijjs = mjs_crs(jj) + mybinctr
            DO ji = 2, nlei_crs
               ijis = mis_crs(ji) + mxbinctr
               p_gphi_crs(ji,jj) = p_gphi(ijis,ijjs)
               p_glam_crs(ji,jj) = p_glam(ijis,ijjs)
            ENDDO
         ENDDO
      CASE ( 'U' )
         DO jj =  nldj_crs, nlej_crs
            ijjs = mjs_crs(jj) + mybinctr
            DO ji = 2, nlei_crs
               ijis = mis_crs(ji)
               p_gphi_crs(ji,jj) = p_gphi(ijis,ijjs)
               p_glam_crs(ji,jj) = p_glam(ijis,ijjs)
            ENDDO
         ENDDO
      CASE ( 'V' )
         DO jj =  nldj_crs, nlej_crs
            ijjs = mjs_crs(jj)
            DO ji = 2, nlei_crs
               ijis = mis_crs(ji) + mxbinctr
               p_gphi_crs(ji,jj) = p_gphi(ijis,ijjs)
               p_glam_crs(ji,jj) = p_glam(ijis,ijjs)
            ENDDO
         ENDDO
      CASE ( 'F' )
         DO jj =  nldj_crs, nlej_crs
            ijjs = mjs_crs(jj)
            DO ji = 2, nlei_crs
               ijis = mis_crs(ji)
               p_gphi_crs(ji,jj) = p_gphi(ijis,ijjs)
               p_glam_crs(ji,jj) = p_glam(ijis,ijjs)
            ENDDO
         ENDDO
      END SELECT

      ! Retroactively add back the boundary halo cells.
      !LOLO: ADD!
      !CALL crs_lbc_lnk( p_gphi_crs, cd_type, 1.0 )
      !CALL crs_lbc_lnk( p_glam_crs, cd_type, 1.0 )

      ! Fill up jrow=1 which is zeroed out or not handled by lbc_lnk and lbc_nfd
      SELECT CASE ( cd_type )
      CASE ( 'T', 'V' )
         DO ji = 2, nlei_crs
            ijis = mis_crs(ji) + mxbinctr
            p_gphi_crs(ji,1) = p_gphi(ijis,1)
            p_glam_crs(ji,1) = p_glam(ijis,1)
         ENDDO
      CASE ( 'U', 'F' )
         DO ji = 2, nlei_crs
            ijis = mis_crs(ji)
            p_gphi_crs(ji,1) = p_gphi(ijis,1)
            p_glam_crs(ji,1) = p_glam(ijis,1)
         ENDDO
      END SELECT
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
      INTEGER :: ijie,ijje,ijrs

      !!----------------------------------------------------------------
      ! Initialize

      DO jk = 1, jpk
         DO ji = 2, nlei_crs
            ijie = mie_crs(ji)
            DO jj = nldj_crs, nlej_crs
               ijje = mje_crs(jj)   ;   ijrs =  mje_crs(jj) - mjs_crs(jj)
               ! Only for a factro 3 coarsening
               SELECT CASE ( cd_type )
               CASE ( 'T' )
                  IF( ijrs == 0 .OR. ijrs == 1 ) THEN
                     ! Si à la frontière sud on a pas assez de maille de la grille mère
                     p_e1_crs(ji,jj) = p_e1(ijie-1,ijje) * nn_factx
                     p_e2_crs(ji,jj) = p_e2(ijie-1,ijje) * nn_facty
                  ELSE
                     p_e1_crs(ji,jj) = p_e1(ijie-1,ijje-1) * nn_factx
                     p_e2_crs(ji,jj) = p_e2(ijie-1,ijje-1) * nn_facty
                  ENDIF
               CASE ( 'U' )
                  IF( ijrs == 0 .OR. ijrs == 1 ) THEN
                     ! Si à la frontière sud on a pas assez de maille de la grille mère
                     p_e1_crs(ji,jj) = p_e1(ijie,ijje) * nn_factx
                     p_e2_crs(ji,jj) = p_e2(ijie,ijje) * nn_facty
                  ELSE
                     p_e1_crs(ji,jj) = p_e1(ijie,ijje-1) * nn_factx
                     p_e2_crs(ji,jj) = p_e2(ijie,ijje-1) * nn_facty
                  ENDIF
               CASE ( 'V' )
                  p_e1_crs(ji,jj) = p_e1(ijie-1,ijje) * nn_factx
                  p_e2_crs(ji,jj) = p_e2(ijie-1,ijje) * nn_facty
               CASE ( 'F' )
                  p_e1_crs(ji,jj) = p_e1(ijie,ijje) * nn_factx
                  p_e2_crs(ji,jj) = p_e2(ijie,ijje) * nn_facty
               END SELECT
            ENDDO
         ENDDO
      ENDDO

      !LOLO: ADD:
      !CALL crs_lbc_lnk( p_e1_crs, cd_type, 1.0, pval=1.0 )
      !CALL crs_lbc_lnk( p_e2_crs, cd_type, 1.0, pval=1.0 )

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
      INTEGER                                 :: ji, jj, jk , ii, ij, je_2

      REAL(wp), DIMENSION(:,:,:), POINTER     :: zvol, zmask
      !!----------------------------------------------------------------

      ALLOCATE ( zvol(jpi,jpj,jpk), zmask(jpi,jpj,jpk) )

      p_fld1_crs(:,:,:) = 0.0
      p_fld2_crs(:,:,:) = 0.0

      DO jk = 1, jpk
         zvol(:,:,jk) =  p_e1(:,:) * p_e2(:,:) * p_e3(:,:,jk)
      ENDDO

      zmask(:,:,:) = 0.0
      IF( cd_type == 'W' ) THEN
         zmask(:,:,1) = p_mask(:,:,1)
         DO jk = 2, jpk
            zmask(:,:,jk) = p_mask(:,:,jk-1)
         ENDDO
      ELSE
         DO jk = 1, jpk
            zmask(:,:,jk) = p_mask(:,:,jk)
         ENDDO
      ENDIF

      IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
         IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
            je_2 = mje_crs(2)
            DO jk = 1, jpk
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                  p_fld1_crs(ii,2,jk) =  zvol(ji,je_2  ,jk) + zvol(ji+1,je_2  ,jk) + zvol(ji+2,je_2  ,jk)  &
                     &                 + zvol(ji,je_2-1,jk) + zvol(ji+1,je_2-1,jk) + zvol(ji+2,je_2-1,jk)
                  !
                  zdAm =  zvol(ji  ,je_2,jk) * zmask(ji  ,je_2,jk)  &
                     &   + zvol(ji+1,je_2,jk) * zmask(ji+1,je_2,jk)  &
                     &   + zvol(ji+2,je_2,jk) * zmask(ji+2,je_2,jk)
                  !
                  p_fld2_crs(ii,2,jk) = zdAm / p_fld1_crs(ii,2,jk)
               ENDDO
            ENDDO
         ENDIF
      ELSE
         je_2 = mjs_crs(2)
         DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               p_fld1_crs(ii,2,jk) =  zvol(ji,je_2  ,jk) + zvol(ji+1,je_2  ,jk) + zvol(ji+2,je_2  ,jk)  &
                  &                + zvol(ji,je_2+1,jk) + zvol(ji+1,je_2+1,jk) + zvol(ji+2,je_2+1,jk)  &
                  &                + zvol(ji,je_2+2,jk) + zvol(ji+1,je_2+2,jk) + zvol(ji+2,je_2+2,jk)
               !
               zdAm = zvol(ji  ,je_2  ,jk) * zmask(ji  ,je_2  ,jk)  &
                  &  + zvol(ji+1,je_2  ,jk) * zmask(ji+1,je_2  ,jk)  &
                  &  + zvol(ji+2,je_2  ,jk) * zmask(ji+2,je_2  ,jk)  &
                  &  + zvol(ji  ,je_2+1,jk) * zmask(ji  ,je_2+1,jk)  &
                  &  + zvol(ji+1,je_2+1,jk) * zmask(ji+1,je_2+1,jk)  &
                  &  + zvol(ji+2,je_2+1,jk) * zmask(ji+2,je_2+1,jk)  &
                  &  + zvol(ji  ,je_2+2,jk) * zmask(ji  ,je_2+2,jk)  &
                  &  + zvol(ji+1,je_2+2,jk) * zmask(ji+1,je_2+2,jk)  &
                  &  + zvol(ji+2,je_2+2,jk) * zmask(ji+2,je_2+2,jk)
               !
               p_fld2_crs(ii,2,jk) = zdAm / p_fld1_crs(ii,2,jk)
            ENDDO
         ENDDO
      ENDIF

      DO jk = 1, jpk
         DO jj  = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
               ii  = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
               ij  = ( jj - njstr ) * rfacty_r + 3
               !
               p_fld1_crs(ii,ij,jk) =  zvol(ji,jj  ,jk) + zvol(ji+1,jj  ,jk) + zvol(ji+2,jj  ,jk)  &
                  &                 + zvol(ji,jj+1,jk) + zvol(ji+1,jj+1,jk) + zvol(ji+2,jj+1,jk)  &
                  &                 + zvol(ji,jj+2,jk) + zvol(ji+1,jj+2,jk) + zvol(ji+2,jj+2,jk)
               !
               zdAm =  zvol(ji  ,jj  ,jk) * zmask(ji  ,jj  ,jk)  &
                  &   + zvol(ji+1,jj  ,jk) * zmask(ji+1,jj  ,jk)  &
                  &   + zvol(ji+2,jj  ,jk) * zmask(ji+2,jj  ,jk)  &
                  &   + zvol(ji  ,jj+1,jk) * zmask(ji  ,jj+1,jk)  &
                  &   + zvol(ji+1,jj+1,jk) * zmask(ji+1,jj+1,jk)  &
                  &   + zvol(ji+2,jj+1,jk) * zmask(ji+2,jj+1,jk)  &
                  &   + zvol(ji  ,jj+2,jk) * zmask(ji  ,jj+2,jk)  &
                  &   + zvol(ji+1,jj+2,jk) * zmask(ji+1,jj+2,jk)  &
                  &   + zvol(ji+2,jj+2,jk) * zmask(ji+2,jj+2,jk)
               !
               p_fld2_crs(ii,ij,jk) = zdAm / p_fld1_crs(ii,ij,jk)
            ENDDO
         ENDDO
      ENDDO
      !                                             !  Retroactively add back the boundary halo cells.
      !LOLO: ADD:
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
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in)           :: p_fld   ! T, U, V or W on parent grid
      CHARACTER(len=3),                         INTENT(in)           :: cd_op    ! Operation SUM, MAX or MIN
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
      INTEGER  :: ii, ij, ijie, ijje, je_2
      REAL(wp) :: zflcrs, zsfcrs
      REAL(wp), DIMENSION(:,:,:), POINTER :: zsurf, zsurfmsk, zmask
      !!----------------------------------------------------------------

      p_fld_crs(:,:,:) = 0.0

      SELECT CASE ( cd_op )

      CASE ( 'VOL' )

         ALLOCATE ( zsurf(jpi,jpj,jpk), zsurfmsk(jpi,jpj,jpk) )

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )
            IF( cd_type == 'T' ) THEN
               DO jk = 1, jpk
                  zsurf   (:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk) *  p_mask(:,:,jk)
                  zsurfmsk(:,:,jk) =  zsurf(:,:,jk)
               ENDDO
            ELSE
               zsurf   (:,:,1) =  p_e12(:,:) * p_e3(:,:,1)
               zsurfmsk(:,:,1) =  zsurf(:,:,1) *  p_mask(:,:,1)
               DO jk = 2, jpk
                  zsurf   (:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk)
                  zsurfmsk(:,:,jk) =  zsurf(:,:,jk) * p_mask(:,:,jk-1)
               ENDDO
            ENDIF

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO jk = 1, jpk
                     DO ji = nistr, niend, nn_factx
                        ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                        zflcrs =  p_fld(ji  ,je_2,jk) * zsurfmsk(ji  ,je_2,jk)   &
                           &     + p_fld(ji+1,je_2,jk) * zsurfmsk(ji+1,je_2,jk)   &
                           &     + p_fld(ji+2,je_2,jk) * zsurfmsk(ji+2,je_2,jk)

                        zsfcrs =  zsurf(ji,je_2,jk) + zsurf(ji+1,je_2,jk) + zsurf(ji+2,je_2,jk)
                        !
                        p_fld_crs(ii,2,jk) = zflcrs
                        IF( zsfcrs /= 0.0 )  p_fld_crs(ii,2,jk) = zflcrs / zsfcrs
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO jk = 1, jpk
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     zflcrs =  p_fld(ji  ,je_2  ,jk) * zsurfmsk(ji  ,je_2  ,jk) &
                        &     + p_fld(ji+1,je_2  ,jk) * zsurfmsk(ji+1,je_2  ,jk) &
                        &     + p_fld(ji+2,je_2  ,jk) * zsurfmsk(ji+2,je_2  ,jk) &
                        &     + p_fld(ji  ,je_2+1,jk) * zsurfmsk(ji  ,je_2+1,jk) &
                        &     + p_fld(ji+1,je_2+1,jk) * zsurfmsk(ji+1,je_2+1,jk) &
                        &     + p_fld(ji+2,je_2+1,jk) * zsurfmsk(ji+2,je_2+1,jk) &
                        &     + p_fld(ji  ,je_2+2,jk) * zsurfmsk(ji  ,je_2+2,jk) &
                        &     + p_fld(ji+1,je_2+2,jk) * zsurfmsk(ji+1,je_2+2,jk) &
                        &     + p_fld(ji+2,je_2+2,jk) * zsurfmsk(ji+2,je_2+2,jk)

                     zsfcrs =  zsurf(ji,je_2  ,jk) + zsurf(ji+1,je_2  ,jk) + zsurf(ji+2,je_2  ,jk) &
                        &     + zsurf(ji,je_2+1,jk) + zsurf(ji+1,je_2+1,jk) + zsurf(ji+2,je_2+1,jk) &
                        &     + zsurf(ji,je_2+2,jk) + zsurf(ji+1,je_2+2,jk) + zsurf(ji+2,je_2+2,jk)
                     !
                     p_fld_crs(ii,2,jk) = zflcrs
                     IF( zsfcrs /= 0.0 )  p_fld_crs(ii,2,jk) = zflcrs / zsfcrs
                  ENDDO
               ENDDO
            ENDIF
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                     ij = ( jj - njstr ) * rfacty_r + 3
                     zflcrs =  p_fld(ji  ,jj  ,jk) * zsurfmsk(ji  ,jj  ,jk) &
                        &     + p_fld(ji+1,jj  ,jk) * zsurfmsk(ji+1,jj  ,jk) &
                        &     + p_fld(ji+2,jj  ,jk) * zsurfmsk(ji+2,jj  ,jk) &
                        &     + p_fld(ji  ,jj+1,jk) * zsurfmsk(ji  ,jj+1,jk) &
                        &     + p_fld(ji+1,jj+1,jk) * zsurfmsk(ji+1,jj+1,jk) &
                        &     + p_fld(ji+2,jj+1,jk) * zsurfmsk(ji+2,jj+1,jk) &
                        &     + p_fld(ji  ,jj+2,jk) * zsurfmsk(ji  ,jj+2,jk) &
                        &     + p_fld(ji+1,jj+2,jk) * zsurfmsk(ji+1,jj+2,jk) &
                        &     + p_fld(ji+2,jj+2,jk) * zsurfmsk(ji+2,jj+2,jk)

                     zsfcrs =  zsurf(ji,jj  ,jk) + zsurf(ji+1,jj  ,jk) + zsurf(ji+2,jj  ,jk) &
                        &     + zsurf(ji,jj+1,jk) + zsurf(ji+1,jj+1,jk) + zsurf(ji+2,jj+1,jk) &
                        &     + zsurf(ji,jj+2,jk) + zsurf(ji+1,jj+2,jk) + zsurf(ji+2,jj+2,jk)
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     IF( zsfcrs /= 0.0 )  p_fld_crs(ii,ij,jk) = zflcrs / zsfcrs
                  ENDDO
               ENDDO
            ENDDO
         CASE DEFAULT
            STOP
         END SELECT

         DEALLOCATE ( zsurf, zsurfmsk )

      CASE ( 'SUM' )

         ALLOCATE ( zsurfmsk(jpi,jpj,jpk) )

         SELECT CASE ( cd_type )
         CASE( 'W' )
            IF( PRESENT( p_e3 ) ) THEN
               zsurfmsk(:,:,1) =  p_e12(:,:) * p_e3(:,:,1) * p_mask(:,:,1)
               DO jk = 2, jpk
                  zsurfmsk(:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk) * p_mask(:,:,jk-1)
               ENDDO
            ELSE
               zsurfmsk(:,:,1) =  p_e12(:,:) * p_mask(:,:,1)
               DO jk = 2, jpk
                  zsurfmsk(:,:,jk) =  p_e12(:,:) * p_mask(:,:,jk-1)
               ENDDO
            ENDIF
         CASE DEFAULT
            IF( PRESENT( p_e3 ) ) THEN
               DO jk = 1, jpk
                  zsurfmsk(:,:,jk) =  p_e12(:,:) * p_e3(:,:,jk) * p_mask(:,:,jk)
               ENDDO
            ELSE
               DO jk = 1, jpk
                  zsurfmsk(:,:,jk) =  p_e12(:,:) * p_mask(:,:,jk)
               ENDDO
            ENDIF
         END SELECT

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO jk = 1, jpk
                     DO ji = nistr, niend, nn_factx
                        ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                        zflcrs  =  p_fld(ji  ,je_2,jk) * zsurfmsk(ji  ,je_2,jk) &
                           &      + p_fld(ji+1,je_2,jk) * zsurfmsk(ji+1,je_2,jk) &
                           &      + p_fld(ji+2,je_2,jk) * zsurfmsk(ji+2,je_2,jk)
                        !
                        p_fld_crs(ii,2,jk) = zflcrs
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO jk = 1, jpk
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     zflcrs  =  p_fld(ji  ,je_2  ,jk) * zsurfmsk(ji  ,je_2  ,jk)  &
                        &      + p_fld(ji+1,je_2  ,jk) * zsurfmsk(ji+1,je_2  ,jk)  &
                        &      + p_fld(ji+2,je_2  ,jk) * zsurfmsk(ji+2,je_2  ,jk)  &
                        &      + p_fld(ji  ,je_2+1,jk) * zsurfmsk(ji  ,je_2+1,jk)  &
                        &      + p_fld(ji+1,je_2+1,jk) * zsurfmsk(ji+1,je_2+1,jk)  &
                        &      + p_fld(ji+2,je_2+1,jk) * zsurfmsk(ji+2,je_2+1,jk)  &
                        &      + p_fld(ji  ,je_2+2,jk) * zsurfmsk(ji  ,je_2+2,jk)  &
                        &      + p_fld(ji+1,je_2+2,jk) * zsurfmsk(ji+1,je_2+2,jk)  &
                        &      + p_fld(ji+2,je_2+2,jk) * zsurfmsk(ji+2,je_2+2,jk)
                     !
                     p_fld_crs(ii,2,jk) = zflcrs
                  ENDDO
               ENDDO
            ENDIF
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii  = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                     ij  = ( jj - njstr ) * rfacty_r + 3
                     zflcrs  =  p_fld(ji  ,jj  ,jk) * zsurfmsk(ji  ,jj  ,jk)  &
                        &      + p_fld(ji+1,jj  ,jk) * zsurfmsk(ji+1,jj  ,jk)  &
                        &      + p_fld(ji+2,jj  ,jk) * zsurfmsk(ji+2,jj  ,jk)  &
                        &      + p_fld(ji  ,jj+1,jk) * zsurfmsk(ji  ,jj+1,jk)  &
                        &      + p_fld(ji+1,jj+1,jk) * zsurfmsk(ji+1,jj+1,jk)  &
                        &      + p_fld(ji+2,jj+1,jk) * zsurfmsk(ji+2,jj+1,jk)  &
                        &      + p_fld(ji  ,jj+2,jk) * zsurfmsk(ji  ,jj+2,jk)  &
                        &      + p_fld(ji+1,jj+2,jk) * zsurfmsk(ji+1,jj+2,jk)  &
                        &      + p_fld(ji+2,jj+2,jk) * zsurfmsk(ji+2,jj+2,jk)
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO

         CASE( 'V' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  ijje = mje_crs(2)
               ENDIF
            ELSE
               ijje = mjs_crs(2)
            ENDIF
            !
            DO jk = 1, jpk
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  zflcrs  =  p_fld(ji  ,ijje,jk) * zsurfmsk(ji  ,ijje,jk) &
                     &      + p_fld(ji+1,ijje,jk) * zsurfmsk(ji+1,ijje,jk) &
                     &      + p_fld(ji+2,ijje,jk) * zsurfmsk(ji+2,ijje,jk)
                  !
                  p_fld_crs(ii,2,jk) = zflcrs
               ENDDO
            ENDDO
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                     ij   = ( jj - njstr ) * rfacty_r + 3
                     ijje = mje_crs(ij)
                     zflcrs  =  p_fld(ji  ,ijje,jk) * zsurfmsk(ji  ,ijje,jk) &
                        &      + p_fld(ji+1,ijje,jk) * zsurfmsk(ji+1,ijje,jk) &
                        &      + p_fld(ji+2,ijje,jk) * zsurfmsk(ji+2,ijje,jk)
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO

         CASE( 'U' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO jk = 1, jpk
                     DO ji = nistr, niend, nn_factx
                        ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                        ijie = mie_crs(ii)
                        zflcrs  =  p_fld(ijie,je_2,jk) * zsurfmsk(ijie,je_2,jk)
                        p_fld_crs(ii,2,jk) = zflcrs
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO jk = 1, jpk
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ijie = mie_crs(ii)
                     zflcrs =  p_fld(ijie,je_2  ,jk) * zsurfmsk(ijie,je_2  ,jk)  &
                        &     + p_fld(ijie,je_2+1,jk) * zsurfmsk(ijie,je_2+1,jk)  &
                        &     + p_fld(ijie,je_2+2,jk) * zsurfmsk(ijie,je_2+2,jk)

                     p_fld_crs(ii,2,jk) = zflcrs
                  ENDDO
               ENDDO
            ENDIF
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ij   = ( jj - njstr ) * rfacty_r + 3
                     ijie = mie_crs(ii)
                     zflcrs =  p_fld(ijie,jj  ,jk) * zsurfmsk(ijie,jj  ,jk)  &
                        &    + p_fld(ijie,jj+1,jk) * zsurfmsk(ijie,jj+1,jk)  &
                        &    + p_fld(ijie,jj+2,jk) * zsurfmsk(ijie,jj+2,jk)
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
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

         SELECT CASE ( cd_type )
         CASE( 'W' )
            zmask(:,:,1) = p_mask(:,:,1)
            DO jk = 2, jpk
               zmask(:,:,jk) = p_mask(:,:,jk-1)
            ENDDO
         CASE ( 'T' )
            DO jk = 1, jpk
               zmask(:,:,jk) = p_mask(:,:,jk)
            ENDDO
         END SELECT

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO jk = 1, jpk
                     DO ji = nistr, niend, nn_factx
                        ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                        zflcrs =  &
                           & MAX( p_fld(ji  ,je_2,jk) * zmask(ji  ,je_2,jk) - ( 1.- zmask(ji  ,je_2,jk) ) * r_inf ,  &
                           &      p_fld(ji+1,je_2,jk) * zmask(ji+1,je_2,jk) - ( 1.- zmask(ji+1,je_2,jk) ) * r_inf ,  &
                           &      p_fld(ji+2,je_2,jk) * zmask(ji+2,je_2,jk) - ( 1.- zmask(ji+2,je_2,jk) ) * r_inf  )
                        !
                        p_fld_crs(ii,2,jk) = zflcrs
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO jk = 1, jpk
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     zflcrs =  &
                        & MAX( p_fld(ji  ,je_2  ,jk) * zmask(ji  ,je_2  ,jk) - ( 1.- zmask(ji  ,je_2  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,je_2  ,jk) * zmask(ji+1,je_2  ,jk) - ( 1.- zmask(ji+1,je_2  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,je_2  ,jk) * zmask(ji+2,je_2  ,jk) - ( 1.- zmask(ji+2,je_2  ,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,je_2+1,jk) * zmask(ji  ,je_2+1,jk) - ( 1.- zmask(ji  ,je_2+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,je_2+1,jk) * zmask(ji+1,je_2+1,jk) - ( 1.- zmask(ji+1,je_2+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,je_2+1,jk) * zmask(ji+2,je_2+1,jk) - ( 1.- zmask(ji+2,je_2+1,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,je_2+2,jk) * zmask(ji  ,je_2+2,jk) - ( 1.- zmask(ji  ,je_2+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,je_2+2,jk) * zmask(ji+1,je_2+2,jk) - ( 1.- zmask(ji+1,je_2+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,je_2+2,jk) * zmask(ji+2,je_2+2,jk) - ( 1.- zmask(ji+2,je_2+2,jk) ) * r_inf   )
                     !
                     p_fld_crs(ii,2,jk) = zflcrs
                  ENDDO
               ENDDO
            ENDIF
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii  = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                     ij  = ( jj - njstr ) * rfacty_r + 3
                     zflcrs =  &
                        & MAX( p_fld(ji  ,jj  ,jk) * zmask(ji  ,jj  ,jk) - ( 1.- zmask(ji  ,jj  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,jj  ,jk) * zmask(ji+1,jj  ,jk) - ( 1.- zmask(ji+1,jj  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,jj  ,jk) * zmask(ji+2,jj  ,jk) - ( 1.- zmask(ji+2,jj  ,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,jj+1,jk) * zmask(ji  ,jj+1,jk) - ( 1.- zmask(ji  ,jj+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,jj+1,jk) * zmask(ji+1,jj+1,jk) - ( 1.- zmask(ji+1,jj+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,jj+1,jk) * zmask(ji+2,jj+1,jk) - ( 1.- zmask(ji+2,jj+1,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,jj+2,jk) * zmask(ji  ,jj+2,jk) - ( 1.- zmask(ji  ,jj+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,jj+2,jk) * zmask(ji+1,jj+2,jk) - ( 1.- zmask(ji+1,jj+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,jj+2,jk) * zmask(ji+2,jj+2,jk) - ( 1.- zmask(ji+2,jj+2,jk) ) * r_inf   )
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO

         CASE( 'V' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  ijje = mje_crs(2)
               ENDIF
            ELSE
               ijje = mjs_crs(2)
            ENDIF

            DO jk = 1, jpk
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  zflcrs = &
                     & MAX( p_fld(ji  ,ijje,jk) * p_mask(ji  ,ijje,jk) - ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                     &      p_fld(ji+1,ijje,jk) * p_mask(ji+1,ijje,jk) - ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                     &      p_fld(ji+2,ijje,jk) * p_mask(ji+2,ijje,jk) - ( 1.- p_mask(ji,ijje,jk) ) * r_inf )
                  !
                  p_fld_crs(ii,2,jk) = zflcrs
               ENDDO
            ENDDO
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii  = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                     ij  = ( jj - njstr ) * rfacty_r + 3
                     ijje = mje_crs(ij)
                     !
                     zflcrs = &
                        & MAX( p_fld(ji  ,ijje,jk) * p_mask(ji  ,ijje,jk) - ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,ijje,jk) * p_mask(ji+1,ijje,jk) - ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,ijje,jk) * p_mask(ji+2,ijje,jk) - ( 1.- p_mask(ji,ijje,jk) ) * r_inf )
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO


         CASE( 'U' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO jk = 1, jpk
                     DO ji = nistr, niend, nn_factx
                        ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                        ijie = mie_crs(ii)
                        zflcrs = p_fld(ijie,je_2,jk) * p_mask(ijie,je_2,jk) - ( 1.- p_mask(ijie,je_2,jk) ) * r_inf
                        !
                        p_fld_crs(ii,2,jk) = zflcrs
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO jk = 1, jpk
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ijie = mie_crs(ii)
                     zflcrs = &
                        & MAX( p_fld(ijie,je_2  ,jk) * p_mask(ijie,je_2  ,jk) - ( 1.- p_mask(ijie,je_2,jk) ) * r_inf ,  &
                        &      p_fld(ijie,je_2+1,jk) * p_mask(ijie,je_2+1,jk) - ( 1.- p_mask(ijie,je_2,jk) ) * r_inf ,  &
                        &      p_fld(ijie,je_2+2,jk) * p_mask(ijie,je_2+2,jk) - ( 1.- p_mask(ijie,je_2,jk) ) * r_inf  )
                     !
                     p_fld_crs(ii,2,jk) = zflcrs
                  ENDDO
               ENDDO
            ENDIF
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ij   = ( jj - njstr ) * rfacty_r + 3
                     ijie = mie_crs(ii)
                     zflcrs =  &
                        & MAX( p_fld(ijie,jj  ,jk) * p_mask(ijie,jj  ,jk) - ( 1.- p_mask(ijie,jj,jk) ) * r_inf ,  &
                        &      p_fld(ijie,jj+1,jk) * p_mask(ijie,jj+1,jk) - ( 1.- p_mask(ijie,jj,jk) ) * r_inf ,  &
                        &      p_fld(ijie,jj+2,jk) * p_mask(ijie,jj+2,jk) - ( 1.- p_mask(ijie,jj,jk) ) * r_inf  )
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO

         END SELECT

         DEALLOCATE ( zmask )

      CASE ( 'MIN' )      !   Search the min of unmasked grid cells

         ALLOCATE ( zmask(jpi,jpj,jpk) )

         SELECT CASE ( cd_type )
         CASE( 'W' )
            zmask(:,:,1) = p_mask(:,:,1)
            DO jk = 2, jpk
               zmask(:,:,jk) = p_mask(:,:,jk-1)
            ENDDO
         CASE ( 'T' )
            DO jk = 1, jpk
               zmask(:,:,jk) = p_mask(:,:,jk)
            ENDDO
         END SELECT

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO jk = 1, jpk
                     DO ji = nistr, niend, nn_factx
                        ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                        zflcrs =  &
                           & MIN( p_fld(ji  ,je_2,jk) * zmask(ji  ,je_2,jk) + ( 1.- zmask(ji  ,je_2,jk) ) * r_inf ,  &
                           &      p_fld(ji+1,je_2,jk) * zmask(ji+1,je_2,jk) + ( 1.- zmask(ji+1,je_2,jk) ) * r_inf ,  &
                           &      p_fld(ji+2,je_2,jk) * zmask(ji+2,je_2,jk) + ( 1.- zmask(ji+2,je_2,jk) ) * r_inf  )
                        !
                        p_fld_crs(ii,2,jk) = zflcrs
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO jk = 1, jpk
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     zflcrs =  &
                        & MIN( p_fld(ji  ,je_2  ,jk) * zmask(ji  ,je_2  ,jk) + ( 1.- zmask(ji  ,je_2  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,je_2  ,jk) * zmask(ji+1,je_2  ,jk) + ( 1.- zmask(ji+1,je_2  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,je_2  ,jk) * zmask(ji+2,je_2  ,jk) + ( 1.- zmask(ji+2,je_2  ,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,je_2+1,jk) * zmask(ji  ,je_2+1,jk) + ( 1.- zmask(ji  ,je_2+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,je_2+1,jk) * zmask(ji+1,je_2+1,jk) + ( 1.- zmask(ji+1,je_2+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,je_2+1,jk) * zmask(ji+2,je_2+1,jk) + ( 1.- zmask(ji+2,je_2+1,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,je_2+2,jk) * zmask(ji  ,je_2+2,jk) + ( 1.- zmask(ji  ,je_2+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,je_2+2,jk) * zmask(ji+1,je_2+2,jk) + ( 1.- zmask(ji+1,je_2+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,je_2+2,jk) * zmask(ji+2,je_2+2,jk) + ( 1.- zmask(ji+2,je_2+2,jk) ) * r_inf   )
                     !
                     p_fld_crs(ii,2,jk) = zflcrs
                  ENDDO
               ENDDO
            ENDIF
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii  = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                     ij  = ( jj - njstr ) * rfacty_r + 3
                     zflcrs =  &
                        & MIN( p_fld(ji  ,jj  ,jk) * zmask(ji  ,jj  ,jk) + ( 1.- zmask(ji  ,jj  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,jj  ,jk) * zmask(ji+1,jj  ,jk) + ( 1.- zmask(ji+1,jj  ,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,jj  ,jk) * zmask(ji+2,jj  ,jk) + ( 1.- zmask(ji+2,jj  ,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,jj+1,jk) * zmask(ji  ,jj+1,jk) + ( 1.- zmask(ji  ,jj+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,jj+1,jk) * zmask(ji+1,jj+1,jk) + ( 1.- zmask(ji+1,jj+1,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,jj+1,jk) * zmask(ji+2,jj+1,jk) + ( 1.- zmask(ji+2,jj+1,jk) ) * r_inf ,  &
                        &      p_fld(ji  ,jj+2,jk) * zmask(ji  ,jj+2,jk) + ( 1.- zmask(ji  ,jj+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,jj+2,jk) * zmask(ji+1,jj+2,jk) + ( 1.- zmask(ji+1,jj+2,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,jj+2,jk) * zmask(ji+2,jj+2,jk) + ( 1.- zmask(ji+2,jj+2,jk) ) * r_inf   )
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO

         CASE( 'V' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  ijje = mje_crs(2)
               ENDIF
            ELSE
               ijje = mjs_crs(2)
            ENDIF

            DO jk = 1, jpk
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  zflcrs = &
                     & MIN( p_fld(ji  ,ijje,jk) * p_mask(ji  ,ijje,jk) + ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                     &      p_fld(ji+1,ijje,jk) * p_mask(ji+1,ijje,jk) + ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                     &      p_fld(ji+2,ijje,jk) * p_mask(ji+2,ijje,jk) + ( 1.- p_mask(ji,ijje,jk) ) * r_inf )
                  !
                  p_fld_crs(ii,2,jk) = zflcrs
               ENDDO
            ENDDO
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii  = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
                     ij  = ( jj - njstr ) * rfacty_r + 3
                     ijje = mje_crs(ij)
                     zflcrs = &
                        & MIN( p_fld(ji  ,ijje,jk) * p_mask(ji  ,ijje,jk) + ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                        &      p_fld(ji+1,ijje,jk) * p_mask(ji+1,ijje,jk) + ( 1.- p_mask(ji,ijje,jk) ) * r_inf ,  &
                        &      p_fld(ji+2,ijje,jk) * p_mask(ji+2,ijje,jk) + ( 1.- p_mask(ji,ijje,jk) ) * r_inf )
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO


         CASE( 'U' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO jk = 1, jpk
                     DO ji = nistr, niend, nn_factx
                        ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                        ijie = mie_crs(ii)
                        zflcrs = p_fld(ijie,je_2,jk) * p_mask(ijie,je_2,jk) + ( 1.- p_mask(ijie,je_2,jk) ) * r_inf
                        !
                        p_fld_crs(ii,2,jk) = zflcrs
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO jk = 1, jpk
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ijie = mie_crs(ii)
                     zflcrs = &
                        & MIN( p_fld(ijie,je_2  ,jk) * p_mask(ijie,je_2  ,jk) + ( 1.- p_mask(ijie,je_2,jk) ) * r_inf ,  &
                        &      p_fld(ijie,je_2+1,jk) * p_mask(ijie,je_2+1,jk) + ( 1.- p_mask(ijie,je_2,jk) ) * r_inf ,  &
                        &      p_fld(ijie,je_2+2,jk) * p_mask(ijie,je_2+2,jk) + ( 1.- p_mask(ijie,je_2,jk) ) * r_inf  )
                     !
                     p_fld_crs(ii,2,jk) = zflcrs
                  ENDDO
               ENDDO
            ENDIF
            !
            DO jk = 1, jpk
               DO jj  = njstr, njend, nn_facty
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ij   = ( jj - njstr ) * rfacty_r + 3
                     ijie = mie_crs(ii)
                     zflcrs = &
                        & MIN( p_fld(ijie,jj  ,jk) * p_mask(ijie,jj  ,jk) + ( 1.- p_mask(ijie,jj,jk) ) * r_inf ,  &
                        &      p_fld(ijie,jj+1,jk) * p_mask(ijie,jj+1,jk) + ( 1.- p_mask(ijie,jj,jk) ) * r_inf ,  &
                        &      p_fld(ijie,jj+2,jk) * p_mask(ijie,jj+2,jk) + ( 1.- p_mask(ijie,jj,jk) ) * r_inf  )
                     !
                     p_fld_crs(ii,ij,jk) = zflcrs
                     !
                  ENDDO
               ENDDO
            ENDDO

         END SELECT
         !
         DEALLOCATE ( zmask )
         !
      END SELECT
      !
      !LOLO: ADD:
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
      CHARACTER(len=3),                         INTENT(in)           :: cd_op    ! Operation SUM, MAX or MIN
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
      INTEGER  :: ijie, ijje, ii, ij, je_2
      REAL(wp) :: zflcrs, zsfcrs
      REAL(wp), DIMENSION(:,:), POINTER :: zsurfmsk

      !!----------------------------------------------------------------

      p_fld_crs(:,:) = 0.0

      SELECT CASE ( cd_op )

      CASE ( 'VOL' )

         ALLOCATE ( zsurfmsk(jpi,jpj) )
         zsurfmsk(:,:) =  p_e12(:,:) * p_e3(:,:,1) * p_mask(:,:,1)

         IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
            IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
               je_2 = mje_crs(2)
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  zflcrs =  p_fld(ji  ,je_2) * zsurfmsk(ji  ,je_2)   &
                     &     + p_fld(ji+1,je_2) * zsurfmsk(ji+1,je_2)   &
                     &     + p_fld(ji+2,je_2) * zsurfmsk(ji+2,je_2)

                  zsfcrs =  zsurfmsk(ji,je_2) + zsurfmsk(ji+1,je_2) + zsurfmsk(ji+2,je_2)
                  !
                  p_fld_crs(ii,2) = zflcrs
                  IF( zsfcrs /= 0.0 )  p_fld_crs(ii,2) = zflcrs / zsfcrs
               ENDDO
            ENDIF
         ELSE
            je_2 = mjs_crs(2)
            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               zflcrs =  p_fld(ji  ,je_2  ) * zsurfmsk(ji  ,je_2  ) &
                  &     + p_fld(ji+1,je_2  ) * zsurfmsk(ji+1,je_2  ) &
                  &     + p_fld(ji+2,je_2  ) * zsurfmsk(ji+2,je_2  ) &
                  &     + p_fld(ji  ,je_2+1) * zsurfmsk(ji  ,je_2+1) &
                  &     + p_fld(ji+1,je_2+1) * zsurfmsk(ji+1,je_2+1) &
                  &     + p_fld(ji+2,je_2+1) * zsurfmsk(ji+2,je_2+1) &
                  &     + p_fld(ji  ,je_2+2) * zsurfmsk(ji  ,je_2+2) &
                  &     + p_fld(ji+1,je_2+2) * zsurfmsk(ji+1,je_2+2) &
                  &     + p_fld(ji+2,je_2+2) * zsurfmsk(ji+2,je_2+2)

               zsfcrs =  zsurfmsk(ji,je_2  ) + zsurfmsk(ji+1,je_2  ) + zsurfmsk(ji+2,je_2  ) &
                  &     + zsurfmsk(ji,je_2+1) + zsurfmsk(ji+1,je_2+1) + zsurfmsk(ji+2,je_2+1) &
                  &     + zsurfmsk(ji,je_2+2) + zsurfmsk(ji+1,je_2+2) + zsurfmsk(ji+2,je_2+2)
               !
               p_fld_crs(ii,2) = zflcrs
               IF( zsfcrs /= 0.0 )  p_fld_crs(ii,2) = zflcrs / zsfcrs
            ENDDO
         ENDIF
         !
         DO jj  = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
               ii  = ( ji - mis_crs(2) ) * rfactx_r + 2                 ! cordinate in parent grid
               ij  = ( jj - njstr ) * rfacty_r + 3
               zflcrs =  p_fld(ji  ,jj  ) * zsurfmsk(ji  ,jj  ) &
                  &     + p_fld(ji+1,jj  ) * zsurfmsk(ji+1,jj  ) &
                  &     + p_fld(ji+2,jj  ) * zsurfmsk(ji+2,jj  ) &
                  &     + p_fld(ji  ,jj+1) * zsurfmsk(ji  ,jj+1) &
                  &     + p_fld(ji+1,jj+1) * zsurfmsk(ji+1,jj+1) &
                  &     + p_fld(ji+2,jj+1) * zsurfmsk(ji+2,jj+1) &
                  &     + p_fld(ji  ,jj+2) * zsurfmsk(ji  ,jj+2) &
                  &     + p_fld(ji+1,jj+2) * zsurfmsk(ji+1,jj+2) &
                  &     + p_fld(ji+2,jj+2) * zsurfmsk(ji+2,jj+2)

               zsfcrs =  zsurfmsk(ji,jj  ) + zsurfmsk(ji+1,jj  ) + zsurfmsk(ji+2,jj  ) &
                  &     + zsurfmsk(ji,jj+1) + zsurfmsk(ji+1,jj+1) + zsurfmsk(ji+2,jj+1) &
                  &     + zsurfmsk(ji,jj+2) + zsurfmsk(ji+1,jj+2) + zsurfmsk(ji+2,jj+2)
               !
               p_fld_crs(ii,ij) = zflcrs
               IF( zsfcrs /= 0.0 )  p_fld_crs(ii,ij) = zflcrs / zsfcrs
            ENDDO
         ENDDO

         DEALLOCATE ( zsurfmsk )

      CASE ( 'SUM' )

         ALLOCATE ( zsurfmsk(jpi,jpj) )
         IF( PRESENT( p_e3 ) ) THEN
            zsurfmsk(:,:) =  p_e12(:,:) * p_e3(:,:,1) * p_mask(:,:,1)
         ELSE
            zsurfmsk(:,:) =  p_e12(:,:) * p_mask(:,:,1)
         ENDIF

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     zflcrs  =  p_fld(ji  ,je_2) * zsurfmsk(ji  ,je_2) &
                        &      + p_fld(ji+1,je_2) * zsurfmsk(ji+1,je_2) &
                        &      + p_fld(ji+2,je_2) * zsurfmsk(ji+2,je_2)
                     !
                     p_fld_crs(ii,2) = zflcrs
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  zflcrs  =  p_fld(ji  ,je_2  ) * zsurfmsk(ji  ,je_2  )  &
                     &      + p_fld(ji+1,je_2  ) * zsurfmsk(ji+1,je_2  )  &
                     &      + p_fld(ji+2,je_2  ) * zsurfmsk(ji+2,je_2  )  &
                     &      + p_fld(ji  ,je_2+1) * zsurfmsk(ji  ,je_2+1)  &
                     &      + p_fld(ji+1,je_2+1) * zsurfmsk(ji+1,je_2+1)  &
                     &      + p_fld(ji+2,je_2+1) * zsurfmsk(ji+2,je_2+1)  &
                     &      + p_fld(ji  ,je_2+2) * zsurfmsk(ji  ,je_2+2)  &
                     &      + p_fld(ji+1,je_2+2) * zsurfmsk(ji+1,je_2+2)  &
                     &      + p_fld(ji+2,je_2+2) * zsurfmsk(ji+2,je_2+2)
                  !
                  p_fld_crs(ii,2) = zflcrs
               ENDDO
            ENDIF
            !
            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  zflcrs  =  p_fld(ji  ,jj  ) * zsurfmsk(ji  ,jj  )  &
                     &      + p_fld(ji+1,jj  ) * zsurfmsk(ji+1,jj  )  &
                     &      + p_fld(ji+2,jj  ) * zsurfmsk(ji+2,jj  )  &
                     &      + p_fld(ji  ,jj+1) * zsurfmsk(ji  ,jj+1)  &
                     &      + p_fld(ji+1,jj+1) * zsurfmsk(ji+1,jj+1)  &
                     &      + p_fld(ji+2,jj+1) * zsurfmsk(ji+2,jj+1)  &
                     &      + p_fld(ji  ,jj+2) * zsurfmsk(ji  ,jj+2)  &
                     &      + p_fld(ji+1,jj+2) * zsurfmsk(ji+1,jj+2)  &
                     &      + p_fld(ji+2,jj+2) * zsurfmsk(ji+2,jj+2)
                  !
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         CASE( 'V' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  ijje = mje_crs(2)
               ENDIF
            ELSE
               ijje = mjs_crs(2)
            ENDIF

            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               zflcrs  =  p_fld(ji  ,ijje) * zsurfmsk(ji  ,ijje) &
                  &      + p_fld(ji+1,ijje) * zsurfmsk(ji+1,ijje) &
                  &      + p_fld(ji+2,ijje) * zsurfmsk(ji+2,ijje)
               !
               p_fld_crs(ii,2) = zflcrs
            ENDDO

            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  ijje = mje_crs(ij)
                  zflcrs  =  p_fld(ji  ,ijje) * zsurfmsk(ji  ,ijje) &
                     &      + p_fld(ji+1,ijje) * zsurfmsk(ji+1,ijje) &
                     &      + p_fld(ji+2,ijje) * zsurfmsk(ji+2,ijje)
                  !
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         CASE( 'U' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ijie = mie_crs(ii)
                     zflcrs  =  p_fld(ijie,je_2) * zsurfmsk(ijie,je_2)
                     p_fld_crs(ii,2) = zflcrs
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ijie = mie_crs(ii)
                  zflcrs =  p_fld(ijie,je_2  ) * zsurfmsk(ijie,je_2  )  &
                     &     + p_fld(ijie,je_2+1) * zsurfmsk(ijie,je_2+1)  &
                     &     + p_fld(ijie,je_2+2) * zsurfmsk(ijie,je_2+2)

                  p_fld_crs(ii,2) = zflcrs
               ENDDO
            ENDIF

            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  ijie = mie_crs(ii)
                  zflcrs =  p_fld(ijie,jj  ) * zsurfmsk(ijie,jj  )  &
                     &    + p_fld(ijie,jj+1) * zsurfmsk(ijie,jj+1)  &
                     &    + p_fld(ijie,jj+2) * zsurfmsk(ijie,jj+2)
                  !
                  p_fld_crs(ii,ij) = zflcrs
                  !
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

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     zflcrs =  &
                        & MAX( p_fld(ji  ,je_2) * p_mask(ji  ,je_2,1) - ( 1.- p_mask(ji  ,je_2,1) ) * r_inf ,  &
                        &       p_fld(ji+1,je_2) * p_mask(ji+1,je_2,1) - ( 1.- p_mask(ji+1,je_2,1) ) * r_inf ,  &
                        &       p_fld(ji+2,je_2) * p_mask(ji+2,je_2,1) - ( 1.- p_mask(ji+2,je_2,1) ) * r_inf  )
                     !
                     p_fld_crs(ii,2) = zflcrs
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               zflcrs =  &
                  &  MAX( p_fld(ji  ,je_2  ) * p_mask(ji  ,je_2  ,1) - ( 1.- p_mask(ji  ,je_2  ,1) ) * r_inf ,  &
                  &       p_fld(ji+1,je_2  ) * p_mask(ji+1,je_2  ,1) - ( 1.- p_mask(ji+1,je_2  ,1) ) * r_inf ,  &
                  &       p_fld(ji+2,je_2  ) * p_mask(ji+2,je_2  ,1) - ( 1.- p_mask(ji+2,je_2  ,1) ) * r_inf ,  &
                  &       p_fld(ji  ,je_2+1) * p_mask(ji  ,je_2+1,1) - ( 1.- p_mask(ji  ,je_2+1,1) ) * r_inf ,  &
                  &       p_fld(ji+1,je_2+1) * p_mask(ji+1,je_2+1,1) - ( 1.- p_mask(ji+1,je_2+1,1) ) * r_inf ,  &
                  &       p_fld(ji+2,je_2+1) * p_mask(ji+2,je_2+1,1) - ( 1.- p_mask(ji+2,je_2+1,1) ) * r_inf ,  &
                  &       p_fld(ji  ,je_2+2) * p_mask(ji  ,je_2+2,1) - ( 1.- p_mask(ji  ,je_2+2,1) ) * r_inf ,  &
                  &       p_fld(ji+1,je_2+2) * p_mask(ji+1,je_2+2,1) - ( 1.- p_mask(ji+1,je_2+2,1) ) * r_inf ,  &
                  &       p_fld(ji+2,je_2+2) * p_mask(ji+2,je_2+2,1) - ( 1.- p_mask(ji+2,je_2+2,1) ) * r_inf   )
               !
               p_fld_crs(ii,2) = zflcrs
            ENDIF

            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  zflcrs = &
                     &  MAX( p_fld(ji  ,jj  ) * p_mask(ji  ,jj  ,1) - ( 1.- p_mask(ji  ,jj  ,1) ) * r_inf ,  &
                     &       p_fld(ji+1,jj  ) * p_mask(ji+1,jj  ,1) - ( 1.- p_mask(ji+1,jj  ,1) ) * r_inf ,  &
                     &       p_fld(ji+2,jj  ) * p_mask(ji+2,jj  ,1) - ( 1.- p_mask(ji+2,jj  ,1) ) * r_inf ,  &
                     &       p_fld(ji  ,jj+1) * p_mask(ji  ,jj+1,1) - ( 1.- p_mask(ji  ,jj+1,1) ) * r_inf ,  &
                     &       p_fld(ji+1,jj+1) * p_mask(ji+1,jj+1,1) - ( 1.- p_mask(ji+1,jj+1,1) ) * r_inf ,  &
                     &       p_fld(ji+2,jj+1) * p_mask(ji+2,jj+1,1) - ( 1.- p_mask(ji+2,jj+1,1) ) * r_inf ,  &
                     &       p_fld(ji  ,jj+2) * p_mask(ji  ,jj+2,1) - ( 1.- p_mask(ji  ,jj+2,1) ) * r_inf ,  &
                     &       p_fld(ji+1,jj+2) * p_mask(ji+1,jj+2,1) - ( 1.- p_mask(ji+1,jj+2,1) ) * r_inf ,  &
                     &       p_fld(ji+2,jj+2) * p_mask(ji+2,jj+2,1) - ( 1.- p_mask(ji+2,jj+2,1) ) * r_inf   )
                  !
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         CASE( 'V' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  ijje = mje_crs(2)
               ENDIF
            ELSE
               ijje = mjs_crs(2)
            ENDIF

            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               zflcrs = MAX( p_fld(ji  ,ijje) * p_mask(ji  ,ijje,1) - ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                  &           p_fld(ji+1,ijje) * p_mask(ji+1,ijje,1) - ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                  &           p_fld(ji+2,ijje) * p_mask(ji+2,ijje,1) - ( 1.- p_mask(ji,ijje,1) ) * r_inf )
               !
               p_fld_crs(ii,2) = zflcrs
            ENDDO
            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  ijje = mje_crs(ij)
                  !
                  zflcrs = MAX( p_fld(ji  ,ijje) * p_mask(ji  ,ijje,1) - ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                     &           p_fld(ji+1,ijje) * p_mask(ji+1,ijje,1) - ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                     &           p_fld(ji+2,ijje) * p_mask(ji+2,ijje,1) - ( 1.- p_mask(ji,ijje,1) ) * r_inf )
                  !
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         CASE( 'U' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ijie = mie_crs(ii)
                     zflcrs  =  p_fld(ijie,je_2) * p_mask(ijie,je_2,1) - ( 1.- p_mask(ijie,je_2,1) ) * r_inf
                     p_fld_crs(ii,2) = zflcrs
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ijie = mie_crs(ii)
                  zflcrs =  &
                     &  MAX( p_fld(ijie,je_2  ) * p_mask(ijie,je_2  ,1) - ( 1.- p_mask(ijie,je_2,1) ) * r_inf ,  &
                     &       p_fld(ijie,je_2+1) * p_mask(ijie,je_2+1,1) - ( 1.- p_mask(ijie,je_2,1) ) * r_inf ,  &
                     &       p_fld(ijie,je_2+2) * p_mask(ijie,je_2+2,1) - ( 1.- p_mask(ijie,je_2,1) ) * r_inf  )
                  p_fld_crs(ii,2) = zflcrs
               ENDDO
            ENDIF
            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  ijie = mie_crs(ii)
                  zflcrs =  &
                     &  MAX( p_fld(ijie,jj  ) * p_mask(ijie,jj  ,1) - ( 1.- p_mask(ijie,jj,1) ) * r_inf ,  &
                     &       p_fld(ijie,jj+1) * p_mask(ijie,jj+1,1) - ( 1.- p_mask(ijie,jj,1) ) * r_inf ,  &
                     &      p_fld(ijie,jj+2) * p_mask(ijie,jj+2,1) - ( 1.- p_mask(ijie,jj,1) ) * r_inf  )
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         END SELECT

      CASE ( 'MIN' )      !   Search the min of unmasked grid cells

         SELECT CASE ( cd_type )

         CASE( 'T', 'W' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     zflcrs =  &
                        & MIN( p_fld(ji  ,je_2) * p_mask(ji  ,je_2,1) + ( 1.- p_mask(ji  ,je_2,1) ) * r_inf ,  &
                        &       p_fld(ji+1,je_2) * p_mask(ji+1,je_2,1) + ( 1.- p_mask(ji+1,je_2,1) ) * r_inf ,  &
                        &       p_fld(ji+2,je_2) * p_mask(ji+2,je_2,1) + ( 1.- p_mask(ji+2,je_2,1) ) * r_inf  )
                     !
                     p_fld_crs(ii,2) = zflcrs
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               zflcrs =  &
                  &  MIN( p_fld(ji  ,je_2  ) * p_mask(ji  ,je_2  ,1) + ( 1.- p_mask(ji  ,je_2  ,1) ) * r_inf ,  &
                  &       p_fld(ji+1,je_2  ) * p_mask(ji+1,je_2  ,1) + ( 1.- p_mask(ji+1,je_2  ,1) ) * r_inf ,  &
                  &       p_fld(ji+2,je_2  ) * p_mask(ji+2,je_2  ,1) + ( 1.- p_mask(ji+2,je_2  ,1) ) * r_inf ,  &
                  &       p_fld(ji  ,je_2+1) * p_mask(ji  ,je_2+1,1) + ( 1.- p_mask(ji  ,je_2+1,1) ) * r_inf ,  &
                  &       p_fld(ji+1,je_2+1) * p_mask(ji+1,je_2+1,1) + ( 1.- p_mask(ji+1,je_2+1,1) ) * r_inf ,  &
                  &       p_fld(ji+2,je_2+1) * p_mask(ji+2,je_2+1,1) + ( 1.- p_mask(ji+2,je_2+1,1) ) * r_inf ,  &
                  &       p_fld(ji  ,je_2+2) * p_mask(ji  ,je_2+2,1) + ( 1.- p_mask(ji  ,je_2+2,1) ) * r_inf ,  &
                  &       p_fld(ji+1,je_2+2) * p_mask(ji+1,je_2+2,1) + ( 1.- p_mask(ji+1,je_2+2,1) ) * r_inf ,  &
                  &       p_fld(ji+2,je_2+2) * p_mask(ji+2,je_2+2,1) + ( 1.- p_mask(ji+2,je_2+2,1) ) * r_inf   )
               !
               p_fld_crs(ii,2) = zflcrs
            ENDIF

            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  zflcrs = &
                     &  MIN( p_fld(ji  ,jj  ) * p_mask(ji  ,jj  ,1) + ( 1.- p_mask(ji  ,jj  ,1) ) * r_inf ,  &
                     &       p_fld(ji+1,jj  ) * p_mask(ji+1,jj  ,1) + ( 1.- p_mask(ji+1,jj  ,1) ) * r_inf ,  &
                     &       p_fld(ji+2,jj  ) * p_mask(ji+2,jj  ,1) + ( 1.- p_mask(ji+2,jj  ,1) ) * r_inf ,  &
                     &       p_fld(ji  ,jj+1) * p_mask(ji  ,jj+1,1) + ( 1.- p_mask(ji  ,jj+1,1) ) * r_inf ,  &
                     &       p_fld(ji+1,jj+1) * p_mask(ji+1,jj+1,1) + ( 1.- p_mask(ji+1,jj+1,1) ) * r_inf ,  &
                     &       p_fld(ji+2,jj+1) * p_mask(ji+2,jj+1,1) + ( 1.- p_mask(ji+2,jj+1,1) ) * r_inf ,  &
                     &       p_fld(ji  ,jj+2) * p_mask(ji  ,jj+2,1) + ( 1.- p_mask(ji  ,jj+2,1) ) * r_inf ,  &
                     &       p_fld(ji+1,jj+2) * p_mask(ji+1,jj+2,1) + ( 1.- p_mask(ji+1,jj+2,1) ) * r_inf ,  &
                     &       p_fld(ji+2,jj+2) * p_mask(ji+2,jj+2,1) + ( 1.- p_mask(ji+2,jj+2,1) ) * r_inf   )
                  !
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         CASE( 'V' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  ijje = mje_crs(2)
               ENDIF
            ELSE
               ijje = mjs_crs(2)
            ENDIF

            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               zflcrs = MIN( p_fld(ji  ,ijje) * p_mask(ji  ,ijje,1) + ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                  &           p_fld(ji+1,ijje) * p_mask(ji+1,ijje,1) + ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                  &           p_fld(ji+2,ijje) * p_mask(ji+2,ijje,1) + ( 1.- p_mask(ji,ijje,1) ) * r_inf )
               !
               p_fld_crs(ii,2) = zflcrs
            ENDDO
            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  ijje = mje_crs(ij)
                  !
                  zflcrs = MIN( p_fld(ji  ,ijje) * p_mask(ji  ,ijje,1) + ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                     &           p_fld(ji+1,ijje) * p_mask(ji+1,ijje,1) + ( 1.- p_mask(ji,ijje,1) ) * r_inf ,  &
                     &           p_fld(ji+2,ijje) * p_mask(ji+2,ijje,1) + ( 1.- p_mask(ji,ijje,1) ) * r_inf )
                  !
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         CASE( 'U' )

            IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
               IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
                  je_2 = mje_crs(2)
                  DO ji = nistr, niend, nn_factx
                     ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                     ijie = mie_crs(ii)
                     zflcrs  =  p_fld(ijie,je_2) * p_mask(ijie,je_2,1) + ( 1.- p_mask(ijie,je_2,1) ) * r_inf

                     p_fld_crs(ii,2) = zflcrs
                  ENDDO
               ENDIF
            ELSE
               je_2 = mjs_crs(2)
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ijie = mie_crs(ii)
                  zflcrs =  &
                     &  MIN( p_fld(ijie,je_2  ) * p_mask(ijie,je_2  ,1) + ( 1.- p_mask(ijie,je_2,1) ) * r_inf ,  &
                     &       p_fld(ijie,je_2+1) * p_mask(ijie,je_2+1,1) + ( 1.- p_mask(ijie,je_2,1) ) * r_inf ,  &
                     &       p_fld(ijie,je_2+2) * p_mask(ijie,je_2+2,1) + ( 1.- p_mask(ijie,je_2,1) ) * r_inf  )
                  p_fld_crs(ii,2) = zflcrs
               ENDDO
            ENDIF
            DO jj = njstr, njend, nn_facty
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ij   = ( jj - njstr ) * rfacty_r + 3
                  ijie = mie_crs(ii)
                  zflcrs =  &
                     &  MIN( p_fld(ijie,jj  ) * p_mask(ijie,jj  ,1) + ( 1.- p_mask(ijie,jj,1) ) * r_inf ,  &
                     &       p_fld(ijie,jj+1) * p_mask(ijie,jj+1,1) + ( 1.- p_mask(ijie,jj,1) ) * r_inf ,  &
                     &      p_fld(ijie,jj+2) * p_mask(ijie,jj+2,1) + ( 1.- p_mask(ijie,jj,1) ) * r_inf  )
                  p_fld_crs(ii,ij) = zflcrs
                  !
               ENDDO
            ENDDO

         END SELECT
         !
      END SELECT
      !
      !LOLO: ADD:
      !CALL crs_lbc_lnk( p_fld_crs, cd_type, psgn )
      !
   END SUBROUTINE crs_dom_ope_2d

   SUBROUTINE crs_dom_e3( p_e1, p_e2, p_e3, p_sfc_crs, cd_type, p_mask, p_e3_crs, p_e3_max_crs)
      !!----------------------------------------------------------------
      !!  Arguments
      CHARACTER(len=1),                         INTENT(in) :: cd_type      ! grid type T, W ( U, V, F)
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in) :: p_mask       ! Parent grid T mask
      REAL(wp), DIMENSION(jpi,jpj)    ,         INTENT(in) :: p_e1, p_e2   ! 2D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi,jpj,jpk),         INTENT(in) :: p_e3         ! 3D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(in) :: p_sfc_crs ! Coarse grid box east or north face quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout) :: p_e3_crs ! Coarse grid box east or north face quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(inout) :: p_e3_max_crs ! Coarse grid box east or north face quantity

      !! Local variables
      INTEGER ::  ji, jj, jk                   ! dummy loop indices
      INTEGER ::  ijie, ijje, ii, ij, je_2
      REAL(wp) :: ze3crs
      REAL(wp), DIMENSION(:,:,:), POINTER :: zmask, zsurf

      !!----------------------------------------------------------------

      p_e3_crs    (:,:,:) = 0.
      p_e3_max_crs(:,:,:) = 1.


      ALLOCATE ( zmask(jpi,jpj,jpk), zsurf(jpi,jpj,jpk) )

      SELECT CASE ( cd_type )
      CASE( 'W' )
         zmask(:,:,1) = p_mask(:,:,1)
         DO jk = 2, jpk
            zmask(:,:,jk) = p_mask(:,:,jk-1)
         ENDDO
      CASE DEFAULT
         DO jk = 1, jpk
            zmask(:,:,jk) = p_mask(:,:,jk)
         ENDDO
      END SELECT

      DO jk = 1, jpk
         zsurf(:,:,jk) = p_e1(:,:) * p_e2(:,:) * p_e3(:,:,jk)
      ENDDO

      IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
         IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
            je_2 = mje_crs(2)
            DO jk = 1 , jpk
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  ze3crs =   zsurf(ji  ,je_2,jk) * zmask(ji  ,je_2,jk)   &
                     &   + zsurf(ji+1,je_2,jk) * zmask(ji+1,je_2,jk)   &
                     &   + zsurf(ji+2,je_2,jk) * zmask(ji+2,je_2,jk)

                  p_e3_crs(ii,2,jk) = ze3crs / p_sfc_crs(ii,ij,jk)
                  !
                  ze3crs = MAX( p_e3(ji  ,je_2,jk) * zmask(ji  ,je_2,jk),  &
                     &          p_e3(ji+1,je_2,jk) * zmask(ji+1,je_2,jk),  &
                     &          p_e3(ji+2,je_2,jk) * zmask(ji+2,je_2,jk)  )
                  !
                  p_e3_max_crs(ii,2,jk) = ze3crs
               ENDDO
            ENDDO
         ENDIF
      ELSE
         je_2 = mjs_crs(2)
         DO jk = 1 , jpk
            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               ze3crs =  zsurf(ji  ,je_2  ,jk) * zmask(ji  ,je_2  ,jk)   &
                  &    + zsurf(ji+1,je_2  ,jk) * zmask(ji+1,je_2  ,jk)   &
                  &    + zsurf(ji+2,je_2  ,jk) * zmask(ji+2,je_2  ,jk)   &
                  &    + zsurf(ji  ,je_2+1,jk) * zmask(ji  ,je_2+1,jk)   &
                  &    + zsurf(ji+1,je_2+1,jk) * zmask(ji+1,je_2+1,jk)   &
                  &    + zsurf(ji+2,je_2+1,jk) * zmask(ji+2,je_2+1,jk)   &
                  &    + zsurf(ji  ,je_2+2,jk) * zmask(ji  ,je_2+2,jk)   &
                  &    + zsurf(ji+1,je_2+2,jk) * zmask(ji+1,je_2+2,jk)   &
                  &    + zsurf(ji+2,je_2+2,jk) * zmask(ji+2,je_2+2,jk)

               p_e3_crs(ii,2,jk) = ze3crs / p_sfc_crs(ii,2,jk)
               !
               ze3crs = MAX( p_e3(ji  ,je_2  ,jk) * zmask(ji  ,je_2  ,jk),  &
                  &          p_e3(ji+1,je_2  ,jk) * zmask(ji+1,je_2  ,jk),  &
                  &          p_e3(ji+2,je_2  ,jk) * zmask(ji+2,je_2  ,jk),  &
                  &          p_e3(ji  ,je_2+1,jk) * zmask(ji  ,je_2+1,jk),  &
                  &          p_e3(ji+1,je_2+1,jk) * zmask(ji+1,je_2+1,jk),  &
                  &          p_e3(ji+2,je_2+1,jk) * zmask(ji+2,je_2+1,jk),  &
                  &          p_e3(ji  ,je_2+2,jk) * zmask(ji  ,je_2+2,jk),  &
                  &          p_e3(ji+1,je_2+2,jk) * zmask(ji+1,je_2+2,jk),  &
                  &          p_e3(ji+2,je_2+2,jk) * zmask(ji+2,je_2+2,jk) )

               p_e3_max_crs(ii,2,jk) = ze3crs
            ENDDO
         ENDDO
      ENDIF
      DO jk = 1 , jpk
         DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               ij   = ( jj - njstr ) * rfacty_r + 3
               ze3crs =   zsurf(ji  ,jj  ,jk) * zmask(ji  ,jj  ,jk)   &
                  &        + zsurf(ji+1,jj  ,jk) * zmask(ji+1,jj  ,jk)   &
                  &        + zsurf(ji+2,jj  ,jk) * zmask(ji+2,jj  ,jk)   &
                  &        + zsurf(ji  ,jj+1,jk) * zmask(ji  ,jj+1,jk)   &
                  &        + zsurf(ji+1,jj+1,jk) * zmask(ji+1,jj+1,jk)   &
                  &        + zsurf(ji+2,jj+1,jk) * zmask(ji+2,jj+1,jk)   &
                  &        + zsurf(ji  ,jj+2,jk) * zmask(ji  ,jj+2,jk)   &
                  &        + zsurf(ji+1,jj+2,jk) * zmask(ji+1,jj+2,jk)   &
                  &        + zsurf(ji+2,jj+2,jk) * zmask(ji+2,jj+2,jk)

               p_e3_crs(ii,ij,jk) = ze3crs / p_sfc_crs(ii,ij,jk)
               !
               ze3crs = MAX( p_e3(ji  ,jj  ,jk) * zmask(ji  ,jj  ,jk),  &
                  &          p_e3(ji+1,jj  ,jk) * zmask(ji+1,jj  ,jk),  &
                  &          p_e3(ji+2,jj  ,jk) * zmask(ji+2,jj  ,jk),  &
                  &          p_e3(ji  ,jj+1,jk) * zmask(ji  ,jj+1,jk),  &
                  &          p_e3(ji+1,jj+1,jk) * zmask(ji+1,jj+1,jk),  &
                  &          p_e3(ji+2,jj+1,jk) * zmask(ji+2,jj+1,jk),  &
                  &          p_e3(ji  ,jj+2,jk) * zmask(ji  ,jj+2,jk),  &
                  &          p_e3(ji+1,jj+2,jk) * zmask(ji+1,jj+2,jk),  &
                  &          p_e3(ji+2,jj+2,jk) * zmask(ji+2,jj+2,jk) )

               p_e3_max_crs(ii,ij,jk) = ze3crs
            ENDDO
         ENDDO
      ENDDO

      !LOLO: ADD:
      !CALL crs_lbc_lnk( p_e3_crs    , cd_type, 1.0, pval=1.0 )
      !CALL crs_lbc_lnk( p_e3_max_crs, cd_type, 1.0, pval=1.0 )
      !
      DEALLOCATE ( zsurf, zmask )
      !
   END SUBROUTINE crs_dom_e3

   SUBROUTINE crs_dom_sfc( p_mask, cd_type, p_surf_crs, p_surf_crs_msk,  p_e1, p_e2, p_e3 )

      !!  Arguments
      CHARACTER(len=1),                         INTENT(in)           :: cd_type      ! grid type T, W ( U, V, F)
      REAL(wp), DIMENSION(jpi,jpj,jpk)        , INTENT(in)           :: p_mask       ! Parent grid T mask
      REAL(wp), DIMENSION(jpi,jpj)            , INTENT(in), OPTIONAL :: p_e1, p_e2         ! 3D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi,jpj,jpk)        , INTENT(in), OPTIONAL :: p_e3         ! 3D tracer T or W on parent grid
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(out)          :: p_surf_crs ! Coarse grid box east or north face quantity
      REAL(wp), DIMENSION(jpi_crs,jpj_crs,jpk), INTENT(out)          :: p_surf_crs_msk ! Coarse grid box east or north face quantity

      !! Local variables
      INTEGER  :: ji, jj, jk                   ! dummy loop indices
      INTEGER  :: ii, ij, je_2
      REAL(wp), DIMENSION(:,:,:), POINTER :: zsurf, zsurfmsk
      !!----------------------------------------------------------------
      ! Initialize


      ALLOCATE ( zsurf(jpi,jpj,jpk), zsurfmsk(jpi,jpj,jpk) )
      !
      SELECT CASE ( cd_type )

      CASE ('W')
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e1(:,:) * p_e2(:,:)
         ENDDO
         zsurfmsk(:,:,1) = zsurf(:,:,1) * p_mask(:,:,1)
         DO jk = 2, jpk
            zsurfmsk(:,:,jk) = zsurf(:,:,jk) * p_mask(:,:,jk-1)
         ENDDO

      CASE ('V')
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e1(:,:) * p_e3(:,:,jk)
         ENDDO
         DO jk = 1, jpk
            zsurfmsk(:,:,jk) = zsurf(:,:,jk) * p_mask(:,:,jk)
         ENDDO

      CASE ('U')
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e2(:,:) * p_e3(:,:,jk)
         ENDDO
         DO jk = 1, jpk
            zsurfmsk(:,:,jk) = zsurf(:,:,jk) * p_mask(:,:,jk)
         ENDDO

      CASE DEFAULT
         DO jk = 1, jpk
            zsurf(:,:,jk) = p_e1(:,:) * p_e2(:,:)
         ENDDO
         DO jk = 1, jpk
            zsurfmsk(:,:,jk) = zsurf(:,:,jk) * p_mask(:,:,jk)
         ENDDO
      END SELECT

      IF( nldj_crs == 1 .AND. ( ( mje_crs(2) - mjs_crs(2) ) < 2 ) ) THEN     !!cc bande du sud style ORCA2
         IF( mje_crs(2) - mjs_crs(2) == 1 ) THEN
            je_2 = mje_crs(2)
            DO jk = 1, jpk
               DO ji = nistr, niend, nn_factx
                  ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
                  !
                  p_surf_crs    (ii,2,jk) =  zsurf(ji,je_2  ,jk) + zsurf(ji+1,je_2  ,jk) + zsurf(ji+2,je_2  ,jk) &
                     &                      + zsurf(ji,je_2-1,jk) + zsurf(ji+1,je_2-1,jk) + zsurf(ji+2,je_2-1,jk)  ! Why ?????
                  !
                  p_surf_crs_msk(ii,2,jk) =  zsurfmsk(ji,je_2,jk) + zsurfmsk(ji+1,je_2,jk) + zsurfmsk(ji+2,je_2,jk)
                  !
               ENDDO
            ENDDO
         ENDIF
      ELSE
         je_2 = mjs_crs(2)
         DO jk = 1, jpk
            DO ji = nistr, niend, nn_factx
               ii   = ( ji - mis_crs(2) ) * rfactx_r + 2
               !
               p_surf_crs    (ii,2,jk) =  zsurf(ji,je_2  ,jk) + zsurf(ji+1,je_2  ,jk) + zsurf(ji+2,je_2  ,jk)  &
                  &                   + zsurf(ji,je_2+1,jk) + zsurf(ji+1,je_2+1,jk) + zsurf(ji+2,je_2+1,jk)  &
                  &                   + zsurf(ji,je_2+2,jk) + zsurf(ji+1,je_2+2,jk) + zsurf(ji+2,je_2+2,jk)

               p_surf_crs_msk(ii,2,jk) =  zsurfmsk(ji,je_2  ,jk) + zsurfmsk(ji+1,je_2  ,jk) + zsurfmsk(ji+2,je_2  ,jk)  &
                  &                   + zsurfmsk(ji,je_2+1,jk) + zsurfmsk(ji+1,je_2+1,jk) + zsurfmsk(ji+2,je_2+1,jk)  &
                  &                   + zsurfmsk(ji,je_2+2,jk) + zsurfmsk(ji+1,je_2+2,jk) + zsurfmsk(ji+2,je_2+2,jk)
            ENDDO
         ENDDO
      ENDIF

      DO jk = 1, jpk
         DO jj = njstr, njend, nn_facty
            DO ji = nistr, niend, nn_factx
               ii = ( ji - mis_crs(2) ) * rfactx_r + 2
               ij = ( jj - njstr ) * rfacty_r + 3
               !
               p_surf_crs    (ii,ij,jk) =  zsurf(ji,jj  ,jk) + zsurf(ji+1,jj  ,jk) + zsurf(ji+2,jj  ,jk)  &
                  &                    + zsurf(ji,jj+1,jk) + zsurf(ji+1,jj+1,jk) + zsurf(ji+2,jj+1,jk)  &
                  &                    + zsurf(ji,jj+2,jk) + zsurf(ji+1,jj+2,jk) + zsurf(ji+2,jj+2,jk)

               p_surf_crs_msk(ii,ij,jk) =  zsurfmsk(ji,jj  ,jk) + zsurfmsk(ji+1,jj  ,jk) + zsurfmsk(ji+2,jj  ,jk)  &
                  &                    + zsurfmsk(ji,jj+1,jk) + zsurfmsk(ji+1,jj+1,jk) + zsurfmsk(ji+2,jj+1,jk)  &
                  &                    + zsurfmsk(ji,jj+2,jk) + zsurfmsk(ji+1,jj+2,jk) + zsurfmsk(ji+2,jj+2,jk)
            ENDDO
         ENDDO
      ENDDO

      !LOLO: ADD:
      !CALL crs_lbc_lnk( p_surf_crs    , cd_type, 1.0, pval=1.0 )
      !CALL crs_lbc_lnk( p_surf_crs_msk, cd_type, 1.0, pval=1.0 )

      DEALLOCATE ( zsurfmsk, zsurf )

   END SUBROUTINE crs_dom_sfc

   SUBROUTINE crs_dom_def
      !!----------------------------------------------------------------
      !!               *** SUBROUTINE crs_dom_def ***
      !! ** Purpose :  Three applications.
      !!               1) Define global domain indice of the croasening grid
      !!               2) Define local domain indice of the croasening grid
      !!               3) Define the processor domain indice for a croasening grid
      !!----------------------------------------------------------------
      !!
      !!  local variables

      INTEGER  :: ji,jj,jk,ijjgloT,ijis,ijie,ijjs,ijje,jn      ! dummy indices
      INTEGER  :: ierr                                ! allocation error status


      ! 1.a. Define global domain indices  : take into account the interior domain only ( removes i/j=1 , i/j=jpiglo/jpjglo ) then add 2/3 grid points
      jpiglo_crs   = INT( (jpiglo - 2) / nn_factx ) + 2
      !    jpjglo_crs   = INT( (jpjglo - 2) / nn_facty ) + 2  ! the -2 removes j=1, j=jpj
      !    jpjglo_crs   = INT( (jpjglo - 2) / nn_facty ) + 3
      jpjglo_crs   = INT( (jpjglo - MOD(jpjglo, nn_facty)) / nn_facty ) + 3
      jpiglo_crsm1 = jpiglo_crs - 1
      jpjglo_crsm1 = jpjglo_crs - 1

      jpi_crs = ( jpiglo_crs   - 2 * jpreci + (jpni-1) ) / jpni + 2 * jpreci
      jpj_crs = ( jpjglo_crsm1 - 2 * jprecj + (jpnj-1) ) / jpnj + 2 * jprecj

      IF( noso < 0 ) jpj_crs = jpj_crs + 1    ! add a local band on southern processors

      jpi_crsm1   = jpi_crs - 1
      jpj_crsm1   = jpj_crs - 1
      nperio_crs  = jperio
      npolj_crs   = npolj

      ierr = crs_dom_alloc()          ! allocate most coarse grid arrays

      !LOLO: skip:
      ! 2.a Define processor domain
      nimpp_crs  = 1
      njmpp_crs  = 1
      nlci_crs   = jpi_crs
      nlcj_crs   = jpj_crs
      nldi_crs   = 1
      nldj_crs   = 1
      nlei_crs   = jpi_crs
      nlej_crs   = jpj_crs



      !                         Save the parent grid information
      !LOLO: not needed...?:
      !jpi_full    = jpi
      !jpj_full    = jpj
      !jpim1_full  = jpim1
      !jpjm1_full  = jpjm1
      !nperio_full = nperio

      !npolj_full  = npolj
      !jpiglo_full = jpiglo
      !jpjglo_full = jpjglo

      !nlcj_full   = nlcj
      !nlci_full   = nlci
      !nldi_full   = nldi
      !nldj_full   = nldj
      !nlei_full   = nlei
      !nlej_full   = nlej
      !nimpp_full  = nimpp
      !njmpp_full  = njmpp

      !nlcit_full(:)  = nlcit(:)
      !nldit_full(:)  = nldit(:)
      !nleit_full(:)  = nleit(:)
      !nimppt_full(:) = nimppt(:)
      !nlcjt_full(:)  = nlcjt(:)
      !nldjt_full(:)  = nldjt(:)
      !nlejt_full(:)  = nlejt(:)
      !njmppt_full(:) = njmppt(:)
      !LOLO.
      
      CALL dom_grid_crs  !swich de grille


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

      CALL dom_grid_glo

      mxbinctr   = INT( nn_factx * 0.5 )
      mybinctr   = INT( nn_facty * 0.5 )

      nrestx = MOD( nn_factx, 2 )   ! check if even- or odd- numbered reduction factor
      nresty = MOD( nn_facty, 2 )

      IF ( nrestx == 0 ) THEN
         mxbinctr = mxbinctr - 1
      ENDIF

      IF ( nresty == 0 ) THEN
         mybinctr = mybinctr - 1
         IF ( jperio == 3 .OR. jperio == 4 )  nperio_crs = jperio + 2
         IF ( jperio == 5 .OR. jperio == 6 )  nperio_crs = jperio - 2

         IF ( npolj == 3 ) npolj_crs = 5
         IF ( npolj == 5 ) npolj_crs = 3
      ENDIF

      rfactxy = nn_factx * nn_facty

      ! 2.b. Set up bins for coarse grid, horizontal only.
      ierr = crs_dom_alloc2()

      mis2_crs(:) = 0      ;      mie2_crs(:) = 0
      mjs2_crs(:) = 0      ;      mje2_crs(:) = 0


      SELECT CASE ( nn_binref )

      CASE ( 0 )

         SELECT CASE ( nperio )


         CASE ( 0, 1, 3, 4 )    !   3, 4 : T-Pivot at North Fold

            DO ji = 2, jpiglo_crsm1
               ijie = ( ji * nn_factx ) - nn_factx   !cc
               ijis = ijie - nn_factx + 1
               mis2_crs(ji) = ijis
               mie2_crs(ji) = ijie
            ENDDO
            IF ( jpiglo - 1 - mie2_crs(jpiglo_crsm1) <= nn_factx ) mie2_crs(jpiglo_crsm1) = jpiglo - 2

            ! Handle first the northernmost bin
            IF ( nn_facty == 2 ) THEN   ;    ijjgloT = jpjglo - 1
            ELSE                        ;    ijjgloT = jpjglo
            ENDIF

            DO jj = 2, jpjglo_crs
               ijje = ijjgloT - nn_facty * ( jj - 3 )
               ijjs = ijje - nn_facty + 1
               mjs2_crs(jpjglo_crs-jj+2) = ijjs
               mje2_crs(jpjglo_crs-jj+2) = ijje
            ENDDO

         CASE ( 2 )
            WRITE(numout,*)  'crs_init, jperio=2 not supported'

         CASE ( 5, 6 )    ! F-pivot at North Fold

            DO ji = 2, jpiglo_crsm1
               ijie = ( ji * nn_factx ) - nn_factx
               ijis = ijie - nn_factx + 1
               mis2_crs(ji) = ijis
               mie2_crs(ji) = ijie
            ENDDO
            IF ( jpiglo - 1 - mie2_crs(jpiglo_crsm1) <= nn_factx ) mie_crs(jpiglo_crsm1)  = jpiglo - 2

            ! Treat the northernmost bin separately.
            jj = 2
            ijje = jpj - nn_facty * ( jj - 2 )
            IF ( nn_facty == 3 ) THEN   ;  ijjs = ijje - 1
            ELSE                        ;  ijjs = ijje - nn_facty + 1
            ENDIF
            mjs2_crs(jpj_crs-jj+1) = ijjs
            mje2_crs(jpj_crs-jj+1) = ijje

            ! Now bin the rest, any remainder at the south is lumped in the southern bin
            DO jj = 3, jpjglo_crsm1
               ijje = jpjglo - nn_facty * ( jj - 2 )
               ijjs = ijje - nn_facty + 1
               IF ( ijjs <= nn_facty )  ijjs = 2
               mjs2_crs(jpj_crs-jj+1)   = ijjs
               mje2_crs(jpj_crs-jj+1)   = ijje
            ENDDO

         CASE DEFAULT
            WRITE(numout,*) 'crs_init. Only jperio = 0, 1, 3, 4, 5, 6 supported'

         END SELECT

      CASE (1 )
         WRITE(numout,*) 'crs_init.  Equator-centered bins option not yet available'

      END SELECT

      ! Pad the boundaries, do not know if it is necessary
      mis2_crs(2) = 1             ;  mis2_crs(jpiglo_crs) = mie2_crs(jpiglo_crs - 1) + 1
      mie2_crs(2) = nn_factx      ;  mie2_crs(jpiglo_crs) = jpiglo
      !
      mjs2_crs(1) = 1
      mje2_crs(1) = 1
      !
      mje2_crs(2) = mjs2_crs(3)-1 ;  mje2_crs(jpjglo_crs) = jpjglo
      mjs2_crs(2) = 1             ;  mjs2_crs(jpjglo_crs) = mje2_crs(jpjglo_crs) - nn_facty + 1


      mis_crs(:) = mis2_crs(:)
      mie_crs(:) = mie2_crs(:)
      mjs_crs(:) = mjs2_crs(:)
      mje_crs(:) = mje2_crs(:)
      !
      nistr = mis_crs(2)  ;   niend = mis_crs(nlci_crs - 1)
      njstr = mjs_crs(3)  ;   njend = mjs_crs(nlcj_crs - 1)
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
      REAL(wp), DIMENSION(:,:)  , POINTER :: zmbk
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

      zmbk(:,:) = 0.0
      zmbk(:,:) = REAL( mbathy_crs(:,:), wp )
      !LOLO: ADD:
      !CALL crs_lbc_lnk(zmbk,'T',1.0)
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
      zmbk(:,:) = REAL( mbku_crs(:,:), wp )
      !LOLO: ADD:
      !CALL crs_lbc_lnk(zmbk,'U',1.0)
      mbku_crs  (:,:) = MAX( INT( zmbk(:,:) ), 1 )
      zmbk(:,:) = REAL( mbkv_crs(:,:), wp )
      !LOLO: ADD:
      !CALL crs_lbc_lnk(zmbk,'V',1.0)
      mbkv_crs  (:,:) = MAX( INT( zmbk(:,:) ), 1 )
      !
      DEALLOCATE ( zmbk )
      !
   END SUBROUTINE crs_dom_bat


END MODULE mod_crs
