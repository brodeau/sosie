MODULE MOD_INIT

   USE mod_conf

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GET_ARGUMENTS, READ_NMLST, REMINDER

   !! Declaration of namelist :
   !! -------------------------

   NAMELIST /ndom_src/ csource, ivect, l_reg_src, cf_src, cv_src, cv_t_src, &
      &                cf_x_src, cv_lon_src, cv_lat_src, cf_lsm_src, cv_lsm_src, &
      &                ewper_src, cf_z_src, cv_z_src, ctype_z_src,    &
      &                cf_bathy_src, cv_bathy_src, ssig_src

   NAMELIST /ndom_trg/ ctarget, l_reg_trg, cf_x_trg, cv_lon_trg, cv_lat_trg, cf_lsm_trg,   &
      &                cv_lsm_trg, ewper_trg,cf_z_trg, cv_z_trg, ctype_z_trg, &
      &                cf_bathy_trg, cv_bathy_trg, ssig_trg

   NAMELIST /ninterp/  cmethod, idrown, l_save_drwn, ismooth, jt1, jt2, jplev, vmax, vmin, ismooth_out, ibx_xtra_sm

   NAMELIST /noutput/  cv_out, cu_out, cln_out, cv_t_out, cd_out, cextra, &
      &                lmout, rmiss_val, lct, t0, t_stp, cv_z_out

   PRIVATE usage


CONTAINS


   SUBROUTINE GET_ARGUMENTS()

      INTEGER            :: iargc, jarg
      CHARACTER(len=400) :: cr
      CHARACTER(LEN=2), DIMENSION(3), PARAMETER :: clist_opt = (/ '-h','-p','-f' /)

      PRINT *, ''

      jarg = 0

      DO WHILE ( jarg < iargc() )

         jarg = jarg + 1
         CALL getarg(jarg,cr)

         SELECT CASE (trim(cr))

         CASE('-h')
            call usage()

         CASE('-f')

            IF ( jarg + 1 > iargc() ) THEN
               PRINT *, 'ERROR: Missing namelist name!' ; call usage()
            ELSE

               jarg = jarg + 1
               CALL getarg(jarg,cr)

               IF ( ANY(clist_opt == trim(cr)) ) THEN
                  PRINT *, 'ERROR: "', trim(cr), '" is definitively not the name of the namelist!'
                  call usage()
               ELSE
                  cf_nml_sosie = trim(cr)
               END IF

            END IF

         CASE DEFAULT
            PRINT *, 'Unknown option: ', trim(cr) ; PRINT *, ''
            CALL usage()

         END SELECT

      END DO

      PRINT *, ''
      PRINT *, ' * namelist = ', trim(cf_nml_sosie)
      PRINT *

   END SUBROUTINE GET_ARGUMENTS


   SUBROUTINE READ_NMLST(iex)

      INTEGER, INTENT(in) :: iex !: 1 => sosie, 2 => corr_vect.x

      LOGICAL :: lexist

      INQUIRE(FILE=trim(cf_nml_sosie), EXIST=lexist )
      IF ( .NOT. lexist ) THEN
         WRITE(*,'("ERROR: Namelist ",a," cannot be found!")') '"'//TRIM(cf_nml_sosie)//'"'
         CALL usage()
      END IF

      PRINT *, ''
      PRINT *, 'Opening namelist "', TRIM(cf_nml_sosie), '"'

      OPEN( UNIT=numnam, FILE=trim(cf_nml_sosie), FORM='FORMATTED', STATUS='OLD' )

      !! Reading source section:
      READ( numnam, ndom_src )

      !! Reading target section:
      REWIND( numnam )
      READ( numnam, ndom_trg )


      idrown%np_penetr = 0 ! defaults
      idrown%nt_smooth = 5
      idrown%l_msk_chg = .false.

      !! Reading interpolation section:
      REWIND( numnam )
      READ( numnam, ninterp )

      !! If 3D interpolation:
      IF (jplev == 0) THEN
         IF ( TRIM(cv_z_out) == '' ) THEN
            cv_z_out = TRIM(cv_z_trg) ; ! keep same name
         ELSE
            PRINT *, '*** Target depth vector in target file will be renamed to: ', TRIM(cv_z_out)
         END IF
      END IF

      !! Reading netcdf section:
      REWIND( numnam )
      READ( numnam, noutput )

      CLOSE(numnam)


      l_drown_src = .FALSE.
      IF ( idrown%np_penetr > 0 ) l_drown_src = .TRUE.

      IF ( (.NOT. l_drown_src).AND.(l_save_drwn) ) THEN
         PRINT *, 'WARNING: forcing "l_save_drwn=.FALSE." because DROWN not called (idrown[1]==0)...'
         PRINT *, '         if you want to "drown" input field, set idrown[1]>0 !'
         PRINT *, ''
         l_save_drwn = .FALSE.
      END IF

      IF ( iex == 1 ) THEN

         !! If this is a vector component, then interpolated files won't be correct
         !! until grid distorsion correction, so changing name of output variable :
         IF ( ivect /= 0 ) THEN
            lmout = .FALSE.         ! never mask for 1st stage of interpolation of the vector
            IF  ( ivect == 1 )  THEN
               cv_out = 'uraw'
            ELSEIF ( ivect == 2 )  THEN
               cv_out = 'vraw'
            ELSE
               PRINT *, 'Error: while interpolating a component of a vector field: specify&
                  & ivect=1 or ivect=2 into the namelist!'
            END IF
         END IF

         !! Output file name:
         !! -----------------

         cf_short= TRIM(cv_out)//'_'//TRIM(csource)//'-'//TRIM(ctarget)//   &
            &    '_'//TRIM(cextra)//'.nc'

         cf_out = trim(cd_out)//'/'//trim(cf_short)

         !! Printing user's settings
         !!-------------------------
         !OPEN(UNIT=6, FORM='FORMATTED', RECL=512)
         !!
         WRITE(6,*)''
         IF ( ivect /= 0 ) THEN
            WRITE(6,*)'The current field is a component of a vector!'
         ELSE
            WRITE(6,*)'The current field is NOT a component of a vector!'
         END IF
         WRITE(6,*)''
         !!
         WRITE(6,*)''
         IF ( l_reg_src ) THEN
            WRITE(6,*)'Source grid is declared as regular'
         ELSE
            WRITE(6,*)'Source grid is declared as irregular'
         END IF
         WRITE(6,*)''
         !!
         !!
         WRITE(6,*)''
         IF ( l_reg_trg ) THEN
            WRITE(6,*)'Target grid is declared as regular'
         ELSE
            WRITE(6,*)'Target grid is declared as irregular'
         END IF
         WRITE(6,*)''
         !!

         IF (TRIM(cv_t_src)=='missing_rec') lct = .TRUE.

         WRITE(6,*)''
         WRITE(6,*)'Source file: ', TRIM(cf_src)
         WRITE(6,*)''
         WRITE(6,*)'Method for 2D interoplation: ', cmethod
         WRITE(6,*)''
         IF ( .NOT. lct ) WRITE(6,*)'Time record name in source file: ', TRIM(cv_t_src)
         WRITE(6,*)''
         WRITE(6,*)'File containing the source grid: ', trim(cf_x_src)
         WRITE(6,*)''
         WRITE(6,*)'Longitude name: ', trim(cv_lon_src)
         WRITE(6,*)''
         WRITE(6,*)'Latitude name: ', trim(cv_lat_src)
         WRITE(6,*)''
         WRITE(6,*)'Source grid mask: ', trim(cf_lsm_src)
         WRITE(6,*)''
         WRITE(6,*)'Source mask variable: ', trim(cv_lsm_src)
         WRITE(6,*)''
         WRITE(6,*)'Variable to be treated:', trim(cv_src)
         WRITE(6,*)''
         WRITE(6,*)'Level to treat:', jplev
         WRITE(6,*)''
         WRITE(6,*)'Name of variable on the output file:', trim(cv_out)
         WRITE(6,*)''
         WRITE(6,*)'Units of variable to be treated:', trim(cu_out)
         WRITE(6,*)''
         WRITE(6,*)'Long name of variable to be treated:', trim(cln_out)
         WRITE(6,*)''
         WRITE(6,*)'Target grid: ', trim(cf_x_trg)
         !!
         IF ( (.NOT. l_reg_trg).and.(trim(cf_x_trg) == 'spheric') ) THEN
            PRINT *, 'Problem! If you set "cf_x_trg" to "spheric", then set "l_reg_trg" to ".TRUE."'
            PRINT *, ''
            STOP
         END IF
         !!
         WRITE(6,*)''
         WRITE(6,*)'Longitude on target grid: ', trim(cv_lon_trg)
         WRITE(6,*)''
         WRITE(6,*)'Latitude on target grid: ', trim(cv_lat_trg)
         WRITE(6,*)''
         WRITE(6,*)'Longitude name in target file: ', trim(cv_lon_trg)
         WRITE(6,*)''
         WRITE(6,*)'Latitude name in target file: ', trim(cv_lat_trg)
         WRITE(6,*)''
         WRITE(6,*)'Target land-sea mask file: ', trim(cf_lsm_trg)
         PRINT *, ''
         WRITE(6,*)'Name of land-sea mask variable on target grid: ', trim(cv_lsm_trg)
         WRITE(6,*)''
         WRITE(6,*)'Output file: ', trim(cf_out)
         WRITE(6,*)''
         WRITE(6,*)'Shall we use DROWN on source fields: ', l_drown_src
         WRITE(6,*)''
         IF ( l_drown_src ) THEN
            WRITE(6,*)' => about DROWN action:'
            WRITE(6,*)'    -> continental penetration in # of pixels:', idrown%np_penetr
            WRITE(6,*)'    -> # of time to apply a smoothing on continents after DROWN:', idrown%nt_smooth
            WRITE(6,*)'    -> is the mask on the source field is changing with time :', idrown%l_msk_chg
            WRITE(6,*)''
         END IF
         WRITE(6,*)'East west periodicity of source grid: ', ewper_src
         WRITE(6,*)''
         WRITE(6,*)'Masking output file: ', lmout
         WRITE(6,*)''
         WRITE(6,*)'Value for masked points in output file: ', rmiss_val
         WRITE(6,*)''
         WRITE(6,*)'East west periodicity of target grid: ', ewper_trg
         WRITE(6,*)''
         IF ( ismooth > 0 ) THEN
            WRITE(6,*)'We are going to smooth field '//TRIM(cv_src)//' prior to interpolation!'
            WRITE(6,*)'     => # smoothing itterations:', ismooth
            WRITE(6,*)''
         END IF

         IF (jplev == 0) THEN
            WRITE(6,*)''
            WRITE(6,*)' Going to perform 3D interpolation (jplev = 0): '
            WRITE(6,*)''
            WRITE(6,*)' => file containing source depth vector: ', trim(cf_z_src)
            WRITE(6,*)' => name of source depth vector: ',         trim(cv_z_src)
            WRITE(6,*)' => file containing target depth vector: ', trim(cf_z_trg)
            WRITE(6,*)' => name of target depth vector: ',         trim(cv_z_trg)
            PRINT *,   ' => Target depth vector in output file will be renamed to: ', trim(cv_z_out)
            WRITE(6,*)''
         END IF


         IF (lct) THEN
            WRITE(6,*)'Time is controled ! '
            WRITE(6,*)'Start at: ', t0
            WRITE(6,*)'Time step is: ', t_stp
            WRITE(6,*)''
         ELSE
            WRITE(6,*)'Time is not controlled and will be the same than in source file! '
         END IF

         WRITE(6,*)''
         WRITE(6,*)'                ---------------------'
         WRITE(6,*)''
         WRITE(6,*) ''
         WRITE(6,*) ''
         !CLOSE(6)


         !! Checking for some "non-sense"
         IF ( (cmethod == 'akima').and.(.NOT. l_reg_src) ) THEN
            PRINT *, 'The Akima method only supports regular source grids!'
            PRINT *, '--> If the grid of the source domain is irregular, '
            PRINT *, '    please change "cmethod" from akima to "bilin" into the namelist!'
            STOP
         END IF

         WRITE(cpat,'(a,"-",a)') trim(csource), trim(ctarget)
         PRINT *, '';  PRINT *, 'Starting interpolation for config "',trim(cpat),'"...'
         PRINT *, '';  PRINT *, ''; PRINT *, ''

      END IF

   END SUBROUTINE READ_NMLST



   SUBROUTINE REMINDER()

      IF ( (TRIM(cmethod) == 'no_xy') .AND. .NOT.(l_itrp_3d) ) THEN
         WRITE(6,*) ''
         WRITE(6,*) ' ERROR: it makes no sense to use "no_xy" as the interpolation method'
         WRITE(6,*) '        if it is not to perform a vertical interpolation (3D)!'
         WRITE(6,*) '        If you really meant "no_xy" then:'
         WRITE(6,*) '           - set "jplev = 0" into the namelist '
         WRITE(6,*) '           - fill the 3D-interpolation related namelist parameters in ndom_src and ndom_trg'
         WRITE(6,*) ''
         STOP
      END IF

      IF ( (jplev == 0).AND.(nk_trg == 1) ) THEN
         WRITE(6,*) ''
         WRITE(6,*) ' ERROR: you want a 3D interpolation (jplev=0) but your target domain'
         WRITE(6,*) '        has less than 2 levels!!! /  nk_trg =', nk_trg
         WRITE(6,*) ''
         STOP
      END IF

      WRITE(6,*) '' ;  WRITE(6,*) ''
      WRITE(6,*) '====================================================='
      WRITE(6,*) '                    Current config:'
      WRITE(6,*) '                    ---------------'  !
      WRITE(6,'(" Source domain dimension          : ",i5,"  x",i5,"  x",i4)') &
         &        ni_src , nj_src, nk_src
      WRITE(6,'(" Output domain dimension         : ",i5,"  x",i5,"  x",i4)') &
         &        ni_trg, nj_trg, nk_trg
      WRITE(6,'(" Number of time steps to proceed : ",i5)') Ntr
      IF (l_itrp_3d) THEN
         WRITE(6,*) 'Type of interpolation           :    3D'
      ELSE
         IF (l3d) THEN
            WRITE(6,'(" Type of interpolation           : 2D on level ",i3," of ",i4)') &
               jplev, nk_src
         ELSE
            WRITE(6,*) 'Type of interpolation           :    2D'
         END IF
      END IF
      WRITE(6,*) '====================================================='
      WRITE(6,*) ''; WRITE(6,*) ''
      !!
   END SUBROUTINE REMINDER



   !! Local routine:

   SUBROUTINE usage()
      !!
      PRINT *,''
      PRINT *,'List of command line options:'
      PRINT *, ''
      PRINT *,' -f  <namelist_file>  => Specify which namelist file to use'
      PRINT *, ''
      PRINT *,' -h                   => Show this message'
      PRINT *, ''
      STOP
      !!
   END SUBROUTINE usage

END MODULE MOD_INIT
