

![SOSIE Logo](https://brodeau.github.io/sosie/sosie_files/sosie_300.svg)


# Getting started with SOSIE

<!-- (Also check the online SOSIE doc for more info on the interpolation method:
https://brodeau.github.io/sosie/) -->

You need a FORTRAN-90 compiler and the netCDF-4 (HDF5) library with f90 support.

For stability, clone branch 3.0 or the older 2.6 corresponding to releases of the same name...

`>> git clone -b '3.0' ...`



Compiling the executables:
--

 * Create `make.macro`, the configurable part of the `Makefile`. This files should be a symbolic
   link or a copy of architecture-related templates found in the `macro` directory,
   this file must specify your compiler and options as well as the path to netcdf

 * It's important to specify the default INTEGER precision to 4 in your compilation flags (switch `-i4` for *ifort* and `-fdefault-integer-4` for *Gfortran*)

 * compile the executables by simply running `make` (gmake)

 * if everything goes well, the `sosie3.x` and `corr_vect.x`
   executables have been created and are ready to be used

 * `sosie3.x` is the executable used to interpolate 2D fields of a scalar and
   "vectors onto regular grids", it requires a *namelist* configuration file,
   the provided template *namelist* `template.namelist` should be documented enough to start now.

 * `corr_vect.x` is used to correct vector components in the case of a distorded
   grid, it needs BOTH components of the vector primarly interpolated with
   `sosie3.x`. By default it only support ORCA grids configuration.
   See the EXAMPLE section for more details.

 You can start to interpolate!
 * Tune the *namelist* file according to your needs...
 * It`s a good idea to use scripts to automatically generate *namelists*

#### Debian-based Gnu-Linux distro like Ubuntu

Using version 10 (or close) of Gfortran:

    >> sudo apt install gfortran-10 libnetcdf-dev libnetcdff-dev
    >> ln -sf macro/make.macro_gfortran-10_Linux make.macro
    >> make



A few working examples / test cases:
--

First download the two archives containing the NetCDF files for the examples:

For generic examples:
* data.tar.gz      &rarr; https://drive.google.com/uc?id=1dpcHXvZyuNTKD9ijeTiLITyFFpnQRC8i&export=download

For ROMS examples:
* data_roms.tar.gz &rarr; https://drive.google.com/uc?authuser=0&id=1ljE3KZu1fqQlKSumm3CtaR0Uf2I9yps3&export=download

Save and extract them into the `examples` sub-directory of SOSIE.

Into this `examples` sub-directory you will find various *namelists* and scripts
to test SOSIE on working configurations.  For each example you will find a
commented and working *namelist* (usually `namelist.exampleX`) explaining the relevant *namelist* tuning.  We encourage you to
use a software like ncview to visualize fields to be interpolated and fields
that have been interpolated.
The common approach to test a given example (# X):

     >> cd examples/
     >> gunzip data/*.gz
     >> sosie3.x -f namelist.exampleX

&nbsp;

#### Example #1: basic 2D scalar field interpolation

Interpolation of Reynolds (2002) Long Term Mean SST onto the ORCA1 grid
(illustrated on Fig. 3-4). Uncompress files `coordinates+tmask_ORCA1.nc.gz` and
`sst.ltm.1971-2000.nc.gz` into `./data/`:

      >> sosie3.x -f namelist.example1

Check `SST_360x180-ORCA1_REYNOLDS_LTM.nc`

&nbsp;

#### Example #2: 3D scalar field interpolation

3D interpolation of Levitus (1998) temperature climatology onto the ORCA1 grid
(only march). Uncompress files `coordinates+tmask_ORCA1.nc.gz` and
`T_levitus_march.nc.gz` into `./data/`

     >> sosie3.x -f namelist.example2
Check `temp_360x180-ORCA1_march.nc`

&nbsp;

#### Example #3: Interpolation from an irregular (ORCA1) to a regular lat-lon 1x1 deg. grid

2D interpolation of a SST snapshot from a random NEMO-ORCA1 simulation onto
lat-lon 1x1 deg. grid using the bilinear algorithm. Uncompress the file
`sst_ORCA1_example.nc.gz` into `./data/`

    >> sosie3.x -f namelist.example3
Check `sst_ORCA1-1x1_test.nc`

&nbsp;

#### Example #4: Interpolation and correction of a 2D vector field from a regular lat-lon 1x1 deg. grid to an irregular grid (ORCA1)

As the ORCA family of grids gets distorded in the northern hemisphere it is
necessary to correct (i.e. rotate) both components of the vector. In this
example the input vector field is the wind at 10 from a few 6-hourly snapshots
of the ERA-INTERIM re-analysis.

      >> cd examples/vector_correction/
Do the "raw" interpolation for the zonal component of the wind:

      >> sosie3.x -f namelist.example4_O1t_x
Do the "raw" interpolation for the meridional component of the wind:

      >> sosie3.x -f namelist.example4_O1t_y
Now that `uraw_1x1-deg-ORCA1_grid_T.nc4` and `uraw_1x1-deg-ORCA1_grid_T.nc4` are created, it's time to *rotate* them onto the T-grid:

      >> corr_vect.x -G T -f namelist.example4_O1t -m ../data/mesh_mask_ORCA1v2_light.nc4
Check `u10_1x1-deg-ORCA1_grid_T.nc4` and `v10_1x1-deg-ORCA1_grid_T.nc4` (vector components on T-points of the grid).

It is possible to do the same correction onto U,V grid points rather than T points:

      >> corr_vect.x -G U -f namelist.example4_O1t -m ../data/mesh_mask_ORCA1v2_light.nc4
Check `u10_1x1-deg-ORCA1_grid_U.nc4` and `v10_1x1-deg-ORCA1_grid_V.nc4`.



&nbsp;

#### Example #5: 2D regular lat-lon to polar stereographic
Interpolation of high-resolution surface 2-meter air temperature from ECMWF onto a polar stereographic projection of the Arctic.
Do the interpolation:

          >> sosie3.x -f namelist.example5

Check `T2M_2560x480-polar-stereo_Arctic.nc`


&nbsp;

#### Example #6: 3D interpolation from ORCA2 to ORCA1 tri-polar grids
Interpolation of a random 3D+time monthly salinity field on the ORCA2 grid to the ORCA1 grid using the bilinear method.
Do the interpolation:

          >> sosie3.x -f namelist.example6

Check `so_ORCA2-ORCA1_test.nc`

&nbsp;

