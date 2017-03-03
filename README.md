

 ---------------------------------------------------------------------------
  Instalation of SOSIE (Sosie is Only a Surface Interpolation Environement)
 ---------------------------------------------------------------------------

 You need a FORTRAN-90 compiler and the netcdf library (with f90 support) installed.


 Compiling the executables :
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~

 - Create 'make.macro', the Makefile configurable part. This files should be a symbolic
   link or a copy of architecture-related templates found in the 'macro' directory, 
   this file must specify your compiler and options as well as the path to netcdf

- It's important to specify the default INTEGER precision to 4 in your compilation flags
   ifort => -i4  Gfortran => -fdefault-integer-4

 - compile the executables by simply running "make" (gmake)

 - if everything goes well, the 'sosie.x' and 'corr_vect.x'
   executables have been created and are ready to be used

 - 'sosie.x' is the executable used to interpolate 2D fields of a scalar and
   "vectors onto regular grids", it requires a 'namelist' configuration file, 
   the provided template namelist should be documented enough to start now.
        --> check 'template_scalar.namelist' for the interpolation of a scalar field
        --> check 'template_U.namelist' and 'template_V.namelist' for a vector

 - 'corr_vect.x' is used to correct vector components in the case of a distorded
   grid, it needs BOTH components of the vector primarly interpolated with 
   'sosie.x'. By default it only support ORCA grids configuration.
   See the EXAMPLE section for more details.

 You can start to interpolate!
 - Tune the 'namelist' file according to your needs...
 - It's a good idea to use scripts to automatically generate namelists 





Getting started with SOSIE, a few working examples:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The same guidelines are available online at http://sosie.sourceforge.net/

Into the "examples" directory you will find namelist and data to test sosie on working configurations.
For each example you will find a commented and working namelist (namelist.exampleX) explaining the relevant namelist tuning.
We encourage you to use a software like ncview to visualize files to be interpolated and interpolated files.
The common approach to test a given example (# X):


  >> cd examples/
  >> gunzip data/*.gz
  >> sosie.x -f namelist.exampleX


Example #1: basic 2D scalar field interpolation
 -  Interpolation of Reynolds (2002) Long Term Mean SST onto the ORCA1 grid (illustrated on Fig. 3-4). Uncompress files coordinates+tmask_ORCA1.nc.gz and sst.ltm.1971-2000.nc.gz in ./data
 -  Do the interpolation: sosie.x -f namelist.example1
 -  Check the newly created SST_360x180-ORCA1_REYNOLDS_LTM.nc


Example #2: 3D scalar field interpolation
 -  3D interpolation of Levitus (1998) temperature climatology onto the ORCA1 grid (only march). Uncompress files coordinates+tmask_ORCA1.nc.gz and T_levitus_march.nc.gz in ./data
 -  Do the interpolation: sosie.x -f namelist.example2
 -  Check the newly created temp_360x180-ORCA1_march.nc


Example #3: Interpolation from an irregular (ORCA1) to a regular lat-lon 1x1 deg. grid
 -  2D interpolation of a SST snapshot from a random NEMO-ORCA1 simulation onto lat-lon 1x1 deg. grid using the bilinear algorithm (illustrated on Fig. 5-6). Uncompress the file sst_ORCA1_example.nc.gz in ./data
 -  Do the interpolation: sosie.x -f namelist.example3
 -  Check the newly created sst_ORCA1-1x1_test.nc


Example #4: Interpolation and correction of a 2D vector field from a regular lat-lon 1x1 deg. grid to an irregular grid (ORCA1) 
As the ORCA family of grids gets distorded in the northern hemisphere it is necessary to correct (i.e. rotate) both components of the vector. In this example the input vector field is the wind at 10 from a few 6-hourly snapshots of the ERA-INTERIM re-analysis.
Do the "raw" interpolation for the zonal component of the wind: >> sosie.x -f namelist.example4_x
Do the "raw" interpolation for the meridional component of the wind: >> sosie.x -f namelist.example4_y
uraw_1x1-deg-ORCA1_vct.nc and vraw_1x1-deg-ORCA1_vct.nc are created, time to correct:
>> corr_vect.x -x u10 -y v10 -f namelist.example4_x -m data/mesh_mask_ORCA1_light.nc
Check the newly created u10_1x1-deg-ORCA1_vct.nc and u10_1x1-deg-ORCA1_vct.nc

Example #5: Interpolation of a vector field from an irregular (ORCA) grid to a regular grid 
Say that you want to interpolate the current field from an NEMO output to a regular grid.
file ORCA1_my_run_grid_U.nc contains the zonal current vozocrtx
file ORCA1_my_run_grid_V.nc contains the meridional current vomecrty
>> corr_vect.x -I -x vozocrtx -y vomecrty -i ORCA1_my_run_grid_U.nc ORCA1_my_run_grid_V.nc -m mesh_mask_ORCA1_light.nc -t time_counter 
files ORCA1_my_run_grid_U_unrotated.nc and ORCA1_my_run_grid_V_unrotated.nc are generated and contain the "unrotated" version of each component of vector field, these fields can be normally interpolated to a regular grid using the bilinear method... 




 Remarks and bug reports to :            brodeau@gmail.com
