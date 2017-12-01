
This example does the interpolation of a surface 10m wind vector field from a
regular grid (ECMWF) to the ORCA2 tri-polar NEMO grid, including the rotation of
the vector to account for the distortion of the ORCA grid in the North.

Two options included in the test (run "./test_corr_vect.sh"):

*** u10, v10 regular => ORCA2: u10 on grid_U and v10 on grid_V ==> setup labeled "O2uv"

*** u10, v10 regular => ORCA2: u10 on grid_T and v10 on grid_T ==> setup labeled "O2t"

