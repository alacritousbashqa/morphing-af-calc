# Morphing Airfoil Estimator
Code meant to provide an initial prediction of the aerodynamic forces on a morphing airfoil. Performs inviscid (panel method) and viscous (Thwaites method) calculations. Wake prediction and turbulent viscous effects are planned to be added.

## Files
- **generateNACA4**: generates the points for a specified 4 digit NACA airfoil with a specified number of panels
- **geo_decomp**: retrieves the midpoints, normals, and tangents of every panel
- **sample_run**: example run file showing the creation of a geometry, necessary environment variables, and running of Thwaites method
- **sample_panel_run**: example run file for just the 2nd order panel method
- **TE_maker**: modifies the trailing edge to end at a specified point
- **Thwaites_panel_1**: performs Thwaites method for laminar boundary layer estimation
- **Thwaites_run**: run file that performs Thwaites method for multiple angles of attack and calculates the force coefficients
- **vortex_2order**: performs inviscid linearly varying vortex panel method
- **getVel**: calculates the induced velocity at a point from the freestream and all the panel strengths
- **getVi**: calculates induced velocity at a point from vortices on a panel
- **wake_run**: run file that calculates the wake using the double wake method
