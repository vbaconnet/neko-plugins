Simple box example case with FST. At the moment this is just a box with periodicity in the z-direction,
normal outflow on the y faces and inflow + outflow in the x direction.

The inflow condition is "user_velocity" which is driven by the user_dirichlet function
in the user file. The boundary conditions are taken from the initial condition field.

This case should run easily on minimum 4 ranks.

Compile using:

```
makeneko ../0* user.f90
``` 
