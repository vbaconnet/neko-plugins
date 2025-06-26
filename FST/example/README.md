Simple channel example case with FST. At the moment this is just a box with periodicity in the z-direction,
no slip on the y faces and inflow + outflow in the x direction.

The inflow condition is "user_velocity" which is driven by the user_dirichlet function
in the user file. 

The outflow condition is "outflow+user" which allows you to set a specified pressure profile.

The two boundary conditions above are taken from the initial condition field.

This case should run easily on minimum 4 ranks.

Compile using:

```
makeneko ../0* user.f90
``` 

# For the swept wing case

**Remember that the number of GLL points = polynomial order + 1**. For example, polynomial order of 7 means `lx1 = 7+1 = 8`.

Checklist (without FST):
- Convert your `re2` mesh to `nmsh` and adjust the boundary labels in `run.case` based on the output of the converter.
- Set your `end_time` and all other time parameters in `run.case`
- Set your tolerances etc in `run.case`
- You shouldn't need to do anything with the user file `user.f90` but you can have a look if you want to.

Checklist (with FST):
- Set the FST parameters in the file `../01_global_params.f90`
- Set some more FST parameters in `run.case` (alpha, start time, linear ramp_up time length) 
