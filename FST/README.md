# Free-Stream turbulence generation

This is a neko implementation of the FST code applied as a boundary condition. The FST is applied with a fringe function in space and in time. 

# Usage

## FST parameters

Set your turbulent length scale, turbulence intensity, free-stream velocity and total wavenumber discretization parameters in `01_global_params.f90`. 

## User file

The initialization, generation, application of the FST is driven by `07_fst_bc_driver.f90`. This file implements three functions that must be used in the user file:
- `fst_bc_driver_initialize()`
- `fst_bc_driver_apply()`
- `fst_bc_driver_finalize()`

An example of usage is given in the user file `example.f90`. Note that to apply the boundary condition we use the `field_dirichlet_update` function which
requires the use of the `user_velocity` boundary condition on the desired boundary (see `example.case`).

## Case file

The driver module uses some parameters that should be given in the case file. Below is the JSON object taken from `example.case` that shows which parameters to use:

```.json
"FST": {
      "enabled": true,   // default is true
      "t_start": 0.0001, // Time at which to start applying FST
      "t_ramp": 0.001,   // Length of the linear ramp in time
      "alpha": 0.2,      // see below for full explanation of what this is
      "ystart": -0.01,  // Lower bound for the fringe function
      "yend": 0.01,     // High bound for the fringe function
      "periodic_z": true // Self-explanatory. If periodic in y add "periodic_y": true
}
```

### Spatial fringe parameters

A smooth fringe function is applied on the 2D inlet plane, which at the moment is assumed to be `(y,z)`.
The shape of this fringe is the one used in SIMSON and by lots of other people:

$$
\lambda (y,z) = \lambda_y(y)\lambda_z(z),
$$

With

$$
\lambda_u = S\left( \frac{u - u_{start}}{\delta_{u,rise}}\right) - S\left( \frac{u - u_{end}}{\delta_{u,fall}}+1\right)
$$

Note that $\lambda_u = 1$ if the direction `u` is set to be periodic.

`_start` and `_end` parameters need to be set by the user, which represent geometrical coordinates. 
By default, and only if the direction is not periodic, `_start` will be set to the minimum coordinate on the boundary (in that direction). 
The same goes for `_end`, it will be by default set to the maximum value.
The quantities $\delta_{u,*}$ are computed as a percentage $\alpha$ of the total boundary length 
in the direction `u`: $\delta_{u,rise} = \delta_{u,fall} = \alpha * L_u$, where 
$L_u$ is the total domain length at the inlet in the direction u. 

The parameters that need to be set are:
- `y_start` and `y_end` if the `y` direction is not periodic
- `z_start` and `z_end` if the `z` direction is not periodic
- `alpha`, which takes a number between 0 and 1 (there is no check if > 1 or < 0).

Note that if any direction is set to periodic, the corresponding parameters _start and _end will not be used.

### Time parameters

It is possible to set a linear ramp in time to gradually apply the FST. Use the following parameters:
- `t_start`, the time at which to start applying the FST 
- `t_ramp`, the length of the linear ramp, where the amplitude of the FST grows linear up to 1 at `t = t_start + t_ramp`
