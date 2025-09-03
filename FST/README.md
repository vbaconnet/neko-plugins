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
\lambda_y = S\left( \frac{y - y_{start}}{\delta_{y,rise}}\right) - S\left( \frac{y - y_{end}}{\delta_{y,fall}}\right) + 1
$$

and the same for z. Note that $\lambda_y(y) = 1$ if `"periodic_y": true`, and the same applies for z.

The parameters that need to be set are:
- `y_start` and `y_end` if the `y` direction is not periodic
- `z_start` and `z_end` if the `z` direction is not periodic

Note that if any direction is set to periodic, the corresponding parameters _start and _end will not be used.

### Time parameters

The parameter `alpha` is used to specify the delta_rise and delta_fall parameters for the firnge function in space. Basically, instead of specifying the rise and fall
lengths for the fringe function, we give a percentage of the total boundary length. This is so that you can avoid applying FST close to the boundaries that might make 
the simulation blow up. If you don't want to do this you can set `alpha` to a very small value or tweak a few things in the `07_fst_bc_driver` module).

You must specify parameters `periodic_y` and `periodic_z` to `true` if your simulation is periodic in any (or both) of those directions. 
