# Free-Stream turbulence generation

This is a neko implementation of the FST code applied as a boundary condition. The FST is applied with a fringe function in space and in time. 

# Usage

## FST parameters

Set your turbulent length scale, turbulence intensity, free-stream velocity and total wavenumber discretization parameters in `01_global_params.f90`. 

## User file

The initialization, generation, application of the FST is driven by `07_fst_bc_driver.f90`. This file implements three functions that must be used in the user file:
- fst_bc_driver_initialize()
- fst_bc_driver_apply()
- fst_bc_driver_finalize()

An example of usage is given in the user file `example.f90`. Note that to apply the boundary condition we use the `field_dirichlet_update` function which
requires the use of the `user_velocity` boundary condition on the desired boundary (see `example.case`).

## Case file

The driver module uses some parameters that should be given in the case file. Below is the JSON object taken from `example.case` that shows which parameters to use:

```.json
"FST": {
      "enabled": true,   // default is true
      "alpha": 0.2,      // see below for full explanation of what this is
      "t_start": 0.0001, // Time at which to start applying FST
      "t_ramp": 0.001,   // Length of the linear ramp in time
      "periodic_z": true 
}
```

The parameter `alpha` is used to specify the delta_rise and delta_fall parameters for the firnge function in space. Basically, instead of specifying the rise and fall
lengths for the fringe function, we give a percentage of the total boundary length. This is so that you can avoid applying FST close to the boundaries that might make 
the simulation blow up. If you don't want to do this you can set `alpha` to a very small value or tweak a few things in the `07_fst_bc_driver` module).

You must specify parameters `periodic_y` and `periodic_z` to `true` if your simulation is periodic in any (or both) of those directions. 
