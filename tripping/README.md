# Tripping

Introduces a volume forcing along a line parallel to the y or z direction.

`trip.f90` implements the `trip` module to be used in the user file. See instructions
in the header of the file.

A user file is also provided as an example.

Compile using:

```
makeneko trip.f90 user.f90
```

Make sure to enable the user source term in the case file as well.
