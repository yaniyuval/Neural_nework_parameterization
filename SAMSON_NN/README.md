# SAMSON NN

## Description

This directory contains the source code from SAMSON (an adaptation of the [SAM](http://rossby.msrc.sunysb.edu/~marat/SAM.html) atmospheric model) with changes made to incorporate a Neural Net parameterisation of convection.

## Building

The following dependencies are required:

- NetCDF and NetCDF-Fortran
- MPI

To build create a symbolic link to your desired `Build` script and `Makefile`.
These may need modifications to run on your system. Guidance can be found in the
[online SAM documentation](http://rossby.msrc.sunysb.edu/~marat/SAM.html), though be
aware that SAMSON is based on an older version of SAM.  
Ensure that you have permissions to run the [`Build`](Build) script and then execute it
with :
```
./Build
```

This should produce `DATA3D`, `OBJ`, and `RESTART` files in your desired location, and
an executable `SAM_RAD_<RADSCHEME>`.

## Running

To run the executable requires setting up a case.
The name of the case in stored in the [`Casename`](Casename) file, with associated
input files in a folder with the case name.  
At minimum a case requires:

- `prm` - namelist for model parameters
- `snd` - initial soundings

It may also contain:

- `grd`
- `lsf`
- `sfc`
- `rad`

as described in the SAM documentation.

The executable can be run using `mpirun`. E,g, to run on 8 processors:
```
mpirun -np 8 SAM_RAD_<RADSCHEME>
```
Results will appear in the locations set up by the Build script and Makefile, linked
from this directory.

### Selecting the NN implementation

There is the option to run using Janni's original implementation of the NN parametersation
or our re-implementation.
For Janni's implementation is contained in the file:

- `nn_convection_flux.f90`

whilst our implementation is contained in:

- `nn_cf_net.f90`
- `nn_convection_flux.f90`
- `nn_interface_SAM.f90`

Default behaviour is to run with Janni's implementation, but this can be changed by
setting `doiccsnn = .true.` in the Case `prm` file.

