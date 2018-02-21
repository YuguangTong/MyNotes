# AstroGK
The source code of AstroGK with both slow mode and Alfven mode antennas is [here](files/agk_160929.tar.bz2). But there are some issues to compile it.
- `antenna.f90` line 419 `use dist_fn_arrays, only: g, vpa` should be commented out. 
- replace `dist_fn.f90` by a [new version](files/dist_fn.f90).
- To compile on TACC/stampede2, add this [makefile](files/Makefile.stampede2) to Makefiles directory. Also make sure to load the correct packages by running:

> export MAKEFLAGS=-IMakefiles
> export GK_SYSTEM='stampede2'

> module purge
> module load intel/17.0.4
> module load impi/17.0.3
> module load netcdf/4.3.3.1
> module load fftw2/2.1.5
> module load hdf5/1.8.16
> module load TACC


