# DiskModel

Initial-condition generator for isolated disk galaxies (dark matter halo, stellar disk,
gas disk, bulge, optional hot gaseous halo and central black hole) in Gadget format.

This is a modified version of **MakeDiskGalaxy** by Volker Springel
(Springel & White 1999; Springel, Di Matteo & Hernquist 2005), extended with
redshift scaling, an extended HI disk, a hot gas halo, and alternative halo structure models.

## Layout

```
src/            C source (+ src/nrsrc/: Numerical Recipes routines)
params/         Example input file for MakeDiskGalaxy
scripts/        Slurm job script and IDL plotting scripts
docs/img/       Example plots
```

## Build

Requires a C compiler and `make`; no external libraries beyond libm.

```bash
make              # builds ./MakeDiskGalaxy
make CC=icc       # different compiler
make clean
```

Compile-time options are set via `OPT` in the `Makefile`
(`REDSHIFT_SCALING`, `V_SCALING`, `MY_STRUCTURE`, `OLD_STRUCTURE`,
`ADIABATIC_CONTRACTION`, `PURE_HALO`).

## Run

```bash
./MakeDiskGalaxy params/input_dwarf_01solar_mass.txt
```

or on a cluster:

```bash
sbatch scripts/script_rusty.sh
```

Outputs are written to the current directory, named after `OutputFile` in the parameter file:

| File     | Contents                                    |
|----------|---------------------------------------------|
| `*.dat`  | Initial conditions (Gadget format)          |
| `*.vc`   | Circular velocity curve and Toomre Q        |
| `*.sd`   | Surface density profiles                    |
| `*.eqs`  | Effective equation of state                 |
| `*.gsd`  | Gas density field                           |

These can get large (tens of GB) and are ignored by git.

## Parameters

The main knobs in the parameter file (see `params/input_dwarf_01solar_mass.txt`):

- **Halo:** `CC` (concentration), `V200`, `LAMBDA` (spin), `REDSHIFT`, `Omega_m0`, `Omega_L0`
- **Mass fractions:** `MD` (disk), `MB` (bulge), `MBH` (black hole), `FBARY` (hot halo)
- **Disk:** `JD`, `GasFraction`, `DiskHeight`, `GasDiskScaleLength`, `HI_GasMassFraction`, ...
- **Resolution:** `N_HALO`, `N_DISK`, `N_GAS`, `N_BULGE`

See the comments in the file for details.

## Citation

If you use this code in published work, please cite the papers it is based on:

- Springel V., White S. D. M., 1999, MNRAS, 307, 162 - *Tidal tails in cold dark matter cosmologies*
  (disk galaxy model and IC method)
- Springel V., Hernquist L., 2003, MNRAS, 339, 289 - *Cosmological smoothed particle hydrodynamics
  simulations: a hybrid multiphase model for star formation* (effective equation of state used for
  the gas disk)
- Springel V., Di Matteo T., Hernquist L., 2005, MNRAS, 361, 776 - *Modelling feedback from stars
  and black holes in galaxy mergers* (gas disks in hydrostatic equilibrium, MakeDiskGalaxy)

## License

GNU General Public License v3, like GADGET-4. See [LICENSE](LICENSE).
The routines in `src/nrsrc/` are from *Numerical Recipes in C* and are not covered by this license.
