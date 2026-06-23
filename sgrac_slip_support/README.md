# sgrac-slip-support v0

`sgrac_slip_support` contains reusable support routines for SGRAC slip-field modules.

It is a support/library package only. It does not build a standalone pipeline executable.

## Routines

- compute triangle cell centers from node coordinates and triangular connectivity
- compute triangle cell areas in m2
- compute scalar seismic moment
- convert between scalar moment and moment magnitude
- rescale a cell slip field to a target scalar moment

All dimensional quantities use S.I. units:

- coordinates in meters
- areas in m2
- slip in meters
- shear modulus in Pa
- scalar moment in N m

## Moment

Scalar seismic moment is computed as:

```text
M0 = mu * sum_i(slip_i * area_i)
```

Magnitude conversion uses:

```text
Mw = (2/3) * (log10(M0) - 9.1)
M0 = 10 ** (1.5 * Mw + 9.1)
```

## Build

```bash
make
```

This compiles `generic.o`, `sgrac_slip_support.o`, and the corresponding module files.
