# sgrac-slip-support v0

`sgrac-slip-support` contains reusable support routines for SGRAC mesh-field and slip-field modules.

It is a support/library package only. It does not build a standalone pipeline executable.

The support layer is split into reusable modules so future field-generation filters, such as `sgrac-pdf`, can reuse generic mesh and VTK utilities without depending on slip-specific physics.

## `sgrac_mesh_field_support`

Generic mesh geometry utilities:

- compute triangle cell centers from node coordinates and triangular connectivity
- compute triangle cell areas in m2

This module contains no slip-specific or earthquake-specific logic.

## `sgrac_slip_support`

Slip-physics utilities:

- compute scalar seismic moment
- convert between scalar moment and moment magnitude
- rescale a cell slip field to a target scalar moment

## `sgrac_vtk_field_io`

Generic VTK legacy ASCII `POLYDATA` field I/O:

- read a triangular mesh and preserve the original text lines
- read arbitrary existing `CELL_DATA` scalar fields
- write or replace arbitrary `CELL_DATA` scalar fields
- preserve unrelated fields and metadata

When adding a new cell scalar, the writer respects VTK section boundaries. If a file has `POINT_DATA`, the new `CELL_DATA` scalar is inserted before it. This module is shared infrastructure for future field-based filters; `sgrac-pdf` is not implemented here.

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

This compiles `generic.o`, `sgrac_mesh_field_support.o`, `sgrac_slip_support.o`, `sgrac_vtk_field_io.o`, and the corresponding module files.
