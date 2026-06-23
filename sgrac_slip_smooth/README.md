# sgrac-slip-smooth v0

`sgrac-slip-smooth` is the first SGRAC Phase B slip-field filter.

It reads a VTK legacy `POLYDATA` triangular rupture mesh, computes a smooth Gaussian slip field from intrinsic mesh distance, and writes the same mesh with a `CELL_DATA` scalar field named `slip`.

## Interface

Pipeline style:

```bash
./sgrac-slip-smooth center_node=1 sigma=4000.0 mw=6.5 < rupture.vtk > slip_raw.vtk
```

File style:

```bash
./sgrac-slip-smooth in=rupture.vtk out=slip_raw.vtk center_node=1 sigma=4000.0 mw=6.5
```

Arguments use the project `forparse` `key=value` convention.

| key | meaning | units |
|---|---|---:|
| `center_node` | 1-based project node index for the Gaussian center | none |
| `sigma` | Gaussian width | m |
| `mw` | target moment magnitude, required if `peak_slip` is absent | none |
| `mu` | shear modulus, physical mode default `3.0e10` | Pa |
| `peak_slip` | optional direct/debug peak slip; overrides `mw` and `mu` | m |
| `in` | optional input VTK file, default stdin | none |
| `out` | optional output VTK file, default stdout | none |

## Model

The Gaussian center is `center_node`, using the project 1-based node-index convention.

`sgrac-slip-smooth` rebuilds the same `trilat-distance` topology used by `sgrac-geometry`, computes nodal geodesic distance from `center_node`, and evaluates the Gaussian from the cell-average distance:

```text
r_i = (dg(n1) + dg(n2) + dg(n3)) / 3
shape_i = exp(-0.5 * r_i^2 / sigma^2)
```

This keeps the smooth slip field intrinsic to the fault surface instead of using straight 3D Euclidean distance.

In physical mode, the Gaussian shape is scaled to the target scalar moment:

```text
M0 = 10 ** (1.5 * Mw + 9.1)
M0 = mu * sum_i(slip_i * area_i)
```

For quick checks, `peak_slip=<m>` may be provided instead. If `peak_slip` is present, `mw` and `mu` are ignored.

The output `slip` field is `CELL_DATA`. Existing `POINT_DATA` and `CELL_DATA` are preserved when possible. If a cell scalar field named `slip` already exists, it is replaced.

## Validation

Using the sample rupture mesh from `sgrac_extract`:

```bash
make
./sgrac-slip-smooth center_node=1 sigma=4000.0 mw=6.5 < ../sgrac_extract/rupture.vtk > slip_raw.vtk
grep -n "CELL_DATA\|SCALARS slip" slip_raw.vtk
```

Expected result:

- output remains VTK legacy `POLYDATA`
- `POINTS` and `POLYGONS` geometry are unchanged
- `slip` is written as `CELL_DATA`
- the number of slip values equals the number of triangles
- slip values are non-negative

## Build

```bash
make
```
