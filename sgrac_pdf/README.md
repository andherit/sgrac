# sgrac-pdf

`sgrac-pdf` builds probability-density fields on an
irregular triangular rupture mesh.

The source package directory is `sgrac_pdf`. The executable is `sgrac-pdf`.

Available modes:

```text
mode=addgaussian
mode=normalize
```

## mode=addgaussian

`mode=addgaussian` adds one intrinsic geodesic Gaussian component to the
temporary `CELL_DATA` scalar field named `pdf_weight`:

```text
d_i = (d_n1 + d_n2 + d_n3) / 3
G_i = exp(-0.5 * (d_i / sigma)**2)
pdf_weight_i = pdf_weight_i + amplitude * G_i
```

The nodal distances are geodesic distances on the mesh from `center_node`,
computed with the same `trilat-distance` routines used by
`sgrac-slip-smooth`. No Euclidean fallback is used.

If the input file has no `pdf_weight` field, the accumulator starts from zero.
If the field already exists, it is read and updated.

No normalization is performed in this mode. `amplitude` is the peak value added
before any future normalization. The Gaussian is not area-weighted.

## mode=normalize

`mode=normalize` converts the temporary accumulator into the final
dimensionless cell-probability `CELL_DATA` scalar field named `pdf`:

```text
I = sum_i(pdf_weight_i * area_i)
pdf_i = pdf_weight_i * area_i / I
sum_i(pdf_i) = 1
```

Triangle areas are recomputed from the mesh geometry. The area correction is
already applied during normalization, so the final `pdf` is a dimensionless
probability mass associated with each triangle, not a probability density.

The normalize output is a clean interface file for K223D. It removes internal
construction fields:

```text
pdf_weight
area
```

It preserves unrelated downstream fields such as:

```text
velocity
surface
rupt_time
slip
pdf
```

If an old `pdf` field exists, it is replaced with the newly normalized one.

## Interface

```bash
sgrac-pdf mode=addgaussian center_node=5702 sigma=50000.0 amplitude=1.0 in=rupture.vtk out=pdf_weight_1.vtk
sgrac-pdf mode=normalize in=pdf_weight_1.vtk out=rupture_pdf.vtk
```

Arguments use the project `forparse` `key=value` convention.

| key | meaning | default |
|---|---|---:|
| `mode` | `addgaussian` or `normalize` | required |
| `center_node` | 1-based node index for `addgaussian` | required for `addgaussian` |
| `sigma` | Gaussian width in metres | required for `addgaussian` |
| `amplitude` | peak value added to `pdf_weight` | `1.0` for `addgaussian` |
| `in` | input VTK file | stdin |
| `out` | output VTK file | stdout |

Diagnostics are written to stderr so stdout remains valid VTK in pipelines.

## Examples

Single Gaussian:

```bash
sgrac-pdf \
    mode=addgaussian \
    center_node=5702 \
    sigma=50000.0 \
    in=rupture.vtk \
    out=pdf_weight_1.vtk
```

Sum of two Gaussians:

```bash
sgrac-pdf \
    mode=addgaussian \
    center_node=5702 \
    sigma=50000.0 \
    amplitude=1.0 \
    in=rupture.vtk |
sgrac-pdf \
    mode=addgaussian \
    center_node=8200 \
    sigma=25000.0 \
    amplitude=0.5 \
    out=pdf_weight_2.vtk
```

The second filter reads and updates the existing `pdf_weight` field.

Complete add-and-normalize pipeline:

```bash
sgrac-pdf \
    mode=addgaussian \
    center_node=97 \
    sigma=50000.0 \
    amplitude=1.0 \
    in=rupture.vtk |
sgrac-pdf \
    mode=addgaussian \
    center_node=3984 \
    sigma=25000.0 \
    amplitude=2.0 |
sgrac-pdf \
    mode=normalize \
    out=rupture_pdf.vtk
```

## Field Convention

`pdf_weight` is an unnormalized intermediate accumulator. It is not the final
normalized `pdf` field.

`pdf` is the final cell-probability field intended for direct compatibility
with the current K223D cumulative sampling implementation. The area correction
has already been applied by `sgrac-pdf`, so K223D can consume `pdf` directly as
a discrete probability attached to each triangle.

## Build

```bash
make
```

The executable is written to:

```text
../bin/sgrac-pdf
```
