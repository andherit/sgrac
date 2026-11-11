# rectangular mesh example

As every seismologist knows, a fault is a flat rectangle with a ratio 2:1.

## step up

All the binary are inside `$HOME_SGRAC/bin` and this directory belongs to your `$PATH`


## slip building

### sgrac-gmsh-support

We use first the `sgrac-gmsh-support` to mesh a rectangle of length `lx` and width `lz` with
a spacing dimension `lc`. It is using the GMSH api.

```bash
sgrac-gmsh-support lx=120000 lz=80000 lc=500 > parent.vtk
```

The center of the "fault" is an embeded node, with `id=0` and coordinates `(0,0,0)`

### sgrac-geometry

The next step is the computation of different fields on the mesh in order to prepare the 
extraction of a submesh associated to a given magnitude earthquake.

The key factor here is the shape of this subfault. In real complex fault surfaces, define a
rectangle is not an easy task. However considering the earthquake rupture surface as a rectangle
has not really a meaning anymore with unstructured meshes.

The first possibility is to consider a anysotropic contour function depending on theta.

<table border="0">
  <tr valign="middle">
    <td align="center" width="33%"><img src="images/area.png"  alt="area" /></td>
    <td align="center" width="33%"><img src="images/theta.png" alt="theta" /></td>
    <td align="center" width="33%"><img src="images/distance.png" alt="distance" /></td>
  </tr>
  <tr valign="middle">
    <td align="center" width="33%">area</td>
    <td align="center" width="33%">theta</td>
    <td align="center" width="33%">distance</td>
  </tr>
</table>