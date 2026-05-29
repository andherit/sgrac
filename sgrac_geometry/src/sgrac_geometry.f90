program sgrac_geometry
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  use forparse
  use LAT_mesh
  use LAT_distance
  use lists
  use distance
  use sgrac_vtk_tri_io
  use sgrac_cell_geometry
  implicit none

  type(mesh) :: amesh
  type(diff) :: adiff
  type(containerc), allocatable :: ntoc(:)
  type(containern), allocatable :: nton(:)
  real(pr), allocatable :: dg(:), area(:), dg_cell(:), centroid(:,:), grad_dg(:,:), theta(:)
  character(len=256) :: infile, outfile
  integer :: ierr, source, source_vtk

  infile = '-'
  outfile = '-'
  source = 1
  source_vtk = -1

  ierr = parse_arg('in', infile)
  ierr = parse_arg('out', outfile)
  ierr = parse_arg('source', source)
  ierr = parse_arg('source_vtk', source_vtk)
  if (ierr == PARSE_OK) source = source_vtk + 1

  call read_tri_vtk(trim(infile), amesh)
  if (source < 1 .or. source > amesh%Nnodes) stop 'sgrac-geometry: source node outside mesh range'

  allocate(dg(amesh%Nnodes))
  allocate(ntoc(amesh%Nnodes), nton(amesh%Nnodes))

  adiff%fast = .true.
  verbose = 0
  call pre_onevsall2d_onvertex(amesh, int(source,pin), dg, ntoc, nton)
  if (.not. associated(ntoc(source)%ptr)) then
     write(error_unit,'(a,i0,a)') 'sgrac-geometry: source node ', source, ' is not used by any triangle'
     write(error_unit,'(a)') 'sgrac-geometry: choose a point on the POLYGONS mesh, or use source_vtk=<0-based VTK point id>'
     stop 1
  endif
  call onevsall2d(amesh, dg, ntoc, nton, adiff)

  allocate(area(amesh%Ncells), dg_cell(amesh%Ncells), theta(amesh%Ncells))
  allocate(centroid(amesh%Ncells,3), grad_dg(amesh%Ncells,3))
  call compute_cell_geometry(amesh, dg, area, dg_cell, centroid, grad_dg, theta)

  call write_geometry_vtk(trim(outfile), amesh, dg, area, dg_cell, centroid, grad_dg, theta)

  call free_nton(nton)
  call free_ntoc(ntoc)
end program sgrac_geometry
