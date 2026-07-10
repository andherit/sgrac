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
  type(sgrac_vtk_metadata_type) :: metadata
  real(pr), allocatable :: dg(:), area(:), dg_cell(:), centroid(:,:), grad_dg(:,:), theta(:)
  real(pr), allocatable :: dg_f1(:), dg_f2(:), dg_f1_cell(:), dg_f2_cell(:), dg_sum_cell(:)
  character(len=256) :: infile, outfile
  character(len=64) :: geometry
  integer :: ierr, ierr_geometry_code, ierr_focus1, ierr_focus2, source, source_vtk
  integer :: geometry_code, focus1_node, focus2_node, focus1_vtk, focus2_vtk

  infile = '-'
  outfile = '-'
  geometry = 'theta_ellipse'
  geometry_code = 2
  source = 1
  source_vtk = -1
  focus1_node = -1
  focus2_node = -1
  focus1_vtk = -1
  focus2_vtk = -1

  ierr = parse_arg('in', infile)
  ierr = parse_arg('out', outfile)
  ierr = parse_arg('geometry', geometry)
  ierr_geometry_code = parse_arg('mask_geometry_type_code', geometry_code)
  ierr = parse_arg('source', source)
  ierr = parse_arg('source_vtk', source_vtk)
  if (ierr == PARSE_OK) source = source_vtk + 1
  ierr_focus1 = parse_arg('focus1_node', focus1_node)
  ierr_focus2 = parse_arg('focus2_node', focus2_node)
  ierr = parse_arg('focus1_vtk', focus1_vtk)
  if (ierr == PARSE_OK) then
     focus1_node = focus1_vtk + 1
     ierr_focus1 = PARSE_OK
  endif
  ierr = parse_arg('focus2_vtk', focus2_vtk)
  if (ierr == PARSE_OK) then
     focus2_node = focus2_vtk + 1
     ierr_focus2 = PARSE_OK
  endif

  if (ierr_geometry_code /= PARSE_OK) then
     select case(trim(geometry))
     case('theta_ellipse','theta')
        geometry_code = 2
     case('foci_ellipse','foci')
        geometry_code = 3
     case default
        write(error_unit,'(a,a)') 'sgrac-geometry: unknown geometry mode: ', trim(geometry)
        stop 1
     end select
  endif

  call read_tri_vtk(trim(infile), amesh)

  verbose = 0
  adiff%fast = .true.

  select case(geometry_code)
  case(2)
     if (source < 1 .or. source > amesh%Nnodes) stop 'sgrac-geometry: source node outside mesh range'

     allocate(dg(amesh%Nnodes))
     call compute_geodesic_from_source(source, dg)

     allocate(area(amesh%Ncells), dg_cell(amesh%Ncells), theta(amesh%Ncells))
     allocate(centroid(amesh%Ncells,3), grad_dg(amesh%Ncells,3))
     call compute_cell_geometry(amesh, dg, area, dg_cell, centroid, grad_dg, theta)

     metadata%present = .true.
     metadata%mask_geometry_type_code = 2
     metadata%distance_source_count = 1
     metadata%source_node_indexes = (/source, 0/)
     metadata%theta_convention_code = 2
     metadata%distance_units_code = 2
     call write_geometry_vtk(trim(outfile), amesh, dg, area, dg_cell, centroid, grad_dg, theta, metadata)
  case(3)
     if (ierr_focus1 /= PARSE_OK .or. ierr_focus2 /= PARSE_OK) then
        write(error_unit,'(a)') 'sgrac-geometry: foci_ellipse requires focus1_node=<id> and focus2_node=<id>'
        stop 1
     endif
     if (focus1_node < 1 .or. focus1_node > amesh%Nnodes) stop 'sgrac-geometry: focus1_node outside mesh range'
     if (focus2_node < 1 .or. focus2_node > amesh%Nnodes) stop 'sgrac-geometry: focus2_node outside mesh range'
     if (focus1_node == focus2_node) stop 'sgrac-geometry: focus nodes must be distinct'

     allocate(dg_f1(amesh%Nnodes), dg_f2(amesh%Nnodes))
     call compute_geodesic_from_source(focus1_node, dg_f1)
     call compute_geodesic_from_source(focus2_node, dg_f2)

     allocate(area(amesh%Ncells), centroid(amesh%Ncells,3))
     allocate(dg_f1_cell(amesh%Ncells), dg_f2_cell(amesh%Ncells), dg_sum_cell(amesh%Ncells))
     ! Foci-ellipse geometry prepares two geodesic distance fields and their
     ! cellwise sum; dg_sum_cell is consumed later by sgrac-mask.
     call compute_foci_cell_geometry(amesh, dg_f1, dg_f2, area, centroid, dg_f1_cell, dg_f2_cell, dg_sum_cell)

     metadata%present = .true.
     metadata%mask_geometry_type_code = 3
     metadata%distance_source_count = 2
     metadata%source_node_indexes = (/focus1_node, focus2_node/)
     metadata%theta_convention_code = 0
     metadata%distance_units_code = 2
     call write_foci_geometry_vtk(trim(outfile), amesh, dg_f1, dg_f2, area, centroid, &
                                  dg_f1_cell, dg_f2_cell, dg_sum_cell, metadata)
  case default
     write(error_unit,'(a,i0)') 'sgrac-geometry: unsupported mask_geometry_type_code = ', geometry_code
     stop 1
  end select
contains

subroutine compute_geodesic_from_source(source_node, distarray)
  integer, intent(in) :: source_node
  real(pr), intent(out) :: distarray(amesh%Nnodes)
  type(containerc), allocatable :: ntoc(:)
  type(containern), allocatable :: nton(:)

  allocate(ntoc(amesh%Nnodes), nton(amesh%Nnodes))
  call pre_onevsall2d_onvertex(amesh, int(source_node,pin), distarray, ntoc, nton)
  if (.not. associated(ntoc(source_node)%ptr)) then
     write(error_unit,'(a,i0,a)') 'sgrac-geometry: source node ', source_node, ' is not used by any triangle'
     write(error_unit,'(a)') 'sgrac-geometry: choose a point on the POLYGONS mesh, or use a *_vtk=<0-based VTK point id>'
     stop 1
  endif
  call onevsall2d(amesh, distarray, ntoc, nton, adiff)
  call free_nton(nton)
  call free_ntoc(ntoc)
end subroutine compute_geodesic_from_source
end program sgrac_geometry
