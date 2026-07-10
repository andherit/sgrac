module sgrac_vtk_tri_io
  use generic
  use LAT_mesh
  implicit none

  character(len=*), parameter :: FIELD_DG = 'dg'
  character(len=*), parameter :: FIELD_DG_CELL = 'dg_cell'
  character(len=*), parameter :: FIELD_THETA = 'theta'
  character(len=*), parameter :: FIELD_GRAD_DG = 'grad_dg'
  character(len=*), parameter :: FIELD_DG_F1 = 'dg_f1'
  character(len=*), parameter :: FIELD_DG_F2 = 'dg_f2'
  character(len=*), parameter :: FIELD_DG_F1_CELL = 'dg_f1_cell'
  character(len=*), parameter :: FIELD_DG_F2_CELL = 'dg_f2_cell'
  character(len=*), parameter :: FIELD_DG_SUM_CELL = 'dg_sum_cell'

  type :: sgrac_vtk_metadata_type
     logical :: present = .false.
     integer :: format_version = 1
     integer :: mask_geometry_type_code = 0
     integer :: distance_source_count = 0
     integer :: source_node_indexes(2) = 0
     integer :: theta_convention_code = 0
     integer :: distance_units_code = 1
  end type sgrac_vtk_metadata_type

contains

subroutine read_tri_vtk(filename, amesh)
  character(*), intent(in) :: filename
  type(mesh), intent(out) :: amesh
  integer :: unit, ios, i, npts, npolys, total_size, nverts
  character(len=256) :: line, key, dtype

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 5
  else
     open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
     if (ios /= 0) stop 'sgrac-geometry: cannot open input VTK file'
  endif

  npts = -1
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) stop 'sgrac-geometry: POINTS section not found'
     read(line,*,iostat=ios) key, npts, dtype
     if (ios == 0 .and. trim(key) == 'POINTS') exit
  enddo

  amesh%Nnodes = npts
  allocate(amesh%px(npts), amesh%py(npts), amesh%pz(npts))
  do i=1,npts
     read(unit,*,iostat=ios) amesh%px(i), amesh%py(i), amesh%pz(i)
     if (ios /= 0) stop 'sgrac-geometry: error while reading POINTS'
  enddo

  npolys = -1
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) stop 'sgrac-geometry: POLYGONS/CELLS section not found'
     read(line,*,iostat=ios) key, npolys, total_size
     if (ios == 0 .and. (trim(key) == 'POLYGONS' .or. trim(key) == 'CELLS')) exit
  enddo

  amesh%Ncells = npolys
  allocate(amesh%cell(npolys,3))
  do i=1,npolys
     read(unit,*,iostat=ios) nverts, amesh%cell(i,1), amesh%cell(i,2), amesh%cell(i,3)
     if (ios /= 0) stop 'sgrac-geometry: error while reading triangle connectivity'
     if (nverts /= 3) stop 'sgrac-geometry: only triangular meshes are supported'
     amesh%cell(i,1:3) = amesh%cell(i,1:3) + 1_pin
  enddo

  if (unit /= 5) close(unit)
end subroutine read_tri_vtk

subroutine write_geometry_vtk(filename, amesh, dg, area, dg_cell, centroid, grad_dg, theta, metadata)
  character(*), intent(in) :: filename
  type(mesh), intent(in) :: amesh
  real(pr), intent(in) :: dg(amesh%Nnodes)
  real(pr), intent(in) :: area(amesh%Ncells), dg_cell(amesh%Ncells), theta(amesh%Ncells)
  real(pr), intent(in) :: centroid(amesh%Ncells,3), grad_dg(amesh%Ncells,3)
  type(sgrac_vtk_metadata_type), intent(in), optional :: metadata
  integer :: unit, ios, i
  real(pr) :: theta_deg(amesh%Ncells)

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-geometry: cannot open output VTK file'
  endif

  theta_deg = theta * 180._pr / acos(-1._pr)

  write(unit,'(a)') '# vtk DataFile Version 3.0'
  write(unit,'(a)') 'SGRAC geometry fields, SI units except theta in degrees'
  write(unit,'(a)') 'ASCII'
  write(unit,'(a)') 'DATASET POLYDATA'
  write(unit,'(a,i0,a)') 'POINTS ', amesh%Nnodes, ' double'
  do i=1,amesh%Nnodes
     write(unit,'(3(es24.16,1x))') amesh%px(i), amesh%py(i), amesh%pz(i)
  enddo
  write(unit,'(a,i0,1x,i0)') 'POLYGONS ', amesh%Ncells, 4*amesh%Ncells
  do i=1,amesh%Ncells
     write(unit,'(i0,3(1x,i0))') 3, amesh%cell(i,1)-1, amesh%cell(i,2)-1, amesh%cell(i,3)-1
  enddo

  if (present(metadata)) call write_sgrac_metadata(unit, metadata)

  write(unit,'(a,i0)') 'POINT_DATA ', amesh%Nnodes
  write(unit,'(a,a,a)') 'SCALARS ', FIELD_DG, ' double 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(unit,'(es24.16)') dg(i)
  enddo

  write(unit,'(a,i0)') 'CELL_DATA ', amesh%Ncells
  call write_cell_scalar(unit, 'area', area, amesh%Ncells)
  call write_cell_scalar(unit, FIELD_DG_CELL, dg_cell, amesh%Ncells)
  call write_cell_scalar(unit, FIELD_THETA, theta_deg, amesh%Ncells)
  call write_cell_vector(unit, 'centroid', centroid, amesh%Ncells)
  call write_cell_vector(unit, FIELD_GRAD_DG, grad_dg, amesh%Ncells)

  if (unit /= 6) close(unit)
end subroutine write_geometry_vtk

subroutine write_foci_geometry_vtk(filename, amesh, dg_f1, dg_f2, area, centroid, &
                                   dg_f1_cell, dg_f2_cell, dg_sum_cell, metadata)
  character(*), intent(in) :: filename
  type(mesh), intent(in) :: amesh
  real(pr), intent(in) :: dg_f1(amesh%Nnodes), dg_f2(amesh%Nnodes)
  real(pr), intent(in) :: area(amesh%Ncells), centroid(amesh%Ncells,3)
  real(pr), intent(in) :: dg_f1_cell(amesh%Ncells), dg_f2_cell(amesh%Ncells), dg_sum_cell(amesh%Ncells)
  type(sgrac_vtk_metadata_type), intent(in), optional :: metadata
  integer :: unit, ios, i

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-geometry: cannot open output VTK file'
  endif

  write(unit,'(a)') '# vtk DataFile Version 3.0'
  write(unit,'(a)') 'SGRAC foci geometry fields, SI units'
  write(unit,'(a)') 'ASCII'
  write(unit,'(a)') 'DATASET POLYDATA'
  write(unit,'(a,i0,a)') 'POINTS ', amesh%Nnodes, ' double'
  do i=1,amesh%Nnodes
     write(unit,'(3(es24.16,1x))') amesh%px(i), amesh%py(i), amesh%pz(i)
  enddo
  write(unit,'(a,i0,1x,i0)') 'POLYGONS ', amesh%Ncells, 4*amesh%Ncells
  do i=1,amesh%Ncells
     write(unit,'(i0,3(1x,i0))') 3, amesh%cell(i,1)-1, amesh%cell(i,2)-1, amesh%cell(i,3)-1
  enddo

  if (present(metadata)) call write_sgrac_metadata(unit, metadata)

  write(unit,'(a,i0)') 'POINT_DATA ', amesh%Nnodes
  call write_point_scalar(unit, FIELD_DG_F1, dg_f1, amesh%Nnodes)
  call write_point_scalar(unit, FIELD_DG_F2, dg_f2, amesh%Nnodes)

  write(unit,'(a,i0)') 'CELL_DATA ', amesh%Ncells
  call write_cell_scalar(unit, 'area', area, amesh%Ncells)
  call write_cell_scalar(unit, FIELD_DG_F1_CELL, dg_f1_cell, amesh%Ncells)
  call write_cell_scalar(unit, FIELD_DG_F2_CELL, dg_f2_cell, amesh%Ncells)
  call write_cell_scalar(unit, FIELD_DG_SUM_CELL, dg_sum_cell, amesh%Ncells)
  call write_cell_vector(unit, 'centroid', centroid, amesh%Ncells)

  if (unit /= 6) close(unit)
end subroutine write_foci_geometry_vtk

subroutine write_sgrac_metadata(unit, metadata)
  integer, intent(in) :: unit
  type(sgrac_vtk_metadata_type), intent(in) :: metadata

  write(unit,'(a)') 'FIELD sgrac_metadata 6'
  write(unit,'(a)') 'sgrac_format_version 1 1 int'
  write(unit,'(i0)') metadata%format_version
  write(unit,'(a)') 'mask_geometry_type_code 1 1 int'
  write(unit,'(i0)') metadata%mask_geometry_type_code
  write(unit,'(a)') 'distance_source_count 1 1 int'
  write(unit,'(i0)') metadata%distance_source_count
  write(unit,'(a)') 'source_node_indexes 2 1 int'
  write(unit,'(2(i0,1x))') metadata%source_node_indexes
  write(unit,'(a)') 'theta_convention_code 1 1 int'
  write(unit,'(i0)') metadata%theta_convention_code
  write(unit,'(a)') 'distance_units_code 1 1 int'
  write(unit,'(i0)') metadata%distance_units_code
end subroutine write_sgrac_metadata

subroutine write_cell_scalar(unit, name, vals, n)
  integer, intent(in) :: unit, n
  character(*), intent(in) :: name
  real(pr), intent(in) :: vals(n)
  integer :: i
  write(unit,'(a,a,a)') 'SCALARS ', trim(name), ' double 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i=1,n
     write(unit,'(es24.16)') vals(i)
  enddo
end subroutine write_cell_scalar

subroutine write_point_scalar(unit, name, vals, n)
  integer, intent(in) :: unit, n
  character(*), intent(in) :: name
  real(pr), intent(in) :: vals(n)
  integer :: i
  write(unit,'(a,a,a)') 'SCALARS ', trim(name), ' double 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i=1,n
     write(unit,'(es24.16)') vals(i)
  enddo
end subroutine write_point_scalar

subroutine write_cell_vector(unit, name, vals, n)
  integer, intent(in) :: unit, n
  character(*), intent(in) :: name
  real(pr), intent(in) :: vals(n,3)
  integer :: i
  write(unit,'(a,a,a)') 'VECTORS ', trim(name), ' double'
  do i=1,n
     write(unit,'(3(es24.16,1x))') vals(i,1), vals(i,2), vals(i,3)
  enddo
end subroutine write_cell_vector

end module sgrac_vtk_tri_io
