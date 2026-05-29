module sgrac_vtk_tri_io
  use generic
  use LAT_mesh
  implicit none
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

subroutine write_geometry_vtk(filename, amesh, dg, area, dg_cell, centroid, grad_dg, theta)
  character(*), intent(in) :: filename
  type(mesh), intent(in) :: amesh
  real(pr), intent(in) :: dg(amesh%Nnodes)
  real(pr), intent(in) :: area(amesh%Ncells), dg_cell(amesh%Ncells), theta(amesh%Ncells)
  real(pr), intent(in) :: centroid(amesh%Ncells,3), grad_dg(amesh%Ncells,3)
  integer :: unit, ios, i

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-geometry: cannot open output VTK file'
  endif

  write(unit,'(a)') '# vtk DataFile Version 3.0'
  write(unit,'(a)') 'SGRAC geometry fields, SI units'
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

  write(unit,'(a,i0)') 'POINT_DATA ', amesh%Nnodes
  write(unit,'(a)') 'SCALARS dg double 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i=1,amesh%Nnodes
     write(unit,'(es24.16)') dg(i)
  enddo

  write(unit,'(a,i0)') 'CELL_DATA ', amesh%Ncells
  call write_cell_scalar(unit, 'area', area, amesh%Ncells)
  call write_cell_scalar(unit, 'dg_cell', dg_cell, amesh%Ncells)
  call write_cell_scalar(unit, 'theta', theta, amesh%Ncells)
  call write_cell_vector(unit, 'centroid', centroid, amesh%Ncells)
  call write_cell_vector(unit, 'grad_dg', grad_dg, amesh%Ncells)

  if (unit /= 6) close(unit)
end subroutine write_geometry_vtk

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
