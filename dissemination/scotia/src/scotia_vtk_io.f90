module scotia_vtk_io
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  use LAT_mesh
  implicit none
contains

subroutine read_scotia_vtk(filename, amesh, slip)
  character(*), intent(in) :: filename
  type(mesh), intent(out) :: amesh
  real(pr), allocatable, intent(out) :: slip(:)
  integer :: unit, ios, i, npts, ncells, total_size, nverts
  character(len=256) :: line, key, name, dtype

  open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'scotia: cannot open input VTK file'

  npts = -1
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) stop 'scotia: POINTS section not found'
     key = ''
     dtype = ''
     read(line,*,iostat=ios) key, npts, dtype
     if (ios == 0 .and. trim(key) == 'POINTS') exit
  enddo

  amesh%Nnodes = npts
  allocate(amesh%px(npts), amesh%py(npts), amesh%pz(npts))
  do i = 1, npts
     read(unit,*,iostat=ios) amesh%px(i), amesh%py(i), amesh%pz(i)
     if (ios /= 0) stop 'scotia: error while reading POINTS'
  enddo

  ncells = -1
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) stop 'scotia: POLYGONS/CELLS section not found'
     key = ''
     read(line,*,iostat=ios) key, ncells, total_size
     if (ios == 0 .and. (trim(key) == 'POLYGONS' .or. trim(key) == 'CELLS')) exit
  enddo

  amesh%Ncells = ncells
  allocate(amesh%cell(ncells,3))
  do i = 1, ncells
     read(unit,*,iostat=ios) nverts, amesh%cell(i,1), amesh%cell(i,2), amesh%cell(i,3)
     if (ios /= 0) stop 'scotia: error while reading triangle connectivity'
     if (nverts /= 3) stop 'scotia: only triangular cells are supported'
     amesh%cell(i,1:3) = amesh%cell(i,1:3) + 1_pin
  enddo

  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) stop 'scotia: CELL_DATA section not found'
     key = ''
     read(line,*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'CELL_DATA') exit
  enddo

  allocate(slip(ncells))
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) stop 'scotia: slip scalar not found'
     key = ''
     name = ''
     dtype = ''
     read(line,*,iostat=ios) key, name, dtype
     if (ios == 0 .and. trim(key) == 'SCALARS' .and. trim(name) == 'slip') exit
  enddo

  read(unit,'(a)',iostat=ios) line
  if (ios /= 0 .or. index(adjustl(line), 'LOOKUP_TABLE') /= 1) then
     stop 'scotia: slip scalar missing LOOKUP_TABLE'
  endif

  do i = 1, ncells
     read(unit,*,iostat=ios) slip(i)
     if (ios /= 0) stop 'scotia: error while reading slip values'
  enddo

  close(unit)
end subroutine read_scotia_vtk

subroutine write_scotia_vtk(filename, amesh, slip, t_const, t_prop, t_inv)
  character(*), intent(in) :: filename
  type(mesh), intent(in) :: amesh
  real(pr), intent(in) :: slip(amesh%Ncells)
  real(pr), intent(in) :: t_const(amesh%Nnodes), t_prop(amesh%Nnodes), t_inv(amesh%Nnodes)
  integer :: unit, ios, i

  open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
  if (ios /= 0) stop 'scotia: cannot open output VTK file'

  write(unit,'(a)') '# vtk DataFile Version 2.0'
  write(unit,'(a)') 'scotia rupture times'
  write(unit,'(a)') 'ASCII'
  write(unit,'(a)') 'DATASET POLYDATA'
  write(unit,'(a,i0,a)') 'POINTS ', amesh%Nnodes, ' double'
  do i = 1, amesh%Nnodes
     write(unit,'(3(es24.16,1x))') amesh%px(i), amesh%py(i), amesh%pz(i)
  enddo

  write(unit,'(a,i0,1x,i0)') 'POLYGONS ', amesh%Ncells, 4*amesh%Ncells
  do i = 1, amesh%Ncells
     write(unit,'(i0,3(1x,i0))') 3, amesh%cell(i,1)-1, amesh%cell(i,2)-1, amesh%cell(i,3)-1
  enddo

  write(unit,'(a,i0)') 'CELL_DATA ', amesh%Ncells
  call write_cell_scalar(unit, 'slip', slip)

  write(unit,'(a,i0)') 'POINT_DATA ', amesh%Nnodes
  call write_node_scalar(unit, 'rupt_time_constant', t_const)
  call write_node_scalar(unit, 'rupt_time_slip_prop', t_prop)
  call write_node_scalar(unit, 'rupt_time_slip_inv', t_inv)

  close(unit)
end subroutine write_scotia_vtk

subroutine write_cell_scalar(unit, name, vals)
  integer, intent(in) :: unit
  character(*), intent(in) :: name
  real(pr), intent(in) :: vals(:)
  integer :: i

  write(unit,'(a,a,a)') 'SCALARS ', trim(name), ' double 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i = 1, size(vals)
     write(unit,'(es24.16)') vals(i)
  enddo
end subroutine write_cell_scalar

subroutine write_node_scalar(unit, name, vals)
  integer, intent(in) :: unit
  character(*), intent(in) :: name
  real(pr), intent(in) :: vals(:)
  integer :: i

  write(unit,'(a,a,a)') 'SCALARS ', trim(name), ' double 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i = 1, size(vals)
     write(unit,'(es24.16)') vals(i)
  enddo
end subroutine write_node_scalar

end module scotia_vtk_io
