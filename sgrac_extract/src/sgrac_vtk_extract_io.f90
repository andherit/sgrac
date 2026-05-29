module sgrac_vtk_extract_io
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  implicit none
  integer, parameter :: line_len = 512
contains

subroutine read_parent_vtk(filename, fieldname, px, py, pz, cell, field_values, npoints, ncells)
  character(*), intent(in) :: filename, fieldname
  real(pr), allocatable, intent(out) :: px(:), py(:), pz(:)
  integer(pin), allocatable, intent(out) :: cell(:,:), field_values(:)
  integer, intent(out) :: npoints, ncells
  character(len=line_len), allocatable :: lines(:)
  integer :: nlines

  call read_text_file(filename, lines, nlines)
  call read_points(lines, nlines, px, py, pz, npoints)
  call read_triangles(lines, nlines, cell, ncells)
  call read_cell_scalar_as_int(lines, nlines, ncells, fieldname, field_values)
end subroutine read_parent_vtk

subroutine write_child_vtk(filename, px, py, pz, cell)
  character(*), intent(in) :: filename
  real(pr), intent(in) :: px(:), py(:), pz(:)
  integer(pin), intent(in) :: cell(:,:)
  integer :: unit, ios, i, npoints, ncells

  npoints = size(px)
  ncells = size(cell, 1)

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-extract: cannot open output file'
  endif

  write(unit,'(a)') '# vtk DataFile Version 3.0'
  write(unit,'(a)') 'SGRAC extracted rupture mesh, SI units'
  write(unit,'(a)') 'ASCII'
  write(unit,'(a)') 'DATASET POLYDATA'
  write(unit,'(a,i0,a)') 'POINTS ', npoints, ' double'
  do i=1,npoints
     write(unit,'(3(es24.16,1x))') px(i), py(i), pz(i)
  enddo
  write(unit,'(a,i0,1x,i0)') 'POLYGONS ', ncells, 4*ncells
  do i=1,ncells
     write(unit,'(i0,3(1x,i0))') 3, cell(i,1)-1, cell(i,2)-1, cell(i,3)-1
  enddo

  if (unit /= 6) close(unit)
end subroutine write_child_vtk

subroutine read_text_file(filename, lines, nlines)
  character(*), intent(in) :: filename
  character(len=line_len), allocatable, intent(out) :: lines(:)
  integer, intent(out) :: nlines
  character(len=line_len), allocatable :: tmp(:)
  character(len=line_len) :: line
  integer :: unit, ios, cap

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 5
  else
     open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
     if (ios /= 0) stop 'sgrac-extract: cannot open input file'
  endif

  cap = 1024
  allocate(lines(cap))
  nlines = 0
  do
     read(unit,'(a)',iostat=ios) line
     if (ios /= 0) exit
     if (nlines == cap) then
        allocate(tmp(cap))
        tmp = lines
        deallocate(lines)
        allocate(lines(2*cap))
        lines(1:cap) = tmp
        deallocate(tmp)
        cap = 2*cap
     endif
     nlines = nlines + 1
     lines(nlines) = line
  enddo

  if (unit /= 5) close(unit)
end subroutine read_text_file

subroutine read_points(lines, nlines, px, py, pz, npoints)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  real(pr), allocatable, intent(out) :: px(:), py(:), pz(:)
  integer, intent(out) :: npoints
  character(len=64) :: key, dtype
  integer :: i, j, ios, point_line

  point_line = -1
  npoints = -1
  do i=1,nlines
     key = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, npoints, dtype
     if (ios == 0 .and. trim(key) == 'POINTS') then
        point_line = i
        exit
     endif
  enddo
  if (point_line < 0) stop 'sgrac-extract: POINTS section not found'
  if (npoints <= 0) stop 'sgrac-extract: POINTS count must be positive'
  if (point_line + npoints > nlines) stop 'sgrac-extract: POINTS section is truncated'

  allocate(px(npoints), py(npoints), pz(npoints))
  do j=1,npoints
     read(lines(point_line+j),*,iostat=ios) px(j), py(j), pz(j)
     if (ios /= 0) stop 'sgrac-extract: error while reading POINTS'
  enddo
end subroutine read_points

subroutine read_triangles(lines, nlines, cell, ncells)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer(pin), allocatable, intent(out) :: cell(:,:)
  integer, intent(out) :: ncells
  character(len=64) :: key
  integer :: i, j, ios, total_size, nverts, poly_line

  poly_line = -1
  ncells = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key, ncells, total_size
     if (ios == 0 .and. trim(key) == 'POLYGONS') then
        poly_line = i
        exit
     endif
     if (ios == 0 .and. trim(key) == 'CELLS') then
        stop 'sgrac-extract: CELLS input is not supported; use POLYDATA POLYGONS'
     endif
  enddo
  if (poly_line < 0) stop 'sgrac-extract: POLYGONS section not found'
  if (ncells <= 0) stop 'sgrac-extract: POLYGONS count must be positive'
  if (poly_line + ncells > nlines) stop 'sgrac-extract: POLYGONS section is truncated'

  allocate(cell(ncells,3))
  do j=1,ncells
     read(lines(poly_line+j),*,iostat=ios) nverts, cell(j,1), cell(j,2), cell(j,3)
     if (ios /= 0) stop 'sgrac-extract: error while reading POLYGONS'
     if (nverts /= 3) stop 'sgrac-extract: only triangular POLYGONS are supported'
     cell(j,1:3) = cell(j,1:3) + 1_pin
  enddo
end subroutine read_triangles

subroutine read_cell_scalar_as_int(lines, nlines, ncells, fieldname, values)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncells
  character(*), intent(in) :: fieldname
  integer(pin), allocatable, intent(out) :: values(:)
  character(len=64) :: key, name, dtype
  integer :: i, j, ios, cell_data_line, lookup_line
  real(pr) :: value_real

  cell_data_line = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'CELL_DATA') then
        cell_data_line = i
        exit
     endif
  enddo
  if (cell_data_line < 0) stop 'sgrac-extract: CELL_DATA section not found'

  allocate(values(ncells))
  i = cell_data_line + 1
  do while (i <= nlines)
     key = ''
     name = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, name, dtype
     if (ios == 0 .and. trim(key) == 'SCALARS' .and. trim(name) == trim(fieldname)) then
        lookup_line = i + 1
        if (lookup_line > nlines) stop 'sgrac-extract: malformed scalar field'
        if (index(adjustl(lines(lookup_line)), 'LOOKUP_TABLE') /= 1) then
           stop 'sgrac-extract: scalar field missing LOOKUP_TABLE line'
        endif
        do j=1,ncells
           if (lookup_line + j > nlines) stop 'sgrac-extract: scalar field is truncated'
           read(lines(lookup_line+j),*,iostat=ios) value_real
           if (ios /= 0) stop 'sgrac-extract: error reading scalar field values'
           values(j) = int(value_real, pin)
        enddo
        return
     endif

     if (ios == 0 .and. trim(key) == 'SCALARS') then
        i = i + 2 + ncells
     else if (ios == 0 .and. trim(key) == 'VECTORS') then
        i = i + 1 + ncells
     else if (ios == 0 .and. trim(key) == 'FIELD') then
        stop 'sgrac-extract: FIELD cell data are not supported in v0'
     else
        i = i + 1
     endif
  enddo

  write(error_unit,'(a,a)') 'sgrac-extract: missing requested CELL_DATA scalar ', trim(fieldname)
  stop 1
end subroutine read_cell_scalar_as_int

end module sgrac_vtk_extract_io
