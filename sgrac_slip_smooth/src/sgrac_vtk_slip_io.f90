module sgrac_vtk_slip_io
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  implicit none
  integer, parameter :: line_len = 512
contains

subroutine read_slip_vtk(filename, lines, nlines, px, py, pz, cell, npoints, ncells)
  character(*), intent(in) :: filename
  character(len=line_len), allocatable, intent(out) :: lines(:)
  integer, intent(out) :: nlines, npoints, ncells
  real(pr), allocatable, intent(out) :: px(:), py(:), pz(:)
  integer(pin), allocatable, intent(out) :: cell(:,:)

  call read_text_file(filename, lines, nlines)
  call read_points(lines, nlines, px, py, pz, npoints)
  call read_triangles(lines, nlines, cell, ncells)
end subroutine read_slip_vtk

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
     if (ios /= 0) stop 'sgrac-slip-smooth: cannot open input file'
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
  if (point_line < 0) stop 'sgrac-slip-smooth: POINTS section not found'
  if (npoints <= 0) stop 'sgrac-slip-smooth: POINTS count must be positive'
  if (point_line + npoints > nlines) stop 'sgrac-slip-smooth: POINTS section is truncated'

  allocate(px(npoints), py(npoints), pz(npoints))
  do j=1,npoints
     read(lines(point_line+j),*,iostat=ios) px(j), py(j), pz(j)
     if (ios /= 0) stop 'sgrac-slip-smooth: error while reading POINTS'
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
        stop 'sgrac-slip-smooth: CELLS input is not supported; use POLYDATA POLYGONS'
     endif
  enddo
  if (poly_line < 0) stop 'sgrac-slip-smooth: POLYGONS section not found'
  if (ncells <= 0) stop 'sgrac-slip-smooth: POLYGONS count must be positive'
  if (poly_line + ncells > nlines) stop 'sgrac-slip-smooth: POLYGONS section is truncated'

  allocate(cell(ncells,3))
  do j=1,ncells
     read(lines(poly_line+j),*,iostat=ios) nverts, cell(j,1), cell(j,2), cell(j,3)
     if (ios /= 0) stop 'sgrac-slip-smooth: error while reading POLYGONS'
     if (nverts /= 3) stop 'sgrac-slip-smooth: only triangular POLYGONS are supported'
     cell(j,1:3) = cell(j,1:3) + 1_pin
  enddo
end subroutine read_triangles

subroutine write_vtk_with_cell_scalar(filename, lines, nlines, fieldname, vals, ncell)
  character(*), intent(in) :: filename, fieldname
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell
  real(pr), intent(in) :: vals(ncell)
  integer :: unit, ios, i, cell_data_line, cell_count
  character(len=64) :: key

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-slip-smooth: cannot open output file'
  endif

  cell_data_line = find_cell_data_line(lines, nlines, cell_count)
  if (cell_data_line > 0 .and. cell_count /= ncell) then
     write(error_unit,'(a,i0,a,i0)') 'sgrac-slip-smooth: CELL_DATA count ', cell_count, &
        ' does not match POLYGONS count ', ncell
     stop 1
  endif

  if (cell_data_line < 0) then
     do i=1,nlines
        write(unit,'(a)') trim(lines(i))
     enddo
     write(unit,'(a,i0)') 'CELL_DATA ', ncell
  else
     i = 1
     do while (i <= nlines)
        key = ''
        if (i > cell_data_line) then
           read(lines(i),*,iostat=ios) key
           if (ios == 0 .and. trim(key) == 'SCALARS') then
              if (is_requested_scalar(lines(i), fieldname)) then
                 if (i + 1 + ncell > nlines) stop 'sgrac-slip-smooth: existing slip field is truncated'
                 i = i + 2 + ncell
                 cycle
              endif
           endif
        endif
        write(unit,'(a)') trim(lines(i))
        i = i + 1
     enddo
  endif

  call write_scalar_real(unit, fieldname, vals, ncell)

  if (unit /= 6) close(unit)
end subroutine write_vtk_with_cell_scalar

integer function find_cell_data_line(lines, nlines, cell_count) result(cell_data_line)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer, intent(out) :: cell_count
  integer :: i, ios
  character(len=64) :: key

  cell_data_line = -1
  cell_count = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key, cell_count
     if (ios == 0 .and. trim(key) == 'CELL_DATA') then
        cell_data_line = i
        return
     endif
  enddo
end function find_cell_data_line

logical function is_requested_scalar(line, fieldname)
  character(len=line_len), intent(in) :: line
  character(*), intent(in) :: fieldname
  integer :: ios
  character(len=64) :: key, name

  key = ''
  name = ''
  read(line,*,iostat=ios) key, name
  is_requested_scalar = ios == 0 .and. trim(key) == 'SCALARS' .and. trim(name) == trim(fieldname)
end function is_requested_scalar

subroutine write_scalar_real(unit, name, vals, n)
  integer, intent(in) :: unit, n
  character(*), intent(in) :: name
  real(pr), intent(in) :: vals(n)
  integer :: i

  write(unit,'(a,a,a)') 'SCALARS ', trim(name), ' double 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i=1,n
     write(unit,'(es24.16)') vals(i)
  enddo
end subroutine write_scalar_real

end module sgrac_vtk_slip_io
