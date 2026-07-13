module sgrac_vtk_field_io
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  implicit none
  integer, parameter :: line_len = 512
  private
  public :: line_len
  public :: read_vtk_mesh_and_lines
  public :: read_cell_scalar
  public :: remove_cell_scalar
  public :: write_or_replace_cell_scalar
  public :: write_replace_and_remove_cell_scalars
contains

subroutine read_vtk_mesh_and_lines(filename, lines, nlines, px, py, pz, cell, npoints, ncells)
  character(*), intent(in) :: filename
  character(len=line_len), allocatable, intent(out) :: lines(:)
  integer, intent(out) :: nlines, npoints, ncells
  real(pr), allocatable, intent(out) :: px(:), py(:), pz(:)
  integer(pin), allocatable, intent(out) :: cell(:,:)

  call read_text_file(filename, lines, nlines)
  call require_polydata(lines, nlines)
  call read_points(lines, nlines, px, py, pz, npoints)
  call read_triangles(lines, nlines, cell, ncells)
end subroutine read_vtk_mesh_and_lines

subroutine read_cell_scalar(lines, nlines, fieldname, ncells, found, vals)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncells
  character(*), intent(in) :: fieldname
  logical, intent(out) :: found
  real(pr), allocatable, intent(out) :: vals(:)
  integer :: cell_data_line, cell_data_end, cell_count, scalar_line, scalar_end
  integer :: i, lookup_line, ios

  allocate(vals(ncells))
  vals = 0._pr
  found = .false.

  call find_cell_data_section(lines, nlines, cell_data_line, cell_data_end, cell_count)
  if (cell_data_line < 0) return
  call validate_cell_data_count(cell_count, ncells)

  call find_cell_scalar_block(lines, nlines, fieldname, ncells, cell_data_line, cell_data_end, &
                              scalar_line, scalar_end)
  if (scalar_line < 0) return

  lookup_line = scalar_line + 1
  if (lookup_line > nlines) stop 'sgrac-vtk-field-io: malformed CELL_DATA scalar'
  do i=1,ncells
     if (lookup_line + i > nlines) stop 'sgrac-vtk-field-io: CELL_DATA scalar is truncated'
     read(lines(lookup_line+i),*,iostat=ios) vals(i)
     if (ios /= 0) stop 'sgrac-vtk-field-io: error reading CELL_DATA scalar values'
  enddo
  found = .true.
end subroutine read_cell_scalar

subroutine write_or_replace_cell_scalar(filename, lines, nlines, fieldname, vals, ncells)
  character(*), intent(in) :: filename, fieldname
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncells
  real(pr), intent(in) :: vals(ncells)
  integer :: unit, ios, i
  integer :: cell_data_line, cell_data_end, cell_count, point_data_line
  integer :: scalar_line, scalar_end, insert_line

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-vtk-field-io: cannot open output file'
  endif

  call find_cell_data_section(lines, nlines, cell_data_line, cell_data_end, cell_count)
  if (cell_data_line > 0) call validate_cell_data_count(cell_count, ncells)

  if (cell_data_line > 0) then
     call find_cell_scalar_block(lines, nlines, fieldname, ncells, cell_data_line, cell_data_end, &
                                 scalar_line, scalar_end)
  else
     scalar_line = -1
     scalar_end = -1
  endif

  point_data_line = find_point_data_line(lines, nlines)

  if (cell_data_line < 0) then
     if (point_data_line > 0) then
        insert_line = point_data_line
     else
        insert_line = nlines + 1
     endif
     do i=1,insert_line-1
        write(unit,'(a)') trim(lines(i))
     enddo
     write(unit,'(a,i0)') 'CELL_DATA ', ncells
     call write_scalar_real(unit, fieldname, vals, ncells)
     do i=insert_line,nlines
        write(unit,'(a)') trim(lines(i))
     enddo
  else if (scalar_line > 0) then
     do i=1,scalar_line-1
        write(unit,'(a)') trim(lines(i))
     enddo
     call write_scalar_real(unit, fieldname, vals, ncells)
     do i=scalar_end+1,nlines
        write(unit,'(a)') trim(lines(i))
     enddo
  else
     insert_line = cell_data_end + 1
     do i=1,insert_line-1
        write(unit,'(a)') trim(lines(i))
     enddo
     call write_scalar_real(unit, fieldname, vals, ncells)
     do i=insert_line,nlines
        write(unit,'(a)') trim(lines(i))
     enddo
  endif

  if (unit /= 6) close(unit)
end subroutine write_or_replace_cell_scalar

subroutine remove_cell_scalar(filename, lines, nlines, fieldname, ncells)
  character(*), intent(in) :: filename, fieldname
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncells
  character(len=64) :: remove_fields(1)
  real(pr) :: dummy_vals(1)

  remove_fields(1) = fieldname
  dummy_vals = 0._pr
  call write_filtered_cell_data(filename, lines, nlines, '', dummy_vals, ncells, &
                                remove_fields, 1, .false.)
end subroutine remove_cell_scalar

subroutine write_replace_and_remove_cell_scalars(filename, lines, nlines, fieldname, vals, ncells, remove_fields, nremove)
  character(*), intent(in) :: filename, fieldname
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncells, nremove
  real(pr), intent(in) :: vals(ncells)
  character(len=*), intent(in) :: remove_fields(nremove)

  call write_filtered_cell_data(filename, lines, nlines, fieldname, vals, ncells, remove_fields, nremove, .true.)
end subroutine write_replace_and_remove_cell_scalars

subroutine write_filtered_cell_data(filename, lines, nlines, fieldname, vals, ncells, remove_fields, nremove, write_new_scalar)
  character(*), intent(in) :: filename, fieldname
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncells, nremove
  real(pr), intent(in) :: vals(:)
  character(len=*), intent(in) :: remove_fields(nremove)
  logical, intent(in) :: write_new_scalar
  integer :: unit, ios, i, next_i
  integer :: cell_data_line, cell_data_end, cell_count, point_data_line, insert_line
  character(len=64) :: key

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-vtk-field-io: cannot open output file'
  endif

  call find_cell_data_section(lines, nlines, cell_data_line, cell_data_end, cell_count)
  if (cell_data_line > 0) call validate_cell_data_count(cell_count, ncells)
  point_data_line = find_point_data_line(lines, nlines)

  if (cell_data_line < 0) then
     if (point_data_line > 0) then
        insert_line = point_data_line
     else
        insert_line = nlines + 1
     endif
     do i=1,insert_line-1
        write(unit,'(a)') trim(lines(i))
     enddo
     if (write_new_scalar) then
        write(unit,'(a,i0)') 'CELL_DATA ', ncells
        call write_scalar_real(unit, fieldname, vals, ncells)
     endif
     do i=insert_line,nlines
        write(unit,'(a)') trim(lines(i))
     enddo
  else
     insert_line = cell_data_end + 1
     i = 1
     do while (i <= nlines)
        if (i == insert_line .and. write_new_scalar) call write_scalar_real(unit, fieldname, vals, ncells)

        if (i > cell_data_line .and. i <= cell_data_end) then
           key = ''
           read(lines(i),*,iostat=ios) key
           if (ios == 0 .and. trim(key) == 'SCALARS' .and. &
               should_remove_scalar(lines(i), fieldname, remove_fields, nremove, write_new_scalar)) then
              next_i = i + 2 + ncells
              if (next_i - 1 > cell_data_end) stop 'sgrac-vtk-field-io: existing CELL_DATA scalar is truncated'
              i = next_i
              cycle
           endif
        endif

        write(unit,'(a)') trim(lines(i))
        i = i + 1
     enddo
     if (insert_line == nlines + 1 .and. write_new_scalar) call write_scalar_real(unit, fieldname, vals, ncells)
  endif

  if (unit /= 6) close(unit)
end subroutine write_filtered_cell_data

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
     if (ios /= 0) stop 'sgrac-vtk-field-io: cannot open input file'
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

subroutine require_polydata(lines, nlines)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer :: i, ios
  character(len=64) :: key, dtype

  do i=1,nlines
     key = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, dtype
     if (ios == 0 .and. trim(key) == 'DATASET') then
        if (trim(dtype) /= 'POLYDATA') stop 'sgrac-vtk-field-io: DATASET must be POLYDATA'
        return
     endif
  enddo
  stop 'sgrac-vtk-field-io: DATASET POLYDATA section not found'
end subroutine require_polydata

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
  if (point_line < 0) stop 'sgrac-vtk-field-io: POINTS section not found'
  if (npoints <= 0) stop 'sgrac-vtk-field-io: POINTS count must be positive'
  if (point_line + npoints > nlines) stop 'sgrac-vtk-field-io: POINTS section is truncated'

  allocate(px(npoints), py(npoints), pz(npoints))
  do j=1,npoints
     read(lines(point_line+j),*,iostat=ios) px(j), py(j), pz(j)
     if (ios /= 0) stop 'sgrac-vtk-field-io: error while reading POINTS'
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
        stop 'sgrac-vtk-field-io: CELLS input is not supported; use POLYDATA POLYGONS'
     endif
  enddo
  if (poly_line < 0) stop 'sgrac-vtk-field-io: POLYGONS section not found'
  if (ncells <= 0) stop 'sgrac-vtk-field-io: POLYGONS count must be positive'
  if (total_size /= 4 * ncells) stop 'sgrac-vtk-field-io: POLYGONS size must match triangular connectivity'
  if (poly_line + ncells > nlines) stop 'sgrac-vtk-field-io: POLYGONS section is truncated'

  allocate(cell(ncells,3))
  do j=1,ncells
     read(lines(poly_line+j),*,iostat=ios) nverts, cell(j,1), cell(j,2), cell(j,3)
     if (ios /= 0) stop 'sgrac-vtk-field-io: error while reading POLYGONS'
     if (nverts /= 3) stop 'sgrac-vtk-field-io: only triangular POLYGONS are supported'
     cell(j,1:3) = cell(j,1:3) + 1_pin
  enddo
end subroutine read_triangles

subroutine find_cell_data_section(lines, nlines, cell_data_line, cell_data_end, cell_count)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer, intent(out) :: cell_data_line, cell_data_end, cell_count
  integer :: i, ios
  character(len=64) :: key

  cell_data_line = -1
  cell_data_end = -1
  cell_count = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key, cell_count
     if (ios == 0 .and. trim(key) == 'CELL_DATA') then
        cell_data_line = i
        exit
     endif
  enddo
  if (cell_data_line < 0) return

  cell_data_end = nlines
  do i=cell_data_line+1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'POINT_DATA') then
        cell_data_end = i - 1
        return
     endif
  enddo
end subroutine find_cell_data_section

subroutine validate_cell_data_count(cell_count, ncells)
  integer, intent(in) :: cell_count, ncells

  if (cell_count /= ncells) then
     write(error_unit,'(a,i0,a,i0)') 'sgrac-vtk-field-io: CELL_DATA count ', cell_count, &
        ' does not match POLYGONS count ', ncells
     stop 1
  endif
end subroutine validate_cell_data_count

subroutine find_cell_scalar_block(lines, nlines, fieldname, ncells, cell_data_line, cell_data_end, &
                                  scalar_line, scalar_end)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncells, cell_data_line, cell_data_end
  character(*), intent(in) :: fieldname
  integer, intent(out) :: scalar_line, scalar_end
  integer :: i, ios
  character(len=64) :: key

  scalar_line = -1
  scalar_end = -1
  i = cell_data_line + 1
  do while (i <= cell_data_end)
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'SCALARS') then
        if (is_requested_scalar(lines(i), fieldname)) then
           if (i + 1 + ncells > cell_data_end) stop 'sgrac-vtk-field-io: existing CELL_DATA scalar is truncated'
           scalar_line = i
           scalar_end = i + 1 + ncells
           return
        endif
        i = i + 2 + ncells
     else if (ios == 0 .and. trim(key) == 'VECTORS') then
        i = i + 1 + ncells
     else if (ios == 0 .and. trim(key) == 'FIELD') then
        i = skip_field_data(lines, nlines, i, cell_data_end)
     else
        i = i + 1
     endif
  enddo
end subroutine find_cell_scalar_block

integer function skip_field_data(lines, nlines, field_line, section_end) result(next_line)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, field_line, section_end
  integer :: i, ifield, nfields, ncomp, ntuple, ios
  character(len=64) :: key, field_name, dtype

  read(lines(field_line),*,iostat=ios) key, field_name, nfields
  if (ios /= 0) stop 'sgrac-vtk-field-io: malformed FIELD data header'
  i = field_line + 1
  do ifield=1,nfields
     if (i > section_end) stop 'sgrac-vtk-field-io: FIELD data is truncated'
     read(lines(i),*,iostat=ios) field_name, ncomp, ntuple, dtype
     if (ios /= 0) stop 'sgrac-vtk-field-io: malformed FIELD data array'
     i = i + 1 + ntuple
     if (i - 1 > nlines) stop 'sgrac-vtk-field-io: FIELD data is truncated'
  enddo
  next_line = i
end function skip_field_data

integer function find_point_data_line(lines, nlines) result(point_data_line)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer :: i, ios
  character(len=64) :: key

  point_data_line = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'POINT_DATA') then
        point_data_line = i
        return
     endif
  enddo
end function find_point_data_line

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

logical function should_remove_scalar(line, fieldname, remove_fields, nremove, remove_fieldname)
  character(len=line_len), intent(in) :: line
  character(*), intent(in) :: fieldname
  character(len=*), intent(in) :: remove_fields(nremove)
  integer, intent(in) :: nremove
  logical, intent(in) :: remove_fieldname
  integer :: i

  should_remove_scalar = .false.
  if (remove_fieldname .and. is_requested_scalar(line, fieldname)) then
     should_remove_scalar = .true.
     return
  endif
  do i=1,nremove
     if (is_requested_scalar(line, trim(remove_fields(i)))) then
        should_remove_scalar = .true.
        return
     endif
  enddo
end function should_remove_scalar

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

end module sgrac_vtk_field_io
