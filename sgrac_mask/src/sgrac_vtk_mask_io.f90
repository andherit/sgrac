module sgrac_vtk_mask_io
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  implicit none
  integer, parameter :: line_len = 512
contains

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
     if (ios /= 0) stop 'sgrac-mask: cannot open input file'
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

subroutine write_text_file_with_mask(filename, lines, nlines, rtheta, phi, mask, ncell)
  character(*), intent(in) :: filename
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell
  real(pr), intent(in) :: rtheta(ncell), phi(ncell)
  integer, intent(in) :: mask(ncell)
  integer :: unit, ios, i, last_line

  if (len_trim(filename) == 0 .or. trim(filename) == '-') then
     unit = 6
  else
     open(newunit=unit, file=trim(filename), status='replace', action='write', iostat=ios)
     if (ios /= 0) stop 'sgrac-mask: cannot open output file'
  endif

  last_line = last_line_before_mask_fields(lines, nlines)
  do i=1,last_line
     write(unit,'(a)') trim(lines(i))
  enddo

  call write_scalar_real(unit, 'Rtheta', rtheta, ncell)
  call write_scalar_real(unit, 'phi', phi, ncell)
  call write_scalar_int(unit, 'mask', mask, ncell)

  if (unit /= 6) close(unit)
end subroutine write_text_file_with_mask

subroutine get_polydata_counts(lines, nlines, npoints, ncells)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer, intent(out) :: npoints, ncells
  integer :: i, ios, total_size
  character(len=64) :: key, dtype

  npoints = -1
  ncells = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios /= 0) cycle

     select case(trim(key))
     case('POINTS')
        dtype = ''
        read(lines(i),*,iostat=ios) key, npoints, dtype
        if (ios /= 0) stop 'sgrac-mask: malformed POINTS header'
     case('POLYGONS','CELLS')
        read(lines(i),*,iostat=ios) key, ncells, total_size
        if (ios /= 0) stop 'sgrac-mask: malformed POLYGONS/CELLS header'
     end select

     if (npoints >= 0 .and. ncells >= 0) exit
  enddo
  if (npoints < 0) stop 'sgrac-mask: POINTS section not found'
  if (ncells < 0) stop 'sgrac-mask: POLYGONS/CELLS section not found'
end subroutine get_polydata_counts

subroutine read_cell_scalar(lines, nlines, ncell, fieldname, vals)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell
  character(*), intent(in) :: fieldname
  real(pr), intent(out) :: vals(ncell)
  integer :: i, j, ios, cell_data_line, lookup_line
  character(len=64) :: key, name, dtype

  cell_data_line = -1
  do i=1,nlines
     read(lines(i),*,iostat=ios) key, j
     if (ios == 0 .and. trim(key) == 'CELL_DATA') then
        cell_data_line = i
        exit
     endif
  enddo
  if (cell_data_line < 0) stop 'sgrac-mask: CELL_DATA section not found'

  i = cell_data_line + 1
  do while (i <= nlines)
     key = ''
     name = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, name, dtype
     if (ios == 0 .and. trim(key) == 'SCALARS' .and. trim(name) == trim(fieldname)) then
        lookup_line = i + 1
        if (lookup_line > nlines) stop 'sgrac-mask: malformed scalar field'
        do j=1,ncell
           if (lookup_line + j > nlines) stop 'sgrac-mask: scalar field is truncated'
           read(lines(lookup_line+j),*,iostat=ios) vals(j)
           if (ios /= 0) stop 'sgrac-mask: error reading scalar values'
        enddo
        return
     endif

     if (ios == 0 .and. trim(key) == 'SCALARS') then
        i = i + 2 + ncell
     else if (ios == 0 .and. trim(key) == 'VECTORS') then
        i = i + 1 + ncell
     else if (ios == 0 .and. trim(key) == 'FIELD') then
        stop 'sgrac-mask: FIELD cell data are not supported in v0'
     else
        i = i + 1
     endif
  enddo

  write(error_unit,*) 'sgrac-mask: missing required CELL_DATA scalar ', trim(fieldname)
  stop 1
end subroutine read_cell_scalar

integer function last_line_before_mask_fields(lines, nlines) result(last_line)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer :: i, ios, cell_data_line
  character(len=64) :: key, name

  last_line = nlines
  cell_data_line = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'CELL_DATA') then
        cell_data_line = i
        exit
     endif
  enddo
  if (cell_data_line < 0) return

  do i=cell_data_line+1,nlines
     key = ''
     name = ''
     read(lines(i),*,iostat=ios) key, name
     if (ios /= 0) cycle
     if (trim(key) == 'SCALARS') then
        select case(trim(name))
        case('Rtheta','phi','mask')
           last_line = i - 1
           return
        end select
     endif
  enddo
end function last_line_before_mask_fields

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

subroutine write_scalar_int(unit, name, vals, n)
  integer, intent(in) :: unit, n
  character(*), intent(in) :: name
  integer, intent(in) :: vals(n)
  integer :: i
  write(unit,'(a,a,a)') 'SCALARS ', trim(name), ' int 1'
  write(unit,'(a)') 'LOOKUP_TABLE default'
  do i=1,n
     write(unit,'(i0)') vals(i)
  enddo
end subroutine write_scalar_int

end module sgrac_vtk_mask_io
