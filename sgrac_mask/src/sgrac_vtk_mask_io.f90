module sgrac_vtk_mask_io
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  implicit none
  integer, parameter :: line_len = 512
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

subroutine read_sgrac_metadata(lines, nlines, metadata)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  type(sgrac_vtk_metadata_type), intent(out) :: metadata
  integer :: i, ios, nfields, ifield, ncomp, ntuple, skip_lines
  character(len=64) :: key, name, dtype, field_name

  metadata = sgrac_vtk_metadata_type()

  do i=1,nlines
     key = ''
     field_name = ''
     read(lines(i),*,iostat=ios) key, field_name, nfields
     if (ios == 0 .and. trim(key) == 'FIELD' .and. trim(field_name) == 'sgrac_metadata') exit
  enddo

  if (i > nlines) return
  metadata%present = .true.

  i = i + 1
  do ifield=1,nfields
     if (i > nlines) exit
     name = ''
     dtype = ''
     read(lines(i),*,iostat=ios) name, ncomp, ntuple, dtype
     if (ios /= 0) exit
     if (i + 1 > nlines) exit

     select case(trim(name))
     case('sgrac_format_version')
        read(lines(i+1),*,iostat=ios) metadata%format_version
     case('mask_geometry_type_code')
        read(lines(i+1),*,iostat=ios) metadata%mask_geometry_type_code
     case('distance_source_count')
        read(lines(i+1),*,iostat=ios) metadata%distance_source_count
     case('source_node_indexes')
        read(lines(i+1),*,iostat=ios) metadata%source_node_indexes(1:min(2,ncomp*ntuple))
     case('theta_convention_code')
        read(lines(i+1),*,iostat=ios) metadata%theta_convention_code
     case('distance_units_code')
        read(lines(i+1),*,iostat=ios) metadata%distance_units_code
     end select

     skip_lines = max(1, ntuple)
     i = i + 1 + skip_lines
  enddo
end subroutine read_sgrac_metadata

subroutine validate_sgrac_geometry_fields(lines, nlines, npoints, ncell, metadata, ok)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, npoints, ncell
  type(sgrac_vtk_metadata_type), intent(in) :: metadata
  logical, intent(out) :: ok

  ok = .true.
  if (.not. metadata%present) return
  call validate_sgrac_geometry_mode(lines, nlines, npoints, ncell, metadata%mask_geometry_type_code, ok)
end subroutine validate_sgrac_geometry_fields

subroutine validate_sgrac_geometry_mode(lines, nlines, npoints, ncell, geometry_code, ok)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, npoints, ncell, geometry_code
  logical, intent(out) :: ok

  ok = .true.

  select case(geometry_code)
  case(0)
     return
  case(1)
     call require_cell_scalar(lines, nlines, ncell, FIELD_DG_CELL, ok)
  case(2)
     call require_cell_scalar(lines, nlines, ncell, FIELD_DG_CELL, ok)
     call require_cell_scalar(lines, nlines, ncell, FIELD_THETA, ok)
  case(3)
     if (.not. has_point_scalar(lines, nlines, npoints, FIELD_DG_F1)) then
        write(error_unit,'(a,a)') 'sgrac-mask: warning: missing POINT_DATA scalar ', FIELD_DG_F1
     endif
     if (.not. has_point_scalar(lines, nlines, npoints, FIELD_DG_F2)) then
        write(error_unit,'(a,a)') 'sgrac-mask: warning: missing POINT_DATA scalar ', FIELD_DG_F2
     endif
     if (.not. has_cell_scalar(lines, nlines, ncell, FIELD_DG_F1_CELL)) then
        write(error_unit,'(a,a)') 'sgrac-mask: warning: missing CELL_DATA scalar ', FIELD_DG_F1_CELL
     endif
     if (.not. has_cell_scalar(lines, nlines, ncell, FIELD_DG_F2_CELL)) then
        write(error_unit,'(a,a)') 'sgrac-mask: warning: missing CELL_DATA scalar ', FIELD_DG_F2_CELL
     endif
     call require_cell_scalar(lines, nlines, ncell, FIELD_DG_SUM_CELL, ok)
  case default
     write(error_unit,'(a,i0)') 'sgrac-mask: warning: unknown mask_geometry_type_code = ', geometry_code
  end select
end subroutine validate_sgrac_geometry_mode

subroutine require_cell_scalar(lines, nlines, ncell, fieldname, ok)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell
  character(*), intent(in) :: fieldname
  logical, intent(inout) :: ok

  if (.not. has_cell_scalar(lines, nlines, ncell, fieldname)) then
     write(error_unit,'(a,a)') 'sgrac-mask: warning: missing CELL_DATA scalar ', trim(fieldname)
     ok = .false.
  endif
end subroutine require_cell_scalar

logical function has_point_scalar(lines, nlines, npoints, fieldname) result(found)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, npoints
  character(*), intent(in) :: fieldname
  integer :: i, ios, point_data_line
  character(len=64) :: key, name, dtype

  found = .false.
  point_data_line = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'POINT_DATA') then
        point_data_line = i
        exit
     endif
  enddo
  if (point_data_line < 0) return

  i = point_data_line + 1
  do while (i <= nlines)
     key = ''
     name = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, name, dtype
     if (ios == 0 .and. trim(key) == 'CELL_DATA') return
     if (ios == 0 .and. trim(key) == 'SCALARS' .and. trim(name) == trim(fieldname)) then
        found = .true.
        return
     endif

     if (ios == 0 .and. trim(key) == 'SCALARS') then
        i = i + 2 + npoints
     else if (ios == 0 .and. trim(key) == 'VECTORS') then
        i = i + 1 + npoints
     else
        i = i + 1
     endif
  enddo
end function has_point_scalar

logical function has_cell_scalar(lines, nlines, ncell, fieldname) result(found)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell
  character(*), intent(in) :: fieldname
  integer :: i, j, ios, cell_data_line
  character(len=64) :: key, name, dtype

  found = .false.
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

  i = cell_data_line + 1
  do while (i <= nlines)
     key = ''
     name = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, name, dtype
     if (ios == 0 .and. trim(key) == 'SCALARS' .and. trim(name) == trim(fieldname)) then
        found = .true.
        return
     endif

     if (ios == 0 .and. trim(key) == 'SCALARS') then
        i = i + 2 + ncell
     else if (ios == 0 .and. trim(key) == 'VECTORS') then
        i = i + 1 + ncell
     else if (ios == 0 .and. trim(key) == 'FIELD') then
        read(lines(i),*,iostat=ios) key, name, j
        i = i + 1
     else
        i = i + 1
     endif
  enddo
end function has_cell_scalar

subroutine read_point_scalar(lines, nlines, npoints, fieldname, vals)
  character(len=line_len), intent(in) :: lines(:)
  integer, intent(in) :: nlines, npoints
  character(*), intent(in) :: fieldname
  real(pr), intent(out) :: vals(npoints)
  integer :: i, j, ios, point_data_line, lookup_line
  character(len=64) :: key, name, dtype

  point_data_line = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key
     if (ios == 0 .and. trim(key) == 'POINT_DATA') then
        point_data_line = i
        exit
     endif
  enddo
  if (point_data_line < 0) stop 'sgrac-mask: POINT_DATA section not found'

  i = point_data_line + 1
  do while (i <= nlines)
     key = ''
     name = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, name, dtype
     if (ios == 0 .and. trim(key) == 'CELL_DATA') exit
     if (ios == 0 .and. trim(key) == 'SCALARS' .and. trim(name) == trim(fieldname)) then
        lookup_line = i + 1
        if (lookup_line > nlines) stop 'sgrac-mask: malformed point scalar field'
        do j=1,npoints
           if (lookup_line + j > nlines) stop 'sgrac-mask: point scalar field is truncated'
           read(lines(lookup_line+j),*,iostat=ios) vals(j)
           if (ios /= 0) stop 'sgrac-mask: error reading point scalar values'
        enddo
        return
     endif

     if (ios == 0 .and. trim(key) == 'SCALARS') then
        i = i + 2 + npoints
     else if (ios == 0 .and. trim(key) == 'VECTORS') then
        i = i + 1 + npoints
     else
        i = i + 1
     endif
  enddo

  write(error_unit,*) 'sgrac-mask: missing required POINT_DATA scalar ', trim(fieldname)
  stop 1
end subroutine read_point_scalar

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
