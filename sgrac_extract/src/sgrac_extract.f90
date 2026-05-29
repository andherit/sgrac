program sgrac_extract
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  use forparse
  use sgrac_vtk_extract_io
  implicit none

  character(len=256) :: infile, outfile, field
  integer :: ierr, value, npoints, ncells, nchild_points, nchild_cells
  integer :: i, j, icell
  real(pr), allocatable :: px(:), py(:), pz(:), child_px(:), child_py(:), child_pz(:)
  integer(pin), allocatable :: cell(:,:), child_cell(:,:), field_values(:), old2new(:)
  logical, allocatable :: selected(:), used(:)

  infile = '-'
  outfile = '-'
  field = 'mask'
  value = 1

  ierr = parse_arg('in', infile)
  ierr = parse_arg('out', outfile)
  ierr = parse_arg('field', field)
  ierr = parse_arg('value', value)
  if (ierr == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-extract: invalid value=<integer>'
     stop 1
  endif

  call read_parent_vtk(trim(infile), trim(field), px, py, pz, cell, field_values, npoints, ncells)

  allocate(selected(ncells), used(npoints), old2new(npoints))
  selected = field_values == int(value, pin)
  nchild_cells = count(selected)
  if (nchild_cells == 0) then
     write(error_unit,'(a,a,a,i0)') 'sgrac-extract: no cells selected for field ', trim(field), ' value ', value
     stop 1
  endif

  used = .false.
  do i=1,ncells
     if (.not. selected(i)) cycle
     do j=1,3
        if (cell(i,j) < 1 .or. cell(i,j) > npoints) then
           write(error_unit,'(a,i0)') 'sgrac-extract: cell references invalid parent node ', cell(i,j)
           stop 1
        endif
        used(cell(i,j)) = .true.
     enddo
  enddo

  old2new = 0_pin
  nchild_points = 0
  if (used(1)) then
     nchild_points = 1
     old2new(1) = 1_pin
     do i=2,npoints
        if (used(i)) then
           nchild_points = nchild_points + 1
           old2new(i) = int(nchild_points, pin)
        endif
     enddo
  else
     do i=1,npoints
        if (used(i)) then
           nchild_points = nchild_points + 1
           old2new(i) = int(nchild_points, pin)
        endif
     enddo
  endif

  allocate(child_px(nchild_points), child_py(nchild_points), child_pz(nchild_points))
  do i=1,npoints
     if (old2new(i) > 0) then
        child_px(old2new(i)) = px(i)
        child_py(old2new(i)) = py(i)
        child_pz(old2new(i)) = pz(i)
     endif
  enddo

  allocate(child_cell(nchild_cells,3))
  icell = 0
  do i=1,ncells
     if (.not. selected(i)) cycle
     icell = icell + 1
     do j=1,3
        child_cell(icell,j) = old2new(cell(i,j))
     enddo
  enddo

  call write_child_vtk(trim(outfile), child_px, child_py, child_pz, child_cell)
end program sgrac_extract
