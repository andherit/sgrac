program sgrac_pdf
  use, intrinsic :: iso_fortran_env, only: error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use generic
  use forparse
  use LAT_mesh
  use LAT_distance
  use lists
  use distance
  use sgrac_mesh_field_support
  use sgrac_vtk_field_io
  implicit none

  character(len=256) :: infile, outfile
  character(len=64) :: mode
  character(len=line_len), allocatable :: lines(:)
  type(mesh) :: amesh
  type(diff) :: adiff
  type(containerc), allocatable :: ntoc(:)
  type(containern), allocatable :: nton(:)
  integer :: ierr_mode, ierr_center, ierr_sigma, ierr_amplitude
  integer :: nlines, npoints, ncells, center_node, i
  real :: sigma_in, amplitude_in
  real(pr) :: sigma, amplitude, r, integral, check_integral, total_area, tol
  real(pr), allocatable :: px(:), py(:), pz(:), dg(:), dg_cell(:)
  real(pr), allocatable :: gaussian(:), pdf_weight(:), previous_weight(:), areas(:), pdf(:)
  integer(pin), allocatable :: cell(:,:)
  logical :: found_pdf_weight
  character(len=64) :: remove_fields(2)

  infile = '-'
  outfile = '-'
  mode = ''
  center_node = -1
  sigma_in = -1.0
  amplitude_in = 1.0

  ierr_mode = parse_arg('mode', mode)
  ierr_center = parse_arg('center_node', center_node)
  ierr_sigma = parse_arg('sigma', sigma_in)
  ierr_amplitude = parse_arg('amplitude', amplitude_in)
  i = parse_arg('in', infile)
  i = parse_arg('out', outfile)

  if (ierr_mode /= PARSE_OK) then
     write(error_unit,'(a)') 'sgrac-pdf: provide mode=<addgaussian|normalize>'
     stop 1
  endif

  select case(trim(mode))
  case('addgaussian')
     call validate_addgaussian_args(ierr_center, ierr_sigma, ierr_amplitude, sigma_in, amplitude_in)
  case('normalize')
     call validate_normalize_args(ierr_center, ierr_sigma, ierr_amplitude)
  case default
     write(error_unit,'(a,a)') 'sgrac-pdf: unsupported mode: ', trim(mode)
     stop 1
  end select

  call read_vtk_mesh_and_lines(trim(infile), lines, nlines, px, py, pz, cell, npoints, ncells)

  select case(trim(mode))
  case('addgaussian')
     call run_addgaussian()
  case('normalize')
     call run_normalize()
  end select

contains

subroutine validate_addgaussian_args(ierr_center, ierr_sigma, ierr_amplitude, sigma_in, amplitude_in)
  integer, intent(in) :: ierr_center, ierr_sigma, ierr_amplitude
  real, intent(in) :: sigma_in, amplitude_in

  if (ierr_center == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-pdf: invalid center_node=<integer>'
     stop 1
  endif
  if (ierr_sigma == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-pdf: invalid sigma=<meters>'
     stop 1
  endif
  if (ierr_amplitude == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-pdf: invalid amplitude=<real>'
     stop 1
  endif
  if (ierr_center /= PARSE_OK) then
     write(error_unit,'(a)') 'sgrac-pdf: provide center_node=<node>'
     stop 1
  endif
  if (ierr_sigma /= PARSE_OK) then
     write(error_unit,'(a)') 'sgrac-pdf: provide sigma=<meters>'
     stop 1
  endif

  sigma = real(sigma_in, pr)
  amplitude = real(amplitude_in, pr)
  if (sigma <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-pdf: sigma must be positive'
     stop 1
  endif
  if (amplitude < 0._pr) then
     write(error_unit,'(a)') 'sgrac-pdf: amplitude must be non-negative'
     stop 1
  endif
end subroutine validate_addgaussian_args

subroutine validate_normalize_args(ierr_center, ierr_sigma, ierr_amplitude)
  integer, intent(in) :: ierr_center, ierr_sigma, ierr_amplitude

  if (ierr_center == PARSE_OK .or. ierr_sigma == PARSE_OK .or. ierr_amplitude == PARSE_OK) then
     write(error_unit,'(a)') 'sgrac-pdf: mode=normalize does not accept center_node, sigma, or amplitude'
     stop 1
  endif
  if (ierr_center == PARSE_TYPE_ERROR .or. ierr_sigma == PARSE_TYPE_ERROR .or. &
      ierr_amplitude == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-pdf: invalid normalize arguments'
     stop 1
  endif
end subroutine validate_normalize_args

subroutine run_addgaussian()
  if (center_node < 1 .or. center_node > npoints) then
     write(error_unit,'(a,i0,a,i0,a)') 'sgrac-pdf: center_node ', center_node, &
        ' outside valid range 1..', npoints, ''
     stop 1
  endif

  call read_cell_scalar(lines, nlines, 'pdf_weight', ncells, found_pdf_weight, previous_weight)

  amesh%Nnodes = npoints
  amesh%Ncells = ncells
  allocate(amesh%px(npoints), amesh%py(npoints), amesh%pz(npoints), amesh%cell(ncells,3))
  amesh%px = px
  amesh%py = py
  amesh%pz = pz
  amesh%cell = cell

  allocate(dg(npoints), dg_cell(ncells), gaussian(ncells), pdf_weight(ncells))
  allocate(ntoc(npoints), nton(npoints))

  adiff%fast = .true.
  verbose = 0
  call pre_onevsall2d_onvertex(amesh, int(center_node, pin), dg, ntoc, nton)
  if (.not. associated(ntoc(center_node)%ptr)) then
     write(error_unit,'(a,i0,a)') 'sgrac-pdf: center_node ', center_node, ' is not used by any triangle'
     stop 1
  endif
  call onevsall2d(amesh, dg, ntoc, nton, adiff)

  do i=1,ncells
     dg_cell(i) = (dg(cell(i,1)) + dg(cell(i,2)) + dg(cell(i,3))) / 3._pr
     r = dg_cell(i)
     gaussian(i) = exp(-0.5_pr * (r / sigma)**2)
  enddo

  pdf_weight = previous_weight + amplitude * gaussian

  write(error_unit,'(a)') 'sgrac-pdf diagnostics:'
  write(error_unit,'(a)') '  mode = addgaussian'
  write(error_unit,'(a,i0)') '  center_node = ', center_node
  write(error_unit,'(a,es24.16)') '  sigma = ', sigma
  write(error_unit,'(a,es24.16)') '  amplitude = ', amplitude
  write(error_unit,'(a,l1)') '  previous pdf_weight found = ', found_pdf_weight
  write(error_unit,'(a,i0)') '  cells = ', ncells
  write(error_unit,'(a,es24.16)') '  gaussian maximum = ', maxval(gaussian)
  write(error_unit,'(a,es24.16)') '  previous weight maximum = ', maxval(previous_weight)
  write(error_unit,'(a,es24.16)') '  final weight maximum = ', maxval(pdf_weight)
  write(error_unit,'(a,es24.16)') '  final weight sum = ', sum(pdf_weight)

  call write_or_replace_cell_scalar(trim(outfile), lines, nlines, 'pdf_weight', pdf_weight, ncells)

  call free_nton(nton)
  call free_ntoc(ntoc)
end subroutine run_addgaussian

subroutine run_normalize()
  call read_cell_scalar(lines, nlines, 'pdf_weight', ncells, found_pdf_weight, pdf_weight)
  if (.not. found_pdf_weight) then
     write(error_unit,'(a)') 'sgrac-pdf: mode=normalize requires CELL_DATA scalar pdf_weight'
     stop 1
  endif

  allocate(areas(ncells), pdf(ncells))
  call compute_triangle_areas(px, py, pz, cell, areas)

  do i=1,ncells
     if (.not. ieee_is_finite(pdf_weight(i))) then
        write(error_unit,'(a,i0)') 'sgrac-pdf: non-finite pdf_weight at cell ', i
        stop 1
     endif
     if (pdf_weight(i) < 0._pr) then
        write(error_unit,'(a,i0)') 'sgrac-pdf: negative pdf_weight at cell ', i
        stop 1
     endif
     if (.not. ieee_is_finite(areas(i))) then
        write(error_unit,'(a,i0)') 'sgrac-pdf: non-finite triangle area at cell ', i
        stop 1
     endif
     if (areas(i) <= 0._pr) then
        write(error_unit,'(a,i0)') 'sgrac-pdf: non-positive triangle area at cell ', i
        stop 1
     endif
  enddo

  integral = sum(pdf_weight * areas)
  total_area = sum(areas)
  if (.not. ieee_is_finite(integral) .or. integral <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-pdf: pdf_weight area integral must be positive'
     stop 1
  endif

  pdf = pdf_weight / integral
  check_integral = sum(pdf * areas)
  tol = 100._pr * epsilon(1._pr) * max(1._pr, abs(check_integral))
  if (abs(check_integral - 1._pr) > tol) then
     write(error_unit,'(a,es24.16)') 'sgrac-pdf: normalized area integral check failed: ', check_integral
     stop 1
  endif

  write(error_unit,'(a)') 'sgrac-pdf diagnostics:'
  write(error_unit,'(a)') '  mode = normalize'
  write(error_unit,'(a,i0)') '  cells = ', ncells
  write(error_unit,'(a,es24.16)') '  minimum pdf_weight = ', minval(pdf_weight)
  write(error_unit,'(a,es24.16)') '  maximum pdf_weight = ', maxval(pdf_weight)
  write(error_unit,'(a,es24.16)') '  weight-area integral = ', integral
  write(error_unit,'(a,es24.16)') '  total rupture area = ', total_area
  write(error_unit,'(a,es24.16)') '  minimum pdf = ', minval(pdf)
  write(error_unit,'(a,es24.16)') '  maximum pdf = ', maxval(pdf)
  write(error_unit,'(a,es24.16)') '  normalized area integral = ', check_integral
  write(error_unit,'(a)') '  removed field = pdf_weight'
  write(error_unit,'(a)') '  removed field = area'
  write(error_unit,'(a)') '  output field = pdf'

  remove_fields(1) = 'pdf_weight'
  remove_fields(2) = 'area'
  call write_replace_and_remove_cell_scalars(trim(outfile), lines, nlines, 'pdf', pdf, ncells, remove_fields, 2)
end subroutine run_normalize

end program sgrac_pdf
