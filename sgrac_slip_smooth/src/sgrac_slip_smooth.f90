program sgrac_slip_smooth
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  use forparse
  use LAT_mesh
  use LAT_distance
  use lists
  use distance
  use sgrac_slip_support
  use sgrac_vtk_slip_io
  implicit none

  type(mesh) :: amesh
  type(diff) :: adiff
  type(containerc), allocatable :: ntoc(:)
  type(containern), allocatable :: nton(:)
  character(len=256) :: infile, outfile
  character(len=line_len), allocatable :: lines(:)
  integer :: ierr_center, ierr_sigma, ierr_peak, ierr_mw, ierr_mu, nlines, npoints, ncells
  integer :: scale_status
  integer :: center_node, i
  real :: sigma_in, peak_slip_in, mw_in, mu_in
  real(pr) :: sigma, peak_slip, mw, mu, target_m0, current_m0, final_m0, r
  real(pr), allocatable :: px(:), py(:), pz(:), dg(:), dg_cell(:), slip(:), slip_scaled(:), areas(:)
  integer(pin), allocatable :: cell(:,:)
  logical :: has_peak, has_mw, has_mu, has_physical_keyword, physical_mode

  infile = '-'
  outfile = '-'
  center_node = -1
  sigma_in = -1.0
  peak_slip_in = -1.0
  mw_in = 0.0
  mu_in = -1.0

  ierr_center = parse_arg('center_node', center_node)
  ierr_sigma = parse_arg('sigma', sigma_in)
  ierr_peak = parse_arg('peak_slip', peak_slip_in)
  ierr_mw = parse_arg('mw', mw_in)
  ierr_mu = parse_arg('mu', mu_in)
  i = parse_arg('in', infile)
  i = parse_arg('out', outfile)

  if (ierr_center == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: invalid center_node=<integer>'
     stop 1
  endif
  if (ierr_sigma == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: invalid sigma=<meters>'
     stop 1
  endif
  if (ierr_peak == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: invalid peak_slip=<meters>'
     stop 1
  endif
  if (ierr_center /= PARSE_OK .or. ierr_sigma /= PARSE_OK) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: provide center_node=<node> sigma=<m> and either mw=<Mw> or peak_slip=<m>'
     stop 1
  endif

  has_peak = ierr_peak == PARSE_OK
  has_mw = ierr_mw == PARSE_OK
  has_mu = ierr_mu == PARSE_OK
  has_physical_keyword = ierr_mw == PARSE_OK .or. ierr_mw == PARSE_TYPE_ERROR .or. &
                         ierr_mu == PARSE_OK .or. ierr_mu == PARSE_TYPE_ERROR
  physical_mode = .not. has_peak

  if (has_peak .and. has_physical_keyword) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: warning: peak_slip is present; mw and mu are ignored'
  endif

  if (physical_mode .and. ierr_mw == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: invalid mw=<moment magnitude>'
     stop 1
  endif
  if (physical_mode .and. ierr_mu == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: invalid mu=<Pa>'
     stop 1
  endif
  if (physical_mode .and. .not. has_mw) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: provide mw=<Mw> or peak_slip=<m>'
     stop 1
  endif

  sigma = real(sigma_in, pr)
  peak_slip = real(peak_slip_in, pr)
  mw = real(mw_in, pr)
  if (has_mu) then
     mu = real(mu_in, pr)
  else
     mu = 3.0e10_pr
  endif
  if (sigma <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: sigma must be positive'
     stop 1
  endif
  if (has_peak .and. peak_slip < 0._pr) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: peak_slip must be non-negative'
     stop 1
  endif
  if (physical_mode .and. mu <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-slip-smooth: mu must be positive'
     stop 1
  endif

  call read_slip_vtk(trim(infile), lines, nlines, px, py, pz, cell, npoints, ncells)
  if (center_node < 1 .or. center_node > npoints) then
     write(error_unit,'(a,i0,a,i0,a)') 'sgrac-slip-smooth: center_node ', center_node, &
        ' outside valid range 1..', npoints, ''
     stop 1
  endif

  amesh%Nnodes = npoints
  amesh%Ncells = ncells
  allocate(amesh%px(npoints), amesh%py(npoints), amesh%pz(npoints), amesh%cell(ncells,3))
  amesh%px = px
  amesh%py = py
  amesh%pz = pz
  amesh%cell = cell

  allocate(dg(npoints), dg_cell(ncells), slip(ncells))
  allocate(ntoc(npoints), nton(npoints))

  adiff%fast = .true.
  verbose = 0
  call pre_onevsall2d_onvertex(amesh, int(center_node, pin), dg, ntoc, nton)
  if (.not. associated(ntoc(center_node)%ptr)) then
     write(error_unit,'(a,i0,a)') 'sgrac-slip-smooth: center_node ', center_node, ' is not used by any triangle'
     stop 1
  endif
  call onevsall2d(amesh, dg, ntoc, nton, adiff)

  do i=1,ncells
     dg_cell(i) = (dg(cell(i,1)) + dg(cell(i,2)) + dg(cell(i,3))) / 3._pr
     r = dg_cell(i)
     slip(i) = exp(-0.5_pr * (r / sigma)**2)
  enddo

  if (physical_mode) then
     allocate(areas(ncells), slip_scaled(ncells))
     call compute_triangle_areas(px, py, pz, cell, areas)
     target_m0 = m0_from_mw(mw)
     current_m0 = scalar_moment(slip, areas, mu)
     call rescale_slip_to_m0(slip, areas, mu, target_m0, slip_scaled, scale_status)
     if (scale_status /= 0) then
        write(error_unit,'(a)') 'sgrac-slip-smooth: cannot scale slip to target moment'
        stop 1
     endif
     slip = slip_scaled
     final_m0 = scalar_moment(slip, areas, mu)

     write(error_unit,'(a)') 'sgrac-slip-smooth diagnostics:'
     write(error_unit,'(a)') '  mode = physical'
     write(error_unit,'(a,es24.16)') '  mw = ', mw
     write(error_unit,'(a,es24.16)') '  mu = ', mu
     write(error_unit,'(a,es24.16)') '  target M0 = ', target_m0
     write(error_unit,'(a,es24.16)') '  unscaled M0 = ', current_m0
     write(error_unit,'(a,es24.16)') '  final M0 = ', final_m0
     write(error_unit,'(a,es24.16)') '  peak slip = ', maxval(slip)
  else
     slip = peak_slip * slip
  endif

  call write_vtk_with_cell_scalar(trim(outfile), lines, nlines, 'slip', slip, ncells)

  call free_nton(nton)
  call free_ntoc(ntoc)
end program sgrac_slip_smooth
