program sgrac_mask
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  use forparse
  use sgrac_vtk_mask_io
  use sgrac_radius_models
  implicit none

  character(len=256) :: infile, outfile, model
  character(len=line_len), allocatable :: lines(:)
  integer :: ierr, ierr_r0, ierr_mw, ierr_stressdrop, ierr_mu, nlines, npoints, ncell, i
  real :: r0_in, anis_in, theta0_in, rmin_in, mw_in, stressdrop_in, mu_in
  real(pr) :: r0, anis, theta0, rmin, mw, stressdrop, mu
  real(pr) :: m0, req, atarget, afinal, relerr, alpha, pi
  real(pr), allocatable :: theta(:), dg_cell(:), area_cell(:), shape(:), rtheta(:), phi(:)
  integer, allocatable :: mask(:)
  logical :: has_r0, has_mw, has_stressdrop, has_mu, has_physical_keyword, physical_mode

  infile = '-'
  outfile = '-'
  model = 'ellipse'
  r0_in = -1.0
  anis_in = 0.0
  theta0_in = 0.0
  rmin_in = 0.0
  mw_in = 0.0
  stressdrop_in = -1.0
  mu_in = -1.0

  ierr = parse_arg('in', infile)
  ierr = parse_arg('out', outfile)
  ierr = parse_arg('model', model)
  ierr_r0 = parse_arg('r0', r0_in)
  ierr = parse_arg('anis', anis_in)
  ierr = parse_arg('theta0', theta0_in)
  ierr = parse_arg('rmin', rmin_in)
  ierr_mw = parse_arg('mw', mw_in)
  ierr_stressdrop = parse_arg('stressdrop', stressdrop_in)
  ierr_mu = parse_arg('mu', mu_in)

  if (ierr_r0 == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid r0 value'
     stop 1
  endif

  has_r0 = ierr_r0 == PARSE_OK
  has_mw = ierr_mw == PARSE_OK
  has_stressdrop = ierr_stressdrop == PARSE_OK
  has_mu = ierr_mu == PARSE_OK
  has_physical_keyword = ierr_mw == PARSE_OK .or. ierr_mw == PARSE_TYPE_ERROR .or. &
                         ierr_stressdrop == PARSE_OK .or. ierr_stressdrop == PARSE_TYPE_ERROR .or. &
                         ierr_mu == PARSE_OK .or. ierr_mu == PARSE_TYPE_ERROR
  physical_mode = .not. has_r0

  if (has_r0 .and. has_physical_keyword) then
     write(error_unit,'(a)') 'sgrac-mask: warning: r0 is present; physical scaling keywords are ignored'
  endif

  if (physical_mode .and. ierr_mw == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid mw value'
     stop 1
  endif
  if (physical_mode .and. ierr_stressdrop == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid stressdrop value'
     stop 1
  endif
  if (physical_mode .and. ierr_mu == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid mu value'
     stop 1
  endif

  if (.not. has_r0 .and. (.not. has_mw .or. .not. has_stressdrop)) then
     write(error_unit,'(a)') 'sgrac-mask: provide either r0=<radius in meters> or mw=<Mw> stressdrop=<Pa>'
     stop 1
  endif

  r0 = real(r0_in, pr)
  anis = real(anis_in, pr)
  theta0 = real(theta0_in, pr)
  rmin = real(rmin_in, pr)
  mw = real(mw_in, pr)
  stressdrop = real(stressdrop_in, pr)
  if (has_mu) then
     mu = real(mu_in, pr)
  else
     mu = 3.0e10_pr
  endif
  pi = acos(-1._pr)

  if (has_r0 .and. r0 <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-mask: r0 must be positive'
     stop 1
  endif
  if (physical_mode .and. stressdrop <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-mask: stressdrop must be positive'
     stop 1
  endif
  if (physical_mode .and. mu <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-mask: mu must be positive'
     stop 1
  endif
  if (.not. physical_mode .and. abs(anis) >= 1._pr) then
     write(error_unit,'(a)') 'sgrac-mask: warning: |anis| >= 1 can produce negative radii before shape check'
  endif

  call read_text_file(trim(infile), lines, nlines)
  call get_polydata_counts(lines, nlines, npoints, ncell)

  allocate(theta(ncell), dg_cell(ncell), area_cell(ncell), shape(ncell), rtheta(ncell), phi(ncell), mask(ncell))
  call read_cell_scalar(lines, nlines, ncell, 'theta', theta)
  call read_cell_scalar(lines, nlines, ncell, 'dg_cell', dg_cell)

  select case(trim(model))
  case('ellipse')
     call ellipse_shape(theta, ncell, anis, theta0, shape)
     if (minval(shape) <= 0._pr) then
        write(error_unit,'(a,es24.16)') 'sgrac-mask: ellipse shape f(theta) is not positive; min f = ', minval(shape)
        stop 1
     endif

     if (physical_mode) then
        call read_cell_scalar(lines, nlines, ncell, 'area', area_cell)

        m0 = 10._pr**(1.5_pr * mw + 9.1_pr)
        req = (7._pr * m0 / (16._pr * stressdrop))**(1._pr / 3._pr)
        atarget = pi * req**2

        call solve_alpha_for_area(shape, dg_cell, area_cell, ncell, atarget, alpha, afinal, relerr)
        call radius_scaled_shape(shape, ncell, alpha, rtheta)
     else
        call radius_ellipse(theta, ncell, r0, anis, theta0, rmin, rtheta)
     endif
  case default
     write(error_unit,'(a,a)') 'sgrac-mask: unknown model: ', trim(model)
     stop 1
  end select

  do i=1,ncell
     phi(i) = dg_cell(i) - rtheta(i)
     if (phi(i) < 0._pr) then
        mask(i) = 1
     else
        mask(i) = 0
     endif
  enddo

  if (physical_mode) then
     afinal = sum(area_cell, mask = mask == 1)
     relerr = abs(afinal - atarget) / atarget
     write(error_unit,'(a)') 'sgrac-mask diagnostics:'
     write(error_unit,'(a)') '  mode = physical'
     write(error_unit,'(a,es24.16)') '  mw = ', mw
     write(error_unit,'(a,es24.16)') '  stressdrop = ', stressdrop
     write(error_unit,'(a,es24.16)') '  mu = ', mu
     write(error_unit,'(a,es24.16)') '  M0 = ', m0
     write(error_unit,'(a,es24.16)') '  req = ', req
     write(error_unit,'(a,es24.16)') '  Atarget = ', atarget
     write(error_unit,'(a,es24.16)') '  Afinal = ', afinal
     write(error_unit,'(a,es24.16)') '  relative area error = ', relerr
     if (afinal > 0._pr) then
        write(error_unit,'(a,es24.16)') '  implied mean slip = ', m0 / (mu * afinal)
     else
        write(error_unit,'(a)') '  implied mean slip = undefined'
     endif
     if (afinal < 0.98_pr * atarget) then
        write(error_unit,'(a)') 'sgrac-mask: warning: final masked area is less than 98% of target area'
     endif
  else
     write(error_unit,'(a)') 'sgrac-mask diagnostics:'
     write(error_unit,'(a)') '  mode = debug'
  endif

  call write_text_file_with_mask(trim(outfile), lines, nlines, rtheta, phi, mask, ncell)
end program sgrac_mask
