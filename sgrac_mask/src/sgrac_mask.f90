program sgrac_mask
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  use forparse
  use sgrac_vtk_mask_io
  use sgrac_radius_models
  use sgrac_mask_border
  implicit none

  character(len=256) :: infile, outfile, model
  character(len=64) :: geometry
  character(len=64) :: border_stop_reason
  character(len=line_len), allocatable :: lines(:)
  type(sgrac_vtk_metadata_type) :: metadata
  integer :: ierr, ierr_r0, ierr_mw, ierr_stressdrop, ierr_mu, ierr_smooth_border
  integer :: ierr_smooth_border_iter_max, ierr_smooth_border_aperture_max, ierr_geometry, ierr_geometry_code, ierr_threshold
  integer :: nlines, npoints, ncell, i, smooth_border_in, smooth_border_iter_max
  integer :: geometry_code
  integer :: border_iter, border_removed, border_added, border_add_candidates
  real :: r0_in, anis_in, theta0_in, rmin_in, mw_in, stressdrop_in, mu_in, threshold_in, smooth_border_aperture_max_in
  real(pr) :: r0, anis, theta0, rmin, mw, stressdrop, mu, threshold, used_threshold, smooth_border_aperture_max
  real(pr) :: border_remove_aperture, border_add_aperture
  real(pr) :: m0, req, atarget, afinal, relerr, alpha, pi, atotal, selected_area
  real(pr) :: d12_from_f1, d12_from_f2, d12, d12_err, tol_threshold, dg_sum_min, dg_sum_max
  real(pr), allocatable :: theta(:), dg_cell(:), dg_sum_cell(:), dg_f1(:), dg_f2(:)
  real(pr), allocatable :: area_cell(:), shape(:), rtheta(:), phi(:)
  integer, allocatable :: mask(:)
  logical :: has_r0, has_mw, has_stressdrop, has_mu, has_physical_keyword, physical_mode
  logical :: smooth_border, metadata_ok, foci_explicit_threshold, foci_magnitude_mode, has_area_cell

  infile = '-'
  outfile = '-'
  model = 'ellipse'
  geometry = ''
  geometry_code = -1
  r0_in = -1.0
  anis_in = 0.0
  theta0_in = 90.0
  rmin_in = 0.0
  mw_in = 0.0
  stressdrop_in = -1.0
  mu_in = -1.0
  threshold_in = -1.0
  smooth_border_in = 0
  smooth_border_iter_max = -1
  smooth_border_aperture_max_in = 60.0
  border_iter = 0
  border_removed = 0
  border_added = 0
  border_add_candidates = 0
  border_stop_reason = 'not_run'
  border_remove_aperture = -1._pr
  border_add_aperture = -1._pr
  used_threshold = -1._pr
  atotal = -1._pr
  selected_area = -1._pr
  has_area_cell = .false.

  ierr = parse_arg('in', infile)
  ierr = parse_arg('out', outfile)
  ierr = parse_arg('model', model)
  ierr_geometry = parse_arg('geometry', geometry)
  ierr_geometry_code = parse_arg('mask_geometry_type_code', geometry_code)
  ierr_r0 = parse_arg('r0', r0_in)
  ierr = parse_arg('anis', anis_in)
  ierr = parse_arg('theta0', theta0_in)
  ierr = parse_arg('rmin', rmin_in)
  ierr_threshold = parse_arg('threshold', threshold_in)
  ierr_mw = parse_arg('mw', mw_in)
  ierr_stressdrop = parse_arg('stressdrop', stressdrop_in)
  ierr_mu = parse_arg('mu', mu_in)
  ierr_smooth_border = parse_arg('smooth_border', smooth_border_in)
  ierr_smooth_border_iter_max = parse_arg('smooth_border_iter_max', smooth_border_iter_max)
  ierr_smooth_border_aperture_max = parse_arg('smooth_border_aperture_max', smooth_border_aperture_max_in)

  if (ierr_r0 == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid r0 value'
     stop 1
  endif
  if (ierr_geometry_code == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid mask_geometry_type_code=<integer>'
     stop 1
  endif
  if (ierr_threshold == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid threshold=<positive distance>'
     stop 1
  endif
  if (ierr_smooth_border == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid smooth_border=<0|1>'
     stop 1
  endif
  if (ierr_smooth_border_iter_max == PARSE_TYPE_ERROR) then
     write(error_unit,'(a)') 'sgrac-mask: invalid smooth_border_iter_max=<positive integer>'
     stop 1
  endif
  if (ierr_smooth_border_aperture_max == PARSE_TYPE_ERROR .or. smooth_border_aperture_max_in <= 0.0) then
     write(error_unit,'(a)') 'sgrac-mask: invalid smooth_border_aperture_max=<positive degrees>'
     stop 1
  endif

  call read_text_file(trim(infile), lines, nlines)
  call get_polydata_counts(lines, nlines, npoints, ncell)
  call read_sgrac_metadata(lines, nlines, metadata)

  if (ierr_geometry_code /= PARSE_OK) then
     if (ierr_geometry == PARSE_OK) then
        select case(trim(geometry))
        case('theta_ellipse','theta','ellipse')
           geometry_code = 2
        case('foci_ellipse','foci')
           geometry_code = 3
        case default
           write(error_unit,'(a,a)') 'sgrac-mask: unknown geometry mode: ', trim(geometry)
           stop 1
        end select
     else if (metadata%present .and. metadata%mask_geometry_type_code /= 0) then
        geometry_code = metadata%mask_geometry_type_code
     else
        geometry_code = 2
     endif
  endif

  call validate_sgrac_geometry_mode(lines, nlines, npoints, ncell, geometry_code, metadata_ok)
  if (.not. metadata_ok) stop 1
  if (metadata%present .and. geometry_code == 3 .and. metadata%distance_source_count /= 2) then
     write(error_unit,'(a,i0)') 'sgrac-mask: foci_ellipse metadata has distance_source_count = ', &
                                metadata%distance_source_count
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
  smooth_border = smooth_border_in /= 0
  foci_explicit_threshold = ierr_threshold == PARSE_OK
  foci_magnitude_mode = geometry_code == 3 .and. .not. foci_explicit_threshold

  if (geometry_code == 2) then
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
  else if (geometry_code == 3) then
     if (ierr_threshold == PARSE_TYPE_ERROR) then
        write(error_unit,'(a)') 'sgrac-mask: invalid threshold=<positive distance>'
        stop 1
     endif
     if (foci_explicit_threshold .and. threshold_in <= 0.0) then
        write(error_unit,'(a)') 'sgrac-mask: foci_ellipse threshold must be positive'
        stop 1
     endif
     if (foci_magnitude_mode .and. ierr_mw == PARSE_TYPE_ERROR) then
        write(error_unit,'(a)') 'sgrac-mask: invalid mw value'
        stop 1
     endif
     if (foci_magnitude_mode .and. ierr_stressdrop == PARSE_TYPE_ERROR) then
        write(error_unit,'(a)') 'sgrac-mask: invalid stressdrop value'
        stop 1
     endif
     if (foci_magnitude_mode .and. (.not. has_mw .or. .not. has_stressdrop)) then
        write(error_unit,'(a)') 'sgrac-mask: foci_ellipse requires either threshold=<distance> or mw=<Mw> stressdrop=<Pa>'
        stop 1
     endif
  else
     write(error_unit,'(a,i0)') 'sgrac-mask: unsupported mask_geometry_type_code = ', geometry_code
     stop 1
  endif

  pi = acos(-1._pr)
  r0 = real(r0_in, pr)
  anis = real(anis_in, pr)
  ! User-facing angles are degrees; radius routines use radians internally.
  theta0 = real(theta0_in, pr) * pi / 180._pr
  rmin = real(rmin_in, pr)
  mw = real(mw_in, pr)
  stressdrop = real(stressdrop_in, pr)
  threshold = real(threshold_in, pr)
  if (has_mu) then
     mu = real(mu_in, pr)
  else
     mu = 3.0e10_pr
  endif

  if (geometry_code == 2 .and. has_r0 .and. r0 <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-mask: r0 must be positive'
     stop 1
  endif
  if (geometry_code == 2 .and. physical_mode .and. stressdrop <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-mask: stressdrop must be positive'
     stop 1
  endif
  if (geometry_code == 3 .and. foci_magnitude_mode .and. stressdrop <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-mask: stressdrop must be positive'
     stop 1
  endif
  if (geometry_code == 2 .and. physical_mode .and. mu <= 0._pr) then
     write(error_unit,'(a)') 'sgrac-mask: mu must be positive'
     stop 1
  endif
  if (geometry_code == 2 .and. .not. physical_mode .and. abs(anis) >= 1._pr) then
     write(error_unit,'(a)') 'sgrac-mask: warning: |anis| >= 1 can produce negative radii before shape check'
  endif

  if (ierr_smooth_border_iter_max /= PARSE_OK) smooth_border_iter_max = ncell
  if (smooth_border_iter_max < 1) then
     write(error_unit,'(a)') 'sgrac-mask: invalid smooth_border_iter_max=<positive integer>'
     stop 1
  endif
  smooth_border_aperture_max = real(smooth_border_aperture_max_in, pr) * pi / 180._pr

  allocate(area_cell(ncell), rtheta(ncell), phi(ncell), mask(ncell))
  area_cell = 0._pr
  if (geometry_code == 2) then
     allocate(theta(ncell), dg_cell(ncell), shape(ncell))
     call read_cell_scalar(lines, nlines, ncell, FIELD_THETA, theta)
     ! The VTK theta field is written in degrees for inspection.
     theta = theta * pi / 180._pr
     call read_cell_scalar(lines, nlines, ncell, FIELD_DG_CELL, dg_cell)

     select case(trim(model))
     case('ellipse')
        call ellipse_shape(theta, ncell, anis, theta0, shape)
        if (minval(shape) <= 0._pr) then
           write(error_unit,'(a,es24.16)') 'sgrac-mask: ellipse shape f(theta) is not positive; min f = ', minval(shape)
           stop 1
        endif

        if (physical_mode) then
           call read_cell_scalar(lines, nlines, ncell, 'area', area_cell)
           has_area_cell = .true.

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
  else
     allocate(dg_sum_cell(ncell))
     call read_cell_scalar(lines, nlines, ncell, FIELD_DG_SUM_CELL, dg_sum_cell)
     dg_sum_min = minval(dg_sum_cell)
     dg_sum_max = maxval(dg_sum_cell)
     has_area_cell = has_cell_scalar(lines, nlines, ncell, 'area')

     if (foci_magnitude_mode .and. .not. has_area_cell) then
        write(error_unit,'(a)') 'sgrac-mask: foci magnitude mode requires CELL_DATA scalar `area`'
        stop 1
     endif
     if (has_area_cell) then
        call read_cell_scalar(lines, nlines, ncell, 'area', area_cell)
        atotal = sum(area_cell)
     endif

     if (foci_explicit_threshold) then
        used_threshold = threshold
     else
        m0 = 10._pr**(1.5_pr * mw + 9.1_pr)
        req = (7._pr * m0 / (16._pr * stressdrop))**(1._pr / 3._pr)
        atarget = pi * req**2

        ! Magnitude gives a circular-crack target area; foci mode selects the
        ! dg_sum_cell threshold by cumulative area on the actual mesh.
        if (atarget > atotal) then
           write(error_unit,'(a)') &
              'sgrac-mask: warning: magnitude-derived target area exceeds available mesh area; selecting whole mesh'
           used_threshold = dg_sum_max
        else
           call select_foci_threshold_by_area(dg_sum_cell, area_cell, ncell, atarget, used_threshold, afinal)
        endif
     endif

     if (metadata%present .and. &
         metadata%source_node_indexes(1) >= 1 .and. metadata%source_node_indexes(1) <= npoints .and. &
         metadata%source_node_indexes(2) >= 1 .and. metadata%source_node_indexes(2) <= npoints .and. &
         has_point_scalar(lines, nlines, npoints, FIELD_DG_F1) .and. &
         has_point_scalar(lines, nlines, npoints, FIELD_DG_F2)) then
        allocate(dg_f1(npoints), dg_f2(npoints))
        call read_point_scalar(lines, nlines, npoints, FIELD_DG_F1, dg_f1)
        call read_point_scalar(lines, nlines, npoints, FIELD_DG_F2, dg_f2)
        d12_from_f1 = dg_f1(metadata%source_node_indexes(2))
        d12_from_f2 = dg_f2(metadata%source_node_indexes(1))
        d12 = 0.5_pr * (d12_from_f1 + d12_from_f2)
        d12_err = abs(d12_from_f1 - d12_from_f2)
        tol_threshold = max(max(1.0e-8_pr * max(1._pr, d12), d12_err), &
                            10._pr * real(epsilon(threshold_in), pr) * max(1._pr, d12))

        if (used_threshold < d12 - tol_threshold) then
           if (foci_explicit_threshold) then
              write(error_unit,'(a)') 'sgrac-mask: threshold is smaller than the geodesic distance between the foci'
              write(error_unit,'(a,es24.16)') '  threshold = ', used_threshold
              write(error_unit,'(a,es24.16)') '  focus distance = ', d12
              write(error_unit,'(a,es24.16)') '  tolerance = ', tol_threshold
              stop 1
           else
              write(error_unit,'(a)') &
                 'sgrac-mask: warning: magnitude-derived threshold is smaller than the focus distance on this discrete mesh'
              write(error_unit,'(a,es24.16)') '  threshold = ', used_threshold
           endif
           write(error_unit,'(a,es24.16)') '  focus distance = ', d12
           write(error_unit,'(a,es24.16)') '  tolerance = ', tol_threshold
        endif
        if (abs(used_threshold - d12) <= tol_threshold) then
           write(error_unit,'(a)') &
              'sgrac-mask: warning: threshold is close to the focus distance; foci ellipse degenerates toward the geodesic segment'
        endif
     endif

     if (used_threshold < dg_sum_min) then
        write(error_unit,'(a)') 'sgrac-mask: warning: threshold is smaller than min(dg_sum_cell); mask will be empty'
     endif
     if (used_threshold >= dg_sum_max) then
        write(error_unit,'(a)') &
           'sgrac-mask: warning: threshold is greater than or equal to max(dg_sum_cell); mask will include the whole mesh'
     endif

     ! In foci mode, Rtheta is not an angular radius law; it stores the
     ! constant foci threshold for compatibility with existing VTK output.
     rtheta = used_threshold
     do i=1,ncell
        phi(i) = dg_sum_cell(i) - used_threshold
        if (phi(i) <= 0._pr) then
           mask(i) = 1
        else
           mask(i) = 0
        endif
     enddo
  endif

  if (smooth_border) then
     call smooth_mask_border(lines, nlines, ncell, phi, mask, smooth_border_iter_max, smooth_border_aperture_max, &
                             border_iter, border_removed, border_added, border_add_candidates, border_stop_reason, &
                             border_remove_aperture, border_add_aperture)
  endif

  if (geometry_code == 3) then
     if (has_area_cell) selected_area = sum(area_cell, mask = mask == 1)
     write(error_unit,'(a)') 'sgrac-mask diagnostics:'
     write(error_unit,'(a)') '  mode = foci_ellipse'
     if (foci_explicit_threshold) then
        write(error_unit,'(a)') '  threshold mode = explicit_threshold'
     else
        write(error_unit,'(a)') '  threshold mode = magnitude'
        write(error_unit,'(a,es24.16)') '  mw = ', mw
        write(error_unit,'(a,es24.16)') '  stressdrop = ', stressdrop
        write(error_unit,'(a,es24.16)') '  M0 = ', m0
        write(error_unit,'(a,es24.16)') '  req = ', req
        write(error_unit,'(a,es24.16)') '  Atarget = ', atarget
     endif
     if (has_area_cell) then
        write(error_unit,'(a,es24.16)') '  Atotal = ', atotal
        write(error_unit,'(a,es24.16)') '  Aselected = ', selected_area
     endif
     write(error_unit,'(a,es24.16)') '  used_threshold = ', used_threshold
     if (metadata%present) then
        write(error_unit,'(a,i0)') '  source node 1 = ', metadata%source_node_indexes(1)
        write(error_unit,'(a,i0)') '  source node 2 = ', metadata%source_node_indexes(2)
     endif
     write(error_unit,'(a,i0)') '  masked cells = ', count(mask == 1)
     if (smooth_border) then
        write(error_unit,'(a,i0)') '  smooth_border iter max = ', smooth_border_iter_max
        write(error_unit,'(a,es24.16)') '  smooth_border aperture max deg = ', smooth_border_aperture_max_in
        write(error_unit,'(a,i0)') '  smooth_border iterations = ', border_iter
        write(error_unit,'(a,i0)') '  smooth_border removed cells = ', border_removed
        write(error_unit,'(a,i0)') '  smooth_border added cells = ', border_added
        write(error_unit,'(a,a)') '  smooth_border stop reason = ', trim(border_stop_reason)
        call write_aperture_diagnostic('  smooth_border final remove aperture deg = ', border_remove_aperture, pi)
        call write_aperture_diagnostic('  smooth_border final add aperture deg = ', border_add_aperture, pi)
        write(error_unit,'(a,i0)') '  smooth_border final add candidates = ', border_add_candidates
     endif
  else if (physical_mode) then
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
     if (smooth_border) then
        write(error_unit,'(a,i0)') '  smooth_border iter max = ', smooth_border_iter_max
        write(error_unit,'(a,es24.16)') '  smooth_border aperture max deg = ', smooth_border_aperture_max_in
        write(error_unit,'(a,i0)') '  smooth_border iterations = ', border_iter
        write(error_unit,'(a,i0)') '  smooth_border removed cells = ', border_removed
        write(error_unit,'(a,i0)') '  smooth_border added cells = ', border_added
        write(error_unit,'(a,a)') '  smooth_border stop reason = ', trim(border_stop_reason)
        call write_aperture_diagnostic('  smooth_border final remove aperture deg = ', border_remove_aperture, pi)
        call write_aperture_diagnostic('  smooth_border final add aperture deg = ', border_add_aperture, pi)
        write(error_unit,'(a,i0)') '  smooth_border final add candidates = ', border_add_candidates
     endif
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
     if (smooth_border) then
        write(error_unit,'(a,i0)') '  smooth_border iter max = ', smooth_border_iter_max
        write(error_unit,'(a,es24.16)') '  smooth_border aperture max deg = ', smooth_border_aperture_max_in
        write(error_unit,'(a,i0)') '  smooth_border iterations = ', border_iter
        write(error_unit,'(a,i0)') '  smooth_border removed cells = ', border_removed
        write(error_unit,'(a,i0)') '  smooth_border added cells = ', border_added
        write(error_unit,'(a,a)') '  smooth_border stop reason = ', trim(border_stop_reason)
        call write_aperture_diagnostic('  smooth_border final remove aperture deg = ', border_remove_aperture, pi)
        call write_aperture_diagnostic('  smooth_border final add aperture deg = ', border_add_aperture, pi)
        write(error_unit,'(a,i0)') '  smooth_border final add candidates = ', border_add_candidates
     endif
  endif

  call write_text_file_with_mask(trim(outfile), lines, nlines, rtheta, phi, mask, ncell)
contains

subroutine write_aperture_diagnostic(label, aperture, pi)
  character(len=*), intent(in) :: label
  real(pr), intent(in) :: aperture, pi

  if (aperture < 0._pr) then
     write(error_unit,'(a)') trim(label)//'none'
  else
     write(error_unit,'(a,es24.16)') label, aperture * 180._pr / pi
  endif
end subroutine write_aperture_diagnostic

subroutine select_foci_threshold_by_area(dg_sum_cell, area_cell, ncell, atarget, threshold, area_selected)
  real(pr), intent(in) :: dg_sum_cell(ncell), area_cell(ncell), atarget
  integer, intent(in) :: ncell
  real(pr), intent(out) :: threshold, area_selected
  integer, allocatable :: order(:)
  integer :: i

  allocate(order(ncell))
  do i=1,ncell
     order(i) = i
  enddo
  call sort_index_by_value(dg_sum_cell, order, 1, ncell)

  area_selected = 0._pr
  threshold = dg_sum_cell(order(1))
  do i=1,ncell
     area_selected = area_selected + area_cell(order(i))
     threshold = dg_sum_cell(order(i))
     if (area_selected >= atarget) exit
  enddo
end subroutine select_foci_threshold_by_area

recursive subroutine sort_index_by_value(vals, idx, left, right)
  real(pr), intent(in) :: vals(:)
  integer, intent(inout) :: idx(:)
  integer, intent(in) :: left, right
  integer :: i, j, tmp
  real(pr) :: pivot

  if (left >= right) return

  i = left
  j = right
  pivot = vals(idx((left + right) / 2))
  do
     do while (vals(idx(i)) < pivot)
        i = i + 1
     enddo
     do while (vals(idx(j)) > pivot)
        j = j - 1
     enddo
     if (i <= j) then
        tmp = idx(i)
        idx(i) = idx(j)
        idx(j) = tmp
        i = i + 1
        j = j - 1
     endif
     if (i > j) exit
  enddo

  if (left < j) call sort_index_by_value(vals, idx, left, j)
  if (i < right) call sort_index_by_value(vals, idx, i, right)
end subroutine sort_index_by_value
end program sgrac_mask
