module sgrac_radius_models
  use generic
  implicit none
contains

subroutine radius_ellipse(theta, ncell, r0, anis, theta0, rmin, rtheta)
  integer, intent(in) :: ncell
  real(pr), intent(in) :: theta(ncell), r0, anis, theta0, rmin
  real(pr), intent(out) :: rtheta(ncell)
  integer :: i

  do i=1,ncell
     rtheta(i) = r0 * (1._pr + anis * cos(2._pr * (theta(i) - theta0)))
     if (rtheta(i) < rmin) rtheta(i) = rmin
  enddo
end subroutine radius_ellipse

subroutine ellipse_shape(theta, ncell, anis, theta0, shape)
  integer, intent(in) :: ncell
  real(pr), intent(in) :: theta(ncell), anis, theta0
  real(pr), intent(out) :: shape(ncell)
  integer :: i

  do i=1,ncell
     shape(i) = 1._pr + anis * cos(2._pr * (theta(i) - theta0))
  enddo
end subroutine ellipse_shape

subroutine radius_scaled_shape(shape, ncell, alpha, rtheta)
  integer, intent(in) :: ncell
  real(pr), intent(in) :: shape(ncell), alpha
  real(pr), intent(out) :: rtheta(ncell)

  rtheta = alpha * shape
end subroutine radius_scaled_shape

subroutine solve_alpha_for_area(shape, dg_cell, area_cell, ncell, atarget, alpha, afinal, relerr)
  integer, intent(in) :: ncell
  real(pr), intent(in) :: shape(ncell), dg_cell(ncell), area_cell(ncell), atarget
  real(pr), intent(out) :: alpha, afinal, relerr
  integer, parameter :: max_iter = 40
  integer :: iter
  real(pr) :: alpha_low, alpha_high, alpha_mid, area_mid, area_high
  real(pr) :: err_mid, best_alpha, best_area, best_err

  alpha_low = 0._pr
  alpha_high = 2._pr * maxval(dg_cell)
  call masked_area_for_alpha(shape, dg_cell, area_cell, ncell, alpha_high, area_high)

  if (area_high < atarget) then
     alpha = alpha_high
     afinal = area_high
     relerr = abs(afinal - atarget) / atarget
     return
  endif

  best_alpha = alpha_high
  best_area = area_high
  best_err = abs(area_high - atarget) / atarget

  do iter=1,max_iter
     alpha_mid = 0.5_pr * (alpha_low + alpha_high)
     call masked_area_for_alpha(shape, dg_cell, area_cell, ncell, alpha_mid, area_mid)
     err_mid = abs(area_mid - atarget) / atarget

     if (err_mid < best_err) then
        best_alpha = alpha_mid
        best_area = area_mid
        best_err = err_mid
     endif
     if (err_mid < 1.e-3_pr) exit

     if (area_mid < atarget) then
        alpha_low = alpha_mid
     else
        alpha_high = alpha_mid
     endif
  enddo

  alpha = best_alpha
  afinal = best_area
  relerr = best_err
end subroutine solve_alpha_for_area

subroutine masked_area_for_alpha(shape, dg_cell, area_cell, ncell, alpha, area_masked)
  integer, intent(in) :: ncell
  real(pr), intent(in) :: shape(ncell), dg_cell(ncell), area_cell(ncell), alpha
  real(pr), intent(out) :: area_masked
  integer :: i

  area_masked = 0._pr
  do i=1,ncell
     if (dg_cell(i) < alpha * shape(i)) area_masked = area_masked + area_cell(i)
  enddo
end subroutine masked_area_for_alpha

end module sgrac_radius_models
