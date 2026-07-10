module sgrac_slip_support
  use generic
  implicit none
contains

pure function scalar_moment(slip, areas, mu) result(m0)
  real(pr), intent(in) :: slip(:), areas(:), mu
  real(pr) :: m0

  m0 = mu * sum(slip * areas)
end function scalar_moment

pure function mw_from_m0(m0) result(mw)
  real(pr), intent(in) :: m0
  real(pr) :: mw

  mw = (2._pr / 3._pr) * (log10(m0) - 9.1_pr)
end function mw_from_m0

pure function m0_from_mw(mw) result(m0)
  real(pr), intent(in) :: mw
  real(pr) :: m0

  m0 = 10._pr**(1.5_pr * mw + 9.1_pr)
end function m0_from_mw

subroutine rescale_slip_to_m0(slip, areas, mu, target_m0, slip_scaled, status)
  real(pr), intent(in) :: slip(:), areas(:), mu, target_m0
  real(pr), intent(out) :: slip_scaled(:)
  integer, intent(out), optional :: status
  real(pr) :: current_m0, scale

  current_m0 = scalar_moment(slip, areas, mu)
  if (current_m0 <= 0._pr .or. mu <= 0._pr .or. target_m0 < 0._pr) then
     slip_scaled = 0._pr
     if (present(status)) status = 1
     return
  endif

  scale = target_m0 / current_m0
  slip_scaled = slip * scale
  if (present(status)) status = 0
end subroutine rescale_slip_to_m0

end module sgrac_slip_support
