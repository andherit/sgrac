module sgrac_cell_geometry
  use generic
  use LAT_mesh
  implicit none
contains

subroutine compute_cell_geometry(amesh, dg, area, dg_cell, centroid, grad_dg, theta)
  type(mesh), intent(in) :: amesh
  real(pr), intent(in) :: dg(amesh%Nnodes)
  real(pr), intent(out) :: area(amesh%Ncells), dg_cell(amesh%Ncells), theta(amesh%Ncells)
  real(pr), intent(out) :: centroid(amesh%Ncells,3), grad_dg(amesh%Ncells,3)
  integer(pin) :: ic, n1, n2, n3
  real(pr), parameter :: tol = 100._pr*epsilon(1._pr)
  real(pr) :: p1(3), p2(3), p3(3), v12(3), v13(3), normal(3), ex(3), ey(3)
  real(pr) :: er(3), edip(3), estrike(3), gvec(3), normal_use(3), ez(3)
  real(pr) :: l12, normn, x3, y3, dfdx, dfdy, f1, f2, f3, tmp

  ez = (/0._pr, 0._pr, 1._pr/)

  do ic=1,amesh%Ncells
     n1 = amesh%cell(ic,1); n2 = amesh%cell(ic,2); n3 = amesh%cell(ic,3)
     p1 = (/amesh%px(n1), amesh%py(n1), amesh%pz(n1)/)
     p2 = (/amesh%px(n2), amesh%py(n2), amesh%pz(n2)/)
     p3 = (/amesh%px(n3), amesh%py(n3), amesh%pz(n3)/)
     v12 = p2-p1
     v13 = p3-p1
     normal = cross3(v12, v13)
     normn = norm3(normal)
     area(ic) = 0.5_pr * normn
     centroid(ic,:) = (p1+p2+p3)/3._pr
     dg_cell(ic) = (dg(n1)+dg(n2)+dg(n3))/3._pr

     if (normn <= epsilon(1._pr)) then
        grad_dg(ic,:) = 0._pr
        theta(ic) = 0._pr
        cycle
     endif

     normal = normal / normn
     l12 = norm3(v12)
     if (l12 <= epsilon(1._pr)) then
        grad_dg(ic,:) = 0._pr
        theta(ic) = 0._pr
        cycle
     endif
     ex = v12 / l12
     ey = cross3(normal, ex)
     x3 = dot_product(v13, ex)
     y3 = dot_product(v13, ey)
     f1 = dg(n1); f2 = dg(n2); f3 = dg(n3)
     if (abs(y3) <= epsilon(1._pr)) then
        grad_dg(ic,:) = 0._pr
     else
        dfdx = (f2-f1)/l12
        dfdy = (f3-f1-dfdx*x3)/y3
        gvec = dfdx*ex + dfdy*ey
        grad_dg(ic,:) = gvec
     endif

     ! Theta is the local bearing of grad_dg in the fault tangent plane.
     ! z is positive downward; theta=0 points down-dip, and positive theta points
     ! toward the strike direction defined by an upward-looking cell normal.
     ! Fallback theta=0 cases below are geometrically undefined angles, not physical zero bearings.
     gvec = grad_dg(ic,:)
     tmp = norm3(gvec)
     if (tmp <= tol) then
        theta(ic) = 0._pr
        cycle
     endif
     er = gvec / tmp

     normal_use = normal
     if (normal_use(3) > 0._pr) normal_use = -normal_use

     edip = ez - dot_product(ez, normal_use)*normal_use
     tmp = norm3(edip)
     if (tmp <= tol) then
        theta(ic) = 0._pr
        cycle
     endif
     edip = edip / tmp

     estrike = cross3(edip, normal_use)
     tmp = norm3(estrike)
     if (tmp <= tol) then
        theta(ic) = 0._pr
        cycle
     endif
     estrike = estrike / tmp

     theta(ic) = atan2(dot_product(er, estrike), dot_product(er, edip))
  enddo
end subroutine compute_cell_geometry

subroutine compute_foci_cell_geometry(amesh, dg_f1, dg_f2, area, centroid, dg_f1_cell, dg_f2_cell, dg_sum_cell)
  type(mesh), intent(in) :: amesh
  real(pr), intent(in) :: dg_f1(amesh%Nnodes), dg_f2(amesh%Nnodes)
  real(pr), intent(out) :: area(amesh%Ncells), centroid(amesh%Ncells,3)
  real(pr), intent(out) :: dg_f1_cell(amesh%Ncells), dg_f2_cell(amesh%Ncells), dg_sum_cell(amesh%Ncells)
  integer(pin) :: ic, n1, n2, n3
  real(pr) :: p1(3), p2(3), p3(3), v12(3), v13(3)

  do ic=1,amesh%Ncells
     n1 = amesh%cell(ic,1); n2 = amesh%cell(ic,2); n3 = amesh%cell(ic,3)
     p1 = (/amesh%px(n1), amesh%py(n1), amesh%pz(n1)/)
     p2 = (/amesh%px(n2), amesh%py(n2), amesh%pz(n2)/)
     p3 = (/amesh%px(n3), amesh%py(n3), amesh%pz(n3)/)
     v12 = p2-p1
     v13 = p3-p1
     area(ic) = 0.5_pr * norm3(cross3(v12, v13))
     centroid(ic,:) = (p1+p2+p3)/3._pr
     dg_f1_cell(ic) = (dg_f1(n1)+dg_f1(n2)+dg_f1(n3))/3._pr
     dg_f2_cell(ic) = (dg_f2(n1)+dg_f2(n2)+dg_f2(n3))/3._pr
     dg_sum_cell(ic) = dg_f1_cell(ic) + dg_f2_cell(ic)
  enddo
end subroutine compute_foci_cell_geometry

pure function cross3(a,b) result(c)
  real(pr), intent(in) :: a(3), b(3)
  real(pr) :: c(3)
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
end function cross3

pure function norm3(a) result(v)
  real(pr), intent(in) :: a(3)
  real(pr) :: v
  v = sqrt(dot_product(a,a))
end function norm3

end module sgrac_cell_geometry
