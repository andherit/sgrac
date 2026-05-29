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
  real(pr) :: p1(3), p2(3), p3(3), v12(3), v13(3), normal(3), ex(3), ey(3)
  real(pr) :: eh(3), eperp(3), gvec(3), globx(3), globy(3)
  real(pr) :: l12, normn, x3, y3, dfdx, dfdy, f1, f2, f3, tmp

  globx = (/1._pr, 0._pr, 0._pr/)
  globy = (/0._pr, 1._pr, 0._pr/)

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

     eh = globx - dot_product(globx, normal)*normal
     tmp = norm3(eh)
     if (tmp <= 100._pr*epsilon(1._pr)) then
        eh = globy - dot_product(globy, normal)*normal
        tmp = norm3(eh)
     endif
     if (tmp <= 100._pr*epsilon(1._pr)) then
        theta(ic) = 0._pr
     else
        eh = eh/tmp
        eperp = cross3(normal, eh)
        theta(ic) = atan2(dot_product(grad_dg(ic,:), eperp), dot_product(grad_dg(ic,:), eh))
     endif
  enddo
end subroutine compute_cell_geometry

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
