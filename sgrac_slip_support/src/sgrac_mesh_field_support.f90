module sgrac_mesh_field_support
  use generic
  implicit none
contains

subroutine compute_triangle_centers(px, py, pz, cell, centers)
  real(pr), intent(in) :: px(:), py(:), pz(:)
  integer(pin), intent(in) :: cell(:,:)
  real(pr), intent(out) :: centers(:,:)
  integer :: i
  integer(pin) :: n1, n2, n3

  do i=1,size(cell,1)
     n1 = cell(i,1)
     n2 = cell(i,2)
     n3 = cell(i,3)
     centers(i,1) = (px(n1) + px(n2) + px(n3)) / 3._pr
     centers(i,2) = (py(n1) + py(n2) + py(n3)) / 3._pr
     centers(i,3) = (pz(n1) + pz(n2) + pz(n3)) / 3._pr
  enddo
end subroutine compute_triangle_centers

subroutine compute_triangle_areas(px, py, pz, cell, areas)
  real(pr), intent(in) :: px(:), py(:), pz(:)
  integer(pin), intent(in) :: cell(:,:)
  real(pr), intent(out) :: areas(:)
  integer :: i
  integer(pin) :: n1, n2, n3
  real(pr) :: v12(3), v13(3), normal(3)

  do i=1,size(cell,1)
     n1 = cell(i,1)
     n2 = cell(i,2)
     n3 = cell(i,3)
     v12 = (/px(n2)-px(n1), py(n2)-py(n1), pz(n2)-pz(n1)/)
     v13 = (/px(n3)-px(n1), py(n3)-py(n1), pz(n3)-pz(n1)/)
     normal = cross3(v12, v13)
     areas(i) = 0.5_pr * norm3(normal)
  enddo
end subroutine compute_triangle_areas

pure function cross3(a, b) result(c)
  real(pr), intent(in) :: a(3), b(3)
  real(pr) :: c(3)

  c(1) = a(2)*b(3) - a(3)*b(2)
  c(2) = a(3)*b(1) - a(1)*b(3)
  c(3) = a(1)*b(2) - a(2)*b(1)
end function cross3

pure function norm3(a) result(v)
  real(pr), intent(in) :: a(3)
  real(pr) :: v

  v = sqrt(dot_product(a, a))
end function norm3

end module sgrac_mesh_field_support
