module sgrac_mask_border
  use generic
  implicit none
contains

subroutine smooth_mask_border(lines, nlines, ncell, phi, mask, iter_max, aperture_max, niter, nremove, nadd, &
                              nadd_candidates, stop_reason, final_remove_aperture, final_add_aperture)
  character(len=*), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell, iter_max
  real(pr), intent(in) :: aperture_max
  real(pr), intent(in) :: phi(ncell)
  integer, intent(inout) :: mask(ncell)
  integer, intent(out) :: niter, nremove, nadd, nadd_candidates
  character(len=*), intent(out) :: stop_reason
  real(pr), intent(out) :: final_remove_aperture, final_add_aperture
  integer(pin), allocatable :: cell(:,:), ea(:), eb(:), ec(:), ee(:)
  integer, allocatable :: border_count(:), linked_count(:)
  integer, allocatable :: mbe_count(:), mbe_edge1(:), mbe_edge2(:)
  real(pr), allocatable :: px(:), py(:), pz(:), aperture(:)
  integer :: i, k, nedge, iter, remove_cell, add_cell
  integer :: npoints
  real(pr) :: remove_aperture, add_aperture

  call read_points(lines, nlines, npoints, px, py, pz)
  call read_triangles(lines, nlines, ncell, cell)
  if (minval(cell) < 1_pin .or. maxval(cell) > int(npoints, pin)) then
     stop 'sgrac-mask: border smoothing connectivity references an invalid point'
  endif

  niter = 0
  nremove = 0
  nadd = 0
  nadd_candidates = 0
  final_remove_aperture = -1._pr
  final_add_aperture = -1._pr
  stop_reason = 'iter_max'
  if (iter_max <= 0) return

  nedge = 3*ncell
  allocate(ea(nedge), eb(nedge), ec(nedge), ee(nedge))
  k = 0
  do i=1,ncell
     call add_edge(cell(i,1), cell(i,2), i, 1, ea, eb, ec, ee, k)
     call add_edge(cell(i,2), cell(i,3), i, 2, ea, eb, ec, ee, k)
     call add_edge(cell(i,3), cell(i,1), i, 3, ea, eb, ec, ee, k)
  enddo

  call sort_edges(ea, eb, ec, ee, 1, nedge)

  allocate(border_count(ncell), linked_count(ncell), mbe_count(ncell), mbe_edge1(ncell), mbe_edge2(ncell))
  allocate(aperture(ncell))

  do iter=1,iter_max
     call compute_border_geometry(nedge, ea, eb, ec, ee, cell, px, py, pz, mask, border_count, linked_count, &
                                  mbe_count, mbe_edge1, mbe_edge2, aperture)
     call select_remove_cell(ncell, phi, mask, border_count, aperture, remove_cell, remove_aperture)
     call select_add_cell(ncell, phi, mask, linked_count, aperture, add_cell, add_aperture)
     if (remove_cell <= 0) then
        stop_reason = 'no_removal_candidate'
        exit
     endif
     if (add_cell <= 0) then
        stop_reason = 'no_addition_candidate'
        exit
     endif
     if (remove_aperture > aperture_max .or. add_aperture > aperture_max) then
        stop_reason = 'aperture_threshold'
        exit
     endif

     mask(remove_cell) = 0
     mask(add_cell) = 1
     niter = niter + 1
     nremove = nremove + 1
     nadd = nadd + 1
  enddo

  call compute_border_geometry(nedge, ea, eb, ec, ee, cell, px, py, pz, mask, border_count, linked_count, &
                               mbe_count, mbe_edge1, mbe_edge2, aperture)
  nadd_candidates = count(mask == 0 .and. linked_count == 2)
  call select_remove_cell(ncell, phi, mask, border_count, aperture, remove_cell, final_remove_aperture)
  call select_add_cell(ncell, phi, mask, linked_count, aperture, add_cell, final_add_aperture)
end subroutine smooth_mask_border

subroutine compute_border_geometry(nedge, ea, eb, ec, ee, cell, px, py, pz, mask, border_count, linked_count, &
                                   mbe_count, mbe_edge1, mbe_edge2, aperture)
  integer, intent(in) :: nedge
  integer(pin), intent(in) :: ea(:), eb(:), ec(:), ee(:), cell(:,:)
  real(pr), intent(in) :: px(:), py(:), pz(:)
  integer, intent(in) :: mask(:)
  integer, intent(out) :: border_count(:), linked_count(:), mbe_count(:), mbe_edge1(:), mbe_edge2(:)
  real(pr), intent(out) :: aperture(:)
  integer :: k, group_start, group_end, group_size, c1, c2
  integer :: i

  border_count = 0
  linked_count = 0
  mbe_count = 0
  mbe_edge1 = 0
  mbe_edge2 = 0
  aperture = huge(1._pr)

  k = 1
  do while (k <= nedge)
     group_start = k
     do
        if (k > nedge) exit
        if (ea(k) /= ea(group_start) .or. eb(k) /= eb(group_start)) exit
        k = k + 1
     enddo
     group_end = k - 1
     group_size = group_end - group_start + 1

     if (group_size == 2) then
        c1 = ec(group_start)
        c2 = ec(group_start+1)
        if (mask(c1) /= mask(c2)) then
           call record_mbe_edge(c1, ee(group_start), mbe_count, mbe_edge1, mbe_edge2)
           call record_mbe_edge(c2, ee(group_start+1), mbe_count, mbe_edge1, mbe_edge2)
           if (mask(c1) == 1) then
              border_count(c1) = border_count(c1) + 1
              linked_count(c2) = linked_count(c2) + 1
           else
              border_count(c2) = border_count(c2) + 1
              linked_count(c1) = linked_count(c1) + 1
           endif
        endif
     endif
  enddo

  do i=1,size(mbe_count)
     if (mbe_count(i) == 2) then
        aperture(i) = cell_aperture(cell, i, mbe_edge1(i), mbe_edge2(i), px, py, pz)
     endif
  enddo
end subroutine compute_border_geometry

subroutine record_mbe_edge(cell_id, edge_id, mbe_count, mbe_edge1, mbe_edge2)
  integer, intent(in) :: cell_id
  integer(pin), intent(in) :: edge_id
  integer, intent(inout) :: mbe_count(:), mbe_edge1(:), mbe_edge2(:)

  mbe_count(cell_id) = mbe_count(cell_id) + 1
  if (mbe_count(cell_id) == 1) then
     mbe_edge1(cell_id) = int(edge_id)
  else if (mbe_count(cell_id) == 2) then
     mbe_edge2(cell_id) = int(edge_id)
  endif
end subroutine record_mbe_edge

subroutine select_remove_cell(ncell, phi, mask, border_count, aperture, cell_id, best_aperture)
  integer, intent(in) :: ncell
  real(pr), intent(in) :: phi(ncell)
  real(pr), intent(in) :: aperture(ncell)
  integer, intent(in) :: mask(ncell), border_count(ncell)
  integer, intent(out) :: cell_id
  real(pr), intent(out) :: best_aperture
  integer :: i
  real(pr) :: best_phi
  real(pr), parameter :: angle_tol = 100._pr*epsilon(1._pr)

  cell_id = 0
  best_phi = -huge(1._pr)
  best_aperture = -1._pr
  do i=1,ncell
     if (mask(i) == 1 .and. border_count(i) == 2) then
        if (cell_id <= 0 .or. aperture(i) < best_aperture - angle_tol .or. &
            (abs(aperture(i) - best_aperture) <= angle_tol .and. phi(i) > best_phi)) then
           best_aperture = aperture(i)
           best_phi = phi(i)
           cell_id = i
        endif
     endif
  enddo
end subroutine select_remove_cell

subroutine select_add_cell(ncell, phi, mask, linked_count, aperture, cell_id, best_aperture)
  integer, intent(in) :: ncell
  real(pr), intent(in) :: phi(ncell)
  real(pr), intent(in) :: aperture(ncell)
  integer, intent(in) :: mask(ncell), linked_count(ncell)
  integer, intent(out) :: cell_id
  real(pr), intent(out) :: best_aperture
  integer :: i
  real(pr) :: best_phi
  real(pr), parameter :: angle_tol = 100._pr*epsilon(1._pr)

  cell_id = 0
  best_phi = huge(1._pr)
  best_aperture = -1._pr
  do i=1,ncell
     if (mask(i) == 0 .and. linked_count(i) == 2) then
        if (cell_id <= 0 .or. aperture(i) < best_aperture - angle_tol .or. &
            (abs(aperture(i) - best_aperture) <= angle_tol .and. phi(i) < best_phi)) then
           best_aperture = aperture(i)
           best_phi = phi(i)
           cell_id = i
        endif
     endif
  enddo
end subroutine select_add_cell

pure function cell_aperture(cell, cell_id, edge1, edge2, px, py, pz) result(angle)
  integer(pin), intent(in) :: cell(:,:)
  integer, intent(in) :: cell_id, edge1, edge2
  real(pr), intent(in) :: px(:), py(:), pz(:)
  real(pr) :: angle
  real(pr) :: u(3), v(3), cosang, nu, nv
  integer :: common_vertex, node_a, node_b, node_c

  common_vertex = common_edge_vertex(edge1, edge2)
  if (common_vertex <= 0) then
     angle = huge(1._pr)
     return
  endif

  select case(common_vertex)
  case(1)
     node_b = cell(cell_id,1)
     node_a = cell(cell_id,2)
     node_c = cell(cell_id,3)
  case(2)
     node_b = cell(cell_id,2)
     node_a = cell(cell_id,1)
     node_c = cell(cell_id,3)
  case default
     node_b = cell(cell_id,3)
     node_a = cell(cell_id,2)
     node_c = cell(cell_id,1)
  end select

  u = (/px(node_a)-px(node_b), py(node_a)-py(node_b), pz(node_a)-pz(node_b)/)
  v = (/px(node_c)-px(node_b), py(node_c)-py(node_b), pz(node_c)-pz(node_b)/)
  nu = norm3(u)
  nv = norm3(v)
  if (nu <= epsilon(1._pr) .or. nv <= epsilon(1._pr)) then
     angle = huge(1._pr)
     return
  endif

  cosang = dot_product(u, v) / (nu * nv)
  cosang = max(-1._pr, min(1._pr, cosang))
  angle = acos(cosang)
end function cell_aperture

pure integer function common_edge_vertex(edge1, edge2) result(vertex_id)
  integer, intent(in) :: edge1, edge2

  vertex_id = 0
  if ((edge1 == 1 .and. edge2 == 2) .or. (edge1 == 2 .and. edge2 == 1)) vertex_id = 2
  if ((edge1 == 2 .and. edge2 == 3) .or. (edge1 == 3 .and. edge2 == 2)) vertex_id = 3
  if ((edge1 == 1 .and. edge2 == 3) .or. (edge1 == 3 .and. edge2 == 1)) vertex_id = 1
end function common_edge_vertex

pure function norm3(a) result(v)
  real(pr), intent(in) :: a(3)
  real(pr) :: v

  v = sqrt(dot_product(a,a))
end function norm3

subroutine read_points(lines, nlines, npoints, px, py, pz)
  character(len=*), intent(in) :: lines(:)
  integer, intent(in) :: nlines
  integer, intent(out) :: npoints
  real(pr), allocatable, intent(out) :: px(:), py(:), pz(:)
  character(len=64) :: key, dtype
  integer :: i, j, ios, point_line

  point_line = -1
  npoints = -1
  do i=1,nlines
     key = ''
     dtype = ''
     read(lines(i),*,iostat=ios) key, npoints, dtype
     if (ios == 0 .and. trim(key) == 'POINTS') then
        point_line = i
        exit
     endif
  enddo
  if (point_line < 0) stop 'sgrac-mask: POINTS section not found for border smoothing'
  if (point_line + npoints > nlines) stop 'sgrac-mask: POINTS section is truncated'

  allocate(px(npoints), py(npoints), pz(npoints))
  do j=1,npoints
     read(lines(point_line+j),*,iostat=ios) px(j), py(j), pz(j)
     if (ios /= 0) stop 'sgrac-mask: error while reading border smoothing coordinates'
  enddo
end subroutine read_points

subroutine read_triangles(lines, nlines, ncell, cell)
  character(len=*), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell
  integer(pin), allocatable, intent(out) :: cell(:,:)
  character(len=64) :: key
  integer :: i, j, ios, total_size, nverts, poly_line, ncells_file

  poly_line = -1
  ncells_file = -1
  do i=1,nlines
     key = ''
     read(lines(i),*,iostat=ios) key, ncells_file, total_size
     if (ios == 0 .and. (trim(key) == 'POLYGONS' .or. trim(key) == 'CELLS')) then
        poly_line = i
        exit
     endif
  enddo
  if (poly_line < 0) stop 'sgrac-mask: POLYGONS/CELLS section not found for border smoothing'
  if (ncells_file /= ncell) stop 'sgrac-mask: cell count mismatch while reading border smoothing connectivity'
  if (poly_line + ncell > nlines) stop 'sgrac-mask: POLYGONS/CELLS section is truncated'

  allocate(cell(ncell,3))
  do j=1,ncell
     read(lines(poly_line+j),*,iostat=ios) nverts, cell(j,1), cell(j,2), cell(j,3)
     if (ios /= 0) stop 'sgrac-mask: error while reading border smoothing connectivity'
     if (nverts /= 3) stop 'sgrac-mask: border smoothing supports triangular cells only'
     cell(j,1:3) = cell(j,1:3) + 1_pin
  enddo
end subroutine read_triangles

subroutine add_edge(n1, n2, icell, iedge, ea, eb, ec, ee, nedge)
  integer(pin), intent(in) :: n1, n2
  integer, intent(in) :: icell, iedge
  integer(pin), intent(inout) :: ea(:), eb(:), ec(:), ee(:)
  integer, intent(inout) :: nedge

  nedge = nedge + 1
  ea(nedge) = min(n1, n2)
  eb(nedge) = max(n1, n2)
  ec(nedge) = int(icell, pin)
  ee(nedge) = int(iedge, pin)
end subroutine add_edge

recursive subroutine sort_edges(ea, eb, ec, ee, left, right)
  integer(pin), intent(inout) :: ea(:), eb(:), ec(:), ee(:)
  integer, intent(in) :: left, right
  integer :: i, j, pivot
  integer(pin) :: pa, pb

  if (left >= right) return
  i = left
  j = right
  pivot = (left + right) / 2
  pa = ea(pivot)
  pb = eb(pivot)

  do
     do while (edge_less(ea(i), eb(i), pa, pb))
        i = i + 1
     enddo
     do while (edge_less(pa, pb, ea(j), eb(j)))
        j = j - 1
     enddo
     if (i <= j) then
        call swap_edge(ea, eb, ec, ee, i, j)
        i = i + 1
        j = j - 1
     endif
     if (i > j) exit
  enddo

  if (left < j) call sort_edges(ea, eb, ec, ee, left, j)
  if (i < right) call sort_edges(ea, eb, ec, ee, i, right)
end subroutine sort_edges

logical function edge_less(a1, b1, a2, b2)
  integer(pin), intent(in) :: a1, b1, a2, b2

  edge_less = a1 < a2 .or. (a1 == a2 .and. b1 < b2)
end function edge_less

subroutine swap_edge(ea, eb, ec, ee, i, j)
  integer(pin), intent(inout) :: ea(:), eb(:), ec(:), ee(:)
  integer, intent(in) :: i, j
  integer(pin) :: ta, tb, tc, te

  if (i == j) return
  ta = ea(i); tb = eb(i); tc = ec(i); te = ee(i)
  ea(i) = ea(j); eb(i) = eb(j); ec(i) = ec(j); ee(i) = ee(j)
  ea(j) = ta; eb(j) = tb; ec(j) = tc; ee(j) = te
end subroutine swap_edge

end module sgrac_mask_border
