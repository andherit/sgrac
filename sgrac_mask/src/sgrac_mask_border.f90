module sgrac_mask_border
  use, intrinsic :: iso_fortran_env, only: error_unit
  use generic
  implicit none
contains

subroutine smooth_mask_border(lines, nlines, ncell, phi, mask, nremove, nadd, nadd_candidates)
  character(len=*), intent(in) :: lines(:)
  integer, intent(in) :: nlines, ncell
  real(pr), intent(in) :: phi(ncell)
  integer, intent(inout) :: mask(ncell)
  integer, intent(out) :: nremove, nadd, nadd_candidates
  integer(pin), allocatable :: cell(:,:), ea(:), eb(:), ec(:), ee(:)
  integer, allocatable :: border_count(:), linked_count(:), add_candidates(:)
  logical, allocatable :: remove_flag(:), used_candidate(:)
  integer :: i, j, k, nedge, group_start, group_end, group_size
  integer :: c1, c2, best_idx, best_cell
  real(pr) :: best_phi

  call read_triangles(lines, nlines, ncell, cell)

  nedge = 3*ncell
  allocate(ea(nedge), eb(nedge), ec(nedge), ee(nedge))
  k = 0
  do i=1,ncell
     call add_edge(cell(i,1), cell(i,2), i, 1, ea, eb, ec, ee, k)
     call add_edge(cell(i,2), cell(i,3), i, 2, ea, eb, ec, ee, k)
     call add_edge(cell(i,3), cell(i,1), i, 3, ea, eb, ec, ee, k)
  enddo

  call sort_edges(ea, eb, ec, ee, 1, nedge)

  allocate(border_count(ncell), linked_count(ncell), remove_flag(ncell))
  border_count = 0
  linked_count = 0
  remove_flag = .false.

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

  do i=1,ncell
     remove_flag(i) = mask(i) == 1 .and. border_count(i) == 2
  enddo
  nremove = count(remove_flag)
  nadd_candidates = count(mask == 0 .and. linked_count == 2)
  nadd = min(nremove, nadd_candidates)

  do i=1,ncell
     if (remove_flag(i)) mask(i) = 0
  enddo

  if (nadd <= 0) return

  allocate(add_candidates(nadd_candidates), used_candidate(nadd_candidates))
  k = 0
  do i=1,ncell
     if (mask(i) == 0 .and. linked_count(i) == 2 .and. .not. remove_flag(i)) then
        k = k + 1
        add_candidates(k) = i
     endif
  enddo
  used_candidate = .false.

  do j=1,nadd
     best_idx = 0
     best_cell = 0
     best_phi = huge(1._pr)
     do i=1,nadd_candidates
        if (used_candidate(i)) cycle
        if (phi(add_candidates(i)) < best_phi) then
           best_phi = phi(add_candidates(i))
           best_idx = i
           best_cell = add_candidates(i)
        endif
     enddo
     if (best_idx <= 0) exit
     mask(best_cell) = 1
     used_candidate(best_idx) = .true.
  enddo
end subroutine smooth_mask_border

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
