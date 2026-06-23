program sgrac_gmsh_support
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use, intrinsic :: iso_c_binding, only: c_int, c_size_t, c_double
  use forparse
  use gmsh
  implicit none

  integer, parameter :: pr = selected_real_kind(12, 100)
  integer, parameter :: TRI3 = 2
  integer, parameter :: GMSH_ALGORITHM = 5

  real :: lx_in, lz_in, lc_in, x0_in, y0_in, z0_in
  real(pr) :: lx, lz, lc, x0, y0, z0
  integer :: ierr
  character(256) :: outname
  logical :: use_file
  integer :: unit_out

  integer(c_int) :: ierr_g
  type(gmsh_t) :: gmsh_api
  integer :: pc, p1, p2, p3, p4, l1, l2, l3, l4, loop_tag, surf_tag
  integer(c_int), allocatable :: element_types(:)
  integer(c_size_t), allocatable :: node_tags(:), elem_tags(:), elem_tags_n(:)
  integer(c_size_t), allocatable :: elem_node_tags(:), elem_node_tags_n(:)
  real(c_double), allocatable :: coord(:), param_coord(:)
  integer, allocatable :: tag2local(:)
  integer :: nnode, ntri, max_tag, itype, off, tri_off, i, center_node
  integer, allocatable :: tri(:,:), tri_used(:,:)
  real(pr), allocatable :: xyz(:,:), xyz_used(:,:)

  if (need_help()) then
     call print_help()
     stop
  endif

  lx_in = 120000.0
  lz_in = 80000.0
  lc_in = 500.0
  x0_in = 0.0
  y0_in = 0.0
  z0_in = 0.0
  outname = ''

  ierr = parse_arg('lx', lx_in)
  ierr = parse_arg('lz', lz_in)
  ierr = parse_arg('lc', lc_in)
  ierr = parse_arg('x0', x0_in)
  ierr = parse_arg('y0', y0_in)
  ierr = parse_arg('z0', z0_in)
  ierr = parse_arg('out', outname)

  lx = real(lx_in, pr)
  lz = real(lz_in, pr)
  lc = real(lc_in, pr)
  x0 = real(x0_in, pr)
  y0 = real(y0_in, pr)
  z0 = real(z0_in, pr)

  if (lx <= 0.0_pr .or. lz <= 0.0_pr .or. lc <= 0.0_pr) then
     write(error_unit,*) 'ERROR: lx, lz and lc must be positive S.I. lengths in meters.'
     stop 1
  endif

  use_file = len_trim(outname) > 0
  if (use_file) then
     open(newunit=unit_out, file=trim(outname), status='replace', action='write', iostat=ierr)
     if (ierr /= 0) then
        write(error_unit,*) 'ERROR: cannot open output file: ', trim(outname)
        stop 1
     endif
  else
     unit_out = output_unit
  endif

  call gmsh_api%initialize(readConfigFiles=.false., run=.false., ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmshInitialize')
  call gmsh_api%option%setNumber('General.Terminal', 0.0_c_double, ierr_g)
  call gmsh_api%option%setNumber('Mesh.Algorithm', real(GMSH_ALGORITHM, c_double), ierr_g)
  call gmsh_api%option%setNumber('Mesh.CharacteristicLengthMin', real(lc, c_double), ierr_g)
  call gmsh_api%option%setNumber('Mesh.CharacteristicLengthMax', real(lc, c_double), ierr_g)
  call gmsh_api%option%setNumber('Mesh.SaveAll', 1.0_c_double, ierr_g)

  call gmsh_api%model%add('sgrac_parent_support', ierr_g)
  call check_gmsh(ierr_g, 'gmshModelAdd')

  pc = gmsh_api%model%geo%addPoint(real(x0, c_double), real(y0, c_double), real(z0, c_double), &
                            real(lc, c_double), ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addPoint pc')
  p1 = gmsh_api%model%geo%addPoint(real(x0 - 0.5_pr*lx, c_double), real(y0, c_double), real(z0 - 0.5_pr*lz, c_double), &
                            real(lc, c_double), ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addPoint p1')
  p2 = gmsh_api%model%geo%addPoint(real(x0 + 0.5_pr*lx, c_double), real(y0, c_double), real(z0 - 0.5_pr*lz, c_double), &
                            real(lc, c_double), ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addPoint p2')
  p3 = gmsh_api%model%geo%addPoint(real(x0 + 0.5_pr*lx, c_double), real(y0, c_double), real(z0 + 0.5_pr*lz, c_double), &
                            real(lc, c_double), ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addPoint p3')
  p4 = gmsh_api%model%geo%addPoint(real(x0 - 0.5_pr*lx, c_double), real(y0, c_double), real(z0 + 0.5_pr*lz, c_double), &
                            real(lc, c_double), ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addPoint p4')

  l1 = gmsh_api%model%geo%addLine(p1, p2, ierr=ierr_g)
  l2 = gmsh_api%model%geo%addLine(p2, p3, ierr=ierr_g)
  l3 = gmsh_api%model%geo%addLine(p3, p4, ierr=ierr_g)
  l4 = gmsh_api%model%geo%addLine(p4, p1, ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addLine')

  loop_tag = gmsh_api%model%geo%addCurveLoop([int(l1,c_int), int(l2,c_int), int(l3,c_int), int(l4,c_int)], ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addCurveLoop')
  surf_tag = gmsh_api%model%geo%addPlaneSurface([int(loop_tag,c_int)], ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmsh_api%model%geo%addPlaneSurface')

  call gmsh_api%model%geo%synchronize(ierr_g)
  call check_gmsh(ierr_g, 'gmshModelGeoSynchronize')
  call gmsh_api%model%mesh%embed(0, [int(pc,c_int)], 2, surf_tag, ierr_g)
  call check_gmsh(ierr_g, 'gmshModelMeshEmbed center')
  call gmsh_api%model%mesh%generate(2, ierr_g)
  call check_gmsh(ierr_g, 'gmshModelMeshGenerate')

  call gmsh_api%model%mesh%getNodes(node_tags, coord, param_coord, dim=-1, tag=-1, &
                             includeBoundary=.true., returnParametricCoord=.false., ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmshModelMeshGetNodes')
  call gmsh_api%model%mesh%getElements(element_types, elem_tags, elem_tags_n, elem_node_tags, elem_node_tags_n, &
                                dim=2, tag=surf_tag, ierr=ierr_g)
  call check_gmsh(ierr_g, 'gmshModelMeshGetElements')

  nnode = size(node_tags)
  if (nnode <= 0) then
     write(error_unit,*) 'ERROR: gmsh returned no nodes.'
     call gmsh_api%finalize(ierr_g)
     stop 1
  endif

  max_tag = int(maxval(node_tags))
  allocate(tag2local(max_tag)); tag2local = 0
  allocate(xyz(3,nnode))
  do i = 1, nnode
     tag2local(int(node_tags(i))) = i
     xyz(1,i) = real(coord(3*i-2), pr)
     xyz(2,i) = real(coord(3*i-1), pr)
     xyz(3,i) = real(coord(3*i  ), pr)
  enddo
  center_node = find_node_by_coordinate(xyz, x0, y0, z0, 1.e-9_pr*max(lc, 1._pr))

  ntri = 0
  do itype = 1, size(element_types)
     if (element_types(itype) == TRI3) ntri = int(elem_node_tags_n(itype)) / 3
  enddo
  if (ntri <= 0) then
     write(error_unit,*) 'ERROR: gmsh returned no 3-node triangles on the support surface.'
     call gmsh_api%finalize(ierr_g)
     stop 1
  endif

  allocate(tri(3,ntri))
  off = 1
  tri_off = 0
  do itype = 1, size(element_types)
     if (element_types(itype) == TRI3) then
        do i = 1, int(elem_node_tags_n(itype)) / 3
           tri_off = tri_off + 1
           tri(1,tri_off) = tag2local(int(elem_node_tags(off + 3*i - 3)))
           tri(2,tri_off) = tag2local(int(elem_node_tags(off + 3*i - 2)))
           tri(3,tri_off) = tag2local(int(elem_node_tags(off + 3*i - 1)))
        enddo
     endif
     off = off + int(elem_node_tags_n(itype))
  enddo

  call compact_used_nodes(xyz, tri, center_node, xyz_used, tri_used)
  call write_vtk_polydata(unit_out, xyz_used, tri_used, lx, lz, lc, GMSH_ALGORITHM)

  if (use_file) close(unit_out)
  call gmsh_api%finalize(ierr_g)

contains

  logical function need_help()
    need_help = helpisneeded('help')
    if (.not. need_help) need_help = helpisneeded('-h')
    if (.not. need_help) need_help = helpisneeded('--help')
  end function need_help

  subroutine check_gmsh(ierr, where)
    integer(c_int), intent(in) :: ierr
    character(len=*), intent(in) :: where
    if (ierr /= 0_c_int) then
       write(error_unit,*) 'ERROR in ', trim(where), ': gmsh ierr=', ierr
       call gmsh_api%finalize()
       stop 1
    endif
  end subroutine check_gmsh

  subroutine write_vtk_polydata(unit, xyz, tri, lx, lz, lc, algorithm)
    integer, intent(in) :: unit
    real(pr), intent(in) :: xyz(:,:)
    integer, intent(in) :: tri(:,:)
    real(pr), intent(in) :: lx, lz, lc
    integer, intent(in) :: algorithm
    integer :: i, nnode, ntri
    nnode = size(xyz,2)
    ntri = size(tri,2)

    write(unit,'(a)') '# vtk DataFile Version 3.0'
    write(unit,'(a)') 'SGRAC gmsh parent support mesh on x-z fault plane; units: SI meters'
    write(unit,'(a)') 'ASCII'
    write(unit,'(a)') 'DATASET POLYDATA'
    write(unit,'(a,1x,i0,1x,a)') 'POINTS', nnode, 'double'
    do i = 1, nnode
       write(unit,'(3(es24.16,1x))') xyz(1,i), xyz(2,i), xyz(3,i)
    enddo
    write(unit,'(a,1x,i0,1x,i0)') 'POLYGONS', ntri, 4*ntri
    do i = 1, ntri
       write(unit,'(i1,1x,i0,1x,i0,1x,i0)') 3, tri(1,i)-1, tri(2,i)-1, tri(3,i)-1
    enddo
    write(unit,'(a,1x,i0)') 'FIELD FieldData', 4
    write(unit,'(a)') 'lx_m 1 1 double'
    write(unit,'(es24.16)') lx
    write(unit,'(a)') 'lz_m 1 1 double'
    write(unit,'(es24.16)') lz
    write(unit,'(a)') 'lc_m 1 1 double'
    write(unit,'(es24.16)') lc
    write(unit,'(a)') 'gmsh_algorithm 1 1 int'
    write(unit,'(i0)') algorithm
  end subroutine write_vtk_polydata

  subroutine compact_used_nodes(xyz_in, tri_in, first_node, xyz_out, tri_out)
    real(pr), intent(in) :: xyz_in(:,:)
    integer, intent(in) :: tri_in(:,:)
    integer, intent(in) :: first_node
    real(pr), allocatable, intent(out) :: xyz_out(:,:)
    integer, allocatable, intent(out) :: tri_out(:,:)
    logical, allocatable :: used(:)
    integer, allocatable :: old2new(:)
    integer :: i, j, inode, nused

    allocate(used(size(xyz_in,2)))
    allocate(old2new(size(xyz_in,2)))
    used = .false.
    old2new = 0

    do i = 1, size(tri_in,2)
       do j = 1, 3
          inode = tri_in(j,i)
          if (inode < 1 .or. inode > size(xyz_in,2)) then
             write(error_unit,*) 'ERROR: triangle references node outside the coordinate array.'
             stop 1
          endif
          used(inode) = .true.
       enddo
    enddo

    nused = count(used)
    if (first_node < 1 .or. first_node > size(xyz_in,2) .or. .not. used(first_node)) then
       write(error_unit,*) 'ERROR: center point is not used by any triangle.'
       stop 1
    endif
    allocate(xyz_out(3,nused))
    allocate(tri_out(3,size(tri_in,2)))

    nused = 1
    old2new(first_node) = 1
    xyz_out(:,1) = xyz_in(:,first_node)
    do i = 1, size(xyz_in,2)
       if (used(i) .and. i /= first_node) then
          nused = nused + 1
          old2new(i) = nused
          xyz_out(:,nused) = xyz_in(:,i)
       endif
    enddo

    do i = 1, size(tri_in,2)
       do j = 1, 3
          tri_out(j,i) = old2new(tri_in(j,i))
       enddo
    enddo
  end subroutine compact_used_nodes

  integer function find_node_by_coordinate(xyz, x, y, z, tol) result(inode)
    real(pr), intent(in) :: xyz(:,:)
    real(pr), intent(in) :: x, y, z, tol
    real(pr) :: d2, d2_min
    integer :: i

    inode = 0
    d2_min = huge(1._pr)
    do i = 1, size(xyz,2)
       d2 = (xyz(1,i)-x)**2 + (xyz(2,i)-y)**2 + (xyz(3,i)-z)**2
       if (d2 < d2_min) then
          d2_min = d2
          inode = i
       endif
    enddo

    if (inode == 0 .or. d2_min > tol**2) then
       write(error_unit,*) 'ERROR: gmsh did not return the embedded center node.'
       stop 1
    endif
  end function find_node_by_coordinate

  subroutine print_help()
    write(output_unit,'(a)') 'sgrac-gmsh-support: generate a triangular parent support mesh with gmsh.'
    write(output_unit,'(a)') ''
    write(output_unit,'(a)') 'Usage:'
    write(output_unit,'(a)') '  sgrac-gmsh-support lx=120000 lz=80000 lc=500 > parent.vtk'
    write(output_unit,'(a)') '  sgrac-gmsh-support lx=120000 lz=80000 lc=500 out=parent.vtk'
    write(output_unit,'(a)') ''
    write(output_unit,'(a)') 'Parameters, all S.I. units:'
    write(output_unit,'(a)') '  lx        support length in meters      default 120000'
    write(output_unit,'(a)') '  lz        support depth span in meters  default 80000'
    write(output_unit,'(a)') '  lc        target gmsh size in meters    default 500'
    write(output_unit,'(a)') '  x0,y0,z0  support center in meters      default 0,0,0'
    write(output_unit,'(a)') '  out       optional output VTK file; stdout by default'
    write(output_unit,'(a)') ''
    write(output_unit,'(a)') 'Output: VTK legacy POLYDATA, triangles only, no CELL_TYPES.'
  end subroutine print_help

end program sgrac_gmsh_support
