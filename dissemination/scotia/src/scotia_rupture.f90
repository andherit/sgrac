program scotia_rupture
  use generic
  use lists
  use LAT_mesh
  use LAT_time
  use time
  use scotia_vtk_io
  implicit none

  type(mesh) :: amesh
  type(diff) :: adiff
  character(len=256) :: input_file, output_file, arg
  integer(pin) :: source_vtk_id, source_node
  real(pr), allocatable :: slip(:)
  real(pr), allocatable :: v_const(:), v_prop(:), v_inv(:)
  real(pr), allocatable :: t_const(:), t_prop(:), t_inv(:)
  real(pr), parameter :: vr_const = 3000._pr
  real(pr), parameter :: vr_min = 1500._pr
  real(pr), parameter :: vr_max = 4500._pr

  input_file = 'input.vtk'
  output_file = 'rupture_times.vtk'
  source_vtk_id = 908_pin

  if (command_argument_count() >= 1) call get_command_argument(1, input_file)
  if (command_argument_count() >= 2) call get_command_argument(2, output_file)
  if (command_argument_count() >= 3) then
     call get_command_argument(3, arg)
     read(arg,*) source_vtk_id
  endif

  call read_scotia_vtk(input_file, amesh, slip)

  source_node = source_vtk_id + 1_pin
  if (source_node < 1_pin .or. source_node > amesh%Nnodes) then
     write(0,*) 'scotia: source VTK id out of range: ', source_vtk_id
     stop 1
  endif

  adiff%fast = .true.

  allocate(v_const(amesh%Ncells), v_prop(amesh%Ncells), v_inv(amesh%Ncells))
  allocate(t_const(amesh%Nnodes), t_prop(amesh%Nnodes), t_inv(amesh%Nnodes))

  v_const = vr_const
  call slip_velocity(slip, vr_min, vr_max, .false., v_prop)
  call slip_velocity(slip, vr_min, vr_max, .true., v_inv)

  write(*,*) 'scotia: source vtk id/node = ', source_vtk_id, source_node
  write(*,*) 'scotia: constant velocity case'
  call compute_case(amesh, adiff, source_node, v_const, t_const)
  write(*,*) 'scotia: slip-proportional capped velocity case'
  call compute_case(amesh, adiff, source_node, v_prop, t_prop)
  write(*,*) 'scotia: slip-inverse capped velocity case'
  call compute_case(amesh, adiff, source_node, v_inv, t_inv)

  call write_scotia_vtk(output_file, amesh, slip, t_const, t_prop, t_inv)
  write(*,*) 'scotia: wrote ', trim(output_file)

contains

subroutine compute_case(amesh, adiff, source_node, v_case, t_out)
  type(mesh), intent(in) :: amesh
  type(diff), intent(inout) :: adiff
  integer(pin), intent(in) :: source_node
  real(pr), intent(in) :: v_case(amesh%Ncells)
  real(pr), intent(out) :: t_out(amesh%Nnodes)
  real(pr), allocatable :: traveltime(:)
  type(containern), allocatable :: nton(:)

  call pre_timeonevsall2d_onvertex(amesh, source_node, traveltime, nton)
  call timeonevsall2d(amesh, v_case, traveltime, nton, adiff)
  t_out = traveltime
  call free_nton(nton)
  deallocate(nton, traveltime)
end subroutine compute_case

subroutine slip_velocity(slip, vmin, vmax, inverse, vel)
  real(pr), intent(in) :: slip(:), vmin, vmax
  logical, intent(in) :: inverse
  real(pr), intent(out) :: vel(size(slip))
  real(pr) :: smin, smax, frac
  integer :: i

  smin = minval(slip)
  smax = maxval(slip)

  do i = 1, size(slip)
     if (smax > smin) then
        frac = (slip(i) - smin) / (smax - smin)
     else
        frac = 0.5_pr
     endif
     if (inverse) frac = 1._pr - frac
     vel(i) = max(vmin, min(vmax, vmin + frac*(vmax - vmin)))
  enddo
end subroutine slip_velocity

end program scotia_rupture
