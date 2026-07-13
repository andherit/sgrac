FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src

OBJS = generic.o sgrac_mesh_field_support.o sgrac_slip_support.o sgrac_vtk_field_io.o

all: $(OBJS)

sgrac_mesh_field_support.o: $(SRC)/sgrac_mesh_field_support.f90 generic.o
	$(FC) $(FFLAGS) -c $<

sgrac_slip_support.o: $(SRC)/sgrac_slip_support.f90 generic.o
	$(FC) $(FFLAGS) -c $<

sgrac_vtk_field_io.o: $(SRC)/sgrac_vtk_field_io.f90 generic.o
	$(FC) $(FFLAGS) -c $<

generic.o: $(SRC)/generic.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod
