FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src
SUPPORT_SRC = ../sgrac_slip_support/src

OBJS = generic.o lists.o LAT_mesh.o LAT_distance.o distance.o forparse.o \
       sgrac_slip_support.o sgrac_vtk_slip_io.o sgrac_slip_smooth.o

all: sgrac-slip-smooth

sgrac-slip-smooth: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

generic.o: $(SRC)/generic.f90
	$(FC) $(FFLAGS) -c $<

forparse.o: $(SRC)/forparse.f90
	$(FC) $(FFLAGS) -c $<

lists.o: $(SRC)/lists.f90 generic.o
	$(FC) $(FFLAGS) -c $<

LAT_mesh.o: $(SRC)/LAT_mesh.f90 generic.o lists.o
	$(FC) $(FFLAGS) -c $<

LAT_distance.o: $(SRC)/LAT_distance.f90 generic.o LAT_mesh.o
	$(FC) $(FFLAGS) -c $<

distance.o: $(SRC)/distance.f90 generic.o lists.o LAT_mesh.o LAT_distance.o
	$(FC) $(FFLAGS) -c $<

sgrac_slip_support.o: $(SUPPORT_SRC)/sgrac_slip_support.f90 generic.o
	$(FC) $(FFLAGS) -c $<

sgrac_vtk_slip_io.o: $(SRC)/sgrac_vtk_slip_io.f90 generic.o
	$(FC) $(FFLAGS) -c $<

sgrac_slip_smooth.o: $(SRC)/sgrac_slip_smooth.f90 generic.o forparse.o LAT_mesh.o LAT_distance.o lists.o distance.o sgrac_slip_support.o sgrac_vtk_slip_io.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod sgrac-slip-smooth
