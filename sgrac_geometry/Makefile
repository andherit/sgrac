FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src

OBJS = generic.o lists.o LAT_mesh.o LAT_distance.o distance.o forparse.o sgrac_vtk_tri_io.o sgrac_cell_geometry.o sgrac_geometry.o

all: sgrac-geometry

sgrac-geometry: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod sgrac-geometry
