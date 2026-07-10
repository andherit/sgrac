FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src
BIN_DIR ?= ../bin
EXE = $(BIN_DIR)/sgrac-geometry

OBJS = generic.o lists.o LAT_mesh.o LAT_distance.o distance.o forparse.o sgrac_vtk_tri_io.o sgrac_cell_geometry.o sgrac_geometry.o

all: $(EXE)

$(EXE): $(OBJS) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXE)
