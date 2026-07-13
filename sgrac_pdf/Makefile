FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src
SLIP_SRC = ../sgrac_slip_smooth/src
SUPPORT_SRC = ../sgrac_slip_support/src
BIN_DIR ?= ../bin
EXE = $(BIN_DIR)/sgrac-pdf

OBJS = generic.o lists.o LAT_mesh.o LAT_distance.o distance.o forparse.o \
       sgrac_mesh_field_support.o sgrac_vtk_field_io.o sgrac_pdf.o

all: $(EXE)

$(EXE): $(OBJS) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

generic.o: $(SLIP_SRC)/generic.f90
	$(FC) $(FFLAGS) -c $<

forparse.o: $(SLIP_SRC)/forparse.f90
	$(FC) $(FFLAGS) -c $<

lists.o: $(SLIP_SRC)/lists.f90 generic.o
	$(FC) $(FFLAGS) -c $<

LAT_mesh.o: $(SLIP_SRC)/LAT_mesh.f90 generic.o lists.o
	$(FC) $(FFLAGS) -c $<

LAT_distance.o: $(SLIP_SRC)/LAT_distance.f90 generic.o LAT_mesh.o
	$(FC) $(FFLAGS) -c $<

distance.o: $(SLIP_SRC)/distance.f90 generic.o lists.o LAT_mesh.o LAT_distance.o
	$(FC) $(FFLAGS) -c $<

sgrac_mesh_field_support.o: $(SUPPORT_SRC)/sgrac_mesh_field_support.f90 generic.o
	$(FC) $(FFLAGS) -c $<

sgrac_vtk_field_io.o: $(SUPPORT_SRC)/sgrac_vtk_field_io.f90 generic.o
	$(FC) $(FFLAGS) -c $<

sgrac_pdf.o: $(SRC)/sgrac_pdf.f90 generic.o forparse.o LAT_mesh.o LAT_distance.o lists.o distance.o sgrac_mesh_field_support.o sgrac_vtk_field_io.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXE)
