FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src
BIN_DIR ?= ../bin
EXE = $(BIN_DIR)/sgrac-mask

OBJS = generic.o forparse.o sgrac_vtk_mask_io.o sgrac_radius_models.o sgrac_mask_border.o sgrac_mask.o

all: $(EXE)

$(EXE): $(OBJS) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

sgrac_mask_border.o: $(SRC)/sgrac_mask_border.f90 generic.o
	$(FC) $(FFLAGS) -c $<

sgrac_mask.o: $(SRC)/sgrac_mask.f90 generic.o forparse.o sgrac_vtk_mask_io.o sgrac_radius_models.o sgrac_mask_border.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXE)
