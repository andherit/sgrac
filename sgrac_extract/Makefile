FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src
BIN_DIR ?= ../bin
EXE = $(BIN_DIR)/sgrac-extract

OBJS = generic.o forparse.o sgrac_vtk_extract_io.o sgrac_extract.o

all: $(EXE)

$(EXE): $(OBJS) | $(BIN_DIR)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXE)
