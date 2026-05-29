FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src

OBJS = generic.o forparse.o sgrac_vtk_extract_io.o sgrac_extract.o

all: sgrac-extract

sgrac-extract: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod sgrac-extract
