FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all
SRC = src

OBJS = generic.o forparse.o sgrac_vtk_mask_io.o sgrac_radius_models.o sgrac_mask.o

all: sgrac-mask

sgrac-mask: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod sgrac-mask
