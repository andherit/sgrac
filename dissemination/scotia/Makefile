FC = gfortran
FFLAGS ?= -O2 -Wall -Wextra -fcheck=all -fbacktrace
SRC = src
TRILAT_SRC ?= /home/herrero/Atelier/mesh/trilat-time/src
OBJ_DIR = build/obj
MOD_DIR = build/mod
MOD_DIR_ABS = $(CURDIR)/$(MOD_DIR)

OBJS = \
	$(OBJ_DIR)/generic.o \
	$(OBJ_DIR)/lists.o \
	$(OBJ_DIR)/LAT_mesh.o \
	$(OBJ_DIR)/LAT_time.o \
	$(OBJ_DIR)/time_debug.o \
	$(OBJ_DIR)/time.o \
	$(OBJ_DIR)/scotia_vtk_io.o \
	$(OBJ_DIR)/scotia_rupture.o

all: scotia-rupture-times

scotia-rupture-times: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

$(OBJ_DIR)/generic.o: $(TRILAT_SRC)/generic.f90 | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $< -o generic.o

$(OBJ_DIR)/lists.o: $(TRILAT_SRC)/lists.f90 $(OBJ_DIR)/generic.o | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $< -o lists.o

$(OBJ_DIR)/LAT_mesh.o: $(TRILAT_SRC)/LAT_mesh.f90 $(OBJ_DIR)/generic.o $(OBJ_DIR)/lists.o | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $< -o LAT_mesh.o

$(OBJ_DIR)/LAT_time.o: $(TRILAT_SRC)/LAT_time.f90 $(OBJ_DIR)/generic.o $(OBJ_DIR)/LAT_mesh.o | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $< -o LAT_time.o

$(OBJ_DIR)/time_debug.o: $(TRILAT_SRC)/time_debug.f90 $(OBJ_DIR)/generic.o | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $< -o time_debug.o

$(OBJ_DIR)/time.o: $(TRILAT_SRC)/time.f90 $(OBJ_DIR)/generic.o $(OBJ_DIR)/lists.o $(OBJ_DIR)/LAT_mesh.o $(OBJ_DIR)/LAT_time.o $(OBJ_DIR)/time_debug.o | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $< -o time.o

$(OBJ_DIR)/scotia_vtk_io.o: $(SRC)/scotia_vtk_io.f90 $(OBJ_DIR)/generic.o $(OBJ_DIR)/LAT_mesh.o | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $(CURDIR)/$< -o scotia_vtk_io.o

$(OBJ_DIR)/scotia_rupture.o: $(SRC)/scotia_rupture.f90 $(OBJ_DIR)/generic.o $(OBJ_DIR)/lists.o $(OBJ_DIR)/LAT_mesh.o $(OBJ_DIR)/LAT_time.o $(OBJ_DIR)/time.o $(OBJ_DIR)/scotia_vtk_io.o | $(OBJ_DIR) $(MOD_DIR)
	cd $(OBJ_DIR) && $(FC) $(FFLAGS) -J$(MOD_DIR_ABS) -I$(MOD_DIR_ABS) -c $(CURDIR)/$< -o scotia_rupture.o

run: scotia-rupture-times
	./scotia-rupture-times input.vtk rupture_times.vtk 908

.PHONY: all run clean

clean:
	rm -rf build scotia-rupture-times
	rm -f *.o *.mod
