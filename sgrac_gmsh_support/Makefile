FC ?= gfortran
FFLAGS ?= -O2 -Wall -Wextra
GMSH_LIB_DIR ?= $(HOME)/software/gmsh-4.15.2-Linux64-sdk/lib
LDFLAGS ?= -L$(GMSH_LIB_DIR) -Wl,-rpath,$(GMSH_LIB_DIR)
LDLIBS ?= -lgmsh

SRC_DIR := src
OBJ_DIR := build/obj
MOD_DIR := build/mod

TARGET := sgrac-gmsh-support
OBJS := $(OBJ_DIR)/forparse.o $(OBJ_DIR)/gmsh.o $(OBJ_DIR)/sgrac_gmsh_support.o

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

$(OBJ_DIR)/forparse.o: $(SRC_DIR)/forparse.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/gmsh.o: $(SRC_DIR)/gmsh.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/sgrac_gmsh_support.o: $(SRC_DIR)/sgrac_gmsh_support.f90 $(OBJ_DIR)/forparse.o $(OBJ_DIR)/gmsh.o | $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

clean:
	rm -rf build $(TARGET)
