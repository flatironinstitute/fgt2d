# Makefile for fgt2d
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
#FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy
# -pg -no-pie is for profiling
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy -pg -no-pie

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 
OMP = OFF

LBLAS = -lblas -llapack

FGT_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FGT_INSTALL_DIR = ${HOME}/lib
endif

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

LIBNAME=libfgt2d
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LFGTLINKLIB = -lfgt2d
LLINKLIB = -lfgt2d


# For your OS, override the above by placing make variables in make.inc
-include make.inc

# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
  FFLAGS += $(OMPFLAGS)
  LIBS += $(OMPLIBS)
  DYLIBS += $(OMPLIBS)
  F2PYDYLIBS += $(OMPLIBS)
endif

LIBS += $(LBLAS) $(LDBLASINC)
DYLIBS += $(LBLAS) $(LDBLASINC)

#
# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/prini_new.o \
	$(COM)/sort.o \
	$(COM)/sparse_reps.o \
	$(COM)/cumsum.o \
	$(COM)/fmmcommon2d.o \
	$(COM)/pts_tree2d2.o \
	$(COM)/tree_routs4.o 

# Gauss transform objects
GAUSS = src
GOBJS = $(GAUSS)/fgt2d.o \
	$(GAUSS)/g2ddirect.o \
	$(GAUSS)/g2dhlall.o \
	$(GAUSS)/g2dsoeall.o \
	$(GAUSS)/g2dsxall.o 


# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

OBJS = $(COMOBJS) $(GOBJS) 



.PHONY: usage lib install test test-dyn python 

default: usage

usage:
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for fgt2d. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests"
	@echo "  make test-dyn - test successful installation by validation tests linked to dynamic library"
	@echo "  make python - compile and test python interfaces using python"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"
	@echo "-------------------------------------------------------------------------"

#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

#
# build the library...
#
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif

$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/

$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(FFLAGS) $(OBJS) -o $(DYNAMICLIB)
	mv $(DYNAMICLIB) lib/

install: $(STATICLIB) $(DYNAMICLIB)
	echo $(FGT_INSTALL_DIR)
	mkdir -p $(FGT_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(FGT_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(FGT_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(FGT_INSTALL_DIR)/
	@echo "Make sure to include " $(FGT_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FGT_INSTALL_DIR)  " -lfgt2d"


#
# testing routines
#
test: $(STATICLIB) $(TOBJS) test/fgt2d 
	(cd test; ./run_fgttest.sh)


test-dyn: $(DYNAMICLIB)  $(TOBJS) test/fgt2d-dyn 
	(cd test; ./run_fgttest.sh)

test/fgt2d:
	$(FC) $(FFLAGS) test/test_fgt2d.f $(TOBJS) -o test/int2-test-fgt2d -Llib-static $(LLINKLIB)  


#
# Linking test files to dynamic libraries
#

test/fgt2d-dyn:
	$(FC) $(FFLAGS) test/test_fgt2d.f $(TOBJS) -o test/int2-test-fgt2d -L$(FGT_INSTALL_DIR)  $(LLINKLIB)  


#
# build the python bindings/interface
#
python: $(STATICLIB)
	cd python && export FGT2D_LIBS='$(LIBS)' && pip install -e . 

#
# housekeeping routines
#
clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f test/int2-*

objclean: 
	rm -f $(OBJS) $(TOBJS)
	rm -f test/*.o 
