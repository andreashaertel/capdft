# SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@physik.uni-freiburg.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Binary name ##################################################################
BNAME = main
# Compiler #####################################################################
CXX=g++
# Libraries to be linked #######################################################
LIBS=-lm -lfftw3 -lgsl -lgslcblas
# Compiler flags ###############################################################
CXXFLAGS=-pg -pedantic -std=c++2a -O2 -Wall
# Directories for binaries, objects and sources ################################
BDIR = bin/
ODIR = obj/
SDIR = src/
# All files from the normal repository that are used in this example ###########
REPOSDIR = ../../src/
REPOFILE =
REPOFILE += banded_matrix
REPOFILE += convergence_criterion
REPOFILE += convergence_criterion_max_dev
REPOFILE += convergence_criterion_steps
REPOFILE += convergence_criterion_nan
REPOFILE += data_frame
REPOFILE += functional
REPOFILE += functional_fmt_planar
REPOFILE += functional_es_mf_planar
REPOFILE += functional_es_delta_planar
REPOFILE += integration
REPOFILE += iterator
REPOFILE += properties
REPOFILE += planar_poisson_solver
REPOFILE += planar_poisson_solver_any
# ParameterHandler files #######################################################
REPOSDIR2 = ../../../parameter_handler/src/
REPOFILE2 = parameter_handler
# Get host name ################################################################
HOST = $(shell hostname)
# Offsets for binaries and objects #############################################
BOFF = $(BDIR)$(HOST).
OOFF = $(ODIR)$(HOST).
# Create empty file collection variables #######################################
SOURCE =
HEADERS =
OBJECTS =
REPOSOURCE =
REPOHEADERS =
REPOOBJECTS =
# Header files and object files from main repo #################################
REPOSOURCE += $(addprefix $(REPOSDIR), $(addsuffix .cpp, $(REPOFILE)))
REPOSOURCE2 = $(addprefix $(REPOSDIR2), $(addsuffix .cpp, $(REPOFILE2)))
REPOHEADERS += $(addprefix $(REPOSDIR), $(addsuffix .hpp, $(REPOFILE)))
REPOHEADERS2 = $(addprefix $(REPOSDIR2), $(addsuffix .hpp, $(REPOFILE2)))
REPOOBJECTS += $(patsubst $(REPOSDIR)%.cpp,$(OOFF)%.o, $(REPOSOURCE))
REPOOBJECTS2 = $(patsubst $(REPOSDIR2)%.cpp,$(OOFF)%.o, $(REPOSOURCE2))
# Standard header files and object files #######################################
HEADERS += $(wildcard $(SDIR)*.h)
HEADERS += $(wildcard $(SDIR)*.hpp)
SOURCE += $(wildcard $(SDIR)*.cpp)
OBJECTS += $(patsubst $(SDIR)%.cpp,$(OOFF)%.o,$(SOURCE))
# Get host name ################################################################
HOST = $(shell hostname)
# Offsets for binaries and objects #############################################
BOFF = $(BDIR)$(HOST).
OOFF = $(ODIR)$(HOST).
# TARGETS ######################################################################

.PHONY: default all init todo checkstyle clean compile

default: init compile

all: checkstyle init compile todo

init:
	@mkdir -p bin
	@mkdir -p obj

todo:
	@echo "================================ TODO ================================"
	@grep -n "TODO[(]*.*[)]*[:]*.*" $(SOURCE) $(HEADERS) $(REPOSOURCE) $(REPOHEADERS)


# Google C++ style checker:
# Install cpplint with
# $pip3 install cpplint
# And add ~/.local/bin to PATH variable.
# Checks all files in the src directory.
checkstyle:
	@cpplint $(SOURCE) $(HEADERS) $(REPOSOURCE) $(REPOHEADERS)

clean:
	rm -f $(OOFF)*.o
	rm -f $(BOFF)*

compile: $(BOFF)$(BNAME)

# $<  first dependency
# $@  targets
# $+  list of all dependencies
# $^  list of all dependencies, no doubles

$(BOFF)$(BNAME): $(OBJECTS) $(REPOOBJECTS) $(REPOOBJECTS2)
	$(CXX) -pg -o $@ $^ $(LIBS)

$(OOFF)%.o: $(SDIR)%.cpp $(HEADERS) $(REPOHEADERS) $(REPOHEADERS2)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
$(REPOOBJECTS): $(REPOSOURCE) $(HEADERS) $(REPOHEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $(REPOSDIR)$(subst .,,$(suffix $(basename $@))).cpp
$(REPOOBJECTS2): $(REPOSOURCE2) $(HEADER2)
	$(CXX) $(CXXFLAGS) -c -o $@ $(REPOSDIR2)$(subst .,,$(suffix $(basename $@))).cpp

