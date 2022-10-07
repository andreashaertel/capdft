# SPDX-FileCopyrightText: 2021 Moritz Bültmann <moritz.bueltmann@physik.uni-freiburg.de>
# SPDX-License-Identifier: LGPL-3.0-or-later
# Binary name ##################################################################
BNAME=main
# Compiler #####################################################################
CXX=g++
# Libraries to be linked #######################################################
LIBS=-lm -lfftw3 -lgsl -lgslcblas
# Compiler flags ###############################################################
CXXFLAGS=-pedantic -std=c++2a -O2 -Wall
# Directories for binaries, objects and sources ################################
BDIR = bin/
ODIR = obj/
SDIR = src/
# All files from the normal repository that are used in this example ###########
REPOSDIR = ../../src/
REPOFILE =
REPOFILE += convergence_criterion
REPOFILE += convergence_criterion_max_dev
REPOFILE += convergence_criterion_steps
REPOFILE += convergence_criterion_nan
REPOFILE += data_frame
REPOFILE += functional
REPOFILE += functional_fmt_cartesian
#REPOFILE += functional_es_mf_cartesian
#REPOFILE += functional_es_delta_cartesian
REPOFILE += iterator
REPOFILE += properties
#REPOFILE += planar_poisson_solver
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
REPOHEADERS += $(addprefix $(REPOSDIR), $(addsuffix .hpp, $(REPOFILE)))
REPOOBJECTS += $(patsubst $(REPOSDIR)%.cpp,$(OOFF)%.o, $(REPOSOURCE))
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

$(BOFF)$(BNAME): $(OBJECTS) $(REPOOBJECTS)
	$(CXX) -o $@ $^ $(LIBS)

$(OOFF)%.o: $(SDIR)%.cpp $(HEADERS) $(REPOHEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
$(REPOOBJECTS): $(REPOSOURCE) $(HEADERS) $(REPOHEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $(REPOSDIR)$(subst .,,$(suffix $(basename $@))).cpp

