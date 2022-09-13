# SPDX-FileCopyrightText: 2021 Andreas Härtel <http://andreashaertel.anno1982.de/>
# SPDX-FileCopyrightText: 2022 Moritz Bültmann <http://andreashaertel.anno1982.de/>
# SPDX-License-Identifier: LGPL-3.0-or-later
################################################################################
#
# Makefile
# for the project capdft 
#
# The project contains the C++ library capdft. 
#
################################################################################
## binary name #################################################################
LIBNAME = capdft
# Compiler #####################################################################
CXX = g++
# Include directories ##########################################################
INCDIRS = -I./ -I./include
# Library directory ############################################################
LIBDIRS = -L./ -L./bin
# Libraries to be linked #######################################################
LIBS = -lm -lpthread -lfftw3 -lgsl -lgslcblas
# Compiler flags ###############################################################
CXXFLAGS = -pedantic -std=c++2a -O2 -Wall -fopenmp $(INCDIRS)
# Linker flags #################################################################
LDFLAGS = -fopenmp $(LIBDIRS)
# Directories for binaries, objects and sources ################################
IDIR = include
BDIR = bin
ODIR = obj
SDIR = src
# Directories for examples #####################################################
EXAMPLEDIR = examples
EXAMPLES = $(wildcard $(EXAMPLEDIR)/*)
EXAMPLEMAKES = $(patsubst %,%/Makefile,$(EXAMPLES))
# Get host name ################################################################
HOST = $(shell hostname)
# Offsets for binaries and objects #############################################
BOFF = $(BDIR)/$(HOST).
OOFF = $(ODIR)/$(HOST).
IOFF = $(IDIR)/$(HOST).
# Header files and object files ################################################
HEADERS = $(wildcard $(SDIR)/*.hpp)
SOURCES = $(wildcard $(SDIR)/*.cpp)
OBJECTS = $(patsubst $(SDIR)/%.cpp,$(OOFF)%.o,$(SOURCES))
# System information files #####################################################
SYSINFO = $(OOFF)systeminfo
# TARGETS ######################################################################

.PHONY: default all init checkstyle compile info clean examples

default: compile bind

all: checkstyle info compile bind examples

# Initialize repository
init:
	@echo " Check for existing directories ... "
	mkdir -p obj
	mkdir -p include
	mkdir -p bin

# Google C++ style checker:
# Install cpplint with
# $pip3 install cpplint
# And add ~/.local/bin to PATH variable.
# Checks all files in the src directory.
checkstyle:
	cpplint src/*

# Compile all objects and write system information
compile: info $(OBJECTS)

# Write system information
info:
	@echo "HOSTNAME=$(shell hostname)" > $(SYSINFO)
	@echo "DATE=$(shell date)" >> $(SYSINFO)
	@echo "SYSTEM=$(shell uname -a)" >> $(SYSINFO)
	@echo "COMPILERVERSION=$(shell $(CXX) --version)" >> $(SYSINFO)
	@echo "CPUINFO=$(shell lshw -class processor 2>/dev/null)" >> $(SYSINFO)
	@echo "DISPLAY=$(shell lspci -vnn | grep "VGA" 2>/dev/null)" >> $(SYSINFO)

# Cleaning all targets (objects, binaries, and final includes)
clean:
	@echo " Cleaning objects ... "
	rm -f $(ODIR)*
	@echo " Cleaning binaries ... "
	rm -f $(BDIR)*
	@echo " Cleaning includes ... "
	rm -f $(IDIR)*

examples: $(EXAMPLES)

# Syntax: 
# $<	the first dependency
# $+	list of all dependencies
# $^	list of all dependencies; repeating entries are droped
# $@	name of the target

################################################################################
# Binding
################################################################################

# Make header and library and programs
bind: $(IDIR)$(LIBNAME).hpp $(BDIR)lib$(HOST).$(LIBNAME).a

# Make header
$(IDIR)$(LIBNAME).hpp: $(HEADERS)
	@echo " Make header $@ ... "
	cat $^ > $@

# Make library
$(BDIR)lib$(HOST).$(LIBNAME).a: $(OBJECTS)
	@echo " Binding $@ ... "
	ar -rc $@ $^

################################################################################
# Compiling
################################################################################

# Compile each source file (depending on its header file)
# 
# If each source only depends on its own header (no cross linking) then we 
# can use:
# $(OOFF)%.o: $(SDIR)%.cpp $(SDIR)%.hpp
# 	@echo " Compiling $< ... "
# 	$(CXX) $(CXXFLAGS) -c -o $@ $<
# However, sometimes we expect to find cross linking. 
# Then, we have to add an explicit rule. 

#################
# Explicit rules: 

$(OOFF)parameter_handler.o: $(SDIR)parameter_handler.cpp $(SDIR)parameter_handler.hpp
	@echo " Compiling $< ... "
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#################
# General rule: 

$(OOFF)%.o: $(SDIR)%.cpp $(SDIR)%.hpp
	@echo " Compiling $< ... "
	$(CXX) $(CXXFLAGS) -c -o $@ $<

################################################################################
# Compiling examples
################################################################################

$(EXAMPLES): FORCE
	@echo "Compiling example \"$(patsubst $(EXAMPLEDIR)/%,%,$@)\" ..."
	cd $@ && make

FORCE: ;
