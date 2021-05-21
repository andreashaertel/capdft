# Binary name ##################################################################
BNAME=main
# Compiler #####################################################################
CXX=g++
# Libraries to be linked #######################################################
LIBS=-lm
# Compiler flags ###############################################################
CXXFLAGS=-pedantic -std=c++2a -O2 -Wall
# Directories for binaries, objects and sources ################################
BDIR = bin/
ODIR = obj/
SDIR = src/
REPOSDIR = ../../src/
REPOFILE = species_properties
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
REPOSOURCE += $(REPOSDIR)$(REPOFILE).cpp
REPOHEADERS += $(REPOSDIR)$(REPOFILE).hpp
REPOOBJECTS += $(patsubst $(REPOSDIR)%.cpp,$(OOFF)%.o,$(REPOSOURCE))
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
default: createdir compile

all: createdir checkstyle compile

createdir:
	mkdir -p bin
	mkdir -p obj

# Google C++ style checker:
# Install cpplint with
# $pip3 install cpplint
# And add ~/.local/bin to PATH variable.
# Checks all files in the src directory.
checkstyle:
	cpplint $(SDIR)*
	cpplint $(REPOSDIR)$(REPOFILE).*

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
$(OOFF)species_properties.o: $(REPOSDIR)$(REPOFILE).cpp $(HEADERS) $(REPOHEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
