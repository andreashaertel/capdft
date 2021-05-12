# Binary name ##################################################################
BNAME=main
# Compiler #####################################################################
CXX=g++
# Include directory ############################################################
INCDIR=
# Library directory ############################################################
LIBDIR=-I./
# Libraries to be linked #######################################################
LIBS=-lm
# Compiler flags ###############################################################
CXXFLAGS=-pedantic -std=c++2a -O2 -Wall $(INCDIR)
# Directories for binaries, objects and sources ################################
BDIR = bin/
ODIR = obj/
SDIR = src/
# Get host name ################################################################
HOST = $(shell hostname)
# Offsets for binaries and objects #############################################
BOFF = $(BDIR)$(HOST).
OOFF = $(ODIR)$(HOST).
# Header files and object files ################################################
HEADERS = $(wildcard $(SDIR)*.h)
HEADERS += $(wildcard $(SDIR)*.hpp)
SOURCE = $(wildcard $(SDIR)*.cpp)
OBJECTS = $(patsubst $(SDIR)%.cpp,$(OOFF)%.o,$(SOURCE))
# Get host name ################################################################
HOST = $(shell hostname)
# Offsets for binaries and objects #############################################
BOFF = $(BDIR)$(HOST).
OOFF = $(ODIR)$(HOST).
# TARGETS ######################################################################
default: compile

all: checkstyle compile

# Google C++ style checker:
# Install cpplint with
# $pip3 install cpplint
# And add ~/.local/bin to PATH variable.
# Checks all files in the src directory.
checkstyle:
	cpplint src/*

clean:
	rm -f $(OOFF)*.o
	rm -f $(BOFF)*

compile: $(BOFF)$(BNAME)

# $<  first dependency
# $@  targets
# $+  list of all dependencies
# $^  list of all dependencies, no doubles

$(BOFF)$(BNAME): $(OBJECTS)
	$(CXX) -o $@ -L$(LIBDIR) $^ $(LIBS)

$(OOFF)%.o: $(SDIR)%.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<
