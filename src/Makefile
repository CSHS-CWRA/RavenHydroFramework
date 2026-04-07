# Alan J. Barton
# Created this file to allow compilation under Ubuntu linux.
#
# Modified from:
#        http://stackoverflow.com/questions/2481269/how-to-make-simple-c-makefile
# Stacksize increase compiler flag found on:
#        http://linuxtoosx.blogspot.ca/2010/10/stack-overflow-increasing-stack-limit.html
#
# Created
# Modified   Feb 2017    Juliane Mai    - add NetCDF support
#            Mar 2017    Juliane Mai    - increase stack size to 256 MB from standard 64MB
#            Jan 2019    Alan J. Barton - added OS name to created exe name
#            Dec 2023    David Huard    - support for python binding
#            Jul 2024    James Craig and Bohao Tang - support for lp_solve
#            Feb 2025    James Craig and Bohao Tang - support for new Mac (Applie Silicon)
###
# appname := raven_rev$(shell svnversion -n . | sed -e 's/:/-/g')_$(shell uname -o).exe
appname := Raven.exe

CXX      := g++
CXXFLAGS := -Wno-deprecated
LDLIBS   :=
LDFLAGS  :=

# OPTION 0) some compilers require the c++11 flag, some may not
CXXFLAGS += -std=c++11 -fPIC

# OPTION 1) include netcdf - uncomment following two commands (assumes netCDF path = /usr/local):
CXXFLAGS += -Dnetcdf
LDLIBS   += -L/usr/local -lnetcdf

# OPTION 1b) include netcdf - for newer MacOS with Apple Silicon (use with option 1 also uncommented)
#CXXFLAGS += -I/opt/homebrew/include
#LDLIBS   += -L/opt/homebrew/lib

# OPTION 2) include lp_solve for water management optimization- uncomment following two commands (assumes liblpsolve55.a in ../lib/lp_solve_unix/ folder):
#CXXFLAGS += -D_LPSOLVE_
#LDLIBS   += -L../lib/lp_solve_unix -llpsolve55

# OPTION 2b) include lp_solve for water management optimization- for newer MacOS with Apple Silicon (use with option 2 also uncommented):
#LDLIBS   += -L/opt/homebrew/lib

# OPTION 3) if you use a OSX/BSD system, uncomment the LDFLAGS line below
# this is to allow for use a 1Gb stack, see http://linuxtoosx.blogspot.ca/2010/10/stack-overflow-increasing-stack-limit.html
#LDFLAGS  += -Wl,-stack_size,0x80000000,-stack_addr,0xf0000000

# OPTION 4) uncomment for M1 and newer MacOS
#CXXFLAGS += -Dfinite=isfinite

srcfiles := $(shell ls *.cpp)
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) -o $(appname) $(objects) $(LDLIBS) $(LDFLAGS)

libraven:
	$(CXX) $(CXXFLAGS) -shared $(LDLIBS) $(shell python3 -m pybind11 --includes) $(LDFLAGS) -I . py/libraven.cpp -o libraven$(shell python3-config --extension-suffix)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend

include .depend
