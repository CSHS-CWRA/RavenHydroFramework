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
###
# appname := raven_rev$(shell svnversion -n . | sed -e 's/:/-/g')_$(shell uname -o).exe
appname := Raven.exe

CXX      := g++
CXXFLAGS := -Wno-deprecated
# some compilers require the c++11 flag, some may not
CXXFLAGS += std=-c++11

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

# include netcdf
CXXFLAGS += -Dnetcdf                # if netcdf is installed use "CXXFLAGS += -Dnetcdf",                    else "CXXFLAGS += "
LDLIBS   := -L/usr/local -lnetcdf   # if netcdf is installed give first path "-L<PATH>"and then "-lnetcdf", else "LDLIBS   := "

# if you use a OSX/BSD system, uncomment the LDFLAGS line below
# this is to allow for use a 1Gb, see http://linuxtoosx.blogspot.ca/2010/10/stack-overflow-increasing-stack-limit.html
# LDFLAGS  := -Wl,-stack_size,0x80000000,-stack_addr,0xf0000000

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) -o $(appname) $(objects) $(LDLIBS) $(LDFLAGS) 

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend

include .depend

