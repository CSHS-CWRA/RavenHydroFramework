# Alan J. Barton
# Created this file to allow compilation under Ubuntu linux.
# Modified from: http://stackoverflow.com/questions/2481269/how-to-make-simple-c-makefile
###
appname := raven_rev.exe #$(shell svnversion -n . | sed -e 's/:/-/g').exe

CXX := g++
CXXFLAGS := -Wno-deprecated

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

# include netcdf
LDLIBS   := -L/usr/local -lnetcdf

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

depend: .depend

.depend: $(srcfiles)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(objects)

dist-clean: clean
	rm -f *~ .depend

include .depend

