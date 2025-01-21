cmake_minimum_required(VERSION 3.20)
# made with help from Zhu Liang-Jun of IGSNRR in Beijing
# modified for GitHub by Trevor James Smith of Ouranos in Montreal
# modified for lpsolve support by Maxim Krassovski of DFO

# optional cmake command line arguments (e.g, "cmake -D COMPILE_LIB=ON" .)
option(COMPILE_LIB "If ON, will create a dynamic lib file (default: OFF)" OFF)
option(COMPILE_EXE "If ON, will create a executable file (default: ON)" ON)
option(PYTHON, "If ON, will create a share library for python (default: OFF)" OFF)
option(LPSOLVE,"If ON, will link to lp_solve optimization library (default: OFF)" OFF)

# Setup Project
PROJECT(Raven CXX)

# Add the cmake directory to the module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find NetCDF
find_package(NetCDF) #may also be 'netCDF'
find_package(netCDF)
# find header & source
file(GLOB HEADER "src/*.h")
file(GLOB SOURCE "src/*.cpp")

# Create library with Python bindings
# To work, install pybind11, making sure it comes with cmake files (conda install -c conda-forge pybind11)
if(PYTHON)
  SET(PYBIND11_NEWPYTHON ON)
  find_package(Python COMPONENTS REQUIRED Interpreter Development)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(libraven MODULE src/py/libraven.cpp)
  include_directories("src")
  target_compile_features(libraven PUBLIC cxx_std_11)
endif()

# creates a shared library - file extension is OS dependent (Linux: .so, Windows: .dll)
if(COMPILE_LIB)
  add_library(ravenbmi SHARED ${SOURCE})
  target_compile_definitions(ravenbmi PUBLIC BMI_LIBRARY)
endif()

# creates an executable - file extension is OS dependent (Linux: none, Windows: .exe)
if(COMPILE_EXE)
  add_executable(Raven
    ${SOURCE}
    ${HEADER}
  )
  target_compile_definitions(Raven PUBLIC STANDALONE)
  set_target_properties(Raven PROPERTIES LINKER_LANGUAGE CXX)
  if(LPSOLVE)
    target_link_directories(Raven PRIVATE lib/lp_solve)  # where liblpsolve55.so is located
    target_link_libraries(Raven lpsolve55)
    add_definitions(-D_LPSOLVE_)
  endif()
endif()

source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

if(NETCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven NetCDF::NetCDF)
elseif(netCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven netcdf)
endif()

# unset cmake variables to avoid polluting the cache
unset(COMPILE_LIB CACHE)
unset(COMPILE_EXE CACHE)
