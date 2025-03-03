cmake_minimum_required(VERSION 3.20)
# made with help from Zhu Liang-Jun of IGSNRR in Beijing
# modified for GitHub and C++11 by Trevor James Smith of Ouranos in Montreal
# modified for lpsolve support by Maxim Krassovski of DFO

# Set the C++ standard to C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF)

# optional cmake command line arguments (e.g, "cmake -D COMPILE_LIB=ON" .)
option(COMPILE_LIB "If ON, will create a dynamic lib file (default: OFF)" OFF)
option(COMPILE_EXE "If ON, will create a executable file (default: ON)" ON)
option(PYTHON, "If ON, will create a share library for python (default: OFF)" OFF)
option(LPSOLVE, "If ON, will link to lp_solve optimization library (default: OFF)" OFF)

# Setup Project
PROJECT(Raven CXX)

# Add the cmake directory to the module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find NetCDF
find_package(NetCDF) #may also be 'netCDF'
if (NOT NetCDF_FOUND)
  message(STATUS "NetCDF not found, trying netCDF")
  find_package(netCDF)
endif()

# Find HDF5
find_package(HDF5)

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
  set_target_properties(Raven PROPERTIES LINKER_LANGUAGE CXX)

  # Remove deprecation warnings for GCC
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(Raven PRIVATE "-Wno-deprecated")
    message(STATUS "Modified compile flags with '-Wno-deprecated'")
  ENDIF()

  if(LPSOLVE)
    target_link_directories(Raven PRIVATE lib/lp_solve)  # where liblpsolve55.so is located
    target_link_libraries(Raven lpsolve55)
    add_definitions(-D_LPSOLVE_)
  endif()
endif()

if(NETCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven NetCDF::NetCDF)
elseif(netCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven netcdf)
endif()

source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

# unset cmake variables to avoid polluting the cache
unset(COMPILE_LIB CACHE)
unset(COMPILE_EXE CACHE)
