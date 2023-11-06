cmake_minimum_required(VERSION 3.20)
# made with help from Zhu Liang-Jun of IGSNRR in Beijing
# modified for GitHub by Trevor James Smith of Ouranos in Montreal

# optional cmake command line arguments (e.g, "cmake -D COMPILE_LIB=ON" .)
option(COMPILE_LIB "If ON, will create a dynamic lib file (default: OFF)" OFF)
option(COMPILE_EXE "If ON, will create a executable file (default: ON)" ON)

# Setup Project
PROJECT(Raven CXX)

# Add the cmake directory to the module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find NetCDF
find_package(NetCDF)
# find header & source
file(GLOB HEADER "src/*.h")
file(GLOB SOURCE "src/*.cpp")

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
endif()
source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

IF(NETCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven NetCDF::NetCDF)
ENDIF()

# unset cmake variables to avoid polluting the cache
unset(COMPILE_LIB CACHE)
unset(COMPILE_EXE CACHE)
