cmake_minimum_required(VERSION 3.20)
# made with help from Zhu Liang-Jun of IGSNRR in Beijing
# modified for GitHub by Trevor James Smith of Ouranos in Montreal

# Setup Project
PROJECT(Raven CXX)

# Add the cmake directory to the module path
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Find NetCDF
find_package(NetCDF)
# find header & source
file(GLOB HEADER "src/*.h")
file(GLOB SOURCE "src/*.cpp")

add_executable(Raven
  ${SOURCE}
  ${HEADER}
)
source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

IF(NETCDF_FOUND)
  add_definitions(-Dnetcdf)
  include_directories(${NetCDF_INCLUDE_DIRS})
  target_link_libraries(Raven ${NetCDF})
ENDIF()

set_target_properties(Raven PROPERTIES LINKER_LANGUAGE CXX)
