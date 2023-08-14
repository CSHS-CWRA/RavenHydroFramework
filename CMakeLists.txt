cmake_minimum_required(VERSION 3.20)
# made with help from Zhu Liang-Jun of IGSNRR in Beijing

# Setup Project
PROJECT(Raven CXX)

# Find NetCDF
find_package(NetCDF)
# find header & source
file(GLOB HEADER "src/*.h")
file(GLOB SOURCE "src/*.cpp")

# creates a shared library - file extension is OS dependent (Linux: .so, Windows: .dll)
add_library(ravenbmi SHARED ${SOURCE})

# creates an executable - file extension is OS dependent (Linux: none, Windows: .exe)
add_executable(Raven
  ${SOURCE}
  ${HEADER}
)
source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

IF(NETCDF_FOUND)
  message("${NETCDF_INCLUDES}")
  add_definitions(-Dnetcdf)
  include_directories(${NETCDF_INCLUDES})
  target_link_libraries(Raven ${NETCDF_LIBRARIES})
ENDIF()

set_target_properties(Raven PROPERTIES LINKER_LANGUAGE CXX)
