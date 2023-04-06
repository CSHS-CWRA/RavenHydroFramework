cmake_minimum_required(VERSION 3.20)

# Setup Project
PROJECT(Raven CXX)

# find header & source
file(GLOB HEADER "*.h")
file(GLOB SOURCE "*.cpp")

add_executable(Raven
  ${SOURCE}
  ${HEADER}
)
source_group("Header Files" FILES ${HEADER})
source_group("Source Files" FILES ${SOURCE})

#set_target_properties(Raven PROPERTIES LINKER_LANGUAGE CXX)
