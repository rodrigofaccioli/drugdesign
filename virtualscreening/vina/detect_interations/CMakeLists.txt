cmake_minimum_required(VERSION 2.8)

# project Information 
project(detect_interactions)
set(PROJECT_VERSION "1.0")


# Set compiler flags
SET( CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )
SET( CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS})
SET( CMAKE_CXX_FLAGS "-lm -fno-stack-protector")

#Set CXX compiler for all files below
set_source_files_properties(src/detect_interactions.c PROPERTIES LANGUAGE CXX )

# set include
#include_directories(include)

# add target
add_executable(detect_interactions src/detect_interactions.c)
#target_link_libraries(detect_interactions)

# install
install(TARGETS detect_interactions DESTINATION bin)