cmake_minimum_required(VERSION 2.8)

# project Information 
project(detect_hbonds)
set(PROJECT_VERSION "1.0")


# Set compiler flags
SET( CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )
SET( CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS})
SET( CMAKE_CXX_FLAGS "-lm -fno-stack-protector")

#Set CXX compiler for all files below
set_source_files_properties(src/detect_hbonds.c PROPERTIES LANGUAGE CXX )

# set include
#include_directories(include)

# add target
add_executable(detect_hbonds src/detect_hbonds.c)

# install
install(TARGETS detect_hbonds DESTINATION bin)
