cmake_minimum_required(VERSION 2.8)

# project Information 
project(virtual_screening)
set(PROJECT_VERSION "1.0")

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# Require MPI for this project:
find_package(MPI REQUIRED)

# Set compiler flags
SET( CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
SET( CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
SET( CMAKE_CXX_FLAGS "-lm -fopenmp")

#Set CXX compiler for all files below
set_source_files_properties(src/vs_main_blocking.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/vs_main_adaptative.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/futil.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/string_owner.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/messages.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/load_parameters.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/docking.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/vina.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/execution_information.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/prepare_ligand.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/vs_main_blocking_seq.c PROPERTIES LANGUAGE CXX )
set_source_files_properties(src/docking_log.c PROPERTIES LANGUAGE CXX )

# set include
include_directories(include)
include_directories(MPI_INCLUDE_PATH)

# add libries
add_library(VS_lib STATIC
src/futil.c
src/string_owner.c
src/messages.c
src/load_parameters.c
src/docking.c
src/vina.c
src/execution_information.c
src/docking_log.c
)#end of VS_lib

# add target
add_executable(vs_main_blocking src/vs_main_blocking.c)
target_link_libraries(vs_main_blocking VS_lib ${MPI_C_LIBRARIES})

add_executable(prepare_ligand src/prepare_ligand.c)
target_link_libraries(prepare_ligand VS_lib ${MPI_C_LIBRARIES})

add_executable(vs_main_blocking_seq src/vs_main_blocking_seq.c)
target_link_libraries(vs_main_blocking_seq VS_lib ${MPI_C_LIBRARIES})

add_executable(vs_main_adaptative src/vs_main_adaptative.c)
target_link_libraries(vs_main_adaptative VS_lib ${MPI_C_LIBRARIES})


# install
install(TARGETS vs_main_blocking DESTINATION bin)
install(TARGETS prepare_ligand DESTINATION bin)
install(TARGETS vs_main_blocking_seq DESTINATION bin)