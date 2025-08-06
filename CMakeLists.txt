cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(Arnoldi LANGUAGES Fortran)

# ---Create directories-----------------------------------------------------------

# lib dir
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# build dir
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

# data .exe output dir
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/data")

# ---Build executables-----------------------------------------------------------

add_library(Ritz SHARED src/arnoldi.f90 src/francis.f90)
target_compile_options(Ritz PRIVATE -O3 -shared -fPIC)

add_executable(eig src/eig.f90)
target_link_libraries(eig Ritz)