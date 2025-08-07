cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(Arnoldi LANGUAGES Fortran)

# ---Build executables-----------------------------------------------------------

add_library(Ritz SHARED src/arnoldi.f90 src/francis.f90)
target_compile_options(Ritz PRIVATE -O3 -shared -fPIC)

add_executable(eigenvalues src/eigenvalues.f90)
target_link_libraries(eigenvalues PRIVATE Ritz)

add_executable(convergence src/convergence.f90)
target_link_libraries(convergence PRIVATE Ritz)