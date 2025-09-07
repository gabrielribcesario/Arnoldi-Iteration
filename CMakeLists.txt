cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(Arnoldi LANGUAGES Fortran)

# ---Build libraries-------------------------------------------------------------

add_library(Ritz SHARED src/arnoldi.f90 src/francis.f90)
target_compile_options(Ritz PRIVATE -O3 -shared -fPIC -msse4.2)

add_library(Random SHARED src/random.f90)
target_compile_options(Random PRIVATE -shared -fPIC)

# ---Build executables-----------------------------------------------------------

add_executable(clock src/clock.f90)
target_link_libraries(clock PRIVATE Ritz Random)

add_executable(convergence src/convergence.f90)
target_link_libraries(convergence PRIVATE Ritz Random)