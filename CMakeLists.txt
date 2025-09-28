cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(ArnoldiIteration LANGUAGES Fortran)

# Use GNU-style directory hierarchy
include(GNUInstallDirs)

# make install output dir
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install")

# List of targets
set(INSTALLABLE_TARGETS "")

# ---Build libraries-------------------------------------------------------------

add_library(Ritz SHARED src/arnoldi.f90 src/francis.f90)
target_compile_options(Ritz PRIVATE -O3 -shared -fPIC -mavx2)
list(APPEND INSTALLABLE_TARGETS Ritz)

add_library(Random SHARED src/random.f90)
target_compile_options(Random PRIVATE -shared -fPIC)
list(APPEND INSTALLABLE_TARGETS Random)

# ---Build executables-----------------------------------------------------------

add_executable(clock src/clock.f90)
target_link_libraries(clock PRIVATE Ritz Random)
list(APPEND INSTALLABLE_TARGETS clock)

add_executable(convergence src/convergence.f90)
target_link_libraries(convergence PRIVATE Ritz Random)
list(APPEND INSTALLABLE_TARGETS convergence)

# ---Install libraries/executables-----------------------------------------------

# Add configuration file
configure_file("${PROJECT_SOURCE_DIR}/cmake/setvars.sh.in"
               "${CMAKE_INSTALL_PREFIX}/setvars.sh" @ONLY)

# Install targets
foreach(target ${INSTALLABLE_TARGETS})
    install(TARGETS "${target}"
            EXPORT "${PROJECT_NAME}Targets"
            LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
            RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endforeach()

# Install the Python script
install(FILES src/convergence.py 
        DESTINATION ${CMAKE_INSTALL_BINDIR} 
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)