# CMake project file for SteamTesting

cmake_minimum_required (VERSION 2.6)
project (SteamTesting)
enable_language (Fortran)

# add the SteamPolynomials project
add_subdirectory ("../SteamPolynomials" ${CMAKE_BINARY_DIR})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
get_filename_component (default_prefix "/usr/local" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the /usr/local directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -funderscoring -march=native -fopenmp -fPIC -Wall -O3 -fno-trapping-math -fno-range-check")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-std=f2008 -ffree-form -funderscoring -march=native -Wall -fPIC -O0 -g3")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort
  set (CMAKE_Fortran_FLAGS_RELEASE "-free -O2 -xHost -mcmodel=large -parallel -openmp -fPIC -shared-intel")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-free -O0 -xHost -mcmodel=large -fPIC -debug full -g3 -shared-intel")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Fortran source files
set (SOURCES
"src/SteamTablesTesting.f"
"src/SteamTesting.f")

# build program
add_executable (steam-testing ${SOURCES} )
include_directories (${CMAKE_BINARY_DIR}/SteamPolynomials)
target_link_libraries (steam-testing ${CMAKE_BINARY_DIR}/SteamPolynomials/libSteamPolynomials.so)

# install executable
install (TARGETS steam-testing RUNTIME DESTINATION "bin")
