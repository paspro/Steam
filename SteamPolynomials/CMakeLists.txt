# CMake project file for SteamPolynomials

cmake_minimum_required (VERSION 2.6)
project (SteamPolynomials)
enable_language (Fortran)

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
  # The GNU Fortran Compiler (gfortran)
  set (CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -funderscoring -march=native -fopenmp -Wall -fPIC -O3 -fno-trapping-math -fno-range-check")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-std=f2008 -ffree-form -funderscoring -march=native -Wall -fPIC -O0 -g3")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # The Intel Fortran Compiler (ifort)
  set (CMAKE_Fortran_FLAGS_RELEASE "-free -O2 -xHost -mcmodel=large -parallel -openmp -shared-intel")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-free -O0 -xHost -mcmodel=large -debug full -g3 -shared-intel")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Fortran source files
set (SOURCES 
"src/SteamVersion.f"
"src/SteamPrecision.f"
"src/SteamUnits.f"
"src/SteamConstants.f"
"src/SteamRegions.f"
"src/SteamRegion1.f"
"src/SteamRegion2.f"
"src/SteamRegion3.f"
"src/SteamRegion4.f"
"src/SteamRegion5.f"
"src/SteamBoundaries.f"
"src/SteamViscosity.f"
"src/SteamConductivity.f"
"src/SteamStates.f"
"src/SteamNumerics.f"
"src/SteamTables.f")

# Fortran module files
set (MODULES 
"${CMAKE_BINARY_DIR}/SteamVersion.mod"
"${CMAKE_BINARY_DIR}/SteamPrecision.mod"
"${CMAKE_BINARY_DIR}/SteamUnits.mod"
"${CMAKE_BINARY_DIR}/SteamConstants.mod"
"${CMAKE_BINARY_DIR}/SteamRegions.mod"
"${CMAKE_BINARY_DIR}/SteamRegion1.mod"
"${CMAKE_BINARY_DIR}/SteamRegion2.mod"
"${CMAKE_BINARY_DIR}/SteamRegion3.mod"
"${CMAKE_BINARY_DIR}/SteamRegion4.mod"
"${CMAKE_BINARY_DIR}/SteamRegion5.mod"
"${CMAKE_BINARY_DIR}/SteamBoundaries.mod"
"${CMAKE_BINARY_DIR}/SteamViscosity.mod"
"${CMAKE_BINARY_DIR}/SteamConductivity.mod"
"${CMAKE_BINARY_DIR}/SteamStates.mod"
"${CMAKE_BINARY_DIR}/SteamNumerics.mod"
"${CMAKE_BINARY_DIR}/SteamTables.mod")

# build library
add_library (SteamPolynomials SHARED ${SOURCES} )

# install library and module files
install (FILES ${MODULES} DESTINATION "include/SteamPolynomials")
install (TARGETS SteamPolynomials LIBRARY DESTINATION "lib") 
