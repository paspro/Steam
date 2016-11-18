# CMake project file for SteamTesting

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (SteamTesting)
ENABLE_LANGUAGE (Fortran)

# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

# default installation
GET_FILENAME_COMPONENT (default_prefix "/usr/local" ABSOLUTE)
SET (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the /usr/local directory."
      FORCE)

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  SET (CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -funderscoring -march=native -fopenmp -fPIC -Wall -O3 -fno-trapping-math -fno-range-check")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-std=f2008 -ffree-form -funderscoring -march=native -Wall -fPIC -O0 -g3")
ELSEIF (Fortran_COMPILER_NAME MATCHES "IFort.*")
  # IFort
  SET (CMAKE_Fortran_FLAGS_RELEASE "-free -O2 -xHost -mcmodel=large -parallel -openmp -fPIC -shared-intel")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-free -O0 -xHost -mcmodel=large -fPIC -debug full -g3 -shared-intel")
ELSE (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  MESSAGE ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  MESSAGE ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  MESSAGE ("No optimized Fortran compiler flags are known, we just try -O2...")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
ENDIF (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Fortran source files
SET (SOURCES
"src/SteamTablesTesting.f"
"src/SteamTesting.f")

# build program
ADD_EXECUTABLE (steam-testing ${SOURCES} )
INCLUDE_DIRECTORIES (${CMAKE_BINARY_DIR}/../SteamPolynomials)
TARGET_LINK_LIBRARIES (steam-testing ${CMAKE_BINARY_DIR}/../SteamPolynomials/libSteamPolynomials.so)

# install executable
INSTALL (TARGETS steam-testing RUNTIME DESTINATION "bin")
