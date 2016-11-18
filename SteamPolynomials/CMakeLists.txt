# CMake project file for SteamPolynomials

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (SteamPolynomials)
ENABLE_LANGUAGE (Fortran)

# make sure that the default is a RELEASE
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

# default INSTALLation
GET_FILENAME_COMPONENT (default_prefix "/usr/local" ABSOLUTE)
SET (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the INSTALLation directory; by default it INSTALLs in the /usr/local directory."
      FORCE)

# FFLAGS depend on the compiler
GET_FILENAME_COMPONENT (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

IF (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # The GNU Fortran Compiler (gfortran)
  SET (CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -ffree-form -funderscoring -march=native -fopenmp -Wall -fPIC -O3 -fno-trapping-math -fno-range-check")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-std=f2008 -ffree-form -funderscoring -march=native -Wall -fPIC -O0 -g3")
ELSEIF (Fortran_COMPILER_NAME MATCHES "IFort.*")
  # The Intel Fortran Compiler (IFort)
  SET (CMAKE_Fortran_FLAGS_RELEASE "-free -O2 -xHost -mcmodel=large -parallel -openmp -shared-intel")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-free -O0 -xHost -mcmodel=large -debug full -g3 -shared-intel")
ELSE (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  MESSAGE ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  MESSAGE ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  MESSAGE ("No optimized Fortran compiler flags are known, we just try -O2...")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  SET (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
ENDIF (Fortran_COMPILER_NAME MATCHES "gfortran.*")

# Fortran source files
SET (SOURCES 
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
SET (MODULES 
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
ADD_LIBRARY (SteamPolynomials SHARED ${SOURCES} )

# Install library and module files
INSTALL (FILES ${MODULES} DESTINATION "include/SteamPolynomials")
INSTALL (TARGETS SteamPolynomials LIBRARY DESTINATION "lib") 
