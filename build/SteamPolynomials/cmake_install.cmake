# Install script for directory: /home/panos/Development/Steam/SteamPolynomials

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/SteamPolynomials" TYPE FILE FILES
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamVersion.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamPrecision.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamUnits.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamConstants.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamRegions.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamRegion1.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamRegion2.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamRegion3.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamRegion4.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamRegion5.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamBoundaries.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamViscosity.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamConductivity.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamStates.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamNumerics.mod"
    "/home/panos/Development/Steam/build/SteamPolynomials/SteamTables.mod"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSteamPolynomials.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSteamPolynomials.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSteamPolynomials.so"
         RPATH "")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/panos/Development/Steam/build/SteamPolynomials/libSteamPolynomials.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSteamPolynomials.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSteamPolynomials.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSteamPolynomials.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/panos/Development/Steam/build/SteamPolynomials/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
