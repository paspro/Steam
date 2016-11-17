#
# \file CMakeLists.txt
#
# \brief     Steam Library
# \author    Dr. Panos Asproulis (panos.asproulis@renuda.com)
# \date      2016
# \copyright Panos Asproulis (2014-2016). All Rights Reserved.
#
cmake_minimum_required (VERSION 2.6)
##
# Common packages
##
SET(SUBDIRS_COMMON
  SteamPolynomials
  SteamTesting
)

SET(SUBDIRS
  ${SUBDIRS_COMMON}
)

FOREACH(dir ${SUBDIRS})
 ADD_SUBDIRECTORY(${dir})
ENDFOREACH(dir ${SUBDIRS})
