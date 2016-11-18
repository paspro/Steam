#
# \file CMakeLists.txt
#
# \brief     Steam Library
# \author    Dr. Panos Asproulis (panos.asproulis@renuda.com)
# \date      2016
# \copyright Panos Asproulis (2014-2016). All Rights Reserved.
#

CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (Steam)

EXECUTE_PROCESS (COMMAND ${CMAKE_COMMAND} 
                 -H${CMAKE_CURRENT_SOURCE_DIR}/SteamPolynomials
                 -B${CMAKE_CURRENT_BINARY_DIR}/SteamPolynomials)

EXECUTE_PROCESS (COMMAND ${CMAKE_COMMAND} --build
                 ${CMAKE_CURRENT_BINARY_DIR}/SteamPolynomials)

EXECUTE_PROCESS (COMMAND ${CMAKE_COMMAND} 
                 -H${CMAKE_CURRENT_SOURCE_DIR}/SteamTesting
                 -B${CMAKE_CURRENT_BINARY_DIR}/SteamTesting)

EXECUTE_PROCESS (COMMAND ${CMAKE_COMMAND} --build
                 ${CMAKE_CURRENT_BINARY_DIR}/SteamTesting)

        