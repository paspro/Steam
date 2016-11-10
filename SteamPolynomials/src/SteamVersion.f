!-----------------------------------------------------------------------
!
!> @file SteamVersion.f
!
!> @details
!> This module defines the version of the Steam Polynomials library
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamVersion

    IMPLICIT NONE
    !
    !> Major version number
    !
    INTEGER, PARAMETER :: major_version = 1
    !
    !> Minor version number
    !
    INTEGER, PARAMETER :: minor_version = 0
    !
    !> Build number
    !
    INTEGER, PARAMETER :: build = 2
    !
    !> Version status
    !
    CHARACTER(LEN=*), PARAMETER :: code_status = "release"

END MODULE SteamVersion
