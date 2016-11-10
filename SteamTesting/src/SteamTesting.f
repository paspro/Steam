!-----------------------------------------------------------------------
!
!> @file SteamTesting.f
!
!> @details
!> This program tests the computation of the thermodynamic properties of
!> steam performed by the SteamPolynomials library.
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
PROGRAM SteamTesting

    USE SteamTablesTesting

    CALL test_region1
    CALL test_region2
    CALL test_region3
    CALL test_region4
    CALL test_boundary23
    CALL test_boundary2bc
    CALL test_boundary2ab
    CALL test_boundary3ab
    CALL test_viscosity
    CALL test_thermal_conductivity
    CALL test_cubic_interpolation
    CALL test_bilinear_interpolation
    CALL test_table_generation

END PROGRAM SteamTesting
