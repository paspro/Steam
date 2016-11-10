!-----------------------------------------------------------------------
!
!> @file SteamConstants.f
!
!> @details
!> This module defines constants required for the computation of the steam
!> polynomials according to the revised release on the IAPWS Industrial
!> Formulation 1997 for the Thermodynamic Properties of Water and Steam,
!> August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamConstants

    USE SteamVersion
    USE SteamPrecision
    USE SteamUnits

    IMPLICIT NONE

    !-----------------------------------------------------------------------
    !
    !> Steam Constants
    !
    !-----------------------------------------------------------------------

    !
    !> Steam constant R in [J/Kg.K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: IAPWS97_R = 461.526D+00
    !
    !> Steam maximum pressure in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: IAPWS97_PMAX = 100.0D+00
    !
    !> Steam minimum temperature in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: IAPWS97_TMIN = 273.15D+00
    !
    !> Steam maximum temperature in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: IAPWS97_TMAX = 2273.15D+00
    !
    !> Steam critical temperature in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: IAPWS97_TCRIT = 647.096D+00
    !
    !> Steam critical pressure in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: IAPWS97_PCRIT = 22.064D+00
    !
    !> Steam critical density in [Kg/m3]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: IAPWS97_RHOCRIT = 322.0D+00

    !-----------------------------------------------------------------------
    !
    !> Constants for Region 1
    !
    !-----------------------------------------------------------------------

    !
    !> Maximum temperature for region 1 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_TMAX = 623.15D+00
    !
    !> The value of pressure at the point of maximum temperature
    !> that intersects the line that represents region 4 in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_4_TMAX_P = 16.5292D+00

    !-----------------------------------------------------------------------
    !
    !> Constants for Region 2
    !
    !-----------------------------------------------------------------------

    !
    !> Steam maximum temperature in region 2 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_TMAX = 1073.15D+00
    !
    !> Temperature at the point where the line that represents region 4
    !> intersects region 2 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_4_T = 863.15D+00
    !
    !> Pressure at the boundary between the subregions 2a and 2b
    !> in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2A_2B_P = 4.0D+00
    !
    !> Minimum pressure at the intersection of regions 2b and 2c in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_2C_PMIN = 6.54670D+00

    !-----------------------------------------------------------------------
    !
    !> Constants for Region 3
    !
    !-----------------------------------------------------------------------

    !
    !> The entropy at the boundary of regions 3a and 3b in [J/Kg.K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_3B_S = 4.41202148223476D+03
    !
    !> The minimum enthalpy at the boundary of regions 3 and 4 in [J/Kg]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_HMIN = 1.670858218D+06
    !
    !> The maximum enthalpy at the boundary of regions 3 and 4 in [J/Kg]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_HMAX = 2.563592004D+06
    !
    !> The minimum entropy at the boundary of regions 3 and 4 in [J/Kg.K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_SMIN = 3.778281340D+06
    !
    !> The maximum entropy at the boundary of regions 3 and 4 in [J/Kg.K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_SMAX = 5.210887825D+06

    !-----------------------------------------------------------------------
    !
    !> Constants for Region 4
    !
    !-----------------------------------------------------------------------

    !
    !> Maximum temperature for region 4 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_4_TMAX = 647.096D+00

    !-----------------------------------------------------------------------
    !
    !> Constants for Region 5
    !
    !-----------------------------------------------------------------------

    !
    !> Maximum pressure in region 5 in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_5_PMAX = 50.0D+00

END MODULE SteamConstants
