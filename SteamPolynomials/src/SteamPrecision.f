!-----------------------------------------------------------------------
!
!> @file SteamPrecision.f
!
!> @details
!> This module defines the numerical precision used for the computation
!> of the steam polynomials according to the revised release on the
!> IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
!> Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamPrecision

    USE ISO_FORTRAN_ENV

    IMPLICIT NONE
    !-----------------------------------------------------------------------
    !
    !> Arithmetic Precision
    !
    !-----------------------------------------------------------------------

    !
    !> Floating point precision, requesting 64-bits
    !
    INTEGER, PARAMETER :: REAL_HIGH = REAL64
    !
    !> Integer precision, requesting 32-bits
    !
    INTEGER, PARAMETER :: INT_HIGH = INT32
    
    !-----------------------------------------------------------------------
    !
    !> Numerical constants in the specified precision
    !
    !-----------------------------------------------------------------------
    
    REAL(KIND=REAL_HIGH), PARAMETER :: &
        zero     = 0.0D+00,            &
        onesixth = 1.0D+00/6.0D+00,    &
        quarter  = 0.25D+00,           &
        onethird = 1.0D+00/3.0D+00,    &
        half     = 0.5D+00,            &
        one      = 1.0D+00,            &
        two      = 2.0D+00,            &
        three    = 3.0D+00,            &
        four     = 4.0D+00,            &
        five     = 5.0D+00,            &
        six      = 6.0D+00,            &
        seven    = 7.0D+00,            &
        eight    = 8.0D+00,            &
        nine     = 9.0D+00,            &
        ten      = 10.0D+00

END MODULE SteamPrecision


