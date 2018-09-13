!-----------------------------------------------------------------------
!
!> @file SteamViscosity.f
!
!> @details
!> This module implements the polynomials which compute the viscocity
!> of steam according to the release on the IAPWS Formulation 2008 for
!> the Viscosity of Ordinary Water Substance (September 2008)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamViscosity

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Viscocity
    !
    PUBLIC viscosity

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function computes the viscosity in the dilute-gas limit.
    !
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The viscosity in the dilute-gas limit.
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION mu0(tau)

        IMPLICIT NONE
        !
        !> Constant coefficients "I"
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:4) :: I = &
                                                             (/0, 1, 2, 3/)
        !
        !> Constant coefficients "H"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:4) :: H = &
                                                           (/1.67752, 2.20462, 0.6366564, -0.241605/)
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: tau
        REAL(KIND=REAL_HIGH) :: mu0
        !
        ! Compute the viscosity in the dilute-gas limit.
        !
        mu0 = 100.0D+00 / (sqrt(tau) * SUM(H * (tau**I)))

    END FUNCTION mu0
    !-------------------------------------------------------------------------
    !
    !> This function computes the contribution to viscosity due to finite
    !> density.
    !
    !> @param[in] del Dimensionless density parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The contribution to viscosity due to finite density.
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION mu1(del, tau)

        IMPLICIT NONE
        !
        !> Constant coefficients "H"
        !
        REAL(KIND=REAL_HIGH), DIMENSION(1:6, 1:7) :: H
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: del, tau
        REAL(KIND=REAL_HIGH) :: mu1
        !
        ! Variables
        !
        INTEGER(KIND=INT_HIGH) :: i, j
        REAL(KIND=REAL_HIGH) :: sum, tau1
        !
        ! Compute the contribution to viscosity due to finite density
        !
        sum = zero
        H(1, :) = (/5.20094D-01, 2.22531D-01, -2.81378D-01, 1.61913D-01, -3.25372D-02, 0.0D+00, 0.0D+00/)
        H(2, :) = (/8.50895D-02, 9.99115D-01, -9.06851D-01, 2.57399D-01, 0.0D+00, 0.0D+00, 0.0D+00/)
        H(3, :) = (/-1.08374D+00, 1.88797D+00, -7.72479D-01, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00/)
        H(4, :) = (/-2.89555D-01, 1.26613D+00, -4.89837D-01, 0.0D+00, 6.98452D-02, 0.0D+00, -4.35673D-03/)
        H(5, :) = (/0.0D+00, 0.0D+00, -2.57040D-01, 0.0D+00, 0.0D+00, 8.72102D-03, 0.0D+00/)
        H(6, :) = (/0.0D+00, 1.20573D-01, 0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00, -5.93264D-04/)

        DO i = 1, 6
            tau1 = (tau - one)**(i - 1)
            DO j = 1, 7
                IF (H(i, j) == zero) CYCLE
                sum = sum + H(i, j) * tau1 * ((del - one)**(j - 1))
            END DO
        END DO

        mu1 = EXP(del * sum)

    END FUNCTION mu1
    !-------------------------------------------------------------------------
    !
    !> This function computes the viscosity as a function of density and
    !> temperature
    !
    !> @param[in] density     The steam density (Kg/m3)
    !> @param[in] temperature The steam temperature (K)
    !
    !> @return The viscosity in [Pa.sec]
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION viscosity(density, temperature)

        IMPLICIT NONE
        !
        !> Star viscosity for the computation of viscosity(density,temperature)
        !> in [Pa.sec]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: VISCOSITY_MUSTAR = 1.0D-06
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: viscosity
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: del, tau
        !
        ! Compute the viscosity
        !
        del = density / IAPWS97_RHOCRIT
        tau = IAPWS97_TCRIT / temperature
        viscosity = VISCOSITY_MUSTAR * mu0(tau) * mu1(del, tau)

    END FUNCTION viscosity

END MODULE SteamViscosity
