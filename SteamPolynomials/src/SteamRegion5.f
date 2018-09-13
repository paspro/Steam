!-----------------------------------------------------------------------
!
!> @file SteamRegion5.f
!
!> @details
!> This module implements the polynomials which compute the thermodynamic
!> properties of steam in region 5 according to the revised release on the
!> IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
!> Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamRegion5

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Thermodynamic properties computed by polynomials
    !
    PUBLIC specific_internal_energy
    PUBLIC specific_volume
    PUBLIC specific_entropy
    PUBLIC specific_enthalpy
    PUBLIC speed_of_sound
    PUBLIC specific_isobaric_heat_capacity
    PUBLIC specific_isochoric_heat_capacity
    PUBLIC ratio_of_specific_heats
    PUBLIC specific_gibbs_free_energy
    !
    !> Constant coefficients "J" for the ideal gas part
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:6) :: J = &
                                                         (/0, 1, -3, -2, -1, 2/)
    !
    !> Constant coefficients "n" for the ideal gas part
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:6) :: n = &
                                                       (/-0.13179983674201D+02, 0.68540841634434D+01, -0.24805148933466D-01, &
                                                         0.36901534980333D+00, -0.31161318213925D+01, -0.32961626538917D+00/)
    !
    !> Constant coefficients "Ir" for the ideal gas part
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:6) :: Ir = &
                                                         (/1, 1, 1, 2, 2, 3/)
    !
    !> Constant coefficients "Jr" for the ideal gas part
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:6) :: Jr = &
                                                         (/1, 2, 3, 3, 9, 7/)
    !
    !> Constant coefficients "nr" for the ideal gas part
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:6) :: nr = &
                                                       (/0.15736404855259D-02, 0.90153761673944D-03, -0.50270077677648D-02, &
                                                         0.22440037409485D-05, -0.41163275453471D-05, 0.37919454822955D-07/)
    !
    !> Star pressure for region 5 in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_5_PSTAR = 1.0D+00
    !
    !> Star temperature for region 5 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_5_TSTAR = 1000.0D+00

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function computes the dimensionless specific Gibbs free energy
    !> for the ideal gas part
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The dimensionless specific Gibbs free energy
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gam0(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gam0
        !
        ! Polynomial expression
        !
        gam0 = LOG(pi) + SUM(n * (tau**J))

    END FUNCTION gam0
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the specific Gibbs free energy
    !> with respect to the parameter tau for the ideal gas part
    !
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d(Gibbs)/d(tau)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gam0tau(tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: tau
        REAL(KIND=REAL_HIGH) :: gam0tau
        !
        ! Polynomial expression
        !
        gam0tau = SUM(n * J * (tau**(J - 1)))

    END FUNCTION gam0tau
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the specific Gibbs
    !> free energy with respect to the parameter tau for the ideal gas part
    !
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d2(Gibbs)/d(tau)2
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gam0tautau(tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: tau
        REAL(KIND=REAL_HIGH) :: gam0tautau
        !
        ! Polynomial expression
        !
        gam0tautau = SUM(n * J * (J - 1) * (tau**(J - 2)))

    END FUNCTION gam0tautau
    !-------------------------------------------------------------------------
    !
    !> This function computes the dimensionless specific Gibbs free energy
    !> for the residual part
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The dimensionless specific Gibbs free energy
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamr(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gamr
        !
        ! Polynomial expression
        !
        gamr = SUM(nr * (pi**Ir) * (tau**Jr))

    END FUNCTION gamr
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the specific Gibbs free energy
    !> with respect to the parameter pi for the residual part
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d(Gibbs)/d(pi)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamrpi(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gamrpi
        !
        ! Polynomial expression
        !
        gamrpi = SUM(nr * Ir * (pi**(Ir - 1)) * (tau**Jr))

    END FUNCTION gamrpi
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the specific Gibbs
    !> free energy with respect to the parameter pi for the residual part
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d2(Gibbs)/d(pi)2
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamrpipi(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gamrpipi
        !
        ! Polynomial expression
        !
        gamrpipi = SUM(nr * Ir * (Ir - 1) * (pi**(Ir - 2)) * (tau**Jr))

    END FUNCTION gamrpipi
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the specific Gibbs free energy
    !> with respect to the parameter tau for the residual part
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d(Gibbs)/d(tau)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamrtau(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gamrtau
        !
        ! Polynomial expression
        !
        gamrtau = SUM(nr * (pi**Ir) * Jr * (tau**(Jr - 1)))

    END FUNCTION gamrtau
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the specific Gibbs
    !> free energy with respect to the parameter tau for the residual part
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d2(Gibbs)/d(tau)2
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamrtautau(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gamrtautau
        !
        ! Polynomial expression
        !
        gamrtautau = SUM(nr * (pi**Ir) * Jr * (Jr - 1) * (tau**(Jr - 2)))

    END FUNCTION gamrtautau
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the specific Gibbs
    !> free energy with respect to the parameters pi and tau for the
    !> residual part
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d2(Gibbs)/d(pi)d(tau)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamrpitau(pi, tau)

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gamrpitau
        !
        ! Polynomial expression
        !
        gamrpitau = SUM(nr * Ir * (pi**(Ir - 1)) * Jr * (tau**(Jr - 1)))

    END FUNCTION gamrpitau
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific internal energy with respect to
    !> pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific internal energy [J/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_internal_energy(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_internal_energy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, tau, gamtau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the specific internal energy
        !
        gamtau = gam0tau(tau) + gamrtau(pi, tau)
        specific_internal_energy = IAPWS97_R * temperature * (tau * gamtau - one - pi * gamrpi(pi, tau))

    END FUNCTION specific_internal_energy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific volume with respect to
    !> pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific volume [m3/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_volume(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_volume
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, tau, res1
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the specific volume
        !
        res1 = one + gamrpi(pi, tau) * pi
        specific_volume = IAPWS97_R * temperature * res1 / (pressure * mega)

    END FUNCTION specific_volume
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific entropy with respect to
    !> pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific entropy [J/kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_entropy(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_entropy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, tau, gam, gamtau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the specific entropy
        !
        gam = gam0(pi, tau) + gamr(pi, tau)
        gamtau = gam0tau(tau) + gamrtau(pi, tau)
        specific_entropy = IAPWS97_R * (tau * gamtau - gam)

    END FUNCTION specific_entropy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific enthalpy with respect to
    !> pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific enthalpy [J/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_enthalpy(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_enthalpy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, tau, gamtau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the specific enthalpy
        !
        gamtau = gam0tau(tau) + gamrtau(pi, tau)
        specific_enthalpy = IAPWS97_R * temperature * tau * gamtau

    END FUNCTION specific_enthalpy
    !-------------------------------------------------------------------------
    !
    !> This function computes the speed of sound with respect
    !> to pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The speed of sound [m/sec]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION speed_of_sound(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: speed_of_sound
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: res1, res2, res3, res4, pi, tau, gp
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the speed of sound
        !
        gp = gamrpi(pi, tau)
        res1 = IAPWS97_R * temperature * (one + (two + pi * gp) * pi * gp)
        res2 = one - pi * pi * gamrpipi(pi, tau)
        res3 = (one + pi * gp - tau * pi * gamrpitau(pi, tau))**2
        res4 = tau * tau * (gam0tautau(tau) + gamrtautau(pi, tau))

        speed_of_sound = sqrt(res1 / (res2 + res3 / res4))

    END FUNCTION speed_of_sound
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific isobaric heat capacity with respect
    !> to pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific isobaric heat capacity [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_isobaric_heat_capacity(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_isobaric_heat_capacity
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the specific isobaric heat capacity
        !
        specific_isobaric_heat_capacity = -IAPWS97_R * tau * tau * (gam0tautau(tau) + gamrtautau(pi, tau))

    END FUNCTION specific_isobaric_heat_capacity
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific isochoric heat capacity with respect
    !> to pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific isochoric heat capacity [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_isochoric_heat_capacity(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_isochoric_heat_capacity
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: res1, res2, res3, pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the specific isochoric heat capacity
        !
        res1 = -tau * tau * (gam0tautau(tau) + gamrtautau(pi, tau))
        res2 = (one + pi * gamrpi(pi, tau) - tau * pi * gamrpitau(pi, tau))**2
        res3 = one - pi * pi * gamrpipi(pi, tau)

        specific_isochoric_heat_capacity = IAPWS97_R * (res1 - res2 / res3)

    END FUNCTION specific_isochoric_heat_capacity
    !-------------------------------------------------------------------------
    !
    !> This function computes the ratio of specific heats with respect
    !> to pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The ratio of specific heats
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION ratio_of_specific_heats(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: ratio_of_specific_heats
        !
        ! Compute the ratio of specific heats
        !
        ratio_of_specific_heats = specific_isobaric_heat_capacity(pressure, temperature) / &
                                  specific_isochoric_heat_capacity(pressure, temperature)

    END FUNCTION ratio_of_specific_heats
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific Gibbs free energy with respect
    !> to pressure and temperature.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific Gibbs free energy [J/kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_gibbs_free_energy(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_gibbs_free_energy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_5_PSTAR
        tau = REGION_5_TSTAR / temperature
        !
        ! Compute the specific Gibbs free energy
        !
        specific_gibbs_free_energy = IAPWS97_R * temperature * (gam0(pi, tau) + gamr(pi, tau))

    END FUNCTION specific_gibbs_free_energy

END MODULE SteamRegion5
