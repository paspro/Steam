!-----------------------------------------------------------------------
!
!> @file SteamRegion1.f
!
!> @details
!> This module implements the polynomials which compute the thermodynamic
!> properties of steam in region 1 according to the revised release on the
!> IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
!> Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamRegion1

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
    ! Backward equation
    !
    PUBLIC temperature_ph
    PUBLIC pressure_hs
    !
    !> Constant coefficients "I"
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:34) :: I = &
                                                          (/0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, &
                                                            2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32/)
    !
    !> Constant coefficients "J"
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:34) :: J = &
                                                          (/-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, &
                                                            1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, &
                                                            -38, -39, -40, -41/)
    !
    !> Constant coefficients "n"
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:34) :: n = &
                                                        (/0.14632971213167D+00, -0.84548187169114D+00, -0.37563603672040D+01, &
                                                          0.33855169168385D+01, -0.95791963387872D+00, 0.15772038513228D+00, &
                                                          -0.16616417199501D-01, 0.81214629983568D-03, 0.28319080123804D-03, &
                                                          -0.60706301565874D-03, -0.18990068218419D-01, -0.32529748770505D-01, &
                                                          -0.21841717175414D-01, -0.52838357969930D-04, -0.47184321073267D-03, &
                                                          -0.30001780793026D-03, 0.47661393906987D-04, -0.44141845330846D-05, &
                                                          -0.72694996297594D-15, -0.31679644845054D-04, -0.28270797985312D-05, &
                                                          -0.85205128120103D-09, -0.22425281908000D-05, -0.65171222895601D-06, &
                                                          -0.14341729937924D-12, -0.40516996860117D-06, -0.12734301741641D-08, &
                                                          -0.17424871230634D-09, -0.68762131295531D-18, 0.14478307828521D-19, &
                                                          0.26335781662795D-22, -0.11947622640071D-22, 0.18228094581404D-23, &
                                                          -0.93537087292458D-25/)
    !
    !> Star pressure for region 1 in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_PSTAR = 16.53D+00
    !
    !> Star temperature for region 1 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_TSTAR = 1386.0D+00

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function computes the dimensionless specific Gibbs free energy.
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The dimensionless specific Gibbs free energy
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gam(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: gam
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pia, taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 7.1D+00, b = 1.222D+00
        !
        ! Polynomial expression
        !
        pia = a - pi
        taub = tau - b
        gam = SUM(n * (pia**I) * (taub**J))

    END FUNCTION gam
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the specific Gibbs free energy
    !> with respect to the parameter pi.
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d(Gibbs)/d(pi)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gampi(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: gampi
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pia, taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 7.1D+00, b = 1.222D+00
        !
        ! Polynomial expression
        !
        pia = a - pi
        taub = tau - b
        gampi = SUM(-n * I * ((pia)**(I - 1)) * (taub**J))

    END FUNCTION gampi
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the specific Gibbs
    !> free energy with respect to the parameter pi.
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d2(Gibbs)/d(pi)2
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gampipi(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: gampipi
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pia, taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 7.1D+00, b = 1.222D+00
        !
        ! Polynomial expression
        !
        pia = a - pi
        taub = tau - b
        gampipi = SUM(n * I * (I - 1) * (pia**(I - 2)) * (taub**J))

    END FUNCTION gampipi
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the specific Gibbs free energy
    !> with respect to the parameter tau.
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d(Gibbs)/d(tau)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamtau(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: gamtau
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pia, taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 7.1D+00, b = 1.222D+00
        !
        ! Polynomial expression
        !
        pia = a - pi
        taub = tau - b
        gamtau = SUM(n * (pia**I) * J * (taub**(J - 1)))

    END FUNCTION gamtau
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the specific Gibbs
    !> free energy with respect to the parameter tau.
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d2(Gibbs)/d(tau)2
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gamtautau(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: gamtautau
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pia, taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 7.1D+00, b = 1.222D+00
        !
        ! Polynomial expression
        !
        pia = a - pi
        taub = tau - b
        gamtautau = SUM(n * (pia**I) * J * (J - 1) * (taub**(J - 2)))

    END FUNCTION gamtautau
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the specific Gibbs
    !> free energy with respect to the parameters pi and tau.
    !
    !> @param[in] pi  Dimensionless pressure parameter
    !> @param[in] tau Dimensionless temperature parameter
    !
    !> @return The derivative d2(Gibbs)/d(pi)d(tau)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION gampitau(pi, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: gampitau
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pia, taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 7.1D+00, b = 1.222D+00
        !
        ! Polynomial expression
        !
        pia = a - pi
        taub = tau - b
        gampitau = SUM(-n * I * (pia**(I - 1)) * J * (taub**(J - 1)))

    END FUNCTION gampitau
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
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: specific_internal_energy
        REAL(KIND=REAL_HIGH) :: pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the specific internal energy
        !
        specific_internal_energy = IAPWS97_R * temperature * (tau * gamtau(pi, tau) - pi * gampi(pi, tau))

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
        ! Variables
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        REAL(KIND=REAL_HIGH) :: specific_volume
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the specific volume
        !
        specific_volume = IAPWS97_R * temperature * pi * gampi(pi, tau) / (pressure * mega)

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
        REAL(KIND=REAL_HIGH) :: pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the specific entropy
        !
        specific_entropy = IAPWS97_R * (tau * gamtau(pi, tau) - gam(pi, tau))

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
        REAL(KIND=REAL_HIGH) :: pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the specific enthalpy
        !
        specific_enthalpy = IAPWS97_R * temperature * tau * gamtau(pi, tau)

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
        REAL(KIND=REAL_HIGH) :: res1, res2, res3, res4
        REAL(KIND=REAL_HIGH) :: pi, tau, gp
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the speed of sound
        !
        gp = gampi(pi, tau)
        res1 = IAPWS97_R * temperature
        res2 = (gp - tau * gampitau(pi, tau))**2
        res3 = tau * tau * gamtautau(pi, tau)
        res4 = res2 / res3
        speed_of_sound = gp * sqrt(res1 / (res4 - gampipi(pi, tau)))

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
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the specific isobaric heat capacity
        !
        specific_isobaric_heat_capacity = -IAPWS97_R * tau * tau * gamtautau(pi, tau)

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
        REAL(KIND=REAL_HIGH) :: res1, res2
        REAL(KIND=REAL_HIGH) :: pi, tau
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the specific isochoric heat capacity
        !
        res1 = -tau * tau * gamtautau(pi, tau)
        res2 = (gampi(pi, tau) - tau * gampitau(pi, tau))**2
        specific_isochoric_heat_capacity = IAPWS97_R * (res1 + res2 / gampipi(pi, tau))

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
        pi = pressure / REGION_1_PSTAR
        tau = REGION_1_TSTAR / temperature
        !
        ! Compute the specific Gibbs free energy
        !
        specific_gibbs_free_energy = IAPWS97_R * temperature * gam(pi, tau)

    END FUNCTION specific_gibbs_free_energy
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and enthalpy (backward equation)
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ph(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:20) :: II = &
                                                              (/0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:20) :: JJ = &
                                                              (/0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, &
                                                                32, 32, 32/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:20) :: nn = &
                                                            (/-238.72489924521D+00, 404.21188637945D+00, 113.49746881718D+00, &
                                                              -5.8457616048039D+00, -1.528548241314D-04, -1.0866707695377D-06, &
                                                              -13.391744872602D+00, 43.211039183559D+00, -54.010067170506D+00, &
                                                              30.535892203916D+00, -6.5964749423638D+00, 9.3965400878363D-03, &
                                                              1.157364750534D-07, -2.5858641282073D-05, -4.0644363084799D-09, &
                                                              6.6456186191635D-08, 8.0670734103027D-11, -9.3477771213947D-13, &
                                                              5.8265442020601D-15, -1.5020185953503D-17/)
        !
        !> Star pressure for the backward temperature(pressure, enthalpy)
        !> equation in region 1 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_TPH_PSTAR = 1.0D+00
        !
        !> Star temperature for the backward temperature(pressure, enthalpy)
        !> equation in region 1 in [K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_TPH_TSTAR = 1.0D+00
        !
        !> Star enthalpy for the backward temperature(pressure, enthalpy)
        !> equation in region 1 in [J/kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_TPH_HSTAR = 2500.0D+03
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: temperature_ph
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, e1
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_1_TPH_PSTAR
        e1 = one + (enthalpy / REGION_1_TPH_HSTAR)
        temperature_ph = REGION_1_TPH_TSTAR * SUM(nn * (pi**II) * (e1**JJ))

    END FUNCTION temperature_ph
    !-------------------------------------------------------------------------
    !
    !> This function computes the pressure with respect to enthalpy
    !> and entropy (backward equation)
    !
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The pressure [MPa]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION pressure_hs(enthalpy, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:19) :: II = &
                                                              (/0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:19) :: JJ = &
                                                              (/0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:19) :: nn = &
                                                         (/-0.691997014660582D+00, -0.183612548787560D+02, -0.928332409297335D+01, &
                                                             0.659639569909906D+02, -0.162060388912024D+02, 0.450620017338667D+03, &
                                                              0.854680678224170D+03, 0.607523214001162D+04, 0.326487682621856D+02, &
                                                           -0.269408844582931D+02, -0.319947848334300D+03, -0.928354307043320D+03, &
                                                            0.303634537455249D+02, -0.650540422444146D+02, -0.430991316516130D+04, &
                                                             -0.747512324096068D+03, 0.730000345529245D+03, 0.114284032569021D+04, &
                                                              -0.436407041874559D+03/)
        !
        !> Star pressure for the backward pressure(enthalpy, entropy)
        !> equation in region 1 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_PHS_PSTAR = 100.0D+00
        !
        !> Star enthalpy for the backward pressure(enthalpy, entropy)
        !> equation in region 1 in [K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_PHS_HSTAR = 3400.0D+03
        !
        !> Star entropy for the backward pressure(enthalpy, entropy)
        !> equation in region 1 in [J/kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_1_PHS_SSTAR = 7.6D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.05D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: enthalpy, entropy
        REAL(KIND=REAL_HIGH) :: pressure_hs
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: eta, sigma
        !
        ! Compute the pressure
        !
        eta = enthalpy / REGION_1_PHS_HSTAR + a
        sigma = entropy / REGION_1_PHS_SSTAR + a
        pressure_hs = REGION_1_PHS_PSTAR * SUM(nn * (eta**II) * (sigma**JJ))

    END FUNCTION pressure_hs

END MODULE SteamRegion1
