!-----------------------------------------------------------------------
!
!> @file SteamRegion2.f
!
!> @details
!> This module implements the polynomials which compute the thermodynamic
!> properties of steam in region 2 according to the revised release on the
!> IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
!> Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamRegion2

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
    ! Backward functions
    !
    PUBLIC temperature_ph_region2a
    PUBLIC temperature_ph_region2b
    PUBLIC temperature_ph_region2c
    PUBLIC temperature_ps_region2a
    PUBLIC temperature_ps_region2b
    PUBLIC temperature_ps_region2c
    PUBLIC pressure_hs_region2a
    PUBLIC pressure_hs_region2b
    PUBLIC pressure_hs_region2c
    !
    !> Constant coefficients "J" for the ideal gas part
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:9) :: J = &
                                                         (/0, 1, -5, -4, -3, -2, -1, 2, 3/)
    !
    !> Constant coefficients "n" for the ideal gas part
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:9) :: n = &
                                                       (/-0.96927686500217D+01, 0.10086655968018D+02, -0.56087911283020D-02, &
                                                         0.71452738081455D-01, -0.40710498223928D+00, 0.14240819171444D+01, &
                                                         -0.43839511319450D+01, -0.28408632460772D+00, 0.21268463753307D-01/)
    !
    !> Constant coefficients "Ir" for the residual part
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:43) :: Ir = &
                                                          (/1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, &
                                                            7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, &
                                                            24, 24, 24/)
    !
    !> Constant coefficients "Jr" for the residual part
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:43) :: Jr = &
                                                          (/0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, &
                                                            35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, &
                                                            53, 39, 26, 40, 58/)
    !
    !> Constant coefficients "nr" for the residual part
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:43) :: nr = &
                                                        (/-0.17731742473213D-02, -0.17834862292358D-01, -0.45996013696365D-01, &
                                                          -0.57581259083432D-01, -0.50325278727930D-01, -0.33032641670203D-04, &
                                                          -0.18948987516315D-03, -0.39392777243355D-02, -0.43797295650573D-01, &
                                                          -0.26674547914087D-04, 0.20481737692309D-07, 0.43870667284435D-06, &
                                                          -0.32277677238570D-04, -0.15033924542148D-02, -0.40668253562649D-01, &
                                                          -0.78847309559367D-09, 0.12790717852285D-07, 0.48225372718507D-06, &
                                                          0.22922076337661D-05, -0.16714766451061D-10, -0.21171472321355D-02, &
                                                          -0.23895741934104D+02, -0.59059564324270D-17, -0.12621808899101D-05, &
                                                          -0.38946842435739D-01, 0.11256211360459D-10, -0.82311340897998D+01, &
                                                          0.19809712802088D-07, 0.10406965210174D-18, -0.10234747095929D-12, &
                                                          -0.10018179379511D-08, -0.80882908646985D-10, 0.10693031879409D+00, &
                                                          -0.33662250574171D+00, 0.89185845355421D-24, 0.30629316876232D-12, &
                                                          -0.42002467698208D-05, -0.59056029685639D-25, 0.37826947613457D-05, &
                                                          -0.12768608934681D-14, 0.73087610595061D-28, 0.55414715350778D-16, &
                                                          -0.94369707241210D-06/)
    !
    !> Star pressure for region 2 in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PSTAR = 1.0D+00
    !
    !> Star temperature for region 2 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_TSTAR = 540.0D+00

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
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 0.5D+00
        !
        ! Polynomial expression
        !
        taub = tau - b
        gamr = SUM(nr * (pi**Ir) * (taub**Jr))

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
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 0.5D+00
        !
        ! Polynomial expression
        !
        taub = tau - b
        gamrpi = SUM(nr * Ir * (pi**(Ir - 1)) * (taub**Jr))

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
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 0.5D+00
        !
        ! Polynomial expression
        !
        taub = tau - b
        gamrpipi = SUM(nr * Ir * (Ir - 1) * (pi**(Ir - 2)) * (taub**Jr))

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
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 0.5D+00
        !
        ! Polynomial expression
        !
        taub = tau - b
        gamrtau = SUM(nr * (pi**Ir) * Jr * (taub**(Jr - 1)))

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
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 0.5D+00
        !
        ! Polynomial expression
        !
        taub = tau - b
        gamrtautau = SUM(nr * (pi**Ir) * Jr * (Jr - 1) * (taub**(Jr - 2)))

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
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pi, tau
        REAL(KIND=REAL_HIGH) :: gamrpitau
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: taub
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 0.5D+00
        !
        ! Polynomial expression
        !
        taub = tau - b
        gamrpitau = SUM(nr * Ir * (pi**(Ir - 1)) * Jr * (taub**(Jr - 1)))

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
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
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
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
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
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
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
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
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
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
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
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
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
        REAL(KIND=REAL_HIGH) :: pi, tau, res1, res2, res3
        !
        ! Compute the dimensionless parameters pi and tau
        !
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
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
        pi = pressure / REGION_2_PSTAR
        tau = REGION_2_TSTAR / temperature
        !
        ! Compute the specific Gibbs free energy
        !
        specific_gibbs_free_energy = IAPWS97_R * temperature * (gam0(pi, tau) + gamr(pi, tau))

    END FUNCTION specific_gibbs_free_energy
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and enthalpy (backward equation) for the sub-region 2a
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ph_region2a(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:34) :: II = &
                                                              (/0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, &
                                                                2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:34) :: JJ = &
                                                              (/0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, &
                                                                36, 38, 40, 42, 44, 24, 44, 12, 32, 44, 32, 36, 42, 34, &
                                                                44, 28/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:34) :: nn = &
                                                            (/1089.8952318288D+00, 849.51654495535D+00, -107.81748091826D+00, &
                                                              33.153654801263D+00, -7.4232016790248D+00, 11.765048724356D+00, &
                                                              1.844574935579D+00, -4.1792700549624D+00, 6.2478196935812D+00, &
                                                              -17.344563108114D+00, -200.58176862096D+00, 271.96065473796D+00, &
                                                              -455.11318285818D+00, 3091.9688604755D+00, 252266.40357872D+00, &
                                                              -6.1707422868339D-03, -0.31078046629583D+00, 11.670873077107D+00, &
                                                              128127984.04046D+00, -985549096.23276D+00, 2822454697.3002D+00, &
                                                              -3594897141.0703D+00, 1722734991.3197D+00, -13551.334240775D+00, &
                                                              12848734.66465D+00, 1.3865724283226D+00, 235988.32556514D+00, &
                                                              -13105236.545054D+00, 7399.9835474766D+00, -551966.9703006D+00, &
                                                              3715408.5996233D+00, 19127.72923966D+00, -415351.64835634D+00, &
                                                              -62.459855192507D+00/)
        !
        !> Star pressure for the backward temperature(pressure, enthalpy)
        !> equation in region 2 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PH_PSTAR = 1.0D+00
        !
        !> Star enthalpy for the backward temperature(pressure, enthalpy)
        !> equation in region 2 in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PH_HSTAR = 2000.0D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 2.1D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: temperature_ph_region2a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_2_PH_PSTAR
        eta = enthalpy / REGION_2_PH_HSTAR - b
        temperature_ph_region2a = SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION temperature_ph_region2a
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and enthalpy (backward equation) for the sub-region 2b
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ph_region2b(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:38) :: II = &
                                                              (/0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, &
                                                                2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 7, 7, 9, 9/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:38) :: JJ = &
                                                              (/0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, &
                                                                2, 8, 18, 40, 1, 2, 12, 24, 2, 12, 18, 24, 28, 40, 18, 24, &
                                                                40, 28, 2, 28, 1, 40/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:38) :: nn = &
                                                            (/1489.5041079516D+00, 743.07798314034D+00, -97.708318797837D+00, &
                                                              2.4742464705674D+00, -0.63281320016026D+00, 1.1385952129658D+00, &
                                                              -0.47811863648625D+00, 8.5208123431544D-03, 0.93747147377932D+00, &
                                                              3.3593118604916D+00, 3.3809355601454D+00, 0.16844539671904D+00, &
                                                              0.73875745236695D+00, -0.47128737436186D+00, 0.15020273139707D+00, &
                                                            -0.002176411421975D+00, -0.021810755324761D+00, -0.10829784403677D+00, &
                                                              -0.046333324635812D+00, 7.1280351959551D-05, 1.1032831789999D-04, &
                                                              1.8955248387902D-04, 3.0891541160537D-03, 1.3555504554949D-03, &
                                                              2.8640237477456D-07, -1.0779857357512D-05, -7.6462712454814D-05, &
                                                              1.4052392818316D-05, -3.1083814331434D-05, -1.0302738212103D-06, &
                                                              2.821728163504D-07, 1.2704902271945D-06, 7.3803353468292D-08, &
                                                              -1.1030139238909D-08, -8.1456365207833D-14, -2.5180545682962D-11, &
                                                              -1.7565233969407D-18, 8.6934156344163D-15/)
        !
        !> Star pressure for the backward temperature(pressure, enthalpy)
        !> equation in region 2 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PH_PSTAR = 1.0D+00
        !
        !> Star enthalpy for the backward temperature(pressure, enthalpy)
        !> equation in region 2 in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PH_HSTAR = 2000.0D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 2.0D+00, b = 2.6D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: temperature_ph_region2b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_2_PH_PSTAR - a
        eta = enthalpy / REGION_2_PH_HSTAR - b
        temperature_ph_region2b = SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION temperature_ph_region2b
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and enthalpy (backward equation) for the sub-region 2c
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ph_region2c(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:23) :: II = &
                                                              (/-7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, &
                                                                6, 6, 6, 6, 6, 6, 6, 6/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:23) :: JJ = &
                                                              (/0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, &
                                                                10, 12, 16, 20, 22/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:23) :: nn = &
                                                            (/-3236839855524.2D+00, 7326335090218.1D+00, 358250899454.47D+00, &
                                                              -583401318515.9D+00, -10783068217.47D+00, 20825544563.171D+00, &
                                                              610747.83564516D+00, 859777.2253558D+00, -25745.72360417D+00, &
                                                              31081.088422714D+00, 1208.2315865936D+00, 482.19755109255D+00, &
                                                              3.7966001272486D+00, -10.842984880077D+00, -0.04536417267666D+00, &
                                                              1.4559115658698D-13, 1.126159740723D-12, -1.7804982240686D-11, &
                                                              1.2324579690832D-07, -1.1606921130984D-06, 2.7846367088554D-05, &
                                                              -5.9270038474176D-04, 1.2918582991878D-03/)
        !
        !> Star pressure for the backward temperature(pressure, enthalpy)
        !> equation in region 2 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PH_PSTAR = 1.0D+00
        !
        !> Star enthalpy for the backward temperature(pressure, enthalpy)
        !> equation in region 2 in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PH_HSTAR = 2000.0D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 25.0D+00, b = 1.8D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: temperature_ph_region2c
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_2_PH_PSTAR + a
        eta = enthalpy / REGION_2_PH_HSTAR - b
        temperature_ph_region2c = SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION temperature_ph_region2c
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and entropy (backward equation) for the sub-region 2a
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] entropy  The steam enthalpy [J/Kg.K]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ps_region2a(pressure, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:46) :: II = &
                                                            (/-1.50D+00, -1.50D+00, -1.50D+00, -1.50D+00, -1.50D+00, -1.50D+00, &
                                                              -1.25D+00, -1.25D+00, -1.25D+00, -1.00D+00, -1.00D+00, -1.00D+00, &
                                                              -1.00D+00, -1.00D+00, -1.00D+00, -0.75D+00, -0.75D+00, -0.50D+00, &
                                                              -0.50D+00, -0.50D+00, -0.50D+00, -0.25D+00, -0.25D+00, -0.25D+00, &
                                                              -0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00, 0.50D+00, &
                                                              0.50D+00, 0.50D+00, 0.50D+00, 0.50D+00, 0.50D+00, 0.50D+00, &
                                                              0.75D+00, 0.75D+00, 0.75D+00, 0.75D+00, 1.00D+00, 1.00D+00, &
                                                              1.25D+00, 1.25D+00, 1.50D+00, 1.50D+00/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:46) :: JJ = &
                                                              (/-24, -23, -19, -13, -11, -10, -19, -15, -6, -26, -21, -17, &
                                                                -16, -9, -8, -15, -14, -26, -13, -9, -7, -27, -25, -11, -6, &
                                                                1, 4, 8, 11, 0, 1, 5, 6, 10, 14, 16, 0, 4, 9, 17, 7, 18, 3, &
                                                                15, 5, 18/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:46) :: nn = &
                                                            (/-0.39235983861984D+06, 0.51526573827270D+06, &
                                                              0.40482443161048D+05, -0.32193790923902D+03, &
                                                              0.96961424218694D+02, -0.22867846371773D+02, &
                                                              -0.44942914124357D+06, -0.50118336020166D+04, &
                                                              0.35684463560015D+00, 0.44235335848190D+05, &
                                                              -0.13673388811708D+05, 0.42163260207864D+06, &
                                                              0.22516925837475D+05, 0.47442144865646D+03, &
                                                              -0.14931130797647D+03, -0.19781126320452D+06, &
                                                              -0.23554399470760D+05, -0.19070616302076D+05, &
                                                              0.55375669883164D+05, 0.38293691437363D+04, &
                                                              -0.60391860580567D+03, 0.19363102620331D+04, &
                                                              0.42660643698610D+04, -0.59780638872718D+04, &
                                                              -0.70401463926862D+03, 0.33836784107553D+03, &
                                                              0.20862786635187D+02, 0.33834172656196D-01, &
                                                              -0.43124428414893D-04, 0.16653791356412D+03, &
                                                              -0.13986292055898D+03, -0.78849547999872D+00, &
                                                              0.72132411753872D-01, -0.59754839398283D-02, &
                                                              -0.12141358953904D-04, 0.23227096733871D-06, &
                                                              -0.10538463566194D+02, 0.20718925496502D+01, &
                                                              -0.72193155260427D-01, 0.20749887081120D-06, &
                                                              -0.18340657911379D-01, 0.29036272348696D-06, &
                                                              0.21037527893619D+00, 0.25681239729999D-03, &
                                                              -0.12799002933781D-01, -0.82198102652018D-05/)
        !
        !> Star pressure for the backward temperature(pressure, entropy)
        !> equation in region 2 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PS_PSTAR = 1.0D+00
        !
        !> Star entropy for the backward temperature(pressure, entropy)
        !> equation in region 2 in [J/Kg.K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PS_SSTAR = 2000.0D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 2.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, entropy
        REAL(KIND=REAL_HIGH) :: temperature_ps_region2a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, sigma
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_2_PS_PSTAR
        sigma = entropy / REGION_2_PS_SSTAR - b
        temperature_ps_region2a = SUM(nn * (pi**II) * (sigma**JJ))

    END FUNCTION temperature_ps_region2a
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and entropy (backward equation) for the sub-region 2b
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] entropy  The steam enthalpy [J/Kg.K]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ps_region2b(pressure, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:44) :: II = &
                                                              (/-6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, &
                                                                -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, &
                                                                2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:44) :: JJ = &
                                                              (/0, 11, 0, 11, 0, 1, 11, 0, 1, 11, 12, 0, 1, 6, 10, 0, 1, &
                                                                5, 8, 9, 0, 1, 2, 4, 5, 6, 9, 0, 1, 2, 3, 7, 8, 0, 1, 5, &
                                                                0, 1, 3, 0, 1, 0, 1, 2/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:44) :: nn = &
                                                            (/0.31687665083497D+06, 0.20864175881858D+02, &
                                                              -0.39859399803599D+06, -0.21816058518877D+02, &
                                                              0.22369785194242D+06, -0.27841703445817D+04, &
                                                              0.99207436071480D+01, -0.75197512299157D+05, &
                                                              0.29708605951158D+04, -0.34406878548526D+01, &
                                                              0.38815564249115D+00, 0.17511295085750D+05, &
                                                              -0.14237112854449D+04, 0.10943803364167D+01, &
                                                              0.89971619308495D+00, -0.33759740098958D+04, &
                                                              0.47162885818355D+03, -0.19188241993679D+01, &
                                                              0.41078580492196D+00, -0.33465378172097D+00, &
                                                              0.13870034777505D+04, -0.40663326195838D+03, &
                                                              0.41727347159610D+02, 0.21932549434532D+01, &
                                                              -0.10320050009077D+01, 0.35882943516703D+00, &
                                                              0.52511453726066D-02, 0.12838916450705D+02, &
                                                              -0.28642437219381D+01, 0.56912683664855D+00, &
                                                              -0.99962954584931D-01, -0.32632037778459D-02, &
                                                              0.23320922576723D-03, -0.15334809857450D+00, &
                                                              0.29072288239902D-01, 0.37534702741167D-03, &
                                                              0.17296691702411D-02, -0.38556050844504D-03, &
                                                              -0.35017712292608D-04, -0.14566393631492D-04, &
                                                              0.56420857267269D-05, 0.41286150074605D-07, &
                                                              -0.20684671118824D-07, 0.16409393674725D-08/)
        !
        !> Star pressure for the backward temperature(pressure, entropy)
        !> equation in region 2 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PS_PSTAR = 1.0D+00
        !
        !> Star entropy for the backward temperature(pressure, entropy)
        !> equation in region 2 in [J/Kg.K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PS_SSTAR = 0.7853D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 10.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, entropy
        REAL(KIND=REAL_HIGH) :: temperature_ps_region2b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, sigma
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_2_PS_PSTAR
        sigma = b - entropy / REGION_2_PS_SSTAR
        temperature_ps_region2b = SUM(nn * (pi**II) * (sigma**JJ))

    END FUNCTION temperature_ps_region2b
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and entropy (backward equation) for the sub-region 2c
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] entropy  The steam enthalpy [J/Kg.K]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ps_region2c(pressure, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:30) :: II = &
                                                              (/-2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, &
                                                                4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:30) :: JJ = &
                                                              (/0, 1, 0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, &
                                                                1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:30) :: nn = &
                                                            (/0.90968501005365D+03, 0.24045667088420D+04, &
                                                              -0.59162326387130D+03, 0.54145404128074D+03, &
                                                              -0.27098308411192D+03, 0.97976525097926D+03, &
                                                              -0.46966772959435D+03, 0.14399274604723D+02, &
                                                              -0.19104204230429D+02, 0.53299167111971D+01, &
                                                              -0.21252975375934D+02, -0.31147334413760D+00, &
                                                              0.60334840894623D+00, -0.42764839702509D-01, &
                                                              0.58185597255259D-02, -0.14597008284753D-01, &
                                                              0.56631175631027D-02, -0.76155864584577D-04, &
                                                              0.22440342919332D-03, -0.12561095013413D-04, &
                                                              0.63323132660934D-06, -0.20541989675375D-05, &
                                                              0.36405370390082D-07, -0.29759897789215D-08, &
                                                              0.10136618529763D-07, 0.59925719692351D-11, &
                                                              -0.20677870105164D-10, -0.20874278181886D-10, &
                                                              0.10162166825089D-09, -0.16429828281347D-09/)
        !
        !> Star pressure for the backward temperature(pressure, entropy)
        !> equation in region 2 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PS_PSTAR = 1.0D+00
        !
        !> Star entropy for the backward temperature(pressure, entropy)
        !> equation in region 2 in [J/Kg.K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_PS_SSTAR = 2.92510D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: b = 2.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, entropy
        REAL(KIND=REAL_HIGH) :: temperature_ps_region2c
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, sigma
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_2_PS_PSTAR
        sigma = b - entropy / REGION_2_PS_SSTAR
        temperature_ps_region2c = SUM(nn * (pi**II) * (sigma**JJ))

    END FUNCTION temperature_ps_region2c
    !-------------------------------------------------------------------------
    !
    !> This function computes the pressure with respect to enthalpy
    !> and entropy (backward equation) for the sub-region 2a
    !
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The pressure [MPa]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION pressure_hs_region2a(enthalpy, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:29) :: II = &
                                                              (/0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, &
                                                                3, 3, 3, 3, 3, 4, 5, 5, 6, 7/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:29) :: JJ = &
                                                              (/1, 3, 6, 16, 20, 22, 0, 1, 2, 3, 5, 6, 10, 16, 20, 22, 3, &
                                                                16, 20, 0, 2, 3, 6, 16, 16, 3, 16, 3, 1/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:29) :: nn = &
                                                            (/-0.182575361923032D-01, -0.125229548799536D+00, &
                                                              0.592290437320145D+00, 0.604769706185122D+01, &
                                                              0.238624965444474D+03, -0.298639090222922D+03, &
                                                              0.512250813040750D-01, -0.437266515606486D+00, &
                                                              0.413336902999504D+00, -0.516468254574773D+01, &
                                                              -0.557014838445711D+01, 0.128555037824478D+02, &
                                                              0.114144108953290D+02, -0.119504225652714D+03, &
                                                              -0.284777985961560D+04, 0.431757846408006D+04, &
                                                              0.112894040802650D+01, 0.197409186206319D+04, &
                                                              0.151612444706319D+04, 0.141324451421235D-01, &
                                                              0.585501282219601D+00, -0.297258075863012D+01, &
                                                              0.594567314847319D+01, -0.623656565798905D+04, &
                                                              0.965986235133332D+04, 0.681500934948134D+01, &
                                                              -0.633207286824489D+04, -0.558919224465760D+01, &
                                                              0.400645798472063D-01/)
        !
        !> Star pressure for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2a in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2A_PHS_PSTAR = 4.0D+00
        !
        !> Star enthalpy for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2a in [K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2A_PHS_HSTAR = 4200.0D+03
        !
        !> Star entropy for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2a in [J/kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2A_PHS_SSTAR = 12.0D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.5D+00, b = 1.2D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: enthalpy, entropy
        REAL(KIND=REAL_HIGH) :: pressure_hs_region2a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: eta, sigma
        !
        ! Compute the pressure
        !
        eta = enthalpy / REGION_2A_PHS_HSTAR - a
        sigma = entropy / REGION_2A_PHS_SSTAR - b
        pressure_hs_region2a = REGION_2A_PHS_PSTAR * SUM(nn * (eta**II) * (sigma**JJ))**4

    END FUNCTION pressure_hs_region2a
    !-------------------------------------------------------------------------
    !
    !> This function computes the pressure with respect to enthalpy
    !> and entropy (backward equation) for the sub-region 2b
    !
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The pressure [MPa]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION pressure_hs_region2b(enthalpy, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:33) :: II = &
                                                              (/0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, &
                                                                5, 5, 6, 6, 6, 7, 7, 8, 8, 8, 8, 12, 14/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:33) :: JJ = &
                                                              (/0, 1, 2, 4, 8, 0, 1, 2, 3, 5, 12, 1, 6, 18, 0, 1, 7, 12, &
                                                                1, 16, 1, 12, 1, 8, 18, 1, 16, 1, 3, 14, 18, 10, 16/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:33) :: nn = &
                                                            (/0.801496989929495D-01, -0.543862807146111D+00, &
                                                              0.337455597421283D+00, 0.890555451157450D+01, &
                                                              0.313840736431485D+03, 0.797367065977789D+00, &
                                                              -0.121616973556240D+01, 0.872803386937477D+01, &
                                                              -0.169769781757602D+02, -0.186552827328416D+03, &
                                                              0.951159274344237D+05, -0.189168510120494D+02, &
                                                              -0.433407037194840D+04, 0.543212633012715D+09, &
                                                              0.144793408386013D+00, 0.128024559637516D+03, &
                                                              -0.672309534071268D+05, 0.336972380095287D+08, &
                                                              -0.586634196762720D+03, -0.221403224769889D+11, &
                                                              0.171606668708389D+04, -0.570817595806302D+09, &
                                                              -0.312109693178482D+04, 0.207841384633010D+07, &
                                                              0.305605946157786D+13, 0.322157004314333D+04, &
                                                              0.326810259797295D+12, -0.144104158934487D+04, &
                                                              0.410694867802691D+03, 0.109077066873024D+12, &
                                                              -0.247964654258893D+14, 0.188801906865134D+10, &
                                                              -0.123651009018773D+15/)
        !
        !> Star pressure for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2b in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_PHS_PSTAR = 100.0D+00
        !
        !> Star enthalpy for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2b in [K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_PHS_HSTAR = 4100.0D+03
        !
        !> Star entropy for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2b in [J/kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_PHS_SSTAR = 7.9D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.6D+00, b = 1.01D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: enthalpy, entropy
        REAL(KIND=REAL_HIGH) :: pressure_hs_region2b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: eta, sigma
        !
        ! Compute the pressure
        !
        eta = enthalpy / REGION_2B_PHS_HSTAR - a
        sigma = entropy / REGION_2B_PHS_SSTAR - b
        pressure_hs_region2b = REGION_2B_PHS_PSTAR * SUM(nn * (eta**II) * (sigma**JJ))**4

    END FUNCTION pressure_hs_region2b
    !-------------------------------------------------------------------------
    !
    !> This function computes the pressure with respect to enthalpy
    !> and entropy (backward equation) for the sub-region 2c
    !
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The pressure [MPa]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION pressure_hs_region2c(enthalpy, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:31) :: II = &
                                                              (/0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, &
                                                                3, 3, 4, 5, 5, 5, 5, 6, 6, 10, 12, 16/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:31) :: JJ = &
                                                              (/0, 1, 2, 3, 4, 8, 0, 2, 5, 8, 14, 2, 3, 7, 10, 18, 0, 5, &
                                                                8, 16, 18, 18, 1, 4, 6, 14, 8, 18, 7, 7, 10/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:31) :: nn = &
                                                            (/0.112225607199012D+00, -0.339005953606712D+01, &
                                                              -0.320503911730094D+02, -0.197597305104900D+03, &
                                                              -0.407693861553446D+03, 0.132943775222331D+05, &
                                                              0.170846839774007D+01, 0.373694198142245D+02, &
                                                              0.358144365815434D+04, 0.423014446424664D+06, &
                                                              -0.751071025760063D+09, 0.523446127607898D+02, &
                                                              -0.228351290812417D+03, -0.960652417056937D+06, &
                                                              -0.807059292526074D+08, 0.162698017225669D+13, &
                                                              0.772465073604171D+00, 0.463929973837746D+05, &
                                                              -0.137317885134128D+08, 0.170470392630512D+13, &
                                                              -0.251104628187308D+14, 0.317748830835520D+14, &
                                                              0.538685623675312D+02, -0.553089094625169D+05, &
                                                              -0.102861522421405D+07, 0.204249418756234D+13, &
                                                              0.273918446626977D+09, -0.263963146312685D+16, &
                                                              -0.107890854108088D+10, -0.296492620980124D+11, &
                                                              -0.111754907323424D+16/)
        !
        !> Star pressure for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2c in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2C_PHS_PSTAR = 100.0D+00
        !
        !> Star enthalpy for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2c in [K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2C_PHS_HSTAR = 3500.0D+03
        !
        !> Star entropy for the backward pressure(enthalpy, entropy)
        !> equation in region 2, subregion 2c in [J/kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2C_PHS_SSTAR = 5.9D+03
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.7D+00, b = 1.1D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: enthalpy, entropy
        REAL(KIND=REAL_HIGH) :: pressure_hs_region2c
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: eta, sigma
        !
        ! Compute the pressure
        !
        eta = enthalpy / REGION_2C_PHS_HSTAR - a
        sigma = entropy / REGION_2C_PHS_SSTAR - b
        pressure_hs_region2c = REGION_2C_PHS_PSTAR * SUM(nn * (eta**II) * (sigma**JJ))**4

    END FUNCTION pressure_hs_region2c

END MODULE SteamRegion2
