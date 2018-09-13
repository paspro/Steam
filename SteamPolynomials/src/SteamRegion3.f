!-----------------------------------------------------------------------
!
!> @file SteamRegion3.f
!
!> @details
!> This module implements the polynomials which compute the thermodynamic
!> properties of steam in region 3 according to the revised release on the
!> IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
!> Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamRegion3

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Thermodynamic properties computed by polynomials
    !
    PUBLIC pressure
    PUBLIC specific_internal_energy
    PUBLIC specific_entropy
    PUBLIC specific_enthalpy
    PUBLIC specific_isobaric_heat_capacity
    PUBLIC specific_isochoric_heat_capacity
    PUBLIC ratio_of_specific_heats
    PUBLIC speed_of_sound
    PUBLIC specific_helmoltz_free_energy
    !
    ! Backward functions
    !
    PUBLIC temperature_ph_region3a
    PUBLIC temperature_ph_region3b
    PUBLIC specific_volume_ph_region3a
    PUBLIC specific_volume_ph_region3b
    PUBLIC temperature_ps_region3a
    PUBLIC temperature_ps_region3b
    PUBLIC specific_volume_ps_region3a
    PUBLIC specific_volume_ps_region3b
    !
    !> Constant coefficients "I"
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:39) :: I = &
                                                          (/0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, &
                                                            3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11/)
    !
    !> Constant coefficients "J"
    !
    INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:39) :: J = &
                                                          (/0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, &
                                                            2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, &
                                                            0, 1, 26/)
    !
    !> Constant coefficients "n"
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:39) :: n = &
                                                        (/-0.15732845290239D+02, 0.20944396974307D+02, -0.76867707878716D+01, &
                                                          0.26185947787954D+01, -0.28080781148620D+01, 0.12053369696517D+01, &
                                                          -0.84566812812502D-02, -0.12654315477714D+01, -0.11524407806681D+01, &
                                                          0.88521043984318D+00, -0.64207765181607D+00, 0.38493460186671D+00, &
                                                          -0.85214708824206D+00, 0.48972281541877D+01, -0.30502617256965D+01, &
                                                          0.39420536879154D-01, 0.12558408424308D+00, -0.27999329698710D+00, &
                                                          0.13899799569460D+01, -0.20189915023570D+01, -0.82147637173963D-02, &
                                                          -0.47596035734923D+00, 0.43984074473500D-01, -0.44476435428739D+00, &
                                                          0.90572070719733D+00, 0.70522450087967D+00, 0.10770512626332D+00, &
                                                          -0.32913623258954D+00, -0.50871062041158D+00, -0.22175400873096D-01, &
                                                          0.94260751665092D-01, 0.16436278447961D+00, -0.13503372241348D-01, &
                                                          -0.14834345352472D-01, 0.57922953628084D-03, 0.32308904703711D-02, &
                                                          0.80964802996215D-04, -0.16557679795037D-03, -0.44923899061815D-04/)
    !
    !> Constant coefficient "n1"
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: n1 = 0.10658070028513D+01
    !
    !> Star density for region 3 in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_RHOSTAR = IAPWS97_RHOCRIT
    !
    !> Star temperature for region 3 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_TSTAR = IAPWS97_TCRIT

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function computes the dimensionless specific Helmholtz free energy.
    !
    !> @param[in] delta Dimensionless density parameter
    !> @param[in] tau   Dimensionless temperature parameter
    !
    !> @return The dimensionless specific Helmholtz free energy
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION phi(delta, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: delta, tau
        REAL(KIND=REAL_HIGH) :: phi
        !
        ! Polynomial expression
        !
        phi = n1 * LOG(delta) + SUM(n * (delta**I) * (tau**J))

    END FUNCTION phi
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the dimensionless specific
    !> Helmholtz free energy with respect to the parameter delta.
    !
    !> @param[in] delta Dimensionless density parameter
    !> @param[in] tau   Dimensionless temperature parameter
    !
    !> @return The derivative d(phi)/d(delta)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION phidelta(delta, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: delta, tau
        REAL(KIND=REAL_HIGH) :: phidelta
        !
        ! Polynomial expression
        !
        phidelta = n1 / delta + SUM(n * I * (delta**(I - 1)) * (tau**J))

    END FUNCTION phidelta
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the dimensionless specific
    !> Helmholtz free energy with respect to the parameter delta.
    !
    !> @param[in] delta Dimensionless density parameter
    !> @param[in] tau   Dimensionless temperature parameter
    !
    !> @return The derivative d2(phi)/d(delta)2
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION phideltadelta(delta, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: delta, tau
        REAL(KIND=REAL_HIGH) :: phideltadelta
        !
        ! Polynomial expression
        !
        phideltadelta = -n1 / (delta * delta) + SUM(n * I * (I - 1) * (delta**(I - 2)) * (tau**J))

    END FUNCTION phideltadelta
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the dimensionless specific
    !> Helmholtz free energy with respect to the parameter tau.
    !
    !> @param[in] delta Dimensionless density parameter
    !> @param[in] tau   Dimensionless temperature parameter
    !
    !> @return The derivative d(phi)/d(tau)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION phitau(delta, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: delta, tau
        REAL(KIND=REAL_HIGH) :: phitau
        !
        ! Polynomial expression
        !
        phitau = SUM(n * (delta**I) * J * (tau**(J - 1)))

    END FUNCTION phitau
    !-------------------------------------------------------------------------
    !
    !> This function computes the second derivative of the dimensionless specific
    !> Helmholtz free energy with respect to the parameter tau.
    !
    !> @param[in] delta Dimensionless density parameter
    !> @param[in] tau   Dimensionless temperature parameter
    !
    !> @return The derivative d2(phi)/d(tau)2
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION phitautau(delta, tau)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: delta, tau
        REAL(KIND=REAL_HIGH) :: phitautau
        !
        ! Polynomial expression
        !
        phitautau = SUM(n * (delta**I) * J * (J - 1) * (tau**(J - 2)))

    END FUNCTION phitautau
    !-------------------------------------------------------------------------
    !
    !> This function computes the derivative of the dimensionless specific
    !> Helmholtz free energy with respect to the parameter delta and tau.
    !
    !> @param[in] delta Dimensionless density parameter
    !> @param[in] tau   Dimensionless temperature parameter
    !
    !> @return The derivative d2(phi)/d(delta)d(tau)
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION phideltatau(delta, tau)

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: delta, tau
        REAL(KIND=REAL_HIGH) :: phideltatau
        !
        ! Polynomial expression
        !
        phideltatau = SUM(n * I * (delta**(I - 1)) * J * (tau**(J - 1)))

    END FUNCTION phideltatau
    !-------------------------------------------------------------------------
    !
    !> This function computes pressure with respect to density and temperature
    !
    !> @param[in] density     The steam pressure [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return Pressure in [MPa]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION pressure(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: pressure
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute pressure
        !
        pressure = IAPWS97_R * density * temperature * delta * phidelta(delta, tau) * micro

    END FUNCTION pressure
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific internal energy with respect to
    !> density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific internal energy [J/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_internal_energy(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: specific_internal_energy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute the specific internal energy
        !
        specific_internal_energy = IAPWS97_R * temperature * tau * phitau(delta, tau)

    END FUNCTION specific_internal_energy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific entropy with respect to
    !> density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific entropy [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_entropy(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: specific_entropy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute the specific entropy
        !
        specific_entropy = IAPWS97_R * (tau * phitau(delta, tau) - phi(delta, tau))

    END FUNCTION specific_entropy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific enthalpy with respect to
    !> density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific enthalpy [J/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_enthalpy(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: specific_enthalpy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute the specific enthalpy
        !
        specific_enthalpy = IAPWS97_R * temperature * (tau * phitau(delta, tau) + delta * phidelta(delta, tau))

    END FUNCTION specific_enthalpy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific isobaric heat capacity with respect
    !> to density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific isobaric heat capacity [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_isobaric_heat_capacity(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: specific_isobaric_heat_capacity
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau, res1, res2, res3
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute the specific isobaric heat capacity
        !
        res1 = -tau * tau * phitautau(delta, tau)
        res2 = delta * (phidelta(delta, tau) - tau * phideltatau(delta, tau))**2
        res3 = two * phidelta(delta, tau) + delta * phideltadelta(delta, tau)
        specific_isobaric_heat_capacity = IAPWS97_R * (res1 + res2 / res3)

    END FUNCTION specific_isobaric_heat_capacity
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific isochoric heat capacity with respect
    !> to density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific isochoric heat capacity [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_isochoric_heat_capacity(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: specific_isochoric_heat_capacity
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute the specific isochoric heat capacity
        !
        specific_isochoric_heat_capacity = -IAPWS97_R * tau * tau * phitautau(delta, tau)

    END FUNCTION specific_isochoric_heat_capacity
    !-------------------------------------------------------------------------
    !
    !> This function computes the ratio of specific heats with respect
    !> to density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The ratio of specific heats
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION ratio_of_specific_heats(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: ratio_of_specific_heats
        !
        ! Compute the ratio of specific heats
        !
        ratio_of_specific_heats = specific_isobaric_heat_capacity(density, temperature) / &
                                  specific_isochoric_heat_capacity(density, temperature)

    END FUNCTION ratio_of_specific_heats
    !-------------------------------------------------------------------------
    !
    !> This function computes the speed of sound with respect
    !> to density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The speed of sound [m/sec]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION speed_of_sound(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: speed_of_sound
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau, res1, res2, res3
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute the speed of sound
        !
        res1 = two * delta * phidelta(delta, tau) + delta * delta * phideltadelta(delta, tau)
        res2 = (delta * phidelta(delta, tau) - delta * tau * phideltatau(delta, tau))**2
        res3 = tau * tau * phitautau(delta, tau)
        speed_of_sound = sqrt(IAPWS97_R * temperature * (res1 - res2 / res3))

    END FUNCTION speed_of_sound
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific Helmoltz free energy with respect
    !> to density and temperature.
    !
    !> @param[in] density     The steam density [Kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The specific Helmoltz free energy [J/kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_helmoltz_free_energy(density, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: specific_helmoltz_free_energy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: delta, tau
        !
        ! Compute the dimensionless parameters delta and tau
        !
        delta = density / REGION_3_RHOSTAR
        tau = REGION_3_TSTAR / temperature
        !
        ! Compute the specific Helmoltz free energy
        !
        specific_helmoltz_free_energy = IAPWS97_R * temperature * phi(delta, tau)

    END FUNCTION specific_helmoltz_free_energy
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and enthalpy (backward equation) for the sub-region 3a
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ph_region3a(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:31) :: II = &
                                                              (/-12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, &
                                                                -8, -8, -8, -8, -5, -3, -2, -2, -2, -1, -1, 0, 0, 1, &
                                                                3, 3, 4, 4, 10, 12/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:31) :: JJ = &
                                                              (/0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, &
                                                                0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1, 0, 3, 4, 5/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:31) :: nn = &
                                                            (/-0.133645667811215D-06, 0.455912656802978D-05, &
                                                              -0.146294640700979D-04, 0.639341312970080D-02, &
                                                              0.372783927268847D+03, -0.718654377460447D+04, &
                                                              0.573494752103400D+06, -0.267569329111439D+07, &
                                                              -0.334066283302614D-04, -0.245479214069597D-01, &
                                                              0.478087847764996D+02, 0.764664131818904D-05, &
                                                              0.128350627676972D-02, 0.171219081377331D-01, &
                                                              -0.851007304583213D+01, -0.136513461629781D-01, &
                                                              -0.384460997596657D-05, 0.337423807911655D-02, &
                                                              -0.551624873066791D+00, 0.729202277107470D+00, &
                                                              -0.992522757376041D-02, -0.119308831407288D+00, &
                                                              0.793929190615421D+00, 0.454270731799386D+00, &
                                                              0.209998591259910D+00, -0.642109823904738D-02, &
                                                              -0.235155868604540D-01, 0.252233108341612D-02, &
                                                              -0.764885133368119D-02, 0.136176427574291D-01, &
                                                              -0.133027883575669D-01/)
        !
        !> Star pressure for the backward temperature(pressure, enthalpy)
        !> equation in sub-region 3a in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PH_PSTAR = 100.0D+00
        !
        !> Star enthalpy for the backward temperature(pressure, enthalpy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PH_HSTAR = 2300.0D+03
        !
        !> Star temperature for the backward temperature(pressure, enthalpy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PH_TSTAR = 760.0D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.240D+00, b = 0.615D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: temperature_ph_region3a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3A_PH_PSTAR + a
        eta = enthalpy / REGION_3A_PH_HSTAR - b
        temperature_ph_region3a = REGION_3A_PH_TSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION temperature_ph_region3a
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and enthalpy (backward equation) for the sub-region 3b
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ph_region3b(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:33) :: II = &
                                                              (/-12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, &
                                                                -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, -1, -1, &
                                                                0, 0, 1, 3, 5, 6, 8/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:33) :: JJ = &
                                                              (/0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, &
                                                                5, 0, 4, 2, 4, 6, 10, 14, 16, 0, 2, 1, 1, 1, 1, 1/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:33) :: nn = &
                                                            (/0.323254573644920D-04, -0.127575556587181D-03, &
                                                              -0.475851877356068D-03, 0.156183014181602D-02, &
                                                              0.105724860113781D+00, -0.858514221132534D+02, &
                                                              0.724140095480911D+03, 0.296475810273257D-02, &
                                                              -0.592721983365988D-02, -0.126305422818666D-01, &
                                                              -0.115716196364853D+00, 0.849000969739595D+02, &
                                                              -0.108602260086615D-01, 0.154304475328851D-01, &
                                                              0.750455441524466D-01, 0.252520973612982D-01, &
                                                              -0.602507901232996D-01, -0.307622221350501D+01, &
                                                              -0.574011959864879D-01, 0.503471360939849D+01, &
                                                              -0.925081888584834D+00, 0.391733882917546D+01, &
                                                              -0.773146007130190D+02, 0.949308762098587D+04, &
                                                              -0.141043719679409D+07, 0.849166230819026D+07, &
                                                              0.861095729446704D+00, 0.323346442811720D+00, &
                                                              0.873281936020439D+00, -0.436653048526683D+00, &
                                                              0.286596714529479D+00, -0.131778331276228D+00, &
                                                              0.676682064330275D-02/)
        !
        !> Star pressure for the backward temperature(pressure, enthalpy)
        !> equation in sub-region 3b in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PH_PSTAR = 100.0D+00
        !
        !> Star enthalpy for the backward temperature(pressure, enthalpy)
        !> equation in sub-region 3b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PH_HSTAR = 2800.0D+03
        !
        !> Star temperature for the backward temperature(pressure, enthalpy)
        !> equation in sub-region 3b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PH_TSTAR = 860.0D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.298D+00, b = 0.720D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: temperature_ph_region3b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3B_PH_PSTAR + a
        eta = enthalpy / REGION_3B_PH_HSTAR - b
        temperature_ph_region3b = REGION_3B_PH_TSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION temperature_ph_region3b
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific volume with respect to pressure
    !> and enthalpy (backward equation) for the sub-region 3a
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The specific volume [m3/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_volume_ph_region3a(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:32) :: II = &
                                                              (/-12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, &
                                                                -4, -4, -3, -2, -2, -1, -1, -1, -1, 0, 0, 1, 1, 1, 2, &
                                                                2, 3, 4, 5, 8/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:32) :: JJ = &
                                                              (/6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, &
                                                                16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0, 2, 0, 2, 2, 2/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:32) :: nn = &
                                                            (/0.529944062966028D-02, -0.170099690234461D+00, &
                                                              0.111323814312927D+02, -0.217898123145125D+04, &
                                                              -0.506061827980875D-03, 0.556495239685324D+00, &
                                                              -0.943672726094016D+01, -0.297856807561527D+00, &
                                                              0.939353943717186D+02, 0.192944939465981D-01, &
                                                              0.421740664704763D+00, -0.368914126282330D+07, &
                                                              -0.737566847600639D-02, -0.354753242424366D+00, &
                                                              -0.199768169338727D+01, 0.115456297059049D+01, &
                                                              0.568366875815960D+04, 0.808169540124668D-02, &
                                                              0.172416341519307D+00, 0.104270175292927D+01, &
                                                              -0.297691372792847D+00, 0.560394465163593D+00, &
                                                              0.275234661176914D+00, -0.148347894866012D+00, &
                                                              -0.651142513478515D-01, -0.292468715386302D+01, &
                                                              0.664876096952665D-01, 0.352335014263844D+01, &
                                                              -0.146340792313332D-01, -0.224503486668184D+01, &
                                                              0.110533464706142D+01, -0.408757344495612D-01/)
        !
        !> Star pressure for the backward specific volume(pressure, enthalpy)
        !> equation in sub-region 3a in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PH_PSTAR = 100.0D+00
        !
        !> Star enthalpy for the backward specific volume(pressure, enthalpy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PH_HSTAR = 2100.0D+03
        !
        !> Star specific volume for the backward specific volume(pressure, enthalpy)
        !> equation in sub-region 3a in [m3/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PH_VSTAR = 0.0028D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.128D+00, b = 0.727D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: specific_volume_ph_region3a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3A_PH_PSTAR + a
        eta = enthalpy / REGION_3A_PH_HSTAR - b
        specific_volume_ph_region3a = REGION_3A_PH_VSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION specific_volume_ph_region3a
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific volume with respect to pressure
    !> and enthalpy (backward equation) for the sub-region 3b
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] enthalpy The steam enthalpy [J/Kg]
    !
    !> @return The specific volume [m3/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_volume_ph_region3b(pressure, enthalpy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:30) :: II = &
                                                              (/-12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, &
                                                                -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, -1, 0, 1, 1, 2, 2/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:30) :: JJ = &
                                                              (/0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, &
                                                                2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:30) :: nn = &
                                                            (/-0.225196934336318D-08, 0.140674363313486D-07, &
                                                              0.233784085280560D-05, -0.331833715229001D-04, &
                                                              0.107956778514318D-02, -0.271382067378863D+00, &
                                                              0.107202262490333D+01, -0.853821329075382D+00, &
                                                              -0.215214194340526D-04, 0.769656088222730D-03, &
                                                              -0.431136580433864D-02, 0.453342167309331D+00, &
                                                              -0.507749535873652D+00, -0.100475154528389D+03, &
                                                              -0.219201924648793D+00, -0.321087965668917D+01, &
                                                              0.607567815637771D+03, 0.557686450685932D-03, &
                                                              0.187499040029550D+00, 0.905368030448107D-02, &
                                                              0.285417173048685D+00, 0.329924030996098D-01, &
                                                              0.239897419685483D+00, 0.482754995951394D+01, &
                                                              -0.118035753702231D+02, 0.169490044091791D+00, &
                                                              -0.179967222507787D-01, 0.371810116332674D-01, &
                                                              -0.536288335065096D-01, 0.160697101092520D+01/)
        !
        !> Star pressure for the backward specific volume(pressure, enthalpy)
        !> equation in sub-region 3a in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PH_PSTAR = 100.0D+00
        !
        !> Star enthalpy for the backward specific volume(pressure, enthalpy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PH_HSTAR = 2800.0D+03
        !
        !> Star specific volume for the backward specific volume(pressure, enthalpy)
        !> equation in sub-region 3a in [m3/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PH_VSTAR = 0.0088D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.0661D+00, b = 0.720D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, enthalpy
        REAL(KIND=REAL_HIGH) :: specific_volume_ph_region3b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3B_PH_PSTAR + a
        eta = enthalpy / REGION_3B_PH_HSTAR - b
        specific_volume_ph_region3b = REGION_3B_PH_VSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION specific_volume_ph_region3b
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and entropy (backward equation) for the sub-region 3a
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ps_region3a(pressure, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:33) :: II = &
                                                              (/-12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, &
                                                                -6, -5, -5, -5, -4, -4, -4, -2, -2, -1, -1, 0, 0, 0, 1, 2, &
                                                                2, 3, 8, 8, 10/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:33) :: JJ = &
                                                              (/28, 32, 4, 10, 12, 14, 5, 7, 8, 28, 2, 6, 32, 0, 14, &
                                                                32, 6, 10, 36, 1, 4, 1, 6, 0, 1, 4, 0, 0, 3, 2, 0, 1, 2/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:33) :: nn = &
                                                            (/0.150042008263875D+10, -0.159397258480424D+12, &
                                                              0.502181140217975D-03, -0.672057767855466D+02, &
                                                              0.145058545404456D+04, -0.823889534888890D+04, &
                                                              -0.154852214233853D+00, 0.112305046746695D+02, &
                                                              -0.297000213482822D+02, 0.438565132635495D-02, &
                                                              0.137837838635464D-02, -0.297478527157462D+01, &
                                                              0.971777947349413D+13, -0.571527767052398D-04, &
                                                              0.288307949778420D+05, -0.744428289262703D+14, &
                                                              0.128017324848921D+02, -0.368275545889071D+03, &
                                                              0.664768904779177D+16, 0.449359251958880D-01, &
                                                              -0.422897836099655D+01, -0.240614376434179D+00, &
                                                              -0.474341365254924D+01, 0.724093999126110D+00, &
                                                              0.923874349695897D+00, 0.399043655281015D+01, &
                                                              0.384066651868009D-01, -0.359344365571848D-02, &
                                                              -0.735196448821653D+00, 0.188367048396131D+00, &
                                                              0.141064266818704D-03, -0.257418501496337D-02, &
                                                              0.123220024851555D-02/)
        !
        !> Star pressure for the backward temperature(pressure, entropy)
        !> equation in sub-region 3a in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PS_PSTAR = 100.0D+00
        !
        !> Star entropy for the backward temperature(pressure, entropy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PS_SSTAR = 4.4D+03
        !
        !> Star temperature for the backward temperature(pressure, entropy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PS_TSTAR = 760.0D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.240D+00, b = 0.703D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, entropy
        REAL(KIND=REAL_HIGH) :: temperature_ps_region3a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3A_PS_PSTAR + a
        eta = entropy / REGION_3A_PS_SSTAR - b
        temperature_ps_region3a = REGION_3A_PS_TSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION temperature_ps_region3a
    !-------------------------------------------------------------------------
    !
    !> This function computes the temperature with respect to pressure
    !> and entropy (backward equation) for the sub-region 3b
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION temperature_ps_region3b(pressure, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:28) :: II = &
                                                              (/-12, -12, -12, -12, -8, -8, -8, -6, -6, -6, -5, -5, -5, &
                                                                -5, -5, -4, -3, -3, -2, 0, 2, 3, 4, 5, 6, 8, 12, 14/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:28) :: JJ = &
                                                              (/1, 3, 4, 7, 0, 1, 3, 0, 2, 4, 0, 1, 2, 4, 6, 12, 1, 6, &
                                                                2, 0, 1, 1, 0, 24, 0, 3, 1, 2/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:28) :: nn = &
                                                            (/0.527111701601660D+00, -0.401317830052742D+02, &
                                                              0.153020073134484D+03, -0.224799398218827D+04, &
                                                              -0.193993484669048D+00, -0.140467557893768D+01, &
                                                              0.426799878114024D+02, 0.752810643416743D+00, &
                                                              0.226657238616417D+02, -0.622873556909932D+03, &
                                                              -0.660823667935396D+00, 0.841267087271658D+00, &
                                                              -0.253717501764397D+02, 0.485708963532948D+03, &
                                                              0.880531517490555D+03, 0.265015592794626D+07, &
                                                              -0.359287150025783D+00, -0.656991567673753D+03, &
                                                              0.241768149185367D+01, 0.856873461222588D+00, &
                                                              0.655143675313458D+00, -0.213535213206406D+00, &
                                                              0.562974957606348D-02, -0.316955725450471D+15, &
                                                              -0.699997000152457D-03, 0.119845803210767D-01, &
                                                              0.193848122022095D-04, -0.215095749182309D-04/)
        !
        !> Star pressure for the backward temperature(pressure, entropy)
        !> equation in sub-region 3b in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PS_PSTAR = 100.0D+00
        !
        !> Star entropy for the backward temperature(pressure, entropy)
        !> equation in sub-region 3b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PS_SSTAR = 5.3D+03
        !
        !> Star temperature for the backward temperature(pressure, entropy)
        !> equation in sub-region 3b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PS_TSTAR = 860.0D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.760D+00, b = 0.818D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, entropy
        REAL(KIND=REAL_HIGH) :: temperature_ps_region3b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3B_PS_PSTAR + a
        eta = entropy / REGION_3B_PS_SSTAR - b
        temperature_ps_region3b = REGION_3B_PS_TSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION temperature_ps_region3b
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific volume with respect to pressure
    !> and entropy (backward equation) for the sub-region 3a
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The specific volume [m3/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_volume_ps_region3a(pressure, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:28) :: II = &
                                                              (/-12, -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, &
                                                                -5, -4, -3, -3, -2, -2, -1, -1, 0, 0, 0, 1, 2, 4, 5, 6/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:28) :: JJ = &
                                                              (/10, 12, 14, 4, 8, 10, 20, 5, 6, 14, 16, 28, 1, 5, 2, 4, &
                                                                3, 8, 1, 2, 0, 1, 3, 0, 0, 2, 2, 0/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:28) :: nn = &
                                                            (/0.795544074093975D+02, -0.238261242984590D+04, &
                                                              0.176813100617787D+05, -0.110524727080379D-02, &
                                                              -0.153213833655326D+02, 0.297544599376982D+03, &
                                                              -0.350315206871242D+08, 0.277513761062119D+00, &
                                                              -0.523964271036888D+00, -0.148011182995403D+06, &
                                                              0.160014899374266D+07, 0.170802322663427D+13, &
                                                              0.246866996006494D-03, 0.165326084797980D+01, &
                                                              -0.118008384666987D+00, 0.253798642355900D+01, &
                                                              0.965127704669424D+00, -0.282172420532826D+02, &
                                                              0.203224612353823D+00, 0.110648186063513D+01, &
                                                              0.526127948451280D+00, 0.277000018736321D+00, &
                                                              0.108153340501132D+01, -0.744127885357893D-01, &
                                                              0.164094443541384D-01, -0.680468275301065D-01, &
                                                              0.257988576101640D-01, -0.145749861944416D-03/)
        !
        !> Star pressure for the backward specific volume(pressure, entropy)
        !> equation in sub-region 3a in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PS_PSTAR = 100.0D+00
        !
        !> Star entropy for the backward specific volume(pressure, entropy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PS_SSTAR = 4.4D+03
        !
        !> Star specific volume for the backward specific volume(pressure, entropy)
        !> equation in sub-region 3a in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_PS_VSTAR = 0.0028D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.187D+00, b = 0.755D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, entropy
        REAL(KIND=REAL_HIGH) :: specific_volume_ps_region3a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3A_PS_PSTAR + a
        eta = entropy / REGION_3A_PS_SSTAR - b
        specific_volume_ps_region3a = REGION_3A_PS_VSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION specific_volume_ps_region3a
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific volume with respect to pressure
    !> and entropy (backward equation) for the sub-region 3b
    !
    !> @param[in] pressure The steam pressure [MPa]
    !> @param[in] entropy  The steam entropy [J/Kg.K]
    !
    !> @return The specific volume [m3/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_volume_ps_region3b(pressure, entropy)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:31) :: II = &
                                                              (/-12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, &
                                                                -5, -5, -5, -4, -4, -4, -4, -3, -2, -2, -2, -2, -2, -2, &
                                                                0, 0, 0, 1, 1, 2/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:31) :: JJ = &
                                                              (/0, 1, 2, 3, 5, 6, 0, 1, 2, 4, 0, 1, 2, 3, 0, 1, 2, 3, 1, &
                                                                0, 1, 2, 3, 4, 12, 0, 1, 2, 0, 2, 2/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:31) :: nn = &
                                                            (/0.591599780322238D-04, -0.185465997137856D-02, &
                                                              0.104190510480013D-01, 0.598647302038590D-02, &
                                                              -0.771391189901699D+00, 0.172549765557036D+01, &
                                                              -0.467076079846526D-03, 0.134533823384439D-01, &
                                                              -0.808094336805495D-01, 0.508139374365767D+00, &
                                                              0.128584643361683D-02, -0.163899353915435D+01, &
                                                              0.586938199318063D+01, -0.292466667918613D+01, &
                                                              -0.614076301499537D-02, 0.576199014049172D+01, &
                                                              -0.121613320606788D+02, 0.167637540957944D+01, &
                                                              -0.744135838773463D+01, 0.378168091437659D-01, &
                                                              0.401432203027688D+01, 0.160279837479185D+02, &
                                                              0.317848779347728D+01, -0.358362310304853D+01, &
                                                              -0.115995260446827D+07, 0.199256573577909D+00, &
                                                              -0.122270624794624D+00, -0.191449143716586D+02, &
                                                              -0.150448002905284D-01, 0.146407900162154D+02, &
                                                              -0.327477787188230D+01/)
        !
        !> Star pressure for the backward specific volume(pressure, entropy)
        !> equation in sub-region 3b in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PS_PSTAR = 100.0D+00
        !
        !> Star entropy for the backward specific volume(pressure, entropy)
        !> equation in sub-region 3b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PS_SSTAR = 5.3D+03
        !
        !> Star specific volume for the backward specific volume(pressure, entropy)
        !> equation in sub-region 3b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3B_PS_VSTAR = 0.0088D+00
        !
        ! Constants
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: a = 0.298D+00, b = 0.816D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, entropy
        REAL(KIND=REAL_HIGH) :: specific_volume_ps_region3b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, eta
        !
        ! Compute the temperature
        !
        pi = pressure / REGION_3B_PS_PSTAR + a
        eta = entropy / REGION_3B_PS_SSTAR - b
        specific_volume_ps_region3b = REGION_3B_PS_VSTAR * SUM(nn * (pi**II) * (eta**JJ))

    END FUNCTION specific_volume_ps_region3b

END MODULE SteamRegion3
