!-----------------------------------------------------------------------
!
!> @file SteamRegion4.f
!
!> @details
!> This module implements the polynomials which compute the thermodynamic
!> properties of steam in region 4 according to the revised release on the
!> IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
!> Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamRegion4

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Thermodynamic properties computed by polynomials
    !
    PUBLIC saturation_pressure
    PUBLIC saturation_pressure_gradient
    PUBLIC saturation_temperature
    PUBLIC specific_internal_energy
    PUBLIC specific_volume
    PUBLIC specific_enthalpy
    PUBLIC specific_entropy
    PUBLIC specific_isobaric_heat_capacity
    PUBLIC specific_isochoric_heat_capacity
    PUBLIC ratio_of_specific_heats
    !
    !> Constant coefficients "n"
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:10) :: n = &
        (/  0.11670521452767D+04, -0.72421316703206D+06,         &
           -0.17073846940092D+02,  0.12020824702470D+05,         &
           -0.32325550322333D+07,  0.14915108613530D+02,         &
           -0.48232657361591D+04,  0.40511340542057D+06,         &
           -0.23855557567849D+00,  0.65017534844798D+03 /)
    !
    !> Constant coefficients "K"
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:6) :: K = &
        (/  1.99274064D+00,  1.09965342D+00, -0.510839303D+00,  &
           -1.75493479D+00, -45.5170352D+00,  -6.74694450D+05 /)
    !
    !> Constant coefficients "L"
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:6) :: L = &
        (/  -2.03150240D+00, -2.68302940D+00, -5.38626492D+00,  &
            -17.2991605D+00, -44.7586581D+00, -63.9201063D+00 /)
    !
    !> Star pressure for region 4 in [MPa]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_4_PSTAR = 1.0D+00
    !
    !> Star temperature for region 4 in [K]
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: REGION_4_TSTAR = 1.0D+00

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function computes the saturation pressure with respect to
    !> temperature
    !
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return Pressure [MPa]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION saturation_pressure(temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature
        REAL(KIND=REAL_HIGH) :: saturation_pressure
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: tau, theta, theta2, a, b, c, expr
        !
        ! Compute the saturation pressure
        !
        tau    = temperature/REGION_4_TSTAR
        theta  = tau + n(9)/(tau-n(10))
        theta2 = theta*theta
        a      = theta2 + n(1)*theta + n(2)
        b      = n(3)*theta2 + n(4)*theta + n(5)
        c      = n(6)*theta2 + n(7)*theta + n(8)
        expr   = two*c/(-b+sqrt(b*b-four*a*c))

        saturation_pressure = REGION_4_PSTAR*(expr**4)

    END FUNCTION saturation_pressure
    !-------------------------------------------------------------------------
    !
    !> This function computes the saturation pressure gradient with respect
    !> to temperature
    !
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return The pressure gradient [MPa/K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION saturation_pressure_gradient(temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature
        REAL(KIND=REAL_HIGH) :: saturation_pressure_gradient
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: beta, theta, xbeta, xtheta, &
            dthetadT, dbetadtheta, dpdbeta
        !
        ! Compute the saturation pressure gradient
        !
        beta        = (saturation_pressure(temperature)/REGION_4_PSTAR)**quarter
        theta       = temperature/REGION_4_TSTAR + n(9) / (temperature/REGION_4_TSTAR - n(10))
        xbeta       = (two*beta + n(3))*theta*theta + &
            (two*beta*n(1) + n(4))*theta + two*n(2)*beta + n(5)
        xtheta      = (two*theta + n(1))*beta*beta +  &
            (two*n(3)*theta + n(4))*beta + two*n(6)*theta + n(7)
        dthetadT    = (one - n(9)/(temperature/REGION_4_TSTAR - n(10)))/REGION_4_TSTAR
        dbetadtheta = -xtheta/xbeta;
        dpdbeta     = four*beta*beta*beta*REGION_4_PSTAR

        saturation_pressure_gradient = dpdbeta * dbetadtheta * dthetadT

    END FUNCTION saturation_pressure_gradient
    !-------------------------------------------------------------------------
    !
    !> This function computes the saturation temperature with respect to
    !> pressure
    !
    !> @param[in] pressure The steam pressure [MPa]
    !
    !> @return Temperature [K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION saturation_temperature(pressure)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure
        REAL(KIND=REAL_HIGH) :: saturation_temperature
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pi, beta, beta2, e, f, g, d
        !
        ! Compute the saturation temperature
        !
        pi    = pressure/REGION_4_PSTAR
        beta  = pi**quarter
        beta2 = beta*beta
        e     = beta2 + n(3)*beta + n(6)
        f     = n(1)*beta2 + n(4)*beta + n(7)
        g     = n(2)*beta2 + n(5)*beta + n(8)
        d     = two*g/(-f-sqrt(f*f-four*e*g))

        saturation_temperature = REGION_4_TSTAR*half*(n(10) + &
                                 d - sqrt((n(10)+d)**2-four*(n(9)+n(10)*d)))

    END FUNCTION saturation_temperature
    !-------------------------------------------------------------------------
    !
    !> This function computes the saturation water density with respect to
    !> temperature
    !
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return Saturation water density [Kg/m3]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION saturation_water_density(temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature
        REAL(KIND=REAL_HIGH) :: saturation_water_density
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: tau, tau_1_3, tau_2_3, tau_5_3, &
            tau_16_3, tau_43_3, tau_110_3, delta
        !
        ! Compute the saturation water density
        !
        tau       = one - temperature/IAPWS97_TCRIT
        tau_1_3   = tau**onethird
        tau_2_3   = tau_1_3*tau_1_3
        tau_5_3   = tau*tau_2_3
        tau_16_3  = tau_5_3*tau_5_3*tau_5_3*tau_1_3
        tau_43_3  = tau_16_3*tau_16_3*tau_5_3*tau_5_3*tau_1_3
        tau_110_3 = tau_43_3*tau_43_3*tau_16_3*tau_5_3*tau
        delta     = one + K(1)*tau_1_3 + K(2)*tau_2_3 + K(3)*tau_5_3 &
                        + K(4)*tau_16_3 + K(5)*tau_43_3 + K(6)*tau_110_3

        saturation_water_density = IAPWS97_RHOCRIT*delta

    END FUNCTION saturation_water_density
    !-------------------------------------------------------------------------
    !
    !> This function computes the saturation steam density with respect to
    !> temperature
    !
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return Saturation steam density [Kg/m3]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION saturation_steam_density(temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature
        REAL(KIND=REAL_HIGH) :: saturation_steam_density
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: tau, tau_1_6, tau_2_6, tau_4_6, tau_8_6, &
            tau_16_6, tau_18_6, tau_37_6, tau_71_6, ln_delta
        !
        ! Compute the saturation steam density
        !
        tau      = one - temperature/IAPWS97_TCRIT
        tau_1_6  = tau**onesixth
        tau_2_6  = tau_1_6*tau_1_6
        tau_4_6  = tau_2_6*tau_2_6
        tau_8_6  = tau_4_6*tau_4_6
        tau_16_6 = tau_8_6*tau_8_6
        tau_18_6 = tau_16_6*tau_2_6
        tau_37_6 = tau_18_6*tau_18_6*tau_1_6
        tau_71_6 = tau_37_6*tau_18_6*tau_16_6
        ln_delta = L(1)*tau_2_6 + L(2)*tau_4_6 + L(3)*tau_8_6 + &
                   L(4)*tau_18_6 + L(5)*tau_37_6 + L(6)*tau_71_6

        saturation_steam_density = IAPWS97_RHOCRIT*EXP(ln_delta)

    END FUNCTION saturation_steam_density
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific internal energy with respect to
    !> temperature and steam quality.
    !
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] quality     The steam quality [%]
    !
    !> @return The specific internal energy [J/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_internal_energy(temperature, quality)

        USE SteamRegion1, ONLY : u1 => specific_internal_energy
        USE SteamRegion2, ONLY : u2 => specific_internal_energy
        USE SteamRegion3, ONLY : u3 => specific_internal_energy

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature, quality
        REAL(KIND=REAL_HIGH) :: specific_internal_energy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: uf, ug, psat, rhof, rhog
        !
        ! Compute the specific internal energy
        !
        IF (temperature < REGION_1_TMAX) THEN
            psat = saturation_pressure(temperature)
            uf   = u1(psat, temperature)
            ug   = u2(psat, temperature)
        ELSE
            rhof = saturation_water_density(temperature)
            rhog = saturation_steam_density(temperature)
            uf   = u3(rhof, temperature)
            ug   = u3(rhog, temperature)
        END IF

        specific_internal_energy = uf + quality*(ug-uf)

    END FUNCTION specific_internal_energy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific volume with respect to
    !> temperature and steam quality.
    !
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] quality     The steam quality [%]
    !
    !> @return The specific volume [m3/kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_volume(temperature, quality)

        USE SteamRegion1, ONLY : v1 => specific_volume
        USE SteamRegion2, ONLY : v2 => specific_volume

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature, quality
        REAL(KIND=REAL_HIGH) :: specific_volume
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: vf, vg, psat
        !
        ! Compute the specific volume
        !
        IF (temperature < REGION_1_TMAX) THEN
            psat = saturation_pressure(temperature)
            vf   = v1(psat, temperature)
            vg   = v2(psat, temperature)
        ELSE
            vf   = one/saturation_water_density(temperature)
            vg   = one/saturation_steam_density(temperature)
        END IF

        specific_volume = vf + quality*(vg-vf)

    END FUNCTION specific_volume
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific enthalpy with respect to
    !> temperature and steam quality.
    !
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] quality     The steam quality [%]
    !
    !> @return The specific enthalpy [J/Kg]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_enthalpy(temperature, quality)

        USE SteamRegion1, ONLY : h1 => specific_enthalpy
        USE SteamRegion2, ONLY : h2 => specific_enthalpy
        USE SteamRegion3, ONLY : h3 => specific_enthalpy

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature, quality
        REAL(KIND=REAL_HIGH) :: specific_enthalpy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: hf, hg, psat, rhof, rhog
        !
        ! Compute the specific enthalpy
        !
        IF (temperature < REGION_1_TMAX) THEN
            psat = saturation_pressure(temperature)
            hf   = h1(psat, temperature)
            hg   = h2(psat, temperature)
        ELSE
            rhof = saturation_water_density(temperature)
            rhog = saturation_steam_density(temperature)
            hf   = h3(rhof, temperature)
            hg   = h3(rhog, temperature)
        END IF

        specific_enthalpy = hf + quality*(hg-hf)

    END FUNCTION specific_enthalpy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific entropy with respect to
    !> temperature and steam quality.
    !
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] quality     The steam quality [%]
    !
    !> @return The specific entropy [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_entropy(temperature, quality)

        USE SteamRegion1, ONLY : s1 => specific_entropy
        USE SteamRegion2, ONLY : s2 => specific_entropy
        USE SteamRegion3, ONLY : s3 => specific_entropy

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature, quality
        REAL(KIND=REAL_HIGH) :: specific_entropy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: sf, sg, psat, rhof, rhog
        !
        ! Compute the specific entropy
        !
        IF (temperature < REGION_1_TMAX) THEN
            psat = saturation_pressure(temperature)
            sf   = s1(psat, temperature)
            sg   = s2(psat, temperature)
        ELSE
            rhof = saturation_water_density(temperature)
            rhog = saturation_steam_density(temperature)
            sf   = s3(rhof, temperature)
            sg   = s3(rhog, temperature)
        END IF

        specific_entropy = sf + quality*(sg-sf)

    END FUNCTION specific_entropy
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific isobaric heat capacity with respect
    !> to temperature and steam quality.
    !
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] quality     The steam quality [%]
    !
    !> @return The specific isobaric heat capacity [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_isobaric_heat_capacity(temperature, quality)

        USE SteamRegion1, ONLY : cp1 => specific_isobaric_heat_capacity
        USE SteamRegion2, ONLY : cp2 => specific_isobaric_heat_capacity
        USE SteamRegion3, ONLY : cp3 => specific_isobaric_heat_capacity

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature, quality
        REAL(KIND=REAL_HIGH) :: specific_isobaric_heat_capacity
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: cpf, cpg, psat, rhof, rhog
        !
        ! Compute the specific isobaric heat capacity
        !
        IF (temperature < REGION_1_TMAX) THEN
            psat = saturation_pressure(temperature)
            cpf  = cp1(psat, temperature)
            cpg  = cp2(psat, temperature)
        ELSE
            rhof = saturation_water_density(temperature)
            rhog = saturation_steam_density(temperature)
            cpf  = cp3(rhof, temperature)
            cpg  = cp3(rhog, temperature)
        END IF

        specific_isobaric_heat_capacity = cpf + quality*(cpg-cpf)

    END FUNCTION specific_isobaric_heat_capacity
    !-------------------------------------------------------------------------
    !
    !> This function computes the specific isochoric heat capacity with respect
    !> to temperature and steam quality.
    !
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] quality     The steam quality [%]
    !
    !> @return The specific isochoric heat capacity [J/Kg.K]
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION specific_isochoric_heat_capacity(temperature, quality)

        USE SteamRegion1, ONLY : cv1 => specific_isochoric_heat_capacity
        USE SteamRegion2, ONLY : cv2 => specific_isochoric_heat_capacity
        USE SteamRegion3, ONLY : cv3 => specific_isochoric_heat_capacity

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature, quality
        REAL(KIND=REAL_HIGH) :: specific_isochoric_heat_capacity
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: cvf, cvg, psat, rhof, rhog
        !
        ! Compute the specific isochoric heat capacity
        !
        IF (temperature < REGION_1_TMAX) THEN
            psat = saturation_pressure(temperature)
            cvf  = cv1(psat, temperature)
            cvg  = cv2(psat, temperature)
        ELSE
            rhof = saturation_water_density(temperature)
            rhog = saturation_steam_density(temperature)
            cvf  = cv3(rhof, temperature)
            cvg  = cv3(rhog, temperature)
        END IF

        specific_isochoric_heat_capacity = cvf + quality*(cvg-cvf)

    END FUNCTION specific_isochoric_heat_capacity
    !-------------------------------------------------------------------------
    !
    !> This function computes the ratio of specific heats with respect
    !> to temperature and steam quality.
    !
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] quality     The steam quality [%]
    !
    !> @return The ratio of specific heats
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION ratio_of_specific_heats(temperature, quality)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature, quality
        REAL(KIND=REAL_HIGH) :: ratio_of_specific_heats
        !
        ! Compute the ratio of specific heats
        !
        ratio_of_specific_heats = specific_isobaric_heat_capacity(temperature, quality) / &
            specific_isochoric_heat_capacity(temperature, quality)

    END FUNCTION ratio_of_specific_heats

END MODULE SteamRegion4
