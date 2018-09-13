!-----------------------------------------------------------------------
!
!> @file SteamRegions.f
!
!> @details
!> This module implements function that are capable of determining if
!> a steam state belongs to a particular region according to the revised
!> release on the IAPWS Industrial Formulation 1997 for the Thermodynamic
!> Properties of Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamRegions

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Region detection functions
    !
    PUBLIC in_region1
    PUBLIC in_region2
    PUBLIC in_region2a
    PUBLIC in_region2b
    PUBLIC in_region2c
    PUBLIC in_region3
    PUBLIC in_region3a
    PUBLIC in_region3b
    PUBLIC in_region4
    PUBLIC in_region5

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 1 or not.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 1
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION in_region1(pressure, temperature)

        USE SteamRegion4, ONLY: psat => saturation_pressure

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        LOGICAL :: in_region1
        !
        ! Variables
        !
        LOGICAL :: pflag, tflag
        !
        ! Determine whether the steam state belongs in region 1
        !
        tflag = (temperature >= IAPWS97_TMIN) .AND. &
                (temperature <= REGION_1_TMAX)

        pflag = (pressure >= psat(temperature)) .AND. &
                (pressure <= IAPWS97_PMAX)

        in_region1 = .FALSE.

        IF (pflag .AND. tflag) THEN
            in_region1 = .TRUE.
        END IF

    END FUNCTION in_region1
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 2 or not.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 2
    !
    !-------------------------------------------------------------------------
    FUNCTION in_region2(pressure, temperature)

        USE SteamRegion4, ONLY: psat => saturation_pressure
        USE SteamBoundaries, ONLY: pbnd => boundary23_pressure

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        LOGICAL :: in_region2
        !
        ! Determine whether the steam state belongs in region 2
        !
        in_region2 = .FALSE.

        IF ((temperature >= IAPWS97_TMIN) .AND. &
            (temperature <= REGION_1_TMAX)) THEN

            IF ((pressure > 0.0) .AND. &
                (pressure <= psat(temperature))) THEN
                in_region2 = .TRUE.
            END IF

        ELSE IF ((temperature > REGION_1_TMAX) .AND. &
                 (temperature <= REGION_2_4_T)) THEN

            IF ((pressure > 0.0) .AND. &
                (pressure <= pbnd(temperature))) THEN
                in_region2 = .TRUE.
            END IF

        ELSE IF ((temperature > REGION_2_4_T) .AND. &
                 (temperature <= REGION_2_TMAX)) THEN

            IF ((pressure > 0.0) .AND. &
                (pressure <= IAPWS97_PMAX)) THEN
                in_region2 = .TRUE.
            END IF

        END IF

    END FUNCTION in_region2
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 2, sub-region 2a.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 2, sub-region 2a.
    !
    !-------------------------------------------------------------------------
    FUNCTION in_region2a(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        LOGICAL :: in_region2a
        !
        ! Determine whether the steam state belongs in in region 2, sub-region 2a
        !
        in_region2a = .FALSE.

        IF (pressure <= REGION_2A_2B_P) THEN
            in_region2a = in_region2(pressure, temperature)
        END IF

    END FUNCTION in_region2a
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 2, sub-region 2b.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 2, sub-region 2b.
    !
    !-------------------------------------------------------------------------
    FUNCTION in_region2b(pressure, temperature)

        USE SteamBoundaries, ONLY: p2bc => boundary2bc_pressure
        USE SteamRegion2, ONLY: h2 => specific_enthalpy

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        LOGICAL :: in_region2b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: h
        !
        ! Determine whether the steam state belongs in region 2, sub-region 2b
        !
        in_region2b = .FALSE.

        h = h2(pressure, temperature)

        IF (pressure > REGION_2A_2B_P) THEN
            IF (pressure <= p2bc(h)) THEN
                in_region2b = in_region2(pressure, temperature)
            END IF
        END IF

    END FUNCTION in_region2b
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 2, sub-region 2c.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 2, sub-region 2c.
    !
    !-------------------------------------------------------------------------
    FUNCTION in_region2c(pressure, temperature)

        USE SteamBoundaries, ONLY: p2bc => boundary2bc_pressure
        USE SteamRegion2, ONLY: h2 => specific_enthalpy

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        LOGICAL :: in_region2c
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: h
        !
        ! Determine whether the steam state belongs in region 2, sub-region 2c
        !
        in_region2c = .FALSE.

        h = h2(pressure, temperature)

        IF (pressure > p2bc(h)) THEN
            in_region2c = in_region2(pressure, temperature)
        END IF

    END FUNCTION in_region2c
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 3 or not.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 3
    !
    !-------------------------------------------------------------------------
    FUNCTION in_region3(pressure, temperature)

        USE SteamBoundaries

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        LOGICAL :: in_region3
        !
        ! Variables
        !
        LOGICAL :: pflag, tflag
        !
        ! Determine whether the steam state belongs in region 3
        !
        tflag = (temperature >= REGION_1_TMAX) .AND. &
                (temperature <= boundary23_temperature(pressure))

        pflag = (pressure >= boundary23_pressure(temperature)) .AND. &
                (pressure <= IAPWS97_PMAX)

        in_region3 = .FALSE.

        IF (pflag .AND. tflag) THEN
            in_region3 = .TRUE.
        END IF

    END FUNCTION in_region3
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 3a or not.
    !
    !> @param[in] density     The steam density [kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 3a
    !
    !-------------------------------------------------------------------------
    FUNCTION in_region3a(density, temperature)

        USE SteamRegion3
        USE SteamBoundaries

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        LOGICAL :: in_region3a
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: h3ab, enthalpy, press
        !
        ! Determine whether the steam state belongs in region 3a
        !
        enthalpy = specific_enthalpy(density, temperature)
        press = pressure(density, temperature)
        h3ab = boundary3ab_enthalpy(press)

        IF (enthalpy < h3ab) THEN
            in_region3a = in_region3(press, temperature)
        ELSE
            in_region3a = .FALSE.
        END IF

    END FUNCTION in_region3a
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 3b or not.
    !
    !> @param[in] density     The steam density [kg/m3]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 3b
    !
    !-------------------------------------------------------------------------
    FUNCTION in_region3b(density, temperature)

        USE SteamRegion3
        USE SteamBoundaries
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        LOGICAL :: in_region3b
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: h3ab, enthalpy, press
        !
        ! Determine whether the steam state belongs in region 3a
        !
        enthalpy = specific_enthalpy(density, temperature)
        press = pressure(density, temperature)
        h3ab = boundary3ab_enthalpy(press)

        IF (enthalpy >= h3ab) THEN
            in_region3b = in_region3(press, temperature)
        ELSE
            in_region3b = .FALSE.
        END IF

    END FUNCTION in_region3b
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 4 or not. Since this region is a line one cannot determine
    !> exactly whether a steam state lies exactly on this line or not so an
    !> error percentage is required for this function to work.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !> @param[in] error       The acceptable error [%]
    !
    !> @return True if the steam state resides in region 4
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION in_region4(pressure, temperature, error)

        USE SteamRegion4, ONLY: tsat => saturation_temperature, &
            psat => saturation_pressure

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature, error
        LOGICAL :: in_region4
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: pressure_estimate, temperature_estimate
        REAL(KIND=REAL_HIGH) :: t_error, p_error
        !
        ! Determine whether the steam state belongs in region 4
        !
        in_region4 = .FALSE.

        IF ((temperature >= IAPWS97_TMIN) .AND. &
            (temperature <= REGION_4_TMAX)) THEN

            pressure_estimate = psat(temperature)
            temperature_estimate = tsat(pressure)

            p_error = ABS(pressure - pressure_estimate) / pressure
            t_error = ABS(temperature - temperature_estimate) / temperature

            IF ((p_error <= error) .AND. (t_error <= error)) THEN
                in_region4 = .TRUE.
            END IF

        END IF

    END FUNCTION in_region4
    !-------------------------------------------------------------------------
    !
    !> This function determines whether the specified steam state resides in
    !> region 5 or not.
    !
    !> @param[in] pressure    The steam pressure [MPa]
    !> @param[in] temperature The steam temperature [K]
    !
    !> @return True if the steam state resides in region 5
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION in_region5(pressure, temperature)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure, temperature
        LOGICAL :: in_region5
        !
        ! Variables
        !
        LOGICAL :: pflag, tflag
        !
        ! Determine whether the steam state belongs in region 5
        !
        tflag = (temperature >= REGION_2_TMAX) .AND. &
                (temperature <= IAPWS97_TMAX)

        pflag = (pressure > 0.0) .AND. &
                (pressure <= REGION_5_PMAX)

        in_region5 = .FALSE.

        IF (pflag .AND. tflag) THEN
            in_region5 = .TRUE.
        END IF

    END FUNCTION in_region5

END MODULE SteamRegions
