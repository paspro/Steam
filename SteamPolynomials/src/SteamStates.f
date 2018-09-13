!-----------------------------------------------------------------------
!
!> @file SteamStates.f
!
!> @details
!> This module implements a steam state class hierarchy
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamStates

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Public functions / subroutines
    !
    PUBLIC set_primary_properties
    PUBLIC get_specific_volume
    PUBLIC get_density
    PUBLIC get_pressure
    PUBLIC get_temperature
    PUBLIC get_specific_internal_energy
    PUBLIC get_specific_enthalpy
    PUBLIC get_specific_entropy
    PUBLIC get_specific_isobaric_heat_capacity
    PUBLIC get_specific_isochoric_heat_capacity
    PUBLIC get_ratio_of_specific_heats
    PUBLIC get_viscosity
    PUBLIC get_thermal_conductivity
    PUBLIC print_state

    !-----------------------------------------------------------------------
    !
    !> @details
    !> This type implements a steam state.
    !
    !-------------------------------------------------------------------------
    TYPE, PUBLIC :: SteamState

        PRIVATE
        !
        !> The region the steam state is located
        !
        INTEGER(KIND=INT_HIGH) :: region
        !
        !> The steam pressure in [MPa]
        !
        REAL(KIND=REAL_HIGH) :: pressure
        !
        !> The steam temperature in [K]
        !
        REAL(KIND=REAL_HIGH) :: temperature
        !
        !> The steam density in [kg/m3]
        !
        REAL(KIND=REAL_HIGH) :: density
        !
        !> The steam quality in [%]
        !
        REAL(KIND=REAL_HIGH) :: quality
        !
        !> The specific volume of the steam state in [m3/kg]
        !
        REAL(KIND=REAL_HIGH) :: v
        !
        !> The specific internal energy of the steam state in [J/kg]
        !
        REAL(KIND=REAL_HIGH) :: u
        !
        !> The specific enthalpy of the steam state in [J/kg]
        !
        REAL(KIND=REAL_HIGH) :: h
        !
        !> The specific entropy of the steam state in [J/Kg.k]
        !
        REAL(KIND=REAL_HIGH) :: s
        !
        !> The specific isobaric heat capacity of the steam state
        !> in [J/Kg.K]
        !
        REAL(KIND=REAL_HIGH) :: cp
        !
        !> The specific isochoric heat capacity of the steam state
        !> in [J/Kg.K]
        !
        REAL(KIND=REAL_HIGH) :: cv
        !
        !> The ratio of specific heats
        !
        REAL(KIND=REAL_HIGH) :: ratio_cp_cv
        !
        !> Viscosity in [Pa.sec]
        !
        REAL(KIND=REAL_HIGH) :: mu
        !
        !> Thermal Conductivity in [W/K.m]
        !
        REAL(KIND=REAL_HIGH) :: tk

    END TYPE SteamState

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> @details
    !> This function is setting the primary properties of the steam state.
    !> These properties depend on the region the steam state belongs to.
    !> For regions 1, 2 and 5 these are pressure and temperature. For
    !> region 3 they are density and temperature and for region 4 they are
    !> temperature and steam quality.
    !
    !> @param[in]  this   The SteamState class object
    !> @param[in]  region The region of this steam state
    !> @param[in]  prop1  The first primary property (as described above)
    !> @param[in]  prop2  The second primary property (as described above)
    !> @param[out] valid  Becomes true if the specified thermodynamic
    !                     properties are valid for the region of this state
    !
    !-------------------------------------------------------------------------
    SUBROUTINE set_primary_properties(this, region, prop1, prop2, valid)

        USE SteamRegions
        USE SteamRegion3, ONLY: p3 => pressure
        USE SteamRegion4, ONLY: p4 => saturation_pressure

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        INTEGER(KIND=INT_HIGH), INTENT(IN) :: region
        REAL(KIND=REAL_HIGH), INTENT(IN) :: prop1, prop2
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: p
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid

        this % region = region

        SELECT CASE (region)

        CASE (1)

            this % pressure = prop1
            this % temperature = prop2

            IF (present(valid)) THEN
                valid = in_region1(prop1, prop2)
            END IF

        CASE (2)

            this % pressure = prop1
            this % temperature = prop2

            IF (present(valid)) THEN
                valid = in_region2(prop1, prop2)
            END IF

        CASE (3)

            this % density = prop1
            this % temperature = prop2

            IF (present(valid)) THEN
                p = p3(prop1, prop2)
                valid = in_region3(p, prop2)
            END IF

        CASE (4)

            this % temperature = prop1
            this % quality = prop2

            IF (present(valid)) THEN
                p = p4(prop1)
                valid = in_region4(p, prop1, 1.0D-04)
            END IF

        CASE (5)

            this % pressure = prop1
            this % temperature = prop2

            IF (present(valid)) THEN
                valid = in_region5(prop1, prop2)
            END IF

        END SELECT
        !
        ! Set the other thermodynamic properties to invalid values
        !
        this % v = -1.0D+00
        this % u = -1.0D+00
        this % h = -1.0D+00
        this % s = -1.0D+00
        this % cp = -1.0D+00
        this % cv = -1.0D+00
        this % ratio_cp_cv = -1.0D+00
        this % mu = -1.0D+00
        this % tk = -1.0D+00

    END SUBROUTINE set_primary_properties
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific volume of the steam state
    !
    !> @return The specific volume [m3/kg]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_specific_volume(this)

        USE SteamRegion1, ONLY: v1 => specific_volume
        USE SteamRegion2, ONLY: v2 => specific_volume
        USE SteamRegion4, ONLY: v4 => specific_volume
        USE SteamRegion5, ONLY: v5 => specific_volume

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_specific_volume
        !
        ! Compute the specific volume if not already computed
        !
        get_specific_volume = zero

        IF (this % v > zero) THEN
            get_specific_volume = this % v
            RETURN
        END IF

        SELECT CASE (this % region)

        CASE (1)
            this % v = v1(this % pressure, this % temperature)

        CASE (2)
            this % v = v2(this % pressure, this % temperature)

        CASE (3)
            this % v = one / this % density

        CASE (4)
            this % v = v4(this % temperature, this % quality)

        CASE (5)
            this % v = v5(this % pressure, this % temperature)

        END SELECT

        get_specific_volume = this % v

    END FUNCTION get_specific_volume
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific volume of the steam state
    !
    !> @return The specific volume [m3/kg]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_pressure(this)

        USE SteamRegion3, ONLY: p3 => pressure
        USE SteamRegion4, ONLY: p4 => saturation_pressure

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_pressure
        !
        ! Compute pressure if not already known
        !
        get_pressure = zero

        SELECT CASE (this % region)

        CASE (1)
            get_pressure = this % pressure

        CASE (2)
            get_pressure = this % pressure

        CASE (3)
            this % pressure = p3(this % density, this % temperature)
            get_pressure = this % pressure

        CASE (4)
            this % pressure = p4(this % temperature)
            get_pressure = this % pressure

        CASE (5)
            get_pressure = this % pressure

        END SELECT

    END FUNCTION get_pressure
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the density of the steam state
    !
    !> @return The density in [Kg/m3]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_density(this)

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_density

        get_density = one / get_specific_volume(this)

    END FUNCTION get_density
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the temperature of the steam state
    !
    !> @return The temperature in [K]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_temperature(this)

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_temperature

        get_temperature = this % temperature

    END FUNCTION get_temperature
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific internal energy of the steam state
    !
    !> @return The specific internal energy [J/kg]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_specific_internal_energy(this)

        USE SteamRegion1, ONLY: u1 => specific_internal_energy
        USE SteamRegion2, ONLY: u2 => specific_internal_energy
        USE SteamRegion3, ONLY: u3 => specific_internal_energy
        USE SteamRegion4, ONLY: u4 => specific_internal_energy
        USE SteamRegion5, ONLY: u5 => specific_internal_energy

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_specific_internal_energy
        !
        ! Compute the specific enthalpy if not already computed
        !
        get_specific_internal_energy = zero

        IF (this % u > zero) THEN
            get_specific_internal_energy = this % u
            RETURN
        END IF

        SELECT CASE (this % region)

        CASE (1)
            this % u = u1(this % pressure, this % temperature)

        CASE (2)
            this % u = u2(this % pressure, this % temperature)

        CASE (3)
            this % u = u3(this % density, this % temperature)

        CASE (4)
            this % u = u4(this % temperature, this % quality)

        CASE (5)
            this % u = u5(this % pressure, this % temperature)

        END SELECT

        get_specific_internal_energy = this % u

    END FUNCTION get_specific_internal_energy
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific enthalpy of the steam state
    !
    !> @return The specific enthalpy [J/kg]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_specific_enthalpy(this)

        USE SteamRegion1, ONLY: h1 => specific_enthalpy
        USE SteamRegion2, ONLY: h2 => specific_enthalpy
        USE SteamRegion3, ONLY: h3 => specific_enthalpy
        USE SteamRegion4, ONLY: h4 => specific_enthalpy
        USE SteamRegion5, ONLY: h5 => specific_enthalpy

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_specific_enthalpy
        !
        ! Compute the specific enthalpy if not already computed
        !
        get_specific_enthalpy = zero

        IF (this % h > zero) THEN
            get_specific_enthalpy = this % h
            RETURN
        END IF

        SELECT CASE (this % region)

        CASE (1)
            this % h = h1(this % pressure, this % temperature)

        CASE (2)
            this % h = h2(this % pressure, this % temperature)

        CASE (3)
            this % h = h3(this % density, this % temperature)

        CASE (4)
            this % h = h4(this % temperature, this % quality)

        CASE (5)
            this % h = h5(this % pressure, this % temperature)

        END SELECT

        get_specific_enthalpy = this % h

    END FUNCTION get_specific_enthalpy
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific entropy of the steam state
    !
    !> @return The specific entropy [J/kg.K]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_specific_entropy(this)

        USE SteamRegion1, ONLY: s1 => specific_entropy
        USE SteamRegion2, ONLY: s2 => specific_entropy
        USE SteamRegion3, ONLY: s3 => specific_entropy
        USE SteamRegion4, ONLY: s4 => specific_entropy
        USE SteamRegion5, ONLY: s5 => specific_entropy

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_specific_entropy
        !
        ! Compute the specific entropy if not already computed
        !
        get_specific_entropy = zero

        IF (this % s > zero) THEN
            get_specific_entropy = this % s
            RETURN
        END IF

        SELECT CASE (this % region)

        CASE (1)
            this % s = s1(this % pressure, this % temperature)

        CASE (2)
            this % s = s2(this % pressure, this % temperature)

        CASE (3)
            this % s = s3(this % density, this % temperature)

        CASE (4)
            this % s = s4(this % temperature, this % quality)

        CASE (5)
            this % s = s5(this % pressure, this % temperature)

        END SELECT

        get_specific_entropy = this % s

    END FUNCTION get_specific_entropy
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific isobaric heat capacity of the
    !> steam state
    !
    !> @return The specific isobaric heat capacity [J/kg.K]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_specific_isobaric_heat_capacity(this)

        USE SteamRegion1, ONLY: cp1 => specific_isobaric_heat_capacity
        USE SteamRegion2, ONLY: cp2 => specific_isobaric_heat_capacity
        USE SteamRegion3, ONLY: cp3 => specific_isobaric_heat_capacity
        USE SteamRegion4, ONLY: cp4 => specific_isobaric_heat_capacity
        USE SteamRegion5, ONLY: cp5 => specific_isobaric_heat_capacity

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_specific_isobaric_heat_capacity
        !
        ! Compute the specific isobaric heat capacity if not already computed
        !
        get_specific_isobaric_heat_capacity = zero

        IF (this % cp > zero) THEN
            get_specific_isobaric_heat_capacity = this % cp
            RETURN
        END IF

        SELECT CASE (this % region)

        CASE (1)
            this % cp = cp1(this % pressure, this % temperature)

        CASE (2)
            this % cp = cp2(this % pressure, this % temperature)

        CASE (3)
            this % cp = cp3(this % density, this % temperature)

        CASE (4)
            this % cp = cp4(this % temperature, this % quality)

        CASE (5)
            this % cp = cp5(this % pressure, this % temperature)

        END SELECT

        get_specific_isobaric_heat_capacity = this % cp

    END FUNCTION get_specific_isobaric_heat_capacity
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific isochoric heat capacity of the
    !> steam state
    !
    !> @return The specific isochoric heat capacity [J/kg.K]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_specific_isochoric_heat_capacity(this)

        USE SteamRegion1, ONLY: cv1 => specific_isochoric_heat_capacity
        USE SteamRegion2, ONLY: cv2 => specific_isochoric_heat_capacity
        USE SteamRegion3, ONLY: cv3 => specific_isochoric_heat_capacity
        USE SteamRegion4, ONLY: cv4 => specific_isochoric_heat_capacity
        USE SteamRegion5, ONLY: cv5 => specific_isochoric_heat_capacity

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_specific_isochoric_heat_capacity
        !
        ! Compute the specific isochoric heat capacity if not already computed
        !
        get_specific_isochoric_heat_capacity = zero

        IF (this % cv > zero) THEN
            get_specific_isochoric_heat_capacity = this % cv
            RETURN
        END IF

        SELECT CASE (this % region)

        CASE (1)
            this % cv = cv1(this % pressure, this % temperature)

        CASE (2)
            this % cv = cv2(this % pressure, this % temperature)

        CASE (3)
            this % cv = cv3(this % density, this % temperature)

        CASE (4)
            this % cv = cv4(this % temperature, this % quality)

        CASE (5)
            this % cv = cv5(this % pressure, this % temperature)

        END SELECT

        get_specific_isochoric_heat_capacity = this % cv

    END FUNCTION get_specific_isochoric_heat_capacity
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the specific isochoric heat capacity of the
    !> steam state
    !
    !> @return The specific isochoric heat capacity [J/kg.K]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_ratio_of_specific_heats(this)

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_ratio_of_specific_heats
        !
        ! Compute the ratio of specific heats if not already computed
        !
        this % ratio_cp_cv = get_specific_isobaric_heat_capacity(this) / &
                             get_specific_isochoric_heat_capacity(this)

        get_ratio_of_specific_heats = this % ratio_cp_cv

    END FUNCTION get_ratio_of_specific_heats
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the viscosity of the steam state
    !
    !> @return The viscosity in [Pa.sec]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_viscosity(this)

        USE SteamViscosity

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_viscosity
        !
        ! Compute viscosity if not already computed
        !
        get_viscosity = zero

        IF (this % mu > zero) THEN
            get_viscosity = this % mu
            RETURN
        END IF

        this % mu = viscosity(get_density(this), this % temperature)
        get_viscosity = this % mu

    END FUNCTION get_viscosity
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function returns the thermal conductivity of the steam state
    !
    !> @return The thermal conductivity in [W/K.m]
    !
    !-------------------------------------------------------------------------
    FUNCTION get_thermal_conductivity(this)

        USE SteamConductivity

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        REAL(KIND=REAL_HIGH) :: get_thermal_conductivity
        !
        ! Compute the thermal conductivity if not already computed
        !
        get_thermal_conductivity = zero

        IF (this % tk > zero) THEN
            get_thermal_conductivity = this % tk
            RETURN
        END IF

        this % tk = thermal_conductivity(get_density(this), this % temperature)
        get_thermal_conductivity = this % tk

    END FUNCTION get_thermal_conductivity
    !-----------------------------------------------------------------------
    !
    !> @details
    !> This function prints the state to the screen
    !
    !-----------------------------------------------------------------------
    SUBROUTINE print_state(this)

        IMPLICIT NONE
        !
        ! Arguments
        !
        CLASS(SteamState), INTENT(INOUT) :: this
        !
        ! Print the state
        !
        WRITE (*, *)
        WRITE (*, *) "Steam State - Thermodynamic Properties"
        WRITE (*, *) "======================================"

        WRITE (*, *)
        WRITE (*, *) "Steam Region                     = ", &
            this % region

        WRITE (*, *) "Pressure                         = ", &
            get_pressure(this), " [MPa]"

        WRITE (*, *) "Density                          = ", &
            get_density(this), " [Kg/m3]"

        WRITE (*, *) "Temperature                      = ", &
            get_temperature(this), " [K]"

        WRITE (*, *) "Specific Volume                  = ", &
            get_specific_volume(this), " [m3/Kg]"

        WRITE (*, *) "Specific Internal Energy         = ", &
            get_specific_internal_energy(this), " [J/Kg]"

        WRITE (*, *) "Specific Enthalpy                = ", &
            get_specific_enthalpy(this), " [J/Kg]"

        WRITE (*, *) "Specific Entropy                 = ", &
            get_specific_entropy(this), " [J/Kg.K]"

        WRITE (*, *) "Specific Isobaric Heat Capacity  = ", &
            get_specific_isobaric_heat_capacity(this), " [J/Kg.K]"

        WRITE (*, *) "Specific Isochoric Heat Capacity = ", &
            get_specific_isochoric_heat_capacity(this), " [J/Kg.K]"

        WRITE (*, *) "Ratio of Specific Heats          = ", &
            get_ratio_of_specific_heats(this)

        WRITE (*, *) "Viscosity                        = ", &
            get_viscosity(this), " [Pa.sec]"

        WRITE (*, *) "Thermal Conductivity             = ", &
            get_thermal_conductivity(this), " [W/K.m]"

        WRITE (*, *)

    END SUBROUTINE print_state

END MODULE SteamStates
