!-----------------------------------------------------------------------
!
!> @file SteamBoundaries.f
!
!> @details
!> This module implements the polynomials which compute the thermodynamic
!> properties of steam on various boundary lines between regions according
!> to the revised release on the IAPWS Industrial Formulation 1997 for the
!> Thermodynamic Properties of Water and Steam, August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamBoundaries

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Region boundary functions
    !
    PUBLIC boundary23_pressure
    PUBLIC boundary23_temperature
    PUBLIC boundary2ab_enthalpy
    PUBLIC boundary2bc_pressure
    PUBLIC boundary2bc_enthalpy
    PUBLIC boundary3ab_enthalpy
    PUBLIC boundary34_pressure_h
    PUBLIC boundary34_pressure_s
    !
    !> Constant coefficients "n"
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:5) :: n = &
                                                       (/0.34805185628969D+03, -0.11671859879975D+01, &
                                                         0.10192970039326D-02, 0.57254459862746D+03, &
                                                         0.13918839778870D+02/)
    !
    !> Constant coefficients "m"
    !
    REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:5) :: m = &
                                                       (/0.90584278514723D+03, -0.67955786399241D+00, &
                                                         0.12809002730136D-03, 0.26526571908428D+04, &
                                                         0.45257578905948D+01/)

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function computes pressure on the boundary line between regions
    !> 2 and 3 with respect to temperature
    !
    !> @param[in]  temperature The steam temperature [K]
    !> @param[out] valid       Becomes true if the temperature supplied as
    !                   an argument is within the valid range for
    !                   boundary 2-3.
    !
    !> @return The pressure [MPa]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary23_pressure(temperature, valid)

        IMPLICIT NONE
        !
        !> Star pressure for the computations on the boundary line between
        !> regions 2 and 3 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_3_PSTAR = 1.0D+00
        !
        !> Star temperature for the computations on the boundary line between
        !> regions 2 and 3 in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_3_TSTAR = 1.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: temperature
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary23_pressure
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: t
        !
        ! Compute the pressure
        !
        IF (present(valid)) THEN
            valid = (temperature >= REGION_1_TMAX) .AND. &
                    (temperature <= REGION_2_4_T)
        END IF

        t = temperature / REGION_2_3_TSTAR
        boundary23_pressure = REGION_2_3_PSTAR * (n(1) + (n(2) + n(3) * t) * t)

    END FUNCTION boundary23_pressure
    !-------------------------------------------------------------------------
    !
    !> This function computes temperature on the boundary line between regions
    !> 2 and 3 with respect to pressure
    !
    !> @param[in]  pressure The steam pressure [MPa]
    !> @param[out] valid    Becomes true if the pressure supplied as
    !                an argument is within the valid range for
    !                boundary 2-3.
    !
    !> @return The temperature [K]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary23_temperature(pressure, valid)

        IMPLICIT NONE
        !
        !> Star pressure for the computations on the boundary line between
        !> regions 2 and 3 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_3_PSTAR = 1.0D+00
        !
        !> Star temperature for the computations on the boundary line between
        !> regions 2 and 3 in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2_3_TSTAR = 1.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary23_temperature
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: p
        !
        ! Compute the pressure
        !
        IF (present(valid)) THEN
            valid = (pressure >= REGION_1_4_TMAX_P) .AND. &
                    (pressure <= IAPWS97_PMAX)
        END IF

        p = pressure / REGION_2_3_PSTAR
        boundary23_temperature = REGION_2_3_TSTAR * (n(4) + sqrt((p - n(5)) / n(3)))

    END FUNCTION boundary23_temperature
    !-------------------------------------------------------------------------
    !
    !> This function computes enthalpy on the boundary line between regions
    !> 2a and 2b with respect to the specific entropy
    !
    !> @param[in]  entropy The steam specific entropy [J/Kg.K]
    !> @param[out] valid   Becomes true if the entropy supplied as
    !               an argument is within the valid range for
    !               boundary 2a and 2b.
    !
    !> @return The pressure [MPa]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary2ab_enthalpy(entropy, valid)

        USE SteamRegion2, ONLY: s2 => specific_entropy
        USE SteamRegion4, ONLY: tsat => saturation_temperature

        IMPLICIT NONE
        !
        !> Constant coefficients "l"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:4) :: l = &
                                                           (/-0.349898083432139D+04, 0.257560716905876D+04, &
                                                             -0.421073558227969D+03, 0.276349063799944D+02/)
        !
        !> Star enthalpy for the computations on the boundary line between
        !> sub-regions 2a and 2b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2A_2B_HSTAR = 1000.0D+00
        !
        !> Star entropy for the computations on the boundary line between
        !> sub-regions 2a and 2b in [J/Kg.K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2A_2B_SSTAR = 1000.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: entropy
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary2ab_enthalpy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: sigma, smin
        !
        ! Compute the pressure
        !
        IF (present(valid)) THEN
            smin = s2(REGION_2A_2B_P, tsat(REGION_2A_2B_P))
            valid = (entropy >= smin) .AND. &
                    (entropy <= s2(REGION_2A_2B_P, REGION_2_TMAX))
        END IF

        sigma = entropy / REGION_2A_2B_SSTAR
        boundary2ab_enthalpy = REGION_2A_2B_HSTAR * (l(1) + (l(2) + (l(3) + l(4) * sigma) * sigma) * sigma)

    END FUNCTION boundary2ab_enthalpy
    !-------------------------------------------------------------------------
    !
    !> This function computes pressure on the boundary line between regions
    !> 2b and 2c with respect to the specific enthalpy
    !
    !> @param[in] enthalpy The steam specific enthalpy [J/Kg]
    !> @param[out] valid    Becomes true if the enthalpy supplied as
    !                an argument is within the valid range for
    !                boundary 2b and 2c.
    !
    !> @return The pressure [MPa]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary2bc_pressure(enthalpy, valid)

        IMPLICIT NONE
        !
        !> Star pressure for the computations on the boundary line between
        !> sub-regions 2b and 2c in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_2C_PSTAR = 1.0D+00
        !
        !> Star enthalpy for the computations on the boundary line between
        !> sub-regions 2b and 2c in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_2C_HSTAR = 1000.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: enthalpy
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary2bc_pressure
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: e
        !
        ! Compute the pressure
        !
        IF (present(valid)) THEN
            valid = (enthalpy >= boundary2bc_enthalpy(REGION_2B_2C_PMIN)) .AND. &
                    (enthalpy <= boundary2bc_enthalpy(IAPWS97_PMAX))
        END IF

        e = enthalpy / REGION_2B_2C_HSTAR
        boundary2bc_pressure = REGION_2B_2C_PSTAR * (m(1) + (m(2) + m(3) * e) * e)

    END FUNCTION boundary2bc_pressure
    !-------------------------------------------------------------------------
    !
    !> This function computes enthalpy on the boundary line between regions
    !> 2b and 2c with respect to pressure
    !
    !> @param[in]  pressure The steam pressure [MPa]
    !> @param[out] valid    Becomes true if the pressure supplied as
    !                an argument is within the valid range for
    !                boundary 2b-2c.
    !
    !> @return The enthalpy in [J/Kg]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary2bc_enthalpy(pressure, valid)

        IMPLICIT NONE
        !
        !> Star pressure for the computations on the boundary line between
        !> sub-regions 2b and 2c in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_2C_PSTAR = 1.0D+00
        !
        !> Star enthalpy for the computations on the boundary line between
        !> sub-regions 2b and 2c in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_2B_2C_HSTAR = 1000.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary2bc_enthalpy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: p
        !
        ! Compute the pressure
        !
        IF (present(valid)) THEN
            valid = (pressure >= REGION_2B_2C_PMIN) .AND. &
                    (pressure <= IAPWS97_PMAX)
        END IF

        p = pressure / REGION_2B_2C_PSTAR
        boundary2bc_enthalpy = REGION_2B_2C_HSTAR * (m(4) + sqrt((p - m(5)) / m(3)))

    END FUNCTION boundary2bc_enthalpy
    !-------------------------------------------------------------------------
    !
    !> This function computes enthalpy on the boundary line between regions
    !> 3a and 3b with respect to pressure
    !
    !> @param[in]  pressure The steam pressure [MPa]
    !> @param[out] valid    Becomes true if the pressure supplied as
    !                       an argument is within the valid range for
    !                       boundary 3a-3b.
    !
    !> @return The enthalpy in [J/Kg]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary3ab_enthalpy(pressure, valid)

        IMPLICIT NONE
        !
        !> Constant coefficients "l"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:4) :: l = &
                                                           (/0.201464004206875D+04, 0.374696550136983D+01, &
                                                             -0.219921901054187D-01, 0.875131686009950D-04/)
        !
        !> Star pressure for the computations on the boundary line between
        !> sub-regions 3a and 3b in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_3B_PSTAR = 1.0D+00
        !
        !> Star enthalpy for the computations on the boundary line between
        !> sub-regions 3a and 3b in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3A_3B_HSTAR = 1000.0D+00
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: pressure
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary3ab_enthalpy
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: p
        !
        ! Compute the enthalpy
        !
        IF (present(valid)) THEN
            valid = (pressure >= IAPWS97_PCRIT) .AND. &
                    (pressure <= IAPWS97_PMAX)
        END IF

        p = pressure / REGION_3A_3B_PSTAR
        boundary3ab_enthalpy = REGION_3A_3B_HSTAR * (l(1) + p * (l(2) + p * (l(3) + l(4) * p)))

    END FUNCTION boundary3ab_enthalpy
    !-------------------------------------------------------------------------
    !
    !> This function computes pressure on the boundary line between regions
    !> 3 and 4 with respect to enthalpy
    !
    !> @param[in]  enthalpy The steam enthalpy [KJ/Kg]
    !> @param[out] valid    Becomes true if the enthalpy supplied as
    !                       an argument is within the valid range for
    !                       boundary 3-4.
    !
    !> @return The pressure in [MPa]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary34_pressure_h(enthalpy, valid)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:14) :: II = &
                                                              (/0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:14) :: JJ = &
                                                              (/0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:14) :: nn = &
                                                            (/0.600073641753024D+00, -0.936203654849857D+01, &
                                                              0.246590798594147D+02, -0.107014222858224D+03, &
                                                              -0.915821315805768D+14, -0.862332011700662D+04, &
                                                              -0.235837344740032D+02, 0.252304969384128D+18, &
                                                              -0.389718771997719D+19, -0.333775713645296D+23, &
                                                              0.356499469636328D+11, -0.148547544720641D+27, &
                                                              0.330611514838798D+19, 0.813641294467829D+38/)
        !
        !> Star pressure for the computations on the boundary line between
        !> regions 3 and 4 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_PSTAR = 22.0D+00
        !
        !> Star enthalpy for the computations on the boundary line between
        !> regions 3 and 4 in [J/Kg]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_HSTAR = 2600.0D+03
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: enthalpy
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary34_pressure_h
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: eta, eta1, eta2
        !
        ! Compute the pressure
        !
        IF (present(valid)) THEN
            valid = (enthalpy >= REGION_3_4_HMIN) .AND. &
                    (enthalpy <= REGION_3_4_HMAX)
        END IF

        eta = enthalpy / REGION_3_4_HSTAR
        eta1 = eta - 1.020D+00
        eta2 = eta - 0.608D+00
        boundary34_pressure_h = REGION_3_4_PSTAR * SUM(nn * (eta1**II) * (eta2**JJ))

    END FUNCTION boundary34_pressure_h
    !-------------------------------------------------------------------------
    !
    !> This function computes pressure on the boundary line between regions
    !> 3 and 4 with respect to entropy
    !
    !> @param[in]  entropy The steam entropy [KJ/Kg.K]
    !> @param[out] valid   Becomes true if the entropy supplied as
    !                      an argument is within the valid range for
    !                      boundary 3-4.
    !
    !> @return The pressure in [MPa]
    !
    !-------------------------------------------------------------------------
    FUNCTION boundary34_pressure_s(entropy, valid)

        IMPLICIT NONE
        !
        ! Constant polynomial cofficients
        !
        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:10) :: II = &
                                                              (/0, 1, 1, 4, 12, 12, 16, 24, 28, 32/)

        INTEGER(KIND=INT_HIGH), PARAMETER, DIMENSION(1:10) :: JJ = &
                                                              (/0, 1, 32, 7, 4, 14, 36, 10, 0, 18/)

        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:10) :: nn = &
                                                            (/0.639767553612785D+00, -0.129727445396014D+02, &
                                                              -0.224595125848403D+16, 0.177466741801846D+07, &
                                                              0.717079349571538D+10, -0.378829107169011D+18, &
                                                              -0.955586736431328D+35, 0.187269814676188D+24, &
                                                              0.119254746466473D+12, 0.110649277244882D+37/)
        !
        !> Star pressure for the computations on the boundary line between
        !> regions 3 and 4 in [MPa]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_PSTAR = 22.0D+00
        !
        !> Star entropy for the computations on the boundary line between
        !> regions 3 and 4 in [J/Kg.K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: REGION_3_4_SSTAR = 5.2D+03
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: entropy
        LOGICAL, OPTIONAL, INTENT(OUT) :: valid
        REAL(KIND=REAL_HIGH) :: boundary34_pressure_s
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: sigma, sigma1, sigma2
        !
        ! Compute the pressure
        !
        IF (present(valid)) THEN
            valid = (entropy >= REGION_3_4_SMIN) .AND. &
                    (entropy <= REGION_3_4_SMAX)
        END IF

        sigma = entropy / REGION_3_4_SSTAR
        sigma1 = sigma - 1.030D+00
        sigma2 = sigma - 0.699D+00
        boundary34_pressure_s = REGION_3_4_PSTAR * SUM(nn * (sigma1**II) * (sigma2**JJ))

    END FUNCTION boundary34_pressure_s

END MODULE SteamBoundaries

