!-----------------------------------------------------------------------
!
!> @file SteamTablesTesting.f
!
!> @details
!> This module tests the computed thermodynamic properties of steam for
!> all regions according to the revised release on the IAPWS Industrial
!> Formulation 1997 for the Thermodynamic Properties of Water and Steam,
!> August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamTablesTesting

    USE SteamConstants
    USE SteamRegions
    USE SteamStates
    USE SteamNumerics

CONTAINS
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for region 1
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_region1

        USE SteamRegion1

        IMPLICIT NONE
        !
        ! Variables
        !
        INTEGER :: i
        REAL(KIND=REAL_HIGH) :: w_REAL(3), h_REAL(3), s_REAL(3), &
                                v_REAL(3), u_REAL(3), cp_REAL(3)
        REAL(KIND=REAL_HIGH) :: t(3), p(3), h, s, value, error
        LOGICAL :: valid_state
        TYPE(SteamState) :: ss
        !
        ! Test for Region 1 according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Region 1 Test - Computational Errors"
        WRITE (*, *) "===================================="
        !
        ! Input data
        !
        p(1) = 3.0D+00
        t(1) = 300.0D+00

        p(2) = 80.0D+00
        t(2) = 300.0D+00

        p(3) = 3.0D+00
        t(3) = 500.0D+00
        !
        ! Expected values
        !
        v_REAL(1) = 0.100215168D-02
        h_REAL(1) = 0.115331273D+06
        u_REAL(1) = 0.112324818D+06
        s_REAL(1) = 0.392294792D+03
        cp_REAL(1) = 0.417301218D+04
        w_REAL(1) = 0.150773921D+04

        v_REAL(2) = 0.971180894D-03
        h_REAL(2) = 0.184142828D+06
        u_REAL(2) = 0.106448356D+06
        s_REAL(2) = 0.368563852D+03
        cp_REAL(2) = 0.401008987D+04
        w_REAL(2) = 0.163469054D+04

        v_REAL(3) = 0.120241800D-02
        h_REAL(3) = 0.975542239D+06
        u_REAL(3) = 0.971934985D+06
        s_REAL(3) = 0.258041912D+04
        cp_REAL(3) = 0.465580682D+04
        w_REAL(3) = 0.124071337D+04

        DO i = 1, 3

            WRITE (*, *) "============"
            WRITE (*, *) "Test Case : ", i
            WRITE (*, *) "============"

            CALL set_primary_properties(ss, 1, p(i), t(i), valid_state)

            WRITE (*, *) "Steam state in region 1 ?                = ", valid_state

            value = speed_of_sound(p(i), t(i))
            error = ABS(value - w_REAL(i)) / ABS(w_REAL(i))
            WRITE (*, *) "Error in speed of sound                  = ", error

            value = get_specific_volume(ss)
            error = ABS(value - v_REAL(i)) / ABS(v_REAL(i))
            WRITE (*, *) "Error in specific volume                 = ", error

            value = get_specific_enthalpy(ss)
            error = ABS(value - h_REAL(i)) / ABS(h_REAL(i))
            WRITE (*, *) "Error in specific enthalpy               = ", error

            value = get_specific_entropy(ss)
            error = ABS(value - s_REAL(i)) / ABS(s_REAL(i))
            WRITE (*, *) "Error in specific entropy                = ", error

            value = get_specific_internal_energy(ss)
            error = ABS(value - u_REAL(i)) / ABS(u_REAL(i))
            WRITE (*, *) "Error in specific internal energy        = ", error

            value = get_specific_isobaric_heat_capacity(ss)
            error = ABS(value - cp_REAL(i)) / ABS(cp_REAL(i))
            WRITE (*, *) "Error in specific isobaric heat capacity = ", error

            h = get_specific_enthalpy(ss)
            s = get_specific_entropy(ss)

            value = temperature_ph(p(i), h_REAL(i))
            error = ABS(value - t(i)) / ABS(t(i))
            WRITE (*, *) "Error in backward temperature(p,h)       = ", error

            value = pressure_hs(h_REAL(i), s_REAL(i))
            error = ABS(value - p(i)) / ABS(p(i))
            WRITE (*, *) "Error in backward pressure(h,s)          = ", error

            WRITE (*, *)
            WRITE (*, *) "Attempting to compute temperature by inverting polynomial h = h(p,t)"
            value = function_inverter_y(specific_enthalpy, h_REAL(i), p(i), &
                                        200.0D+00, 600.0D+00, 1.0D-05, &
                                        100, 3, .TRUE.)

            WRITE (*, *)
            WRITE (*, *) "Attempting to compute enthalpy by inverting polynomial p = p(h,s)"
            value = function_inverter_x(pressure_hs, p(i), s_REAL(i), &
                                        0.10D+06, 1.0D+07, 1.0D-05, &
                                        100, 3, .TRUE.)
        END DO

    END SUBROUTINE test_region1
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for region 2
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_region2

        USE SteamRegion2

        IMPLICIT NONE
        !
        ! Variables
        !
        INTEGER :: i
        LOGICAL :: in2a, in2b, in2c, valid_state
        REAL(KIND=REAL_HIGH) :: w_REAL(3), h_REAL(3), s_REAL(3), &
                                v_REAL(3), u_REAL(3), cp_REAL(3)
        REAL(KIND=REAL_HIGH) :: t(3), p(3), t_ph, t_ps, p_hs, value, error
        TYPE(SteamState) :: ss
        !
        ! Test for Region 1 according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Region 2 Test - Computational Errors"
        WRITE (*, *) "===================================="
        !
        ! Input data
        !
        p(1) = 0.0035D+00
        t(1) = 300.0D+00

        p(2) = 0.0035D+00
        t(2) = 700.0D+00

        p(3) = 30.0D+00
        t(3) = 700.0D+00

        !
        ! Expected values
        !
        v_REAL(1) = 0.394913866D+02
        h_REAL(1) = 0.254991145D+07
        u_REAL(1) = 0.241169160D+07
        s_REAL(1) = 0.852238967D+04
        cp_REAL(1) = 0.191399162D+04
        w_REAL(1) = 0.427920172D+03

        v_REAL(2) = 0.923015898D+02
        h_REAL(2) = 0.333568375D+07
        u_REAL(2) = 0.301262819D+07
        s_REAL(2) = 0.101749996D+05
        cp_REAL(2) = 0.208141274D+04
        w_REAL(2) = 0.644289068D+03

        v_REAL(3) = 0.542946619D-02
        h_REAL(3) = 0.263149474D+07
        u_REAL(3) = 0.246861076D+07
        s_REAL(3) = 0.517540298D+04
        cp_REAL(3) = 0.103505092D+05
        w_REAL(3) = 0.480386523D+03

        DO i = 1, 3

            WRITE (*, *) "============"
            WRITE (*, *) "Test Case : ", i
            WRITE (*, *) "============"

            CALL set_primary_properties(ss, 2, p(i), t(i), valid_state)

            in2a = in_region2a(p(i), t(i))
            in2b = in_region2b(p(i), t(i))
            in2c = in_region2c(p(i), t(i))

            WRITE (*, *) "Steam state in region 2 ?                = ", valid_state
            WRITE (*, *) "Steam state in region 2a ?               = ", in2a
            WRITE (*, *) "Steam state in region 2b ?               = ", in2b
            WRITE (*, *) "Steam state in region 2c ?               = ", in2c

            value = speed_of_sound(p(i), t(i))
            error = ABS(value - w_REAL(i)) / ABS(w_REAL(i))
            WRITE (*, *) "Error in speed of sound                  = ", error

            value = get_specific_volume(ss)
            error = ABS(value - v_REAL(i)) / ABS(v_REAL(i))
            WRITE (*, *) "Error in specific volume                 = ", error

            value = get_specific_enthalpy(ss)
            error = ABS(value - h_REAL(i)) / ABS(h_REAL(i))
            WRITE (*, *) "Error in specific enthalpy               = ", error

            value = get_specific_entropy(ss)
            error = ABS(value - s_REAL(i)) / ABS(s_REAL(i))
            WRITE (*, *) "Error in specific entropy                = ", error

            value = get_specific_internal_energy(ss)
            error = ABS(value - u_REAL(i)) / ABS(u_REAL(i))
            WRITE (*, *) "Error in specific internal energy        = ", error

            value = get_specific_isobaric_heat_capacity(ss)
            error = ABS(value - cp_REAL(i)) / ABS(cp_REAL(i))
            WRITE (*, *) "Error in specific isobaric heat capacity = ", error

            IF (in2a) THEN
                t_ph = temperature_ph_region2a(p(i), h_REAL(i))
                t_ps = temperature_ps_region2a(p(i), s_REAL(i))
                p_hs = pressure_hs_region2a(h_REAL(i), s_REAL(i))
            ELSE IF (in2b) THEN
                t_ph = temperature_ph_region2b(p(i), h_REAL(i))
                t_ps = temperature_ps_region2b(p(i), s_REAL(i))
                p_hs = pressure_hs_region2b(h_REAL(i), s_REAL(i))
            ELSE
                t_ph = temperature_ph_region2c(p(i), h_REAL(i))
                t_ps = temperature_ps_region2c(p(i), s_REAL(i))
                p_hs = pressure_hs_region2c(h_REAL(i), s_REAL(i))
            END IF

            error = ABS(t(i) - t_ph) / ABS(t_ph)
            WRITE (*, *) "Error in backward temperature(p,h)       = ", error

            error = ABS(t(i) - t_ps) / ABS(t_ps)
            WRITE (*, *) "Error in backward temperature(p,s)       = ", error

            error = ABS(p(i) - p_hs) / ABS(p_hs)
            WRITE (*, *) "Error in backward pressure(h,s)          = ", error

        END DO

    END SUBROUTINE test_region2
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for region 3
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_region3

        USE SteamRegion3

        IMPLICIT NONE
        !
        ! Variables
        !
        INTEGER :: i
        REAL(KIND=REAL_HIGH) :: w_REAL(3), h_REAL(3), s_REAL(3), &
                                p_REAL(3), u_REAL(3), cp_REAL(3)
        REAL(KIND=REAL_HIGH) :: rho(3), t(3), t_ph, v_ph, t_ps, v_ps, &
                                value, error
        LOGICAL :: valid_state, in3a, in3b
        TYPE(SteamState) :: ss
        !
        ! Test for Region 1 according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Region 3 Test - Computational Errors"
        WRITE (*, *) "===================================="
        !
        ! Input data
        !
        rho(1) = 500.0D+00
        t(1) = 650.0D+00

        rho(2) = 200.0D+00
        t(2) = 650.0D+00

        rho(3) = 500.0D+00
        t(3) = 750.0D+00
        !
        ! Expected values
        !
        p_REAL(1) = 0.255837018D+02
        h_REAL(1) = 0.186343019D+07
        u_REAL(1) = 0.181226279D+07
        s_REAL(1) = 0.405427273D+04
        cp_REAL(1) = 0.138935717D+05
        w_REAL(1) = 0.502005554D+03

        p_REAL(2) = 0.222930643D+02
        h_REAL(2) = 0.237512401D+07
        u_REAL(2) = 0.226365868D+07
        s_REAL(2) = 0.485438792D+04
        cp_REAL(2) = 0.446579342D+05
        w_REAL(2) = 0.383444594D+03

        p_REAL(3) = 0.783095639D+02
        h_REAL(3) = 0.225868845D+07
        u_REAL(3) = 0.210206932D+07
        s_REAL(3) = 0.446971906D+04
        cp_REAL(3) = 0.634165359D+04
        w_REAL(3) = 0.760696041D+03

        DO i = 1, 3

            WRITE (*, *) "============"
            WRITE (*, *) "Test Case : ", i
            WRITE (*, *) "============"

            CALL set_primary_properties(ss, 3, rho(i), t(i), valid_state)

            WRITE (*, *) "Steam state in region 3 ?                = ", valid_state

            in3a = in_region3a(rho(i), t(i))
            in3b = in_region3b(rho(i), t(i))

            WRITE (*, *) "Steam state in region 3a ?               = ", in3a
            WRITE (*, *) "Steam state in region 3b ?               = ", in3b

            value = speed_of_sound(rho(i), t(i))
            error = ABS(value - w_REAL(i)) / ABS(w_REAL(i))
            WRITE (*, *) "Error in speed of sound                  = ", error

            value = get_pressure(ss)
            error = ABS(value - p_REAL(i)) / ABS(p_REAL(i))
            WRITE (*, *) "Error in pressure                        = ", error

            value = get_specific_enthalpy(ss)
            error = ABS(value - h_REAL(i)) / ABS(h_REAL(i))
            WRITE (*, *) "Error in specific enthalpy               = ", error

            value = get_specific_entropy(ss)
            error = ABS(value - s_REAL(i)) / ABS(s_REAL(i))
            WRITE (*, *) "Error in specific entropy                = ", error

            value = get_specific_internal_energy(ss)
            error = ABS(value - u_REAL(i)) / ABS(u_REAL(i))
            WRITE (*, *) "Error in specific internal energy        = ", error

            value = get_specific_isobaric_heat_capacity(ss)
            error = ABS(value - cp_REAL(i)) / ABS(cp_REAL(i))
            WRITE (*, *) "Error in specific isobaric heat capacity = ", error

            IF (in3a) THEN
                t_ph = temperature_ph_region3a(p_REAL(i), h_REAL(i))
                t_ps = temperature_ps_region3a(p_REAL(i), s_REAL(i))
                v_ph = specific_volume_ph_region3a(p_REAL(i), h_REAL(i))
                v_ps = specific_volume_ps_region3a(p_REAL(i), s_REAL(i))
            ELSE
                t_ph = temperature_ph_region3b(p_REAL(i), h_REAL(i))
                t_ps = temperature_ps_region3b(p_REAL(i), s_REAL(i))
                v_ph = specific_volume_ph_region3b(p_REAL(i), h_REAL(i))
                v_ps = specific_volume_ps_region3b(p_REAL(i), s_REAL(i))
            END IF

            error = ABS(t(i) - t_ph) / ABS(t_ph)
            WRITE (*, *) "Error in backward temperature(p,h)       = ", error

            error = ABS(t(i) - t_ps) / ABS(t_ps)
            WRITE (*, *) "Error in backward temperature(p,s)       = ", error

            error = ABS(1.0D+00 / rho(i) - v_ph) / ABS(v_ph)
            WRITE (*, *) "Error in backward specific volume(p,h)   = ", error

            error = ABS(1.0D+00 / rho(i) - v_ps) / ABS(v_ps)
            WRITE (*, *) "Error in backward specific volume(p,s)   = ", error

            WRITE (*, *)
            WRITE (*, *) "Attempting to compute temperature by inverting polynomial h = h(rho,t)"
            value = function_inverter_y(specific_enthalpy, h_REAL(i), rho(i), &
                                        500.0D+00, 800.0D+00, 1.0D-05, &
                                        100, 3, .TRUE.)
            WRITE (*, *)
            WRITE (*, *) "Attempting to compute density by inverting polynomial p = p(rho,t)"
            value = function_inverter_x(pressure, p_REAL(i), t(i), &
                                        100.0D+00, 600.0D+00, 1.0D-05, &
                                        100, 3, .TRUE.)

        END DO

    END SUBROUTINE test_region3
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for region 4
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_region4

        USE SteamRegion4

        IMPLICIT NONE
        !
        ! Variables
        !
        INTEGER :: i
        LOGICAL :: valid_state
        REAL(KIND=REAL_HIGH) :: p(3), t(3), psat, tsat, error
        TYPE(SteamState) :: ss
        !
        ! Test for Region 1 according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Region 4 Test - Computational Errors"
        WRITE (*, *) "===================================="

        p(1) = 0.353658941D-02
        t(1) = 300.0D+00

        p(2) = 0.263889776D+01
        t(2) = 500.0D+00

        p(3) = 0.123443146D+02
        t(3) = 600.0D+00

        DO i = 1, 3

            WRITE (*, *) "============="
            WRITE (*, *) "Test Case A : ", i
            WRITE (*, *) "============="

            CALL set_primary_properties(ss, 4, t(i), 0.5D+00, valid_state)

            WRITE (*, *) "Steam state in region 4 ?                = ", valid_state

            psat = saturation_pressure(t(i))
            error = ABS(psat - p(i)) / p(i)
            WRITE (*, *) "Error in saturation pressure             = ", error

        END DO

        p(1) = 0.1D+00
        t(1) = 0.372755919D+03

        p(2) = 1.0D+00
        t(2) = 0.453035632D+03

        p(3) = 10.0D+00
        t(3) = 0.584149488D+03

        DO i = 1, 3

            WRITE (*, *) "============="
            WRITE (*, *) "Test Case B : ", i
            WRITE (*, *) "============="

            CALL set_primary_properties(ss, 4, t(i), 0.5D+00, valid_state)

            WRITE (*, *) "Steam state in region 4 ?                = ", valid_state

            tsat = saturation_temperature(p(i))
            error = ABS(tsat - t(i)) / t(i)
            WRITE (*, *) "Error in saturation temperature          = ", error

        END DO

    END SUBROUTINE test_region4
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for the boundary line
    !> for regions 2 and 3
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_boundary23

        USE SteamBoundaries

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: t, p, t_real, p_real, error
        !
        ! Test for the boundary 2-3 according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Boundary 2-3 Test - Computational Errors"
        WRITE (*, *) "========================================"

        p_real = 0.165291643D+02
        t_real = 0.623150000D+03

        p = boundary23_pressure(t_real)
        error = ABS(p - p_real) / p_real
        WRITE (*, *) "Error in pressure                          = ", error

        t = boundary23_temperature(p_real)
        error = ABS(t - t_real) / t_real
        WRITE (*, *) "Error in temperature                       = ", error

    END SUBROUTINE test_boundary23
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for the boundary line
    !> for regions 2b and 2c
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_boundary2bc

        USE SteamBoundaries

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: h, p, h_real, p_real, error
        !
        ! Test for the boundary 2-3 according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Boundary 2b-2c Test - Computational Errors"
        WRITE (*, *) "=========================================="

        p_real = 0.1000000000D+03
        h_real = 0.3516004323D+07

        p = boundary2bc_pressure(h_real)
        error = ABS(p - p_real) / p_real
        WRITE (*, *) "Error in pressure                          = ", error

        h = boundary2bc_enthalpy(p_real)
        error = ABS(h - h_real) / h_real
        WRITE (*, *) "Error in specific enthalpy                 = ", error

    END SUBROUTINE test_boundary2bc
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for the boundary line
    !> for regions 2a and 2b
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_boundary2ab

        USE SteamBoundaries

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: s_real, h_real, h, error
        !
        ! Test for the boundary 2a-2b according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Boundary 2a-2b Test - Computational Errors"
        WRITE (*, *) "=========================================="

        s_real = 7000D+00
        h_real = 3376437.884D+00

        h = boundary2ab_enthalpy(s_real)
        error = ABS(h - h_real) / h_real
        WRITE (*, *) "Error in specific enthalpy                 = ", error

    END SUBROUTINE test_boundary2ab
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the steam polynomials for the boundary line
    !> for regions 3a and 3b
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_boundary3ab

        USE SteamBoundaries

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: p_real, h_real, h, error
        !
        ! Test for the boundary 3a-3b according to the IAPWS-IF97 document
        !
        WRITE (*, *)
        WRITE (*, *) "Boundary 3a-3b Test - Computational Errors"
        WRITE (*, *) "=========================================="

        p_real = 25D+00
        h_real = 2.095936454D+06

        h = boundary3ab_enthalpy(p_real)
        error = ABS(h - h_real) / h_real
        WRITE (*, *) "Error in specific enthalpy                 = ", error

    END SUBROUTINE test_boundary3ab
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the piecewise cubic interpolation function
    !> for computing the saturation pressure
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_cubic_interpolation

        USE SteamRegion4
        USE SteamTables

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), DIMENSION(:), ALLOCATABLE :: psat, pprime, tsat
        REAL(KIND=REAL_HIGH) :: rf, dt, tmin, tmax, t, error, p_int, p_real
        REAL(KIND=REAL_HIGH), DIMENSION(1:4) :: ptab, ttab
        INTEGER(KIND=INT_HIGH) :: seed, i, imax, num, iloc, iloc_searched, indx1, indx2

        imax = 1000
        tmin = 273.15D+00
        tmax = 647.096D+00
        dt = (tmax - tmin) / REAL(imax - 1, kind=REAL_HIGH)
        !
        ! Test for the Cubic interpolation (1D)
        !
        WRITE (*, *)
        WRITE (*, *) "Cubic Interpolation - Computational Errors"
        WRITE (*, *) "=========================================="
        !
        ! Compute the table data
        !
        ALLOCATE (psat(1:imax), tsat(1:imax), pprime(1:imax))

        FORALL (i=1:imax) &
            tsat(i) = tmin + dt * REAL(i - 1, kind=REAL_HIGH)

        FORALL (i=1:imax) &
            psat(i) = saturation_pressure(tsat(i))

        FORALL (i=1:imax) &
            pprime(i) = saturation_pressure_gradient(tsat(i))
        !
        ! Produce interpolated data at random locations
        !
        num = 10

        CALL RANDOM_SEED(seed)

        DO i = 1, num

            CALL RANDOM_NUMBER(rf)
            t = tmin + (tmax - tmin) * rf

            iloc = 1 + FLOOR((t - tsat(1)) / dt)
            iloc_searched = binary_search_vector(tsat, t, 1, imax)

            p_real = saturation_pressure(t)
            p_int = interpolate_cubic(tsat(iloc:iloc + 1), psat(iloc:iloc + 1), &
                                      pprime(iloc:iloc + 1), t)

            error = 100.0D+00 * ABS(p_int - p_real) / p_real

            WRITE (*, *) "Polynomial Pressure    = ", p_real, " [MPa]"

            WRITE (*, *) "Cubic with Exact Gradients:"
            WRITE (*, *) "Interpolated Pressure  = ", p_int, " [MPa]"
            WRITE (*, *) "Interpolation Error    = ", error, " [%]"

            ptab(2) = psat(iloc)
            ptab(3) = psat(iloc + 1)
            ttab(2) = tsat(iloc)
            ttab(3) = tsat(iloc + 1)

            IF (iloc == 1) THEN
                ptab(1) = psat(iloc)
                ttab(1) = tsat(iloc)
            ELSE
                ptab(1) = psat(iloc - 1)
                ttab(1) = tsat(iloc - 1)
            END IF

            IF (iloc == imax - 1) THEN
                ptab(4) = psat(iloc + 1)
                ttab(4) = tsat(iloc + 1)
            ELSE
                ptab(4) = psat(iloc + 2)
                ttab(4) = tsat(iloc + 2)
            END IF

            p_int = interpolate_cubic_numerical(ttab, ptab, t)

            error = 100.0D+00 * ABS(p_int - p_real) / p_real

            WRITE (*, *) "Cubic with Numerical Gradients:"
            WRITE (*, *) "Interpolated Pressure  = ", p_int, " [MPa]"
            WRITE (*, *) "Interpolation Error    = ", error, " [%]"

            WRITE (*, '(1x,a,i5,i5)') &
                "Interpolation index    = ", iloc, iloc_searched
            WRITE (*, *) "-------------------------"

        END DO

        DEALLOCATE (psat, tsat, pprime)

    END SUBROUTINE test_cubic_interpolation
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the piecewise bilinear interpolation function
    !> for computing the specific enthalpy in region 5 as a function of
    !> pressure and temperature
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_bilinear_interpolation

        USE SteamRegion5
        USE SteamTables

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), DIMENSION(:, :), ALLOCATABLE :: h
        REAL(KIND=REAL_HIGH), DIMENSION(:), ALLOCATABLE :: p, t
        REAL(KIND=REAL_HIGH) :: dt, tmin, tmax, dp, pmin, pmax, pp, tt, &
                                h_int, h_real, error, rf
        INTEGER(KIND=INT_HIGH) :: seed, i, j, imax, jmax, num, &
                                  iloc, jloc, iloc_searched, jloc_searched

        imax = 100
        jmax = 100

        pmin = 0.1D-02
        pmax = REGION_5_PMAX
        dp = (pmax - pmin) / REAL(imax - 1, kind=REAL_HIGH)

        tmin = REGION_2_TMAX
        tmax = IAPWS97_TMAX
        dt = (tmax - tmin) / REAL(jmax - 1, kind=REAL_HIGH)
        !
        ! Test for the bilinear interpolation (2D)
        !
        WRITE (*, *)
        WRITE (*, *) "Bilinear Interpolation - Computational Errors"
        WRITE (*, *) "============================================="
        !
        ! Compute the table data
        !
        ALLOCATE (p(1:imax), t(1:jmax), h(1:imax, 1:jmax))

        FORALL (i=1:imax) &
            p(i) = pmin + dp * REAL(i - 1, kind=REAL_HIGH)

        FORALL (j=1:jmax) &
            t(j) = tmin + dt * REAL(j - 1, kind=REAL_HIGH)

        FORALL (i=1:imax, j=1:jmax) &
            h(i, j) = specific_enthalpy(p(i), t(j))
        !
        ! Produce interpolated data at random locations
        !
        num = 10

        CALL RANDOM_SEED(seed)

        DO i = 1, num

            CALL RANDOM_NUMBER(rf)
            tt = tmin + (tmax - tmin) * rf

            CALL RANDOM_NUMBER(rf)
            pp = pmin + (pmax - pmin) * rf

            iloc = 1 + FLOOR((pp - p(1)) / dp)
            jloc = 1 + FLOOR((tt - t(1)) / dt)

            iloc_searched = binary_search_vector(p, pp, 1, imax)
            jloc_searched = binary_search_vector(t, tt, 1, jmax)

            h_int = interpolate_bilinear(p(iloc:iloc + 1), t(jloc:jloc + 1), &
                                         h(iloc:iloc + 1, jloc:jloc + 1), pp, tt)

            h_real = specific_enthalpy(pp, tt)
            error = ABS(h_int - h_real) / ABS(h_real)

            WRITE (*, *) "Polynomial Specific Enthalpy    = ", h_real, " [J/Kg]"
            WRITE (*, *) "Interpolation Specific Enthalpy = ", h_int, " [J/Kg]"
            WRITE (*, *) "Estimation Error                = ", error, " [%]"
            WRITE (*, '(1x,4(a,i5),a)') &
                "Interpolation Indices           = (", &
                iloc, ",", jloc, ") - (", iloc_searched, ",", &
                jloc_searched, ")"
            WRITE (*, *) "----------------------------------"

        END DO

        DEALLOCATE (p, t, h)

    END SUBROUTINE test_bilinear_interpolation
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the generation of a table for the computation
    !> of the specific enthalpy in region 5
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_table_generation

        USE SteamRegion5
        USE SteamTables

        IMPLICIT NONE
        !
        ! Variables
        !
        LOGICAL :: passed
        INTEGER(KIND=INT_HIGH) :: i, seed, iloc, jloc, num_errors, &
                                  imax, jmax
        REAL(KIND=REAL_HIGH) :: error, rf, xmin, xmax, ymin, ymax, &
                                er, h_real, h_int, xp, yp
        REAL(KIND=REAL_HIGH), DIMENSION(:), ALLOCATABLE :: x, y
        REAL(KIND=REAL_HIGH), DIMENSION(:, :), ALLOCATABLE :: h
        !
        ! Test for the bicubic interpolation (1D)
        !
        WRITE (*, *)
        WRITE (*, *) "Table Generation Test : Specific Enthalpy in R5"
        WRITE (*, *) "==============================================="

        error = 1.0D-05

        xmin = 0.1D-02
        xmax = REGION_5_PMAX
        ymin = REGION_2_TMAX
        ymax = IAPWS97_TMAX

        CALL generate_table(specific_enthalpy, xmin, xmax, ymin, ymax, &
                            0.025D+00, 0.50D+00, error, x, y, h, .TRUE.)
        !
        ! Verify the accuracy of the table
        !
        CALL RANDOM_SEED(seed)

        WRITE (*, *)
        WRITE (*, *) "Table Verification - Interpolation Errors"
        WRITE (*, *) "========================================="
        WRITE (*, *) "Testing 10E+06 random combinations of properties"

        passed = .TRUE.
        imax = size(x)
        jmax = size(y)
        num_errors = 0

!$OMP         parallel do schedule(static), &
!$OMP         firstprivate(imax,jmax,xmin,ymin,xmax,ymax,error), &
!$OMP         private(i,iloc,jloc,rf,xp,yp,h_real,h_int,er), &
!$OMP         shared(x,y,h,num_errors,passed)
        DO i = 1, 10000000

            CALL RANDOM_NUMBER(rf)
            xp = xmin + (xmax - xmin) * rf

            CALL RANDOM_NUMBER(rf)
            yp = ymin + (ymax - ymin) * rf

            h_real = specific_enthalpy(xp, yp)

            iloc = binary_search_vector(x, xp, 1, imax)
            jloc = binary_search_vector(y, yp, 1, jmax)

            h_int = interpolate_bilinear(x(iloc:iloc + 1), y(jloc:jloc + 1), &
                                         h(iloc:iloc + 1, jloc:jloc + 1), xp, yp)

            er = ABS(h_int - h_real) / ABS(h_real)

            IF (er > error) THEN
!$OMP                 critical
                passed = .FALSE.
                num_errors = num_errors + 1
!$OMP                 end critical
            END IF

        END DO
!$OMP         end parallel do

        IF (passed) THEN
            WRITE (*, *) "All computed errors are within the specified limit."
            WRITE (*, *)
        ELSE
            WRITE (*, *) num_errors, " error(s) have been detected."
            WRITE (*, *)
        END IF

    END SUBROUTINE test_table_generation
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the computation of viscosity
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_viscosity

        USE SteamViscosity

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), DIMENSION(1:11), PARAMETER :: t = &
                                                            (/298.15D+00, 298.15D+00, 373.15D+00, 433.15D+00, 433.15D+00, &
                                                              873.15D+00, 873.15D+00, 873.15D+00, 1173.15D+00, &
                                                              1173.15D+00, 1173.15D+00/)

        REAL(KIND=REAL_HIGH), DIMENSION(1:11), PARAMETER :: rho = &
                                                            (/998.0D+00, 1200.0D+00, 1000.0D+00, 1.0D+00, 1000.0D+00, 1.0D+00, &
                                                              100.0D+00, 600.0D+00, 1.0D+00, 100.0D+00, 400.0D+00/)

        REAL(KIND=REAL_HIGH), DIMENSION(1:11), PARAMETER :: mu = &
                                                            (/889.735100D-06, 1437.649467D-06, 307.883622D-06, &
                                                              14.538324D-06, 217.685358D-06, 32.619287D-06, &
                                                              35.802262D-06, 77.430195D-06, 44.217245D-06, &
                                                              47.640433D-06, 64.154608D-06/)

        INTEGER(KIND=INT_HIGH) :: i
        REAL(KIND=REAL_HIGH) :: visc, error
        !
        ! Test the computation of viscosity
        !
        WRITE (*, *)
        WRITE (*, *) "Viscosity Test - Computational Errors"
        WRITE (*, *) "=========================================="

        DO i = 1, 11
            visc = viscosity(rho(i), t(i))
            error = ABS(visc - mu(i)) / mu(i)
            WRITE (*, *) "Error in viscosity                         = ", error
        END DO

    END SUBROUTINE test_viscosity
    !-----------------------------------------------------------------------
    !
    !> This subroutine tests the computation of thermal conductivity
    !
    !-----------------------------------------------------------------------
    SUBROUTINE test_thermal_conductivity

        USE SteamConductivity

        IMPLICIT NONE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), DIMENSION(1:4), PARAMETER :: t = &
                                                           (/298.15D+00, 298.15D+00, 298.15D+00, 873.15D+00/)

        REAL(KIND=REAL_HIGH), DIMENSION(1:4), PARAMETER :: rho = &
                                                           (/0.0D+00, 998.0D+00, 1200.0D+00, 0.0D+00/)

        REAL(KIND=REAL_HIGH), DIMENSION(1:4), PARAMETER :: k = &
                                                           (/18.4341883D-03, 607.712868D-03, 799.038144D-03, 79.1034659D-03/)

        INTEGER(KIND=INT_HIGH) :: i
        REAL(KIND=REAL_HIGH) :: hk, error
        !
        ! Test the computation of thermal conductivity
        !
        WRITE (*, *)
        WRITE (*, *) "Thermal Conductivity Test - Computational Errors"
        WRITE (*, *) "================================================"

        DO i = 1, 4
            hk = thermal_conductivity(rho(i), t(i))
            error = ABS(hk - k(i)) / k(i)
            WRITE (*, *) "Error in thermal conductivity              = ", error
        END DO

    END SUBROUTINE test_thermal_conductivity

END MODULE SteamTablesTesting
