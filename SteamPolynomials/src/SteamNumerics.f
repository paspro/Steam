!-----------------------------------------------------------------------
!
!> @file SteamNumerics.f
!
!> @details
!> This module contains numerical algorithms for various purposes.
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamNumerics

    USE SteamPrecision

    IMPLICIT NONE
    PRIVATE
    !
    ! Public functions
    !
    PUBLIC function_inverter
    PUBLIC function_inverter_x
    PUBLIC function_inverter_y
    PUBLIC interpolate_cubic
    PUBLIC interpolate_cubic_numerical
    PUBLIC interpolate_bilinear
    PUBLIC interpolate_bicubic
    PUBLIC binary_search_vector
    PUBLIC binary_search_array
    PUBLIC quicksort

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function reverses a generic equation f(x) by computing the value
    !> of "x0" when the value of the equation f(x0) is known with the use of
    !> an accelerated Secant root-finder algorithm operating on the function
    !> g(x) = f(x) - f0, where f0 = f(x0) -> g(x0) = 0.
    !
    !> @param[in] f         The function to invert
    !> @param[in] f0        The value of f0 = f(x0,y0)
    !> @param[in] guess1    An guess for the value of x0
    !> @param[in] guess2    A second guess for the value of x0
    !> @param[in] tolerance The tolerance to utilise for the numerical method
    !> @param[in] maxiter   The maximum number of allowed iterations
    !> @param[in] n_order   The interpolation order used by the method
    !> @param[in] logging   If true then write to the screen logging data
    !
    !> @return The value of "x0" when f(x0) is known
    !
    !-------------------------------------------------------------------------
    FUNCTION function_inverter(f, f0, guess1, guess2, tolerance, &
                               maxiter, n_order, logging)

        IMPLICIT NONE
        !
        ! Arguments
        !
        LOGICAL :: logging
        REAL(KIND=REAL_HIGH) :: function_inverter
        INTEGER(KIND=INT_HIGH), INTENT(IN) :: maxiter, n_order
        REAL(KIND=REAL_HIGH), INTENT(IN) :: f0, guess1, guess2, &
                                            tolerance

        INTERFACE
            PURE FUNCTION f(x)
                USE SteamPrecision
                IMPLICIT NONE
                REAL(KIND=REAL_HIGH) :: f
                REAL(KIND=REAL_HIGH), INTENT(IN) :: x
            END FUNCTION f
        END INTERFACE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: error, s, s_old
        INTEGER(KIND=INT_HIGH) :: iter
        !
        ! Accelerated Secant Algorithm
        !
        IF (n_order < 0) THEN
            WRITE (*, *)
            WRITE (*, *) "Error in Function: function_inverter"
            WRITE (*, *) "The order of interpolation must be a positive integer."
            WRITE (*, *)
        END IF

        IF (logging) THEN
            WRITE (*, *)
            WRITE (*, *) "========================================="
            WRITE (*, *) "Function f(x) Inverter Algorithm for x"
            WRITE (*, *) "========================================="
        END IF

        s_old = guess2

        DO iter = n_order + 1, maxiter
            !
            ! Compute the new value
            !
            s = solution(iter, n_order)
            error = ABS(s_old - s) / ABS(s)
            s_old = s

            IF (logging) THEN
                WRITE (*, "(' Iteration = ', i5, ' X = ', d12.5, &
                        & ' Error = ', d12.5)") &
                      iter, s, error
            END IF
            !
            ! Check for convergence
            !
            IF (error <= tolerance) THEN
                EXIT
            END IF
        END DO

        IF (logging) THEN
            WRITE (*, *) "========================================="
            WRITE (*, *)
        END IF

        function_inverter = s

    CONTAINS
        !---------------
        !> imax
        !---------------
        FUNCTION imax(p)

            IMPLICIT NONE
            !
            ! Arguments
            !
            INTEGER(KIND=INT_HIGH) :: imax
            INTEGER(KIND=INT_HIGH), INTENT(IN) :: p

            IF ((p == -1) .OR. (p == 0)) THEN
                imax = 0
            ELSE IF (p <= n_order) THEN
                imax = p - 1
            ELSE
                imax = n_order
            END IF

        END FUNCTION imax
        !--------------------------------
        !> Solution
        !--------------------------------
        RECURSIVE FUNCTION solution(p, i) RESULT(res)

            IMPLICIT NONE
            !
            ! Arguments
            !
            REAL(KIND=REAL_HIGH) :: res
            INTEGER(KIND=INT_HIGH), INTENT(IN) :: p, i
            !
            ! Variables
            !
            REAL(KIND=REAL_HIGH) :: x1, x2, x3, x4, f1, f2

            IF (p == -1) THEN
                res = guess1
            ELSE IF (p == 0) THEN
                res = guess2
            ELSE
                IF (i == 0) THEN
                    !
                    ! Secant step
                    !
                    x1 = solution(p - 1, imax(p - 1))
                    x2 = solution(p - 2, imax(p - 2))
                    f1 = f(x1) - f0
                    f2 = f(x2) - f0
                    res = x1 - f1 * (x1 - x2) / (f1 - f2)
                ELSE
                    !
                    ! i-th order approximation
                    !
                    x1 = solution(p - 1, i - 1)
                    x2 = solution(p - 1, imax(p - 1))
                    x3 = solution(p, i - 1)
                    x4 = solution(p - i - 2, imax(p - i - 2))
                    res = x3 + (x2 - x3) * (x1 - x3) / (x1 + x2 - x3 - x4)
                END IF
            END IF

        END FUNCTION solution

    END FUNCTION function_inverter
    !-------------------------------------------------------------------------
    !
    !> This function reverses a generic equation f(x,y) by computing the value
    !> of "x0" when the value of the equation f(x0,y0) and y0 are known using
    !> an accelerated Secant root-finder algorithm operating on the function
    !> g(x,y0) = f(x,y0) - f0, where f0 = f(x0,y0) -> g(x0,y0) = 0.
    !
    !> @param[in] f         The function to invert
    !> @param[in] f0        The value of f0 = f(x0,y0)
    !> @param[in] y0        The value y0.
    !> @param[in] guess1    An guess for the value of x0
    !> @param[in] guess2    A second guess for the value of x0
    !> @param[in] tolerance The tolerance to utilise for the numerical method
    !> @param[in] maxiter   The maximum number of allowed iterations
    !> @param[in] n_order   The interpolation order used by the method
    !> @param[in] logging   If true then write to the screen logging data
    !
    !> @return The value of "x0" when f(x0,y0) is known
    !
    !-------------------------------------------------------------------------
    FUNCTION function_inverter_x(f, f0, y0, guess1, guess2, tolerance, &
                                 maxiter, n_order, logging)

        IMPLICIT NONE
        !
        ! Arguments
        !
        LOGICAL :: logging
        REAL(KIND=REAL_HIGH) :: function_inverter_x
        INTEGER(KIND=INT_HIGH), INTENT(IN) :: maxiter, n_order
        REAL(KIND=REAL_HIGH), INTENT(IN) :: f0, y0, guess1, guess2, &
                                            tolerance

        INTERFACE
            PURE FUNCTION f(x, y)
                USE SteamPrecision
                IMPLICIT NONE
                REAL(KIND=REAL_HIGH) :: f
                REAL(KIND=REAL_HIGH), INTENT(IN) :: x, y
            END FUNCTION f
        END INTERFACE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: error, s, s_old
        INTEGER(KIND=INT_HIGH) :: iter
        !
        ! Accelerated Secant Algorithm
        !
        IF (n_order < 0) THEN
            WRITE (*, *)
            WRITE (*, *) "Error in Function: function_inverter_x"
            WRITE (*, *) "The order of interpolation must be a positive integer."
            WRITE (*, *)
        END IF

        IF (logging) THEN
            WRITE (*, *)
            WRITE (*, *) "========================================="
            WRITE (*, *) "Function f(x,y0) Inverter Algorithm for x"
            WRITE (*, *) "========================================="
        END IF

        s_old = guess2

        DO iter = n_order + 1, maxiter
            !
            ! Compute the new value
            !
            s = solution(iter, n_order)
            error = ABS(s_old - s) / ABS(s)
            s_old = s

            IF (logging) THEN
                WRITE (*, "(' Iteration = ', i5, ' X = ', d12.5, &
                        & ' Error = ', d12.5)") &
                      iter, s, error
            END IF
            !
            ! Check for convergence
            !
            IF (error <= tolerance) THEN
                EXIT
            END IF
        END DO

        IF (logging) THEN
            WRITE (*, *) "========================================="
            WRITE (*, *)
        END IF

        function_inverter_x = s

    CONTAINS
        !---------------
        !> imax
        !---------------
        FUNCTION imax(p)

            IMPLICIT NONE
            !
            ! Arguments
            !
            INTEGER(KIND=INT_HIGH) :: imax
            INTEGER(KIND=INT_HIGH), INTENT(IN) :: p

            IF ((p == -1) .OR. (p == 0)) THEN
                imax = 0
            ELSE IF (p <= n_order) THEN
                imax = p - 1
            ELSE
                imax = n_order
            END IF

        END FUNCTION imax
        !--------------------------------
        !> Solution
        !--------------------------------
        RECURSIVE FUNCTION solution(p, i) RESULT(res)

            IMPLICIT NONE
            !
            ! Arguments
            !
            REAL(KIND=REAL_HIGH) :: res
            INTEGER(KIND=INT_HIGH), INTENT(IN) :: p, i
            !
            ! Variables
            !
            REAL(KIND=REAL_HIGH) :: x1, x2, x3, x4, f1, f2

            IF (p == -1) THEN
                res = guess1
            ELSE IF (p == 0) THEN
                res = guess2
            ELSE
                IF (i == 0) THEN
                    !
                    ! Secant step
                    !
                    x1 = solution(p - 1, imax(p - 1))
                    x2 = solution(p - 2, imax(p - 2))
                    f1 = f(x1, y0) - f0
                    f2 = f(x2, y0) - f0
                    res = x1 - f1 * (x1 - x2) / (f1 - f2)
                ELSE
                    !
                    ! i-th order approximation
                    !
                    x1 = solution(p - 1, i - 1)
                    x2 = solution(p - 1, imax(p - 1))
                    x3 = solution(p, i - 1)
                    x4 = solution(p - i - 2, imax(p - i - 2))
                    res = x3 + (x2 - x3) * (x1 - x3) / (x1 + x2 - x3 - x4)
                END IF
            END IF

        END FUNCTION solution

    END FUNCTION function_inverter_x
    !-------------------------------------------------------------------------
    !
    !> This function reverses a generic equation f(x,y) by computing the value
    !> of "y0" when the value of the equation f(x0,y0) and x0 are known using
    !> the numerical rootfinder algorithm of Newton-Raphson operating on the
    !> function g(x0,y) = f(x0,y) - f0, where f0 = f(x0,y0) -> g(x0,y0) = 0.
    !
    !> @param[in] f         The function to invert
    !> @param[in] f0        The value of f0 = f(x0,y0)
    !> @param[in] x0        The value x0.
    !> @param[in] guess1    A guess for the value of x0
    !> @param[in] guess2    A second guess for the value of x0
    !> @param[in] tolerance The tolerance to utilise for the numerical method
    !> @param[in] maxiter   The maximum number of allowed iterations
    !> @param[in] n_order   The interpolation order used by the method
    !> @param[in] logging   If true then write to the screen logging data
    !
    !> @return The value of "y0" when f(x0,y0) is known
    !
    !-------------------------------------------------------------------------
    FUNCTION function_inverter_y(f, f0, x0, guess1, guess2, tolerance, &
                                 maxiter, n_order, logging)

        IMPLICIT NONE
        !
        ! Arguments
        !
        LOGICAL :: logging
        REAL(KIND=REAL_HIGH) :: function_inverter_y
        INTEGER(KIND=INT_HIGH), INTENT(IN) :: maxiter, n_order
        REAL(KIND=REAL_HIGH), INTENT(IN) :: f0, x0, guess1, guess2, &
                                            tolerance

        INTERFACE
            PURE FUNCTION f(x, y)
                USE SteamPrecision
                IMPLICIT NONE
                REAL(KIND=REAL_HIGH) :: f
                REAL(KIND=REAL_HIGH), INTENT(IN) :: x, y
            END FUNCTION f
        END INTERFACE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: error, s, s_old
        INTEGER(KIND=INT_HIGH) :: iter
        !
        ! Accelerated Secant Algorithm
        !
        IF (n_order < 0) THEN
            WRITE (*, *)
            WRITE (*, *) "Error in Function: function_inverter_y"
            WRITE (*, *) "The order of interpolation must be a positive integer."
            WRITE (*, *)
        END IF

        IF (logging) THEN
            WRITE (*, *)
            WRITE (*, *) "========================================="
            WRITE (*, *) "Function f(x0,y) Inverter Algorithm for y"
            WRITE (*, *) "========================================="
        END IF

        s_old = guess2

        DO iter = n_order + 1, maxiter
            !
            ! Compute the new value
            !
            s = solution(iter, n_order)
            error = ABS(s_old - s) / ABS(s)
            s_old = s

            IF (logging) THEN
                WRITE (*, "(' Iteration = ', i5, ' Y = ', d12.5, &
                        & ' Error = ', d12.5)") &
                      iter, s, error
            END IF
            !
            ! Check for convergence
            !
            IF (error <= tolerance) THEN
                EXIT
            END IF
        END DO

        IF (logging) THEN
            WRITE (*, *) "========================================="
            WRITE (*, *)
        END IF

        function_inverter_y = s

    CONTAINS
        !---------------
        !> imax
        !---------------
        FUNCTION imax(p)

            IMPLICIT NONE
            !
            ! Arguments
            !
            INTEGER(KIND=INT_HIGH) :: imax
            INTEGER(KIND=INT_HIGH), INTENT(IN) :: p

            IF ((p == -1) .OR. (p == 0)) THEN
                imax = 0
            ELSE IF (p <= n_order) THEN
                imax = p - 1
            ELSE
                imax = n_order
            END IF

        END FUNCTION imax
        !--------------------------------
        !> Solution
        !--------------------------------
        RECURSIVE FUNCTION solution(p, i) RESULT(res)

            IMPLICIT NONE
            !
            ! Arguments
            !
            REAL(KIND=REAL_HIGH) :: res
            INTEGER(KIND=INT_HIGH), INTENT(IN) :: p, i
            !
            ! Variables
            !
            REAL(KIND=REAL_HIGH) :: y1, y2, y3, y4, f1, f2

            IF (p == -1) THEN
                res = guess1
            ELSE IF (p == 0) THEN
                res = guess2
            ELSE
                IF (i == 0) THEN
                    !
                    ! Secant step
                    !
                    y1 = solution(p - 1, imax(p - 1))
                    y2 = solution(p - 2, imax(p - 2))
                    f1 = f(x0, y1) - f0
                    f2 = f(x0, y2) - f0
                    res = y1 - f1 * (y1 - y2) / (f1 - f2)
                ELSE
                    !
                    ! i-th order approximation
                    !
                    y1 = solution(p - 1, i - 1)
                    y2 = solution(p - 1, imax(p - 1))
                    y3 = solution(p, i - 1)
                    y4 = solution(p - i - 2, imax(p - i - 2))
                    res = y3 + (y2 - y3) * (y1 - y3) / (y1 + y2 - y3 - y4)
                END IF
            END IF

        END FUNCTION solution

    END FUNCTION function_inverter_y
    !-------------------------------------------------------------------------
    !
    !> This function performs a piecewise cubic interpolation in order to
    !> compute the value of a property. The interpolation function has the
    !> form:
    !>
    !> s(x) = SUM(Cmj*(x-xj)^m), where the sum is over m = 0, 3
    !>        and Cmj are the interpolation coefficients computed based on
    !>        the supplied data in the form of vectors x, u(x), du(x)/dx
    !
    !> @param[in] x      The vector with the x-values (size=2)
    !> @param[in] u      The vector with the property values (size=2)
    !> @param[in] uprime The vector with the gradient of the property
    !> @param[in] xi     The x-value to be used for computing the interpolated
    !>                   property value which lies between x(1) and x(2)
    !
    !> @return The interpolated value at location xi
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION interpolate_cubic(x, u, uprime, xi)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: interpolate_cubic
        REAL(KIND=REAL_HIGH), INTENT(IN) :: xi
        REAL(KIND=REAL_HIGH), DIMENSION(2), INTENT(IN) :: x, u, uprime
        !
        ! Constants
        !
        INTEGER(KIND=INT_HIGH), DIMENSION(0:3), PARAMETER :: m = (/0, 1, 2, 3/)
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: dx, dx2, dx3, du
        REAL(KIND=REAL_HIGH), DIMENSION(0:3) :: coeff
        !
        ! Compute the interpolation coefficients
        !
        dx = x(2) - x(1)
        du = u(2) - u(1)
        dx2 = dx * dx
        dx3 = dx2 * dx

        coeff(0) = u(1)
        coeff(1) = uprime(1)
        coeff(2) = three * du / dx2 - (uprime(2) + two * uprime(1)) / dx
        coeff(3) = -two * du / dx3 + (uprime(2) + uprime(1)) / dx2
        !
        ! Perform the interpolation
        !
        interpolate_cubic = SUM(coeff * ((xi - x(1))**m))

    END FUNCTION interpolate_cubic
    !-------------------------------------------------------------------------
    !
    !> This function performs a cubic interpolation in order to compute
    !> the value of a property.
    !
    !> @param[in] x  The vector with the x-values (size=4)
    !> @param[in] u  The vector with the property values (size=4)
    !> @param[in] xi The x-value to be used for computing the interpolated
    !>               property value which lies between x(2) and x(3)
    !
    !> @return The interpolated value at location xi
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION interpolate_cubic_numerical(x, u, xi)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: interpolate_cubic_numerical
        REAL(KIND=REAL_HIGH), INTENT(IN) :: xi
        REAL(KIND=REAL_HIGH), DIMENSION(1:4), INTENT(IN) :: x, u
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), DIMENSION(1:4) :: dx, coeff, div
        !
        ! Compute the interpolation coefficients
        !
        dx(1:4) = (xi - x(1:4))
        div(1) = (x(1) - x(2)) * (x(1) - x(3)) * (x(1) - x(4))
        coeff(1) = dx(2) * dx(3) * dx(4) / div(1)
        div(2) = (x(2) - x(1)) * (x(2) - x(3)) * (x(2) - x(4))
        coeff(2) = dx(1) * dx(3) * dx(4) / div(2)
        div(3) = (x(3) - x(1)) * (x(3) - x(2)) * (x(3) - x(4))
        coeff(3) = dx(1) * dx(2) * dx(4) / div(3)
        div(4) = (x(4) - x(1)) * (x(4) - x(2)) * (x(4) - x(3))
        coeff(4) = dx(1) * dx(2) * dx(3) / div(4)
        !
        ! Perform the interpolation
        !
        interpolate_cubic_numerical = SUM(u * coeff)

    END FUNCTION interpolate_cubic_numerical
    !-------------------------------------------------------------------------
    !
    !> This function performs a piecewise bilinear interpolation in order to
    !> compute the value of a property that is a function of two parameters.
    !> The interpolation function has the form:
    !>
    !> s(x,y) = SUMn(SUMm(Cmnjk*(x-xj)^m*(y-yk)^n), where SUMm is over m = 0, 1,
    !>         SUMn is over n = 0, 1 and Cmnjk are the interpolation coefficients
    !>         computed based on the supplied data in the form of the table x, y,
    !>         u(x,y)
    !
    !> @param[in] x  The vector with the x-values (dim=2)
    !> @param[in] y  The vector with the y-values (dim=2)
    !> @param[in] u  The array with the property values (dim=2x2)
    !> @param[in] xi The x-value to be used for computing the interpolated
    !>               property value.
    !> @param[in] yi The y-value to be used for computing the interpolated
    !>               property value.
    !
    !> @return The interpolated value at location xi, yi
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION interpolate_bilinear(x, y, u, xi, yi)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: interpolate_bilinear
        REAL(KIND=REAL_HIGH), INTENT(IN) :: xi, yi
        REAL(KIND=REAL_HIGH), DIMENSION(1:2), INTENT(IN) :: x, y
        REAL(KIND=REAL_HIGH), DIMENSION(1:2, 1:2), INTENT(IN) :: u
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH) :: dx21, dy21, dx2i, dy2i, dxi1, dyi1, u1d, u2d
        !
        ! Compute the interpolation coefficients
        !
        dx21 = x(2) - x(1)
        dy21 = y(2) - y(1)
        dx2i = x(2) - xi
        dy2i = y(2) - yi
        dxi1 = xi - x(1)
        dyi1 = yi - y(1)
        !
        ! Perform the interpolation
        !
        u1d = (u(1, 1) * dx2i + u(2, 1) * dxi1) * dy2i
        u2d = (u(1, 2) * dx2i + u(2, 2) * dxi1) * dyi1
        interpolate_bilinear = (u1d+u2d) / (dx21 * dy21)

    END FUNCTION interpolate_bilinear
    !-------------------------------------------------------------------------
    !
    !> This function performs a bicubic interpolation in order to
    !> compute the value of a property that is a function of two parameters.
    !
    !> @param[in] x  The vector with the x-values (dim=4)
    !> @param[in] y  The vector with the y-values (dim=4)
    !> @param[in] u  The array with the property values (dim=4x4)
    !> @param[in] xi The x-value to be used for computing the interpolated
    !>               property value.
    !> @param[in] yi The y-value to be used for computing the interpolated
    !>               property value.
    !
    !> @return The interpolated value at location xi, yi
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION interpolate_bicubic(x, y, u, xi, yi)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH) :: interpolate_bicubic
        REAL(KIND=REAL_HIGH), INTENT(IN) :: xi, yi
        REAL(KIND=REAL_HIGH), DIMENSION(1:4), INTENT(IN) :: x, y
        REAL(KIND=REAL_HIGH), DIMENSION(1:4, 1:4), INTENT(IN) :: u
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), DIMENSION(1:4) :: uv
        !
        ! Perform a cubic interpolation for each column
        !
        uv(1) = interpolate_cubic_numerical(y, u(1, 1:4), yi)
        uv(2) = interpolate_cubic_numerical(y, u(2, 1:4), yi)
        uv(3) = interpolate_cubic_numerical(y, u(3, 1:4), yi)
        uv(4) = interpolate_cubic_numerical(y, u(4, 1:4), yi)
        !
        ! Perform a horizontal interpolation
        !
        interpolate_bicubic = interpolate_cubic_numerical(x, uv, xi)

        RETURN
    END FUNCTION interpolate_bicubic
    !-------------------------------------------------------------------------
    !
    !> This function implements a binary search algorithm in order to locate
    !> a value in a vector. If the exact value does not exist it will return
    !> the closest smaller value available in the table. Note that the
    !> supplied vector must be sorted in ascending order.
    !
    !> @param[in]  vec   The vector to use for searching
    !> @param[out] value The value to search for
    !> @param[in]  imin  The minimum vector index to use for searching
    !> @param[in]  imax  The maximum vector index to use for searching
    !
    !> @return The index with the location of the result
    !
    !-------------------------------------------------------------------------
    PURE FUNCTION binary_search_vector(vec, value, imin, imax)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), DIMENSION(:), INTENT(IN) :: vec
        INTEGER(KIND=INT_HIGH), INTENT(IN) :: imin, imax
        REAL(KIND=REAL_HIGH), INTENT(IN) :: value
        INTEGER(KIND=INT_HIGH) :: binary_search_vector
        !
        ! Variables
        !
        INTEGER(KIND=INT_HIGH) :: high, low, med
        REAL(KIND=REAL_HIGH) :: vmed
        !
        ! Search
        !
        high = imax
        low = imin
        med = imin + (imax - imin) / 2

        DO WHILE (high .NE. (low + 1))

            vmed = vec(med)

            IF (vmed <= value) THEN
                low = med
            ELSE
                high = med
            END IF

            med = low + (high - low) / 2

        END DO

        binary_search_vector = med
        RETURN

    END FUNCTION binary_search_vector
    !-------------------------------------------------------------------------
    !
    !> This subroutine implements a binary search algorithm in order to locate
    !> a value in a two dimensional array. If the exact value does not exist
    !> it will return the closest smaller value available in the table.
    !
    !> @param[in]  mat   The matrix (2D) to use for searching
    !> @param[out] value The value to search for
    !> @param[in]  imin  The minimum array i-index to use for searching
    !> @param[in]  imax  The maximum array i-index to use for searching
    !> @param[in]  jmin  The minimum array j-index to use for searching
    !> @param[in]  jmax  The maximum array j-index to use for searching
    !> @param[out] i     The i-index for the location of the value
    !> @param[out] j     The j-index for the location of the value
    !
    !-------------------------------------------------------------------------
    PURE SUBROUTINE binary_search_array(mat, value, imin, imax, jmin, jmax, i, j)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), DIMENSION(:, :), INTENT(IN) :: mat
        REAL(KIND=REAL_HIGH), INTENT(IN) :: value
        INTEGER(KIND=INT_HIGH), INTENT(IN) :: imin, imax, jmin, jmax
        INTEGER(KIND=INT_HIGH), INTENT(OUT) :: i, j
        !
        ! Variables
        !
        INTEGER(KIND=INT_HIGH) :: i_index, j_index
        REAL(KIND=REAL_HIGH) :: value_found, error_new, error_previous
        !
        ! Search
        !
        j_index = binary_search_vector(mat(imin, :), value, jmin, jmax)
        value_found = mat(imin, j_index)
        error_previous = (value - value_found) / ABS(value_found)

        i = imin
        j = j_index

        DO i_index = imin + 1, imax

            j_index = binary_search_vector(mat(i, :), value, jmin, jmax)
            value_found = mat(i, j)
            error_new = (value - value_found) / ABS(value_found)

            IF ((error_new >= zero) .AND. (error_new < error_previous)) THEN
                i = i_index
                j = j_index
            END IF

            error_new = error_previous

        END DO

    END SUBROUTINE binary_search_array
    !-------------------------------------------------------------------------
    !
    !> This subroutine implements the QuickSort algorithm for sorting the
    !> elements of a vector to ascending order.
    !
    !> @param[inout] vec  The vector to sort
    !> @param[out]   indx A vector of indices that represents the modification
    !>                   of the original vector. It can be used in order to
    !>                   apply similar modifications to other associated
    !>                   vectors.
    !
    !-------------------------------------------------------------------------
    RECURSIVE SUBROUTINE quicksort(vec, indx)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), DIMENSION(:), INTENT(INOUT) :: vec
        INTEGER(KIND=INT_HIGH), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: indx
        INTEGER(KIND=INT_HIGH) :: iq
        !
        ! Sorting
        !
        IF (size(vec) > 1) THEN

            IF (present(indx)) THEN
                CALL partition(vec, indx, iq)
                CALL quicksort(vec(:iq - 1), indx(:iq - 1))
                CALL quicksort(vec(iq:), indx(iq:))
            ELSE
                CALL partition(vec, marker=iq)
                CALL quicksort(vec(:iq - 1))
                CALL quicksort(vec(iq:))
            END IF

        END IF

    END SUBROUTINE quicksort
    !-------------------------------------------------------------------------
    !
    !> This subroutine partitions a vector for sorting purposes
    !
    !> @param[inout] vec    The vector to partition
    !> @param[inout] indx   A vector of indices to be partitioned in the
    !>                     same manner
    !> @param[out]   marker The partition marker
    !
    !-------------------------------------------------------------------------
    PURE SUBROUTINE partition(vec, indx, marker)

        IMPLICIT NONE
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), DIMENSION(:), INTENT(INOUT) :: vec
        INTEGER(KIND=INT_HIGH), DIMENSION(:), OPTIONAL, INTENT(INOUT) :: indx
        INTEGER(KIND=INT_HIGH), INTENT(OUT) :: marker
        !
        ! Variables
        !
        INTEGER(KIND=INT_HIGH) :: i, j, itemp
        REAL(KIND=REAL_HIGH) :: rtemp
        REAL(KIND=REAL_HIGH) :: x

        x = vec(1)
        i = 0
        j = size(vec) + 1

        DO
            j = j - 1
            DO
                IF (vec(j) <= x) EXIT
                j = j - 1
            END DO

            i = i + 1
            DO
                IF (vec(i) >= x) EXIT
                i = i + 1
            END DO

            IF (i < j) THEN

                rtemp = vec(i)
                vec(i) = vec(j)
                vec(j) = rtemp

                IF (present(indx)) THEN
                    itemp = indx(i)
                    indx(i) = indx(j)
                    indx(j) = itemp
                END IF

            ELSE IF (i == j) THEN
                marker = i + 1
                RETURN
            ELSE
                marker = i
                RETURN
            END IF
        END DO

    END SUBROUTINE partition

END MODULE SteamNumerics
