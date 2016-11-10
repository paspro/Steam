!-----------------------------------------------------------------------
!
!> @file SteamTables.f
!
!> @details
!> This module contains functions capable of generating steam tables
!> for thermodynamic properties at a specific level of accuracy.
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamTables

    USE SteamPrecision
    USE SteamNumerics

    IMPLICIT NONE
    PRIVATE
    !
    ! Public Procedures
    !
    PUBLIC generate_table

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This subroutine creates a table suitable for using in order to obtain
    !> interpolated properties that depend on two independent variables. The
    !> table is constructed such as it results in interpolatied values of a
    !> specific maximum allowed error.
    !
    !> @param[in] poly      A pointer to the polynomial function that
    !>                     generates the accurate, reference values to be
    !>                     used for building the table. This should be of
    !>                     the form: property = poly(x,y) where x, y are
    !>                     the independent variables
    !> @param[in] xmin      The minimum value for the variable x to use
    !> @param[in] xmax      The maximum value for the variable x to use
    !> @param[in] ymin      The minimum value for the variable y to use
    !> @param[in] ymax      The maximum value for the variable y to use
    !> @param[in] dx        The initial fine step to use for variable x
    !> @param[in] dy        The initial fine step to use for variable y
    !> @param[in] error     The maximum allowed interpolation error (%).
    !> @param[out] x        The 2D array with the values of the variable x
    !>                     forming the constructed table
    !> @param[out] y        The 2D array with the values of the variable y
    !>                     forming the constructed table
    !> @param[out] val      The 2D array with the values of the property the
    !>                     supplied polynomial function computes
    !> @param[in] logging   If true then information will be printed during
    !>                     the execution of the subroutine
    !-------------------------------------------------------------------------
    SUBROUTINE generate_table(poly, xmin, xmax, ymin, ymax, dx, dy, &
                              error, x, y, val, logging)

        IMPLICIT NONE
        !
        ! Arguments
        !
        LOGICAL :: logging
        REAL(KIND=REAL_HIGH), INTENT(IN) :: xmin, xmax, ymin, ymax, dx, dy, error
        REAL(KIND=REAL_HIGH), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: val
        REAL(KIND=REAL_HIGH), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: x, y

        INTERFACE
            PURE FUNCTION poly(x,y) RESULT(res)
                USE SteamPrecision
                REAL(KIND=REAL_HIGH), INTENT(IN) :: x, y
                REAL(KIND=REAL_HIGH) :: res
            END FUNCTION poly
        END INTERFACE
        !
        ! Variables
        !
        REAL(KIND=REAL_HIGH), DIMENSION(1:2) :: xp, yp
        REAL(KIND=REAL_HIGH), DIMENSION(1:2,1:2) :: up
        LOGICAL :: check_i, check_j, can_remove_line
        LOGICAL, DIMENSION(:), ALLOCATABLE :: i_used, j_used
        REAL(KIND=REAL_HIGH), DIMENSION(:), ALLOCATABLE :: x_fine, y_fine
        REAL(KIND=REAL_HIGH), DIMENSION(:,:), ALLOCATABLE :: val_fine
        REAL(KIND=REAL_HIGH) :: xk, yk, error_estimate, prop_real, prop_int

        INTEGER(KIND=INT_HIGH) :: i, j, ik, jk, imax, jmax, imax_new, jmax_new, i_p1, j_p1
        INTEGER(KIND=INT_HIGH) :: icheck, jcheck, icheck_p1, icheck_m1, &
                                  jcheck_p1, jcheck_m1, i_new, j_new,   &
                                  points, points_new, inum, jnum
        !
        ! Compute the number of points to use for the table given the initial dxi, dyi values
        !
        imax   = CEILING((xmax-xmin)/dx, kind=INT_HIGH) + 1
        jmax   = CEILING((ymax-ymin)/dy, kind=INT_HIGH) + 1
        points = imax*jmax

        IF (logging) THEN
            WRITE(*,*)
            WRITE(*,*) "===================================="
            WRITE(*,*) "Dynamic Table Generation Algorithm"
            WRITE(*,*) "===================================="
            WRITE(*,*) "Generating the initial table."
            WRITE(*,*) "Number of columns       = ", imax
            WRITE(*,*) "Number of rows          = ", jmax
            WRITE(*,*) "Total table data points = ", imax*jmax
            WRITE(*,*) "===================================="
            WRITE(*,*)
        END IF
        !
        ! Generate the fine table
        !
        ALLOCATE(x_fine(1:imax), y_fine(1:jmax), val_fine(1:imax, 1:jmax))
        ALLOCATE(i_used(1:imax), j_used(1:jmax))
        !
        ! Fill the fine table with values
        !
        i_used = .TRUE.
        j_used = .TRUE.

        !$omp parallel sections

        !$omp section
        FORALL (i = 1:imax)
            x_fine(i) = xmin + real(i-1, kind=REAL_HIGH)*dx
        END FORALL

        !$omp section
        FORALL (j = 1:jmax)
            y_fine(j) = ymin + real(j-1, kind=REAL_HIGH)*dy
        END FORALL
        !$omp end parallel sections

        FORALL (i = 1:imax, j = 1:jmax)
            val_fine(i,j) = poly(x_fine(i), y_fine(j))
        END FORALL
        !
        ! Check the initial table to verify that the interpolation
        ! errors for this table are bounded by the maximum allowed
        ! interpolation error. If this is not the case then the
        ! optimisation process makes no sense.
        !
        WRITE (*,*) "Checking the interpolation error for the initial table..."

        DO i = 1, imax-1
            xp(1) = x_fine(i)
            xp(2) = x_fine(i+1)

            DO j = 1, jmax-1
                yp(1) = y_fine(j)
                yp(2) = y_fine(j+1)

                up(1,1) = val_fine(i,   j)
                up(1,2) = val_fine(i,   j+1)
                up(2,1) = val_fine(i+1, j)
                up(2,2) = val_fine(i+1, j+1)

                DO ik = 0, 2
                    xk = xp(1) + half*real(ik, kind=8)*dx
                    DO jk = 0, 2
                        yk = yp(1) + half*real(jk, kind=8)*dy
                        !
                        ! Compute the accurate value of the property
                        !
                        prop_real = poly(xk, yk)
                        !
                        ! Compute the interpolated value of the property
                        !
                        prop_int = interpolate_bilinear(xp, yp, up, xk, yk)
                        !
                        ! Estimate the error and check it
                        !
                        error_estimate = 100.0D+00*ABS(prop_int-prop_real)/ &
                                                   ABS(prop_real)

                        IF (error_estimate > error) THEN
                            WRITE (*,*) "Encountered unacceptable error:"
                            WRITE (*,*) "Error = ", error_estimate
                            WRITE (*,*) "Location (i,j) = ", i, j
                            WRITE (*,*) "and (ik,jk)    = ", ik, jk
                            WRITE (*,*) "Properties have the values:"
                            WRITE (*,*) "xp = ", xp(1), xp(2)
                            WRITE (*,*) "yp = ", yp(1), yp(2)
                            WRITE (*,*) "up = ", up(1,1), up(1,2), up(2,1), up(2,2)
                            WRITE (*,*) "xk, yk = ", xk, yk
                            WRITE (*,*) "Polynomial property   = ", prop_real
                            WRITE (*,*) "Interpolated property = ", prop_int
                            WRITE (*,*) "Program terminates"
                            STOP
                        END IF

                    END DO
                END DO
            END DO
        END DO

        WRITE (*,*) "Interpolation errors are bounded."
        !
        ! Check the i and j lines alternatively to see
        ! if they can be removed and still have an interpolation
        ! error within the acceptable limit
        !
        IF (logging) THEN
            WRITE(*,*) "Performing table optimisation..."
        END IF
        !
        ! Start with the initial i and j lines that
        ! can be checked for removal
        !
        icheck  = 2
        jcheck  = 2
        check_i = .TRUE.
        check_j = .TRUE.
        !
        ! Check the i and j lines for removal
        !
        DO WHILE (check_i .OR. check_j)
            !
            ! Check the i-line for removal
            !
            IF (check_i) THEN
                !
                ! Find the previous non-deleted i-line
                !
                icheck_m1 = icheck - 1
                DO WHILE (.NOT. i_used(icheck_m1) )
                    icheck_m1 = icheck_m1 - 1
                END DO
                !
                ! Find the next non-deleted i-line
                !
                icheck_p1 = icheck + 1
                DO WHILE (.NOT. i_used(icheck_p1) )
                    icheck_p1 = icheck_p1 + 1
                END DO
                !
                ! Set the xp values
                !
                xp(1) = x_fine(icheck_m1)
                xp(2) = x_fine(icheck_p1)
                inum  = (xp(2)-xp(1))/dx
                !
                ! Iterate the j-lines and find pairs of
                ! non-deleted successive lines
                !
                j_p1 = 1
                can_remove_line = .TRUE.

                RowCheck: DO WHILE (can_remove_line .AND. (j_p1 < jmax))
                    j = j_p1
                    !
                    ! If this j-line has been deleted move to the next one
                    !
                    IF (.NOT. j_used(j)) THEN
                        CYCLE
                    END IF
                    !
                    ! Find the next non-deleted j-line
                    !
                    j_p1 = j + 1
                    DO WHILE (.NOT. j_used(j_p1) )
                        j_p1 = j_p1 + 1
                    END DO
                    !
                    ! Check the interpolation errors for the
                    ! interval (icheck-1, j) -> (icheck+1, j+1)
                    !
                    yp(1) = y_fine(j)
                    yp(2) = y_fine(j_p1)
                    jnum  = (yp(2)-yp(1))/dy

                    up(1,1) = val_fine(icheck_m1, j)
                    up(1,2) = val_fine(icheck_m1, j_p1)
                    up(2,1) = val_fine(icheck_p1, j)
                    up(2,2) = val_fine(icheck_p1, j_p1)

                    DO i = 0, inum
                        xk = xp(1) + REAL(i, KIND=REAL_HIGH)*dx
                        DO j = 0, jnum
                            yk = yp(1) + REAL(j, KIND=REAL_HIGH)*dy
                            !
                            ! Compute the accurate value of the property
                            !
                            prop_real = poly(xk,yk)
                            !
                            ! Compute the interpolated value of the property
                            !
                            prop_int = interpolate_bilinear(xp, yp, up, xk, yk)
                            !
                            ! Estimate the error and check it
                            !
                            error_estimate = 100.0D+00*ABS(prop_int-prop_real)/ABS(prop_real)

                            IF (error_estimate > error) THEN
                                can_remove_line = .FALSE.
                                EXIT RowCheck
                            END IF
                        END DO
                    END DO

                END DO RowCheck
                !
                ! If the line can be removed then mark
                ! it as such
                !
                IF (can_remove_line) THEN
                    i_used(icheck) = .FALSE.
                END IF

            END IF
            !
            ! Check the j-line for removal
            !
            IF (check_j) THEN
                !
                ! Find the previous non-deleted j-line
                !
                jcheck_m1 = jcheck - 1
                DO WHILE (.NOT. j_used(jcheck_m1) )
                    jcheck_m1 = jcheck_m1 - 1
                END DO
                !
                ! Find the next non-deleted j-line
                !
                jcheck_p1 = jcheck + 1
                DO WHILE (.NOT. j_used(jcheck_p1) )
                    jcheck_p1 = jcheck_p1 + 1
                END DO
                !
                ! Set the yp values
                !
                yp(1) = y_fine(jcheck_m1)
                yp(2) = y_fine(jcheck_p1)
                jnum  = (yp(2)-yp(1))/dy
                !
                ! Iterate the i-lines and find pairs of
                ! non-deleted successive lines
                !
                i_p1 = 1
                can_remove_line = .TRUE.

                ColumnCheck: DO WHILE (can_remove_line .AND. (i_p1 < imax))
                    i = i_p1
                    !
                    ! If this i-line has been deleted move to the next one
                    !
                    IF (.NOT. i_used(i)) THEN
                        CYCLE
                    END IF
                    !
                    ! Find the next non-deleted i-line
                    !
                    i_p1 = i + 1
                    DO WHILE (.NOT. i_used(i_p1) )
                        i_p1 = i_p1 + 1
                    END DO
                    !
                    ! Check the interpolation errors for the
                    ! interval (i, jcheck-1) -> (i+1, jcheck+1)
                    !
                    xp(1) = x_fine(i)
                    xp(2) = x_fine(i_p1)
                    inum  = (xp(2)-xp(1))/dx

                    up(1,1) = val_fine(i,    jcheck_m1)
                    up(1,2) = val_fine(i,    jcheck_p1)
                    up(2,1) = val_fine(i_p1, jcheck_m1)
                    up(2,2) = val_fine(i_p1, jcheck_p1)

                    DO i = 0, inum
                        xk = xp(1) + REAL(i, KIND=REAL_HIGH)*dx
                        DO j = 0, jnum
                            yk = yp(1) + REAL(j, KIND=REAL_HIGH)*dy
                            !
                            ! Compute the accurate value of the property
                            !
                            prop_real = poly(xk,yk)
                            !
                            ! Compute the interpolated value of the property
                            !
                            prop_int = interpolate_bilinear(xp, yp, up, xk, yk)
                            !
                            ! Estimate the error and check it
                            !
                            error_estimate = 100.0D+00*ABS(prop_int-prop_real)/ABS(prop_real)

                            IF (error_estimate > error) THEN
                                can_remove_line = .FALSE.
                                EXIT ColumnCheck
                            END IF
                        END DO
                    END DO

                END DO ColumnCheck
                !
                ! If the line can be removed then mark
                ! it as such
                !
                IF (can_remove_line) THEN
                    j_used(jcheck) = .FALSE.
                END IF

            END IF
            !
            ! Move on to the next i and j lines to check for removal
            ! that have not been already removed
            !
            check_i = .FALSE.
            DO i = icheck+1, imax-1
                IF (i_used(i)) THEN
                    icheck = i
                    check_i = .TRUE.
                    EXIT
                END IF
            END DO

            check_j = .FALSE.
            DO j = jcheck+1, jmax-1
                IF (j_used(j)) THEN
                    jcheck = j
                    check_j = .TRUE.
                    EXIT
                END IF
            END DO

        END DO

        imax_new   = COUNT(i_used)
        jmax_new   = COUNT(j_used)
        points_new = imax_new*jmax_new

        IF (logging) THEN
            WRITE(*,*)
            WRITE(*,*) "===================================="
            WRITE(*,*) "Table Optimisation Results"
            WRITE(*,*) "===================================="
            WRITE(*,*) "Eliminated ", imax-imax_new, " columns"
            WRITE(*,*) "Eliminated ", jmax-jmax_new, " rows"
            WRITE(*,*) "------------------------------------"
            WRITE(*,*) "Optimised table has columns = ", imax_new
            WRITE(*,*) "Optimised table has rows    = ", jmax_new
            WRITE(*,*) "Total table data points     = ", imax_new*jmax_new
            WRITE(*,*) "Table size reduction        = ", &
                        100.0D+00*REAL(points-points_new, KIND=REAL_HIGH) / &
                                  REAL(points, KIND=REAL_HIGH), " %"
            WRITE(*,*) "Maximum allowed error       = ", error, " %"
            WRITE(*,*) "===================================="
            WRITE(*,*)
        END IF
        !
        ! Construct the table data
        !
        ALLOCATE(x(1:imax_new), y(1:jmax_new), val(1:imax_new, 1:jmax_new))

        i_new = 0
        DO i = 1, imax
            IF (i_used(i)) THEN
                i_new    = i_new + 1
                x(i_new) = x_fine(i)
            END IF
        END DO

        j_new = 0
        DO j = 1, jmax
            IF (j_used(j)) THEN
                j_new    = j_new + 1
                y(j_new) = y_fine(j)
            END IF
        END DO

        i_new = 0
        DO i = 1, imax
            IF (i_used(i)) THEN
                i_new = i_new + 1
                j_new = 0
                DO j = 1, jmax
                    IF (j_used(j)) THEN
                        j_new = j_new + 1
                        val(i_new,j_new) = val_fine(i,j)
                    END IF
                END DO
            END IF
        END DO

        DEALLOCATE(x_fine, y_fine, val_fine)

    END SUBROUTINE generate_table

END MODULE SteamTables
