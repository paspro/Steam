!-----------------------------------------------------------------------
!
!> @file SteamConductivity.f
!
!> @details
!> This module implements the polynomials which compute the thermal
!> conductivity of steam according to the release on the IAPWS Formulation 
!> 2011 for the Thermal Conductivity of Ordinary Water Substance
!> (September 2011)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamConductivity

    USE SteamConstants

    IMPLICIT NONE
    PRIVATE
    !
    ! Thermal Conductivity
    !
    PUBLIC thermal_conductivity

CONTAINS
    !-------------------------------------------------------------------------
    !
    !> This function computes the thermal conductivity of steam as a function
    !> of density and temperature.
    !
    !> @param[in] density     The steam density (Kg/m3)
    !> @param[in] temperature The steam temperature (K)
    !
    !> @return The thermal conductivity in [W/K.m]
    !
    !-------------------------------------------------------------------------
    ELEMENTAL FUNCTION thermal_conductivity(density, temperature)

        IMPLICIT NONE
        !
        !> Star temperature for the computation of thermal 
        !> conductivity(density,temperature) in [K]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: THCOND_TSTAR = 647.26D+00
        !
        !> Star density for the computation of thermal 
        !> conductivity(density,temperature) in [Kg/m3]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: THCOND_RHOSTAR = 317.7D+00
        !
        !> Star thermal conductivity for the computation of thermal 
        !> conductivity(density,temperature) in [W/K.m]
        !
        REAL(KIND=REAL_HIGH), PARAMETER :: THCOND_KSTAR = 1.0D+00
        !
        !> Constant coefficients "a"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:4) :: a = &
            (/ 0.0102811D+00, 0.0299621D+00, 0.0156146D+00, -0.00422464D+00 /)
        !
        !> Constant coefficients "b"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:3) :: b = &
            (/ -0.397070D+00, 0.400302D+00, 1.060000D+00 /)
        !
        !> Constant coefficients "bb"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:2) :: bb = &
            (/ -0.171587D+00, 2.392190D+00 /)            
        !
        !> Constant coefficients "c"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:6) :: c = &
            (/ 0.64285700D+00, -4.1171700D+00, -6.17937D+00,        &
               0.00308976D+00,  0.0822994D+00,  10.0932D+00 /)
        !
        !> Constant coefficients "d"
        !
        REAL(KIND=REAL_HIGH), PARAMETER, DIMENSION(1:4) :: d = &
            (/ 0.0701309D+00, 0.0118520D+00, 0.00169937D+00, -1.0200D+00 /)
        !
        ! Arguments
        !
        REAL(KIND=REAL_HIGH), INTENT(IN) :: density, temperature
        REAL(KIND=REAL_HIGH) :: thermal_conductivity
        !
        ! Variables
        !
        INTEGER(KIND=INT_HIGH) :: k
        REAL(KIND=REAL_HIGH) :: rhobar, tbar, troot, tpow, lam, tmp
        REAL(KIND=REAL_HIGH) :: dtbar, dtbarpow, q, s, rhobar18, rhobarq
        !
        ! Compute the heat conductivity
        !
        tbar = temperature / THCOND_TSTAR
        rhobar = density / THCOND_RHOSTAR

        troot = sqrt(tbar)
        tpow  = troot
        lam   = zero

        DO k = 1, 4
            lam  = lam + a(k)*tpow
            tpow = tpow*tbar
        END DO

        tmp = (rhobar+bb(2))*(rhobar+bb(2))
        lam = lam + b(1) + b(2)*rhobar + b(3)*EXP(bb(1)*tmp)

        dtbar = ABS(tbar-one) + c(4)
        dtbarpow = dtbar**(three/five)
        q = two + c(5)/dtbarpow

        IF (tbar >= one) THEN
            s = one/dtbar
        ELSE
            s = c(6)/dtbarpow
        END IF

        rhobar18 = rhobar**1.8D+00
        rhobarq  = rhobar**q

        lam = lam + (d(1)/(tbar**10) + d(2))*rhobar18*         &
                    EXP(c(1)*(one-rhobar*rhobar18))            &
                  + d(3)*s*rhobarq*                            &
                    EXP((q/(one+q))*(one-rhobar*rhobarq))      &
                  + d(4)*EXP(c(2)*(troot**3)+c(3)/(rhobar**5))

        thermal_conductivity = THCOND_KSTAR*lam

    END FUNCTION thermal_conductivity

END MODULE SteamConductivity    
