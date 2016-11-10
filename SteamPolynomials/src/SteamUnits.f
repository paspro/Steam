!-----------------------------------------------------------------------
!
!> @file SteamUnits.f
!
!> @details
!> This module defines the units used for the computation of the steam
!> polynomials according to the revised release on the IAPWS Industrial
!> Formulation 1997 for the Thermodynamic Properties of Water and Steam,
!> August 2007 (IAPWS-IF97)
!
!> @author Panos Asproulis
!> @date 2014
!> @copyright Panos Asproulis (2014-2016). All Rights Reserved.
!
!-----------------------------------------------------------------------
MODULE SteamUnits

    USE SteamPrecision

    IMPLICIT NONE
    !
    !> Unit Multipliers
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: &
        tera  = 1.0D+12, &
        giga  = 1.0D+09, &
        mega  = 1.0D+06, &
        kilo  = 1.0D+03, &
        hecta = 1.0D+02, &
        deca  = 1.0D+01, &
        deci  = 1.0D-01, &
        centi = 1.0D-02, &
        milli = 1.0D-03, &
        micro = 1.0D-06
    !
    !> Base units used by the steam polynomials
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: &
        Kelvin = 1.0D+00, &
        Kg     = 1.0D+00, &
        meter  = 1.0D+00, &
        second = 1.0D+00
    !
    !> Composite units
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: &
        cm     = centi*meter,                 &
        mm     = milli*meter,                 &
        Km     = kilo*meter,                  &
        m2     = meter*meter,                 &
        cm2    = cm*cm,                       &
        mm2    = mm*mm,                       &
        m3     = meter*mm2,                   &
        cm3    = cm*cm2,                      &
        mm3    = mm*mm2,                      &
        minute = 60.0D+00*second,             &
        hour   = 60.0D+00*minute,             &
        day    = 12.0D+00*hour,               &
        week   = 7.0D+00*day,                 &
        Hertz  = 1.0D+00/second,              &
        Newton = Kg*meter/(second*second),    &
        Pascal = Newton/m2,                   &
        bar    = hecta*kilo*Pascal,           &
        KPa    = kilo*Pascal,                 &
        MPa    = mega*Pascal,                 &
        Joule  = Newton*meter,                &
        KJ     = kilo*Joule,                  &
        Watt   = Joule/second,                &
        gr     = milli*Kg,                    &
        J_Kg   = Joule/Kg,                    &
        KJ_Kg  = KJ/Kg,                       &
        KJ_KgK = KJ/Kg/Kelvin,                &
        Kg_m3  = Kg/m3,                       &
        m3_Kg  = m3/Kg
    !
    !> Imperial Units
    !
    REAL(KIND=REAL_HIGH), PARAMETER :: &
        Btu    = 1055.05585262D+00*Joule,     &
        Rankin = 0.556D+00*Kelvin,            &
        RPM    = 1.0D+00/minute,              &
        yard   = 0.9144D+00*meter,            &
        foot   = yard/3.0D+00,                &
        inch   = foot/12.0D+00,               &
        mile   = 1760.0D+00*yard,             &
        lbm    = 0.45359237D+00*Kg

END MODULE SteamUnits
