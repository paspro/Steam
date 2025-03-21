// -------------------------------------------------------------------------------------------------
//
//  This library implements the computation of the steam polynomials according to the revised
//  release on the IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water
//  and Steam, August 2007 (IAPWS-IF97).
//
//  Copyright (c) 2014-2025 by Dr. Panos Asproulis (p.asproulis@icloud.com).
//  All Rights Reserved.
//
// -------------------------------------------------------------------------------------------------

//! Steam Regions Module
//!
//! This module implements functions that are capable of determining if
//! a steam state belongs to a particular region according to the revised
//! release on the IAPWS Industrial Formulation 1997 for the Thermodynamic
//! Properties of Water and Steam, August 2007 (IAPWS-IF97)

use crate::steam_boundaries;
use crate::steam_constants::*;
use crate::steam_region2;
use crate::steam_region3;
use crate::steam_region4;

///
/// This function determines whether the specified steam state resides in
/// region 1 or not.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 1.
///
pub fn in_region1(pressure: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 1
    //
    let tflag = (IAPWS97_TMIN..=REGION_1_TMAX).contains(&temperature);
    let pflag =
        (pressure >= steam_region4::saturation_pressure(temperature)) && (pressure <= IAPWS97_PMAX);

    pflag && tflag
}

///
/// This function determines whether the specified steam state resides in
/// region 2 or not.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 2.
///
pub fn in_region2(pressure: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 2.
    //
    if (IAPWS97_TMIN..=REGION_1_TMAX).contains(&temperature) {
        if (pressure > 0.0) && (pressure <= steam_region4::saturation_pressure(temperature)) {
            return true;
        }
    } else if (temperature > REGION_1_TMAX) && (temperature <= REGION_2_4_T) {
        if (pressure > 0.0) && (pressure <= steam_boundaries::boundary23_pressure(temperature)) {
            return true;
        }
    } else if (temperature > REGION_2_4_T)
        && (temperature <= REGION_2_TMAX)
        && (pressure > 0.0)
        && (pressure <= IAPWS97_PMAX)
    {
        return true;
    }

    false
}

///
/// This function determines whether the specified steam state resides in
/// region 2, sub-region 2a.
///
/// - Arguments:
///   - `pressure` - The steam pressure (MPa).
///   - `temperature` - The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 2, sub-region 2a.
///
pub fn in_region2a(pressure: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 2, sub-region 2a.
    //
    if pressure <= REGION_2A_2B_P {
        return in_region2(pressure, temperature);
    }

    false
}

///
/// This function determines whether the specified steam state resides in
/// region 2, sub-region 2b.
///
/// - Arguments:
///   - `pressure` - The steam pressure (MPa).
///   - `temperature` - The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 2, sub-region 2b.
///
pub fn in_region2b(pressure: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 2, sub-region 2b.
    //
    let h = steam_region2::specific_enthalpy(pressure, temperature);

    if pressure > REGION_2A_2B_P && pressure <= steam_boundaries::boundary2bc_pressure(h) {
        return in_region2(pressure, temperature);
    }

    false
}

///
/// This function determines whether the specified steam state resides in
/// region 2, sub-region 2c.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 2, sub-region 2c.
///
pub fn in_region2c(pressure: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 2, sub-region 2c.
    //
    let h = steam_region2::specific_enthalpy(pressure, temperature);

    if pressure > steam_boundaries::boundary2bc_pressure(h) {
        return in_region2(pressure, temperature);
    }

    false
}

///
/// This function determines whether the specified steam state resides in
/// region 3 or not.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 3.
///
pub fn in_region3(pressure: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 3.
    //
    let tflag = (temperature >= REGION_1_TMAX)
        && (temperature <= steam_boundaries::boundary23_temperature(pressure));

    let pflag = (pressure >= steam_boundaries::boundary23_pressure(temperature))
        && (pressure <= IAPWS97_PMAX);

    pflag && tflag
}

///
/// This function determines whether the specified steam state resides in
/// region 3a or not.
///
/// - Arguments:
///   - `density`: The steam density (kg/m3).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 3a.
///
pub fn in_region3a(density: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 3a.
    //
    let enthalpy = steam_region3::specific_enthalpy(density, temperature);
    let press = steam_region3::pressure(density, temperature);
    let h3ab = steam_boundaries::boundary3ab_enthalpy(press);

    if enthalpy < h3ab {
        return in_region3(press, temperature);
    }

    false
}

///
/// This function determines whether the specified steam state resides in
/// region 3b or not.
///
/// - Arguments:
///   - `density`: The steam density (kg/m3).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 3b.
///
pub fn in_region3b(density: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 3b.
    //
    let enthalpy = steam_region3::specific_enthalpy(density, temperature);
    let press = steam_region3::pressure(density, temperature);
    let h3ab = steam_boundaries::boundary3ab_enthalpy(press);

    if enthalpy >= h3ab {
        return in_region3(press, temperature);
    }

    false
}

///
/// This function determines whether the specified steam state resides in
/// region 4 or not. Since this region is a line one cannot determine
/// exactly whether a steam state lies exactly on this line or not so an
/// error percentage is required for this function to work.
///
/// - Arguments:
///   - `pressure` - The steam pressure (MPa).
///   - `temperature` - The steam temperature (K).
///   - `error` - The acceptable error (%).
///
/// - Returns:
///   - True if the steam state resides in region 4.
///
pub fn in_region4(pressure: f64, temperature: f64, error: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 4.
    //
    if (IAPWS97_TMIN..=REGION_4_TMAX).contains(&temperature) {
        let pressure_estimate = steam_region4::saturation_pressure(temperature);
        let temperature_estimate = steam_region4::saturation_temperature(pressure);

        let p_error = (pressure - pressure_estimate).abs() / pressure;
        let t_error = (temperature - temperature_estimate).abs() / temperature;

        if (p_error <= error) && (t_error <= error) {
            return true;
        }
    }

    false
}

///
/// This function determines whether the specified steam state resides in
/// region 5 or not.
///
/// - Arguments:
///   - `pressure` - The steam pressure (MPa).
///   - `temperature` - The steam temperature (K).
///
/// - Returns:
///   - True if the steam state resides in region 5.
///
pub fn in_region5(pressure: f64, temperature: f64) -> bool {
    //
    // Determine whether the steam state belongs in region 5.
    //
    let tflag = (REGION_2_TMAX..=IAPWS97_TMAX).contains(&temperature);
    let pflag = (pressure > 0.0) && (pressure <= REGION_5_PMAX);

    pflag && tflag
}
