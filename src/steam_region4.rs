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

//! Steam Region 4 Module
//!
//! This module implements the polynomials which compute the thermodynamic
//! properties of steam in region 4 according to the revised release on the
//! IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
//! Water and Steam, August 2007 (IAPWS-IF97)

use crate::steam_constants::*;

//
// Import functions from other regions with aliased names.
//
use crate::steam_region1::{
    specific_enthalpy as h1, specific_entropy as s1, specific_internal_energy as u1,
    specific_isobaric_heat_capacity as cp1, specific_isochoric_heat_capacity as cv1,
    specific_volume as v1,
};

use crate::steam_region2::{
    specific_enthalpy as h2, specific_entropy as s2, specific_internal_energy as u2,
    specific_isobaric_heat_capacity as cp2, specific_isochoric_heat_capacity as cv2,
    specific_volume as v2,
};

use crate::steam_region3::{
    specific_enthalpy as h3, specific_entropy as s3, specific_internal_energy as u3,
    specific_isobaric_heat_capacity as cp3, specific_isochoric_heat_capacity as cv3,
};

///
/// Constant coefficients "N".
///
const N: [f64; 10] = [
    0.11670521452767e4,
    -0.72421316703206e6,
    -0.17073846940092e2,
    0.12020824702470e5,
    -0.32325550322333e7,
    0.14915108613530e2,
    -0.48232657361591e4,
    0.40511340542057e6,
    -0.23855557567849e0,
    0.65017534844798e3,
];

///
/// Constant coefficients "K".
///
const K: [f64; 6] = [
    1.99274064,
    1.09965342,
    -0.510839303,
    -1.75493479,
    -45.5170352,
    -6.74694450e5,
];

///
/// Constant coefficients "L".
///
const L: [f64; 6] = [
    -2.03150240,
    -2.68302940,
    -5.38626492,
    -17.2991605,
    -44.7586581,
    -63.9201063,
];

///
/// Star pressure for region 4 in (MPa).
///
const REGION_4_PSTAR: f64 = 1.0;

///
/// Star temperature for region 4 in (K).
///
const REGION_4_TSTAR: f64 = 1.0;

///
/// This function computes the saturation pressure with respect to
/// temperature.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The pressure (MPa).
///
pub fn saturation_pressure(temperature: f64) -> f64 {
    //
    // Compute the saturation pressure.
    //
    let tau = temperature / REGION_4_TSTAR;
    let theta = tau + N[8] / (tau - N[9]);
    let theta2 = theta * theta;
    let a = theta2 + N[0] * theta + N[1];
    let b = N[2] * theta2 + N[3] * theta + N[4];
    let c = N[5] * theta2 + N[6] * theta + N[7];
    let expr = TWO * c / (-b + (b * b - FOUR * a * c).sqrt());

    REGION_4_PSTAR * expr.powi(4)
}

///
/// This function computes the saturation pressure gradient with respect
/// to temperature.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The pressure gradient (MPa/K).
///
pub fn saturation_pressure_gradient(temperature: f64) -> f64 {
    //
    // Compute the saturation pressure gradient.
    //
    let beta = (saturation_pressure(temperature) / REGION_4_PSTAR).powf(QUARTER);
    let theta = temperature / REGION_4_TSTAR + N[8] / (temperature / REGION_4_TSTAR - N[9]);
    let xbeta = (TWO * beta + N[2]) * theta * theta
        + (TWO * beta * N[0] + N[3]) * theta
        + TWO * N[1] * beta
        + N[4];
    let xtheta = (TWO * theta + N[0]) * beta * beta
        + (TWO * N[2] * theta + N[3]) * beta
        + TWO * N[5] * theta
        + N[6];
    let dthetadt = (ONE - N[8] / (temperature / REGION_4_TSTAR - N[9]).powi(2)) / REGION_4_TSTAR;
    let dbetadtheta = -xtheta / xbeta;
    let dpdbeta = FOUR * beta.powi(3) * REGION_4_PSTAR;

    dpdbeta * dbetadtheta * dthetadt
}

///
/// This function computes the saturation temperature with respect to
/// pressure.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///
/// - Returns:
///   - The temperature (K).
///
pub fn saturation_temperature(pressure: f64) -> f64 {
    //
    // Compute the saturation temperature.
    //
    let pi = pressure / REGION_4_PSTAR;
    let beta = pi.powf(QUARTER);
    let beta2 = beta * beta;
    let e = beta2 + N[2] * beta + N[5];
    let f = N[0] * beta2 + N[3] * beta + N[6];
    let g = N[1] * beta2 + N[4] * beta + N[7];
    let d = TWO * g / (-f - (f * f - FOUR * e * g).sqrt());

    REGION_4_TSTAR * HALF * (N[9] + d - ((N[9] + d).powi(2) - FOUR * (N[8] + N[9] * d)).sqrt())
}

///
/// This function computes the saturation water density with respect to
/// temperature.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The saturation water density (Kg/m3).
///
pub fn saturation_water_density(temperature: f64) -> f64 {
    //
    // Compute the saturation water density.
    //
    let tau = ONE - temperature / IAPWS97_TCRIT;
    let tau_1_3 = tau.powf(ONETHIRD);
    let tau_2_3 = tau_1_3 * tau_1_3;
    let tau_5_3 = tau * tau_2_3;
    let tau_16_3 = tau_5_3 * tau_5_3 * tau_5_3 * tau_1_3;
    let tau_43_3 = tau_16_3 * tau_16_3 * tau_5_3 * tau_5_3 * tau_1_3;
    let tau_110_3 = tau_43_3 * tau_43_3 * tau_16_3 * tau_5_3 * tau;
    let delta = ONE
        + K[0] * tau_1_3
        + K[1] * tau_2_3
        + K[2] * tau_5_3
        + K[3] * tau_16_3
        + K[4] * tau_43_3
        + K[5] * tau_110_3;

    IAPWS97_RHOCRIT * delta
}

///
/// This function computes the saturation steam density with respect to
/// temperature.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The saturation steam density (Kg/m3).
///
pub fn saturation_steam_density(temperature: f64) -> f64 {
    //
    // Compute the saturation steam density.
    //
    let tau = ONE - temperature / IAPWS97_TCRIT;
    let tau_1_6 = tau.powf(ONESIXTH);
    let tau_2_6 = tau_1_6 * tau_1_6;
    let tau_4_6 = tau_2_6 * tau_2_6;
    let tau_8_6 = tau_4_6 * tau_4_6;
    let tau_16_6 = tau_8_6 * tau_8_6;
    let tau_18_6 = tau_16_6 * tau_2_6;
    let tau_37_6 = tau_18_6 * tau_18_6 * tau_1_6;
    let tau_71_6 = tau_37_6 * tau_18_6 * tau_16_6;
    let ln_delta = L[0] * tau_2_6
        + L[1] * tau_4_6
        + L[2] * tau_8_6
        + L[3] * tau_18_6
        + L[4] * tau_37_6
        + L[5] * tau_71_6;

    IAPWS97_RHOCRIT * (ln_delta.exp())
}

///
/// This function computes the specific internal energy with respect to
/// temperature and steam quality.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///   - `quality`: The steam quality (%).
///
/// - Returns:
///   - The specific internal energy (J/Kg).
///
pub fn specific_internal_energy(temperature: f64, quality: f64) -> f64 {
    //
    // Compute the specific internal energy.
    //
    let (uf, ug) = if temperature < REGION_1_TMAX {
        let psat = saturation_pressure(temperature);
        (u1(psat, temperature), u2(psat, temperature))
    } else {
        let rhof = saturation_water_density(temperature);
        let rhog = saturation_steam_density(temperature);
        (u3(rhof, temperature), u3(rhog, temperature))
    };

    uf + quality * (ug - uf)
}

///
/// This function computes the specific volume with respect to
/// temperature and steam quality.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///   - `quality`: The steam quality (%).
///
/// - Returns:
///   - The specific volume (m3/kg).
///
pub fn specific_volume(temperature: f64, quality: f64) -> f64 {
    //
    // Compute the specific volume.
    //
    let (vf, vg) = if temperature < REGION_1_TMAX {
        let psat = saturation_pressure(temperature);
        (v1(psat, temperature), v2(psat, temperature))
    } else {
        (
            ONE / saturation_water_density(temperature),
            ONE / saturation_steam_density(temperature),
        )
    };

    vf + quality * (vg - vf)
}

///
/// This function computes the specific enthalpy with respect to
/// temperature and steam quality.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///   - `quality`: The steam quality (%).
///
/// - Returns:
///   - The specific enthalpy (J/Kg).
///
pub fn specific_enthalpy(temperature: f64, quality: f64) -> f64 {
    //
    // Compute the specific enthalpy.
    //
    let (hf, hg) = if temperature < REGION_1_TMAX {
        let psat = saturation_pressure(temperature);
        (h1(psat, temperature), h2(psat, temperature))
    } else {
        let rhof = saturation_water_density(temperature);
        let rhog = saturation_steam_density(temperature);
        (h3(rhof, temperature), h3(rhog, temperature))
    };

    hf + quality * (hg - hf)
}

///
/// This function computes the specific entropy with respect to
/// temperature and steam quality.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///   - `quality`: The steam quality (%).
///
/// - Returns:
///   - The specific entropy (J/Kg.K).
///
pub fn specific_entropy(temperature: f64, quality: f64) -> f64 {
    //
    // Compute the specific entropy.
    //
    let (sf, sg) = if temperature < REGION_1_TMAX {
        let psat = saturation_pressure(temperature);
        (s1(psat, temperature), s2(psat, temperature))
    } else {
        let rhof = saturation_water_density(temperature);
        let rhog = saturation_steam_density(temperature);
        (s3(rhof, temperature), s3(rhog, temperature))
    };

    sf + quality * (sg - sf)
}

///
/// This function computes the specific isobaric heat capacity with respect
/// to temperature and steam quality.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///   - `quality`: The steam quality (%).
///
/// - Returns:
///   - The specific isobaric heat capacity (J/Kg.K).
///
pub fn specific_isobaric_heat_capacity(temperature: f64, quality: f64) -> f64 {
    //
    // Compute the specific isobaric heat capacity.
    //
    let (cpf, cpg) = if temperature < REGION_1_TMAX {
        let psat = saturation_pressure(temperature);
        (cp1(psat, temperature), cp2(psat, temperature))
    } else {
        let rhof = saturation_water_density(temperature);
        let rhog = saturation_steam_density(temperature);
        (cp3(rhof, temperature), cp3(rhog, temperature))
    };

    cpf + quality * (cpg - cpf)
}

///
/// This function computes the specific isochoric heat capacity with respect
/// to temperature and steam quality.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///   - `quality`: The steam quality (%).
///
/// - Returns:
///   - The specific isochoric heat capacity (J/Kg.K).
///
pub fn specific_isochoric_heat_capacity(temperature: f64, quality: f64) -> f64 {
    //
    // Compute the specific isochoric heat capacity.
    //
    let (cvf, cvg) = if temperature < REGION_1_TMAX {
        let psat = saturation_pressure(temperature);
        (cv1(psat, temperature), cv2(psat, temperature))
    } else {
        let rhof = saturation_water_density(temperature);
        let rhog = saturation_steam_density(temperature);
        (cv3(rhof, temperature), cv3(rhog, temperature))
    };

    cvf + quality * (cvg - cvf)
}

///
/// This function computes the ratio of specific heats with respect
/// to temperature and steam quality.
///
/// - Arguments:
///   - `temperature`: The steam temperature (K).
///   - `quality`: The steam quality (%).
///
/// - Returns:
///   - The ratio of specific heats.
///
pub fn ratio_of_specific_heats(temperature: f64, quality: f64) -> f64 {
    //
    // Compute the ratio of specific heats.
    //
    specific_isobaric_heat_capacity(temperature, quality)
        / specific_isochoric_heat_capacity(temperature, quality)
}
