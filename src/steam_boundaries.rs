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

//! Steam Boundaries Module
//!
//! This module implements the polynomials which compute the thermodynamic
//! properties of steam on various boundary lines between regions according
//! to the revised release on the IAPWS Industrial Formulation 1997 for the
//! Thermodynamic Properties of Water and Steam, August 2007 (IAPWS-IF97)

use crate::steam_constants::*;
use crate::steam_region2;
use crate::steam_region4;

///
/// Constant coefficients "n".
///
const N: [f64; 5] = [
    0.34805185628969e3,
    -0.11671859879975e1,
    0.10192970039326e-2,
    0.57254459862746e3,
    0.13918839778870e2,
];

///
/// Constant coefficients "m".
///
const M: [f64; 5] = [
    0.90584278514723e3,
    -0.67955786399241e0,
    0.12809002730136e-3,
    0.26526571908428e4,
    0.45257578905948e1,
];

///
/// This function computes pressure on the boundary line between regions
/// 2 and 3 with respect to temperature.
///
/// - Arguments:
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The pressure [MPa]
///
pub fn boundary23_pressure(temperature: f64) -> f64 {
    //
    // Star pressure for the computations on the boundary line between
    // regions 2 and 3 in [MPa].
    //
    const REGION_2_3_PSTAR: f64 = 1.0;
    //
    // Star temperature for the computations on the boundary line between
    // regions 2 and 3 in [J/Kg].
    //
    const REGION_2_3_TSTAR: f64 = 1.0;
    //
    // Compute the pressure.
    //
    let t = temperature / REGION_2_3_TSTAR;
    REGION_2_3_PSTAR * (N[0] + (N[1] + N[2] * t) * t)
}

///
/// Checks if the temperature is valid for boundary23_pressure.
///
/// - Arguments:
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - `bool`: True if the temperature is within valid range.
///
pub fn is_valid_boundary23_pressure(temperature: f64) -> bool {
    (temperature >= REGION_1_TMAX) && (temperature <= REGION_2_4_T)
}

///
/// This function computes temperature on the boundary line between regions
/// 2 and 3 with respect to pressure
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///
/// - Returns:
///   - The temperature [K].
///
pub fn boundary23_temperature(pressure: f64) -> f64 {
    //
    // Star pressure for the computations on the boundary line between
    // regions 2 and 3 in [MPa].
    //
    const REGION_2_3_PSTAR: f64 = 1.0;
    //
    // Star temperature for the computations on the boundary line between
    // regions 2 and 3 in [J/Kg].
    //
    const REGION_2_3_TSTAR: f64 = 1.0;
    //
    // Compute the temperature.
    //
    let p = pressure / REGION_2_3_PSTAR;
    REGION_2_3_TSTAR * (N[3] + ((p - N[4]) / N[2]).sqrt())
}

///
/// Checks if the pressure is valid for boundary23_temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///
/// - Returns:
///   - `bool`: True if the pressure is within valid range.
///
pub fn is_valid_boundary23_temperature(pressure: f64) -> bool {
    (pressure >= REGION_1_4_TMAX_P) && (pressure <= IAPWS97_PMAX)
}

///
/// This function computes enthalpy on the boundary line between regions
/// 2a and 2b with respect to the specific entropy.
///
/// - Arguments:
///   - `entropy`: The steam specific entropy [J/Kg.K].
///
/// - Returns:
///   - The enthalpy [J/Kg].
///
pub fn boundary2ab_enthalpy(entropy: f64) -> f64 {
    //
    // Constant coefficients "l"
    //
    const L: [f64; 4] = [
        -0.349898083432139e4,
        0.257560716905876e4,
        -0.421073558227969e3,
        0.276349063799944e2,
    ];
    //
    // Star enthalpy for the computations on the boundary line between
    // sub-regions 2a and 2b in [J/Kg].
    //
    const REGION_2A_2B_HSTAR: f64 = 1000.0;
    //
    // Star entropy for the computations on the boundary line between
    // sub-regions 2a and 2b in [J/Kg.K].
    //
    const REGION_2A_2B_SSTAR: f64 = 1000.0;
    //
    // Compute the enthalpy.
    //
    let sigma = entropy / REGION_2A_2B_SSTAR;
    REGION_2A_2B_HSTAR * (L[0] + (L[1] + (L[2] + L[3] * sigma) * sigma) * sigma)
}

///
/// Checks if the entropy is valid for boundary2ab_enthalpy.
///
/// - Arguments:
///   - `entropy`: The steam specific entropy [J/Kg.K].
///
/// - Returns:
///   - `bool`: True if the entropy is within valid range.
///
pub fn is_valid_boundary2ab_enthalpy(entropy: f64) -> bool {
    let smin = steam_region2::specific_entropy(
        REGION_2A_2B_P,
        steam_region4::saturation_temperature(REGION_2A_2B_P),
    );

    (entropy >= smin) && (entropy <= steam_region2::specific_entropy(REGION_2A_2B_P, REGION_2_TMAX))
}

///
/// This function computes pressure on the boundary line between regions
/// 2b and 2c with respect to the specific enthalpy.
///
/// - Arguments:
///   - `enthalpy`: The steam specific enthalpy [J/Kg].
///
/// - Returns:
///   - The pressure [MPa].
///
pub fn boundary2bc_pressure(enthalpy: f64) -> f64 {
    //
    // Star pressure for the computations on the boundary line between
    // sub-regions 2b and 2c in [MPa].
    //
    const REGION_2B_2C_PSTAR: f64 = 1.0;
    //
    // Star enthalpy for the computations on the boundary line between
    // sub-regions 2b and 2c in [J/Kg].
    //
    const REGION_2B_2C_HSTAR: f64 = 1000.0;
    //
    // Compute the pressure.
    //
    let e = enthalpy / REGION_2B_2C_HSTAR;
    REGION_2B_2C_PSTAR * (M[0] + (M[1] + M[2] * e) * e)
}

///
/// Checks if the enthalpy is valid for boundary2bc_pressure.
///
/// - Arguments:
///   - `enthalpy`: The steam specific enthalpy [J/Kg].
///
/// - Returns:
///   - `bool`: True if the enthalpy is within valid range.
///
pub fn is_valid_boundary2bc_pressure(enthalpy: f64) -> bool {
    (enthalpy >= boundary2bc_enthalpy(REGION_2B_2C_PMIN))
        && (enthalpy <= boundary2bc_enthalpy(IAPWS97_PMAX))
}

///
/// This function computes enthalpy on the boundary line between regions
/// 2b and 2c with respect to pressure.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///
/// - Returns:
///   - The enthalpy [J/Kg].
///
pub fn boundary2bc_enthalpy(pressure: f64) -> f64 {
    //
    // Star pressure for the computations on the boundary line between
    // sub-regions 2b and 2c in [MPa].
    //
    const REGION_2B_2C_PSTAR: f64 = 1.0;
    //
    // Star enthalpy for the computations on the boundary line between
    // sub-regions 2b and 2c in [J/Kg].
    //
    const REGION_2B_2C_HSTAR: f64 = 1000.0;
    //
    // Compute the enthalpy.
    //
    let p = pressure / REGION_2B_2C_PSTAR;
    REGION_2B_2C_HSTAR * (M[3] + ((p - M[4]) / M[2]).sqrt())
}

///
/// Checks if the pressure is valid for boundary2bc_enthalpy.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///
/// - Returns:
///   - `bool`: True if the pressure is within valid range.
///
pub fn is_valid_boundary2bc_enthalpy(pressure: f64) -> bool {
    (pressure >= REGION_2B_2C_PMIN) && (pressure <= IAPWS97_PMAX)
}

///
/// This function computes enthalpy on the boundary line between regions
/// 3a and 3b with respect to pressure.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///
/// - Returns:
///   - The enthalpy [J/Kg].
///
pub fn boundary3ab_enthalpy(pressure: f64) -> f64 {
    //
    // Constant coefficients "l".
    //
    const L: [f64; 4] = [
        0.201464004206875e4,
        0.374696550136983e1,
        -0.219921901054187e-1,
        0.875131686009950e-4,
    ];
    //
    // Star pressure for the computations on the boundary line between
    // sub-regions 3a and 3b in [MPa].
    //
    const REGION_3A_3B_PSTAR: f64 = 1.0;
    //
    // Star enthalpy for the computations on the boundary line between
    // sub-regions 3a and 3b in [J/Kg].
    //
    const REGION_3A_3B_HSTAR: f64 = 1000.0;
    //
    // Compute the enthalpy.
    //
    let p = pressure / REGION_3A_3B_PSTAR;
    REGION_3A_3B_HSTAR * (L[0] + p * (L[1] + p * (L[2] + L[3] * p)))
}

///
/// Checks if the pressure is valid for boundary3ab_enthalpy.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///
/// - Returns:
///   - `bool`: True if the pressure is within valid range.
///
pub fn is_valid_boundary3ab_enthalpy(pressure: f64) -> bool {
    (pressure >= IAPWS97_PCRIT) && (pressure <= IAPWS97_PMAX)
}

///
/// This function computes pressure on the boundary line between regions
/// 3 and 4 with respect to enthalpy.
///
/// - Arguments:
///   - `enthalpy`: The steam enthalpy [J/Kg].
///
/// - Returns:
///   - The pressure [MPa].
///
pub fn boundary34_pressure_h(enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [usize; 14] = [0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36];
    const JJ: [usize; 14] = [0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24];
    const NN: [f64; 14] = [
        0.600073641753024e0,
        -0.936203654849857e1,
        0.246590798594147e2,
        -0.107014222858224e3,
        -0.915821315805768e14,
        -0.862332011700662e4,
        -0.235837344740032e2,
        0.252304969384128e18,
        -0.389718771997719e19,
        -0.333775713645296e23,
        0.356499469636328e11,
        -0.148547544720641e27,
        0.330611514838798e19,
        0.813641294467829e38,
    ];
    //
    // Star pressure for the computations on the boundary line between
    // regions 3 and 4 in [MPa].
    //
    const REGION_3_4_PSTAR: f64 = 22.0;
    //
    // Star enthalpy for the computations on the boundary line between
    // regions 3 and 4 in [J/Kg].
    //
    const REGION_3_4_HSTAR: f64 = 2600.0e3;
    //
    // Compute the pressure.
    //
    let eta = enthalpy / REGION_3_4_HSTAR;
    let eta1 = eta - 1.020;
    let eta2 = eta - 0.608;

    let mut sum = 0.0;
    for i in 0..14 {
        sum += NN[i] * eta1.powi(II[i] as i32) * eta2.powi(JJ[i] as i32);
    }

    REGION_3_4_PSTAR * sum
}

///
/// Checks if the enthalpy is valid for boundary34_pressure_h.
///
/// - Arguments:
///   - `enthalpy`: The steam enthalpy [J/Kg].
///
/// - Returns:
///   - `bool`: True if the enthalpy is within valid range.
///
pub fn is_valid_boundary34_pressure_h(enthalpy: f64) -> bool {
    (enthalpy >= REGION_3_4_HMIN) && (enthalpy <= REGION_3_4_HMAX)
}

///
/// This function computes pressure on the boundary line between regions
/// 3 and 4 with respect to entropy.
///
/// - Arguments:
///   - `entropy`: The steam entropy [J/Kg.K].
///
/// - Returns:
///   - The pressure [MPa].
///
pub fn boundary34_pressure_s(entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [usize; 10] = [0, 1, 1, 4, 12, 12, 16, 24, 28, 32];
    const JJ: [usize; 10] = [0, 1, 32, 7, 4, 14, 36, 10, 0, 18];
    const NN: [f64; 10] = [
        0.639767553612785e0,
        -0.129727445396014e2,
        -0.224595125848403e16,
        0.177466741801846e7,
        0.717079349571538e10,
        -0.378829107169011e18,
        -0.955586736431328e35,
        0.187269814676188e24,
        0.119254746466473e12,
        0.110649277244882e37,
    ];
    //
    // Star pressure for the computations on the boundary line between
    // regions 3 and 4 in [MPa].
    //
    const REGION_3_4_PSTAR: f64 = 22.0;
    //
    // Star entropy for the computations on the boundary line between
    // regions 3 and 4 in [J/Kg.K].
    //
    const REGION_3_4_SSTAR: f64 = 5.2e3;
    //
    // Compute the pressure.
    //
    let sigma = entropy / REGION_3_4_SSTAR;
    let sigma1 = sigma - 1.030;
    let sigma2 = sigma - 0.699;

    let mut sum = 0.0;
    for i in 0..10 {
        sum += NN[i] * sigma1.powi(II[i] as i32) * sigma2.powi(JJ[i] as i32);
    }

    REGION_3_4_PSTAR * sum
}

///
/// Checks if the entropy is valid for boundary34_pressure_s.
///
/// - Arguments:
///   - `entropy`: The steam entropy [J/Kg.K].
///
/// - Returns:
///   - `bool`: True if the entropy is within valid range.
///
pub fn is_valid_boundary34_pressure_s(entropy: f64) -> bool {
    (entropy >= REGION_3_4_SMIN) && (entropy <= REGION_3_4_SMAX)
}
