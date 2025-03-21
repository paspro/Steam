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

//! Steam Region 3 Module
//!
//! This module implements the polynomials which compute the thermodynamic
//! properties of steam in region 3 according to the revised release on the
//! IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
//! Water and Steam, August 2007 (IAPWS-IF97)

use crate::steam_constants::*;
use std::f64;

///
/// Constant coefficients "I".
///
const I: [i32; 39] = [
    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6,
    7, 8, 9, 9, 10, 10, 11,
];

///
/// Constant coefficients "J".
///
const J: [i32; 39] = [
    0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3,
    26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26,
];

///
/// Constant coefficients "n".
///
const N: [f64; 39] = [
    -0.15732845290239e2,
    0.20944396974307e2,
    -0.76867707878716e1,
    0.26185947787954e1,
    -0.28080781148620e1,
    0.12053369696517e1,
    -0.84566812812502e-2,
    -0.12654315477714e1,
    -0.11524407806681e1,
    0.88521043984318e0,
    -0.64207765181607e0,
    0.38493460186671e0,
    -0.85214708824206e0,
    0.48972281541877e1,
    -0.30502617256965e1,
    0.39420536879154e-1,
    0.12558408424308e0,
    -0.27999329698710e0,
    0.13899799569460e1,
    -0.20189915023570e1,
    -0.82147637173963e-2,
    -0.47596035734923e0,
    0.43984074473500e-1,
    -0.44476435428739e0,
    0.90572070719733e0,
    0.70522450087967e0,
    0.10770512626332e0,
    -0.32913623258954e0,
    -0.50871062041158e0,
    -0.22175400873096e-1,
    0.94260751665092e-1,
    0.16436278447961e0,
    -0.13503372241348e-1,
    -0.14834345352472e-1,
    0.57922953628084e-3,
    0.32308904703711e-2,
    0.80964802996215e-4,
    -0.16557679795037e-3,
    -0.44923899061815e-4,
];

///
/// Constant coefficient "n1".
///
const N1: f64 = 0.10658070028513e1;

///
/// Star density for region 3 in [MPa].
///
const REGION_3_RHOSTAR: f64 = IAPWS97_RHOCRIT;

///
/// Star temperature for region 3 in [K].
///
const REGION_3_TSTAR: f64 = IAPWS97_TCRIT;

///
/// This function computes the dimensionless specific Helmholtz free energy.
///
/// - Arguments:
///   - `delta`: Dimensionless density parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// # Returns:
///   - The dimensionless specific Helmholtz free energy.
///
fn phi(delta: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let mut sum = 0.0;
    for i in 0..39 {
        sum += N[i] * delta.powi(I[i]) * tau.powi(J[i]);
    }

    N1 * delta.ln() + sum
}

///
/// This function computes the derivative of the dimensionless specific
/// Helmholtz free energy with respect to the parameter delta.
///
/// - Arguments:
///   - `delta`: Dimensionless density parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d(phi)/d(delta).
///
fn phidelta(delta: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let mut sum = 0.0;
    for i in 0..39 {
        if I[i] > 0 {
            sum += N[i] * I[i] as f64 * delta.powi(I[i] - 1) * tau.powi(J[i]);
        }
    }

    N1 / delta + sum
}

///
/// This function computes the second derivative of the dimensionless specific
/// Helmholtz free energy with respect to the parameter delta.
///
/// - Arguments:
///   - `delta`: Dimensionless density parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(phi)/d(delta)2.
///
fn phideltadelta(delta: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let mut sum = 0.0;
    for i in 0..39 {
        if I[i] > 1 {
            sum += N[i] * I[i] as f64 * (I[i] - 1) as f64 * delta.powi(I[i] - 2) * tau.powi(J[i]);
        }
    }

    -N1 / (delta * delta) + sum
}

///
/// This function computes the derivative of the dimensionless specific
/// Helmholtz free energy with respect to the parameter tau.
///
/// - Arguments:
///   - `delta`: Dimensionless density parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d(phi)/d(tau).
///
fn phitau(delta: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let mut sum = 0.0;
    for i in 0..39 {
        if J[i] != 0 {
            sum += N[i] * delta.powi(I[i]) * J[i] as f64 * tau.powi(J[i] - 1);
        }
    }

    sum
}

///
/// This function computes the second derivative of the dimensionless specific
/// Helmholtz free energy with respect to the parameter tau.
///
/// - Arguments:
///   - `delta`: Dimensionless density parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(phi)/d(tau)2.
///
fn phitautau(delta: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let mut sum = 0.0;
    for i in 0..39 {
        if J[i] > 1 || J[i] < 0 {
            sum += N[i] * delta.powi(I[i]) * J[i] as f64 * (J[i] - 1) as f64 * tau.powi(J[i] - 2);
        }
    }

    sum
}

///
/// This function computes the derivative of the dimensionless specific
/// Helmholtz free energy with respect to the parameter delta and tau.
///
/// - Arguments:
///   - `delta`: Dimensionless density parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(phi)/d(delta)d(tau).
///
fn phideltatau(delta: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let mut sum = 0.0;
    for i in 0..39 {
        if I[i] > 0 && J[i] != 0 {
            sum += N[i] * I[i] as f64 * delta.powi(I[i] - 1) * J[i] as f64 * tau.powi(J[i] - 1);
        }
    }

    sum
}

///
/// This function computes pressure with respect to density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - Pressure in [MPa].
///
pub fn pressure(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;

    //
    // Compute pressure.
    //
    IAPWS97_R * density * temperature * delta * phidelta(delta, tau) * 1.0e-6
}

///
/// This function computes the specific internal energy with respect to
/// density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific internal energy [J/Kg].
/// 
pub fn specific_internal_energy(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;
    //
    // Compute the specific internal energy.
    //
    IAPWS97_R * temperature * tau * phitau(delta, tau)
}

///
/// This function computes the specific entropy with respect to
/// density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific entropy [J/Kg.K].
/// 
pub fn specific_entropy(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;
    //
    // Compute the specific entropy.
    //
    IAPWS97_R * (tau * phitau(delta, tau) - phi(delta, tau))
}

///
/// This function computes the specific enthalpy with respect to
/// density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific enthalpy [J/Kg].
/// 
pub fn specific_enthalpy(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;
    //
    // Compute the specific enthalpy.
    //
    IAPWS97_R * temperature * (tau * phitau(delta, tau) + delta * phidelta(delta, tau))
}

///
/// This function computes the specific isobaric heat capacity with respect
/// to density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific isobaric heat capacity [J/Kg.K].
/// 
pub fn specific_isobaric_heat_capacity(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;
    //
    // Compute the specific isobaric heat capacity.
    //
    let res1 = -tau * tau * phitautau(delta, tau);
    let res2 = delta * (phidelta(delta, tau) - tau * phideltatau(delta, tau)).powi(2);
    let res3 = 2.0 * phidelta(delta, tau) + delta * phideltadelta(delta, tau);

    IAPWS97_R * (res1 + res2 / res3)
}

///
/// This function computes the specific isochoric heat capacity with respect
/// to density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific isochoric heat capacity [J/Kg.K].
/// 
pub fn specific_isochoric_heat_capacity(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;
    //
    // Compute the specific isochoric heat capacity.
    //
    -IAPWS97_R * tau * tau * phitautau(delta, tau)
}

///
/// This function computes the ratio of specific heats with respect
/// to density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The ratio of specific heats.
/// 
pub fn ratio_of_specific_heats(density: f64, temperature: f64) -> f64 {
    //
    // Compute the ratio of specific heats.
    //
    specific_isobaric_heat_capacity(density, temperature)
        / specific_isochoric_heat_capacity(density, temperature)
}

///
/// This function computes the speed of sound with respect
/// to density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The speed of sound [m/sec].
/// 
pub fn speed_of_sound(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;
    //
    // Compute the speed of sound.
    //
    let res1 = 2.0 * delta * phidelta(delta, tau) + delta * delta * phideltadelta(delta, tau);
    let res2 = (delta * phidelta(delta, tau) - delta * tau * phideltatau(delta, tau)).powi(2);
    let res3 = tau * tau * phitautau(delta, tau);

    (IAPWS97_R * temperature * (res1 - res2 / res3)).sqrt()
}

///
/// This function computes the specific Helmoltz free energy with respect
/// to density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific Helmoltz free energy [J/kg].
/// 
pub fn specific_helmoltz_free_energy(density: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters delta and tau.
    //
    let delta = density / REGION_3_RHOSTAR;
    let tau = REGION_3_TSTAR / temperature;
    //
    // Compute the specific Helmoltz free energy.
    //
    IAPWS97_R * temperature * phi(delta, tau)
}

///
/// This function computes the temperature with respect to pressure
/// and enthalpy (backward equation) for the sub-region 3a.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `enthalpy`: The steam enthalpy [J/Kg].
///
/// - Returns:
///   - The temperature [K].
/// 
pub fn temperature_ph_region3a(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 31] = [
        -12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -5, -3, -2, -2, -2,
        -1, -1, 0, 0, 1, 3, 3, 4, 4, 10, 12,
    ];

    const JJ: [i32; 31] = [
        0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, 0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1, 0,
        3, 4, 5,
    ];

    const NN: [f64; 31] = [
        -0.133645667811215e-6,
        0.455912656802978e-5,
        -0.146294640700979e-4,
        0.639341312970080e-2,
        0.372783927268847e3,
        -0.718654377460447e4,
        0.573494752103400e6,
        -0.267569329111439e7,
        -0.334066283302614e-4,
        -0.245479214069597e-1,
        0.478087847764996e2,
        0.764664131818904e-5,
        0.128350627676972e-2,
        0.171219081377331e-1,
        -0.851007304583213e1,
        -0.136513461629781e-1,
        -0.384460997596657e-5,
        0.337423807911655e-2,
        -0.551624873066791e0,
        0.729202277107470e0,
        -0.992522757376041e-2,
        -0.119308831407288e0,
        0.793929190615421e0,
        0.454270731799386e0,
        0.209998591259910e0,
        -0.642109823904738e-2,
        -0.235155868604540e-1,
        0.252233108341612e-2,
        -0.764885133368119e-2,
        0.136176427574291e-1,
        -0.133027883575669e-1,
    ];
    //
    // Star pressure for the backward temperature(pressure, enthalpy)
    // equation in sub-region 3a in [MPa].
    //
    const REGION_3A_PH_PSTAR: f64 = 100.0;
    //
    // Star enthalpy for the backward temperature(pressure, enthalpy)
    // equation in sub-region 3a in [J/Kg].
    //
    const REGION_3A_PH_HSTAR: f64 = 2300.0e3;
    //
    // Star temperature for the backward temperature(pressure, enthalpy)
    // equation in sub-region 3a in [J/Kg].
    //
    const REGION_3A_PH_TSTAR: f64 = 760.0;
    //
    // Constants.
    //
    const A: f64 = 0.240;
    const B: f64 = 0.615;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_3A_PH_PSTAR + A;
    let eta = enthalpy / REGION_3A_PH_HSTAR - B;

    let mut sum = 0.0;
    for i in 0..31 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3A_PH_TSTAR * sum
}

///
/// This function computes the temperature with respect to pressure
/// and enthalpy (backward equation) for the sub-region 3b.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `enthalpy`: The steam enthalpy [J/Kg].
///
/// - Returns:
///   - The temperature [K].
/// 
pub fn temperature_ph_region3b(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 33] = [
        -12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1,
        -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8,
    ];

    const JJ: [i32; 33] = [
        0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, 5, 0, 4, 2, 4, 6, 10, 14, 16, 0, 2,
        1, 1, 1, 1, 1,
    ];

    const NN: [f64; 33] = [
        0.323254573644920e-4,
        -0.127575556587181e-3,
        -0.475851877356068e-3,
        0.156183014181602e-2,
        0.105724860113781e0,
        -0.858514221132534e2,
        0.724140095480911e3,
        0.296475810273257e-2,
        -0.592721983365988e-2,
        -0.126305422818666e-1,
        -0.115716196364853e0,
        0.849000969739595e2,
        -0.108602260086615e-1,
        0.154304475328851e-1,
        0.750455441524466e-1,
        0.252520973612982e-1,
        -0.602507901232996e-1,
        -0.307622221350501e1,
        -0.574011959864879e-1,
        0.503471360939849e1,
        -0.925081888584834e0,
        0.391733882917546e1,
        -0.773146007130190e2,
        0.949308762098587e4,
        -0.141043719679409e7,
        0.849166230819026e7,
        0.861095729446704e0,
        0.323346442811720e0,
        0.873281936020439e0,
        -0.436653048526683e0,
        0.286596714529479e0,
        -0.131778331276228e0,
        0.676682064330275e-2,
    ];
    //
    // Star pressure for the backward temperature(pressure, enthalpy)
    // equation in sub-region 3b in [MPa].
    //
    const REGION_3B_PH_PSTAR: f64 = 100.0;
    //
    // Star enthalpy for the backward temperature(pressure, enthalpy)
    // equation in sub-region 3b in [J/Kg].
    //
    const REGION_3B_PH_HSTAR: f64 = 2800.0e3;
    //
    // Star temperature for the backward temperature(pressure, enthalpy)
    // equation in sub-region 3b in [J/Kg].
    //
    const REGION_3B_PH_TSTAR: f64 = 860.0;
    //
    // Constants.
    //
    const A: f64 = 0.298;
    const B: f64 = 0.720;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_3B_PH_PSTAR + A;
    let eta = enthalpy / REGION_3B_PH_HSTAR - B;

    let mut sum = 0.0;
    for i in 0..33 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3B_PH_TSTAR * sum
}

///
/// This function computes the specific volume with respect to pressure
/// and enthalpy (backward equation) for the sub-region 3a.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `enthalpy`: The steam enthalpy [J/Kg].
///
/// - Returns:
///   - The specific volume [m3/Kg].
/// 
pub fn specific_volume_ph_region3a(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 32] = [
        -12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1,
        0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8,
    ];

    const JJ: [i32; 32] = [
        6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, 16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0, 2,
        0, 2, 2, 2,
    ];

    const NN: [f64; 32] = [
        0.529944062966028e-2,
        -0.170099690234461e0,
        0.111323814312927e2,
        -0.217898123145125e4,
        -0.506061827980875e-3,
        0.556495239685324e0,
        -0.943672726094016e1,
        -0.297856807561527e0,
        0.939353943717186e2,
        0.192944939465981e-1,
        0.421740664704763e0,
        -0.368914126282330e7,
        -0.737566847600639e-2,
        -0.354753242424366e0,
        -0.199768169338727e1,
        0.115456297059049e1,
        0.568366875815960e4,
        0.808169540124668e-2,
        0.172416341519307e0,
        0.104270175292927e1,
        -0.297691372792847e0,
        0.560394465163593e0,
        0.275234661176914e0,
        -0.148347894866012e0,
        -0.651142513478515e-1,
        -0.292468715386302e1,
        0.664876096952665e-1,
        0.352335014263844e1,
        -0.146340792313332e-1,
        -0.224503486668184e1,
        0.110533464706142e1,
        -0.408757344495612e-1,
    ];
    //
    // Star pressure for the backward specific volume(pressure, enthalpy)
    // equation in sub-region 3a in [MPa].
    //
    const REGION_3A_PH_PSTAR: f64 = 100.0;
    //
    // Star enthalpy for the backward specific volume(pressure, enthalpy)
    // equation in sub-region 3a in [J/Kg].
    //
    const REGION_3A_PH_HSTAR: f64 = 2100.0e3;
    //
    // Star specific volume for the backward specific volume(pressure, enthalpy)
    // equation in sub-region 3a in [m3/Kg].
    //
    const REGION_3A_PH_VSTAR: f64 = 0.0028;
    //
    // Constants.
    //
    const A: f64 = 0.128;
    const B: f64 = 0.727;
    //
    // Compute the specific volume.
    //
    let pi = pressure / REGION_3A_PH_PSTAR + A;
    let eta = enthalpy / REGION_3A_PH_HSTAR - B;

    let mut sum = 0.0;
    for i in 0..32 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3A_PH_VSTAR * sum
}

///
/// This function computes the specific volume with respect to pressure
/// and enthalpy (backward equation) for the sub-region 3b.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `enthalpy`: The steam enthalpy [J/Kg].
///
/// - Returns:
///   - The specific volume [m3/Kg].
/// 
pub fn specific_volume_ph_region3b(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 30] = [
        -12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2, -1,
        -1, -1, -1, 0, 1, 1, 2, 2,
    ];

    const JJ: [i32; 30] = [
        0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6,
    ];

    const NN: [f64; 30] = [
        -0.225196934336318e-8,
        0.140674363313486e-7,
        0.233784085280560e-5,
        -0.331833715229001e-4,
        0.107956778514318e-2,
        -0.271382067378863e0,
        0.107202262490333e1,
        -0.853821329075382e0,
        -0.215214194340526e-4,
        0.769656088222730e-3,
        -0.431136580433864e-2,
        0.453342167309331e0,
        -0.507749535873652e0,
        -0.100475154528389e3,
        -0.219201924648793e0,
        -0.321087965668917e1,
        0.607567815637771e3,
        0.557686450685932e-3,
        0.187499040029550e0,
        0.905368030448107e-2,
        0.285417173048685e0,
        0.329924030996098e-1,
        0.239897419685483e0,
        0.482754995951394e1,
        -0.118035753702231e2,
        0.169490044091791e0,
        -0.179967222507787e-1,
        0.371810116332674e-1,
        -0.536288335065096e-1,
        0.160697101092520e1,
    ];
    //
    // Star pressure for the backward specific volume(pressure, enthalpy)
    // equation in sub-region 3b in [MPa].
    //
    const REGION_3B_PH_PSTAR: f64 = 100.0;
    //
    // Star enthalpy for the backward specific volume(pressure, enthalpy)
    // equation in sub-region 3b in [J/Kg].
    //
    const REGION_3B_PH_HSTAR: f64 = 2800.0e3;
    //
    // Star specific volume for the backward specific volume(pressure, enthalpy)
    // equation in sub-region 3b in [m3/Kg].
    //
    const REGION_3B_PH_VSTAR: f64 = 0.0088;
    //
    // Constants.
    //
    const A: f64 = 0.0661;
    const B: f64 = 0.720;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_3B_PH_PSTAR + A;
    let eta = enthalpy / REGION_3B_PH_HSTAR - B;

    let mut sum = 0.0;
    for i in 0..30 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3B_PH_VSTAR * sum
}

///
/// This function computes the temperature with respect to pressure
/// and entropy (backward equation) for the sub-region 3a.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `entropy`: The steam entropy [J/Kg.K].
///
/// - Returns:
///   - The temperature [K].
/// 
pub fn temperature_ps_region3a(pressure: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 33] = [
        -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -6, -6, -5, -5, -5, -4, -4, -4, -2, -2,
        -1, -1, 0, 0, 0, 1, 2, 2, 3, 8, 8, 10,
    ];

    const JJ: [i32; 33] = [
        28, 32, 4, 10, 12, 14, 5, 7, 8, 28, 2, 6, 32, 0, 14, 32, 6, 10, 36, 1, 4, 1, 6, 0, 1, 4, 0,
        0, 3, 2, 0, 1, 2,
    ];

    const NN: [f64; 33] = [
        0.150042008263875e10,
        -0.159397258480424e12,
        0.502181140217975e-3,
        -0.672057767855466e2,
        0.145058545404456e4,
        -0.823889534888890e4,
        -0.154852214233853e0,
        0.112305046746695e2,
        -0.297000213482822e2,
        0.438565132635495e-2,
        0.137837838635464e-2,
        -0.297478527157462e1,
        0.971777947349413e13,
        -0.571527767052398e-4,
        0.288307949778420e5,
        -0.744428289262703e14,
        0.128017324848921e2,
        -0.368275545889071e3,
        0.664768904779177e16,
        0.449359251958880e-1,
        -0.422897836099655e1,
        -0.240614376434179e0,
        -0.474341365254924e1,
        0.724093999126110e0,
        0.923874349695897e0,
        0.399043655281015e1,
        0.384066651868009e-1,
        -0.359344365571848e-2,
        -0.735196448821653e0,
        0.188367048396131e0,
        0.141064266818704e-3,
        -0.257418501496337e-2,
        0.123220024851555e-2,
    ];
    //
    // Star pressure for the backward temperature(pressure, entropy)
    // equation in sub-region 3a in [MPa].
    //
    const REGION_3A_PS_PSTAR: f64 = 100.0;
    //
    // Star entropy for the backward temperature(pressure, entropy)
    // equation in sub-region 3a in [J/Kg].
    //
    const REGION_3A_PS_SSTAR: f64 = 4.4e3;
    //
    // Star temperature for the backward temperature(pressure, entropy)
    // equation in sub-region 3a in [J/Kg].
    //
    const REGION_3A_PS_TSTAR: f64 = 760.0;
    //
    // Constants.
    //
    const A: f64 = 0.240;
    const B: f64 = 0.703;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_3A_PS_PSTAR + A;
    let eta = entropy / REGION_3A_PS_SSTAR - B;

    let mut sum = 0.0;
    for i in 0..33 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3A_PS_TSTAR * sum
}

///
/// This function computes the temperature with respect to pressure
/// and entropy (backward equation) for the sub-region 3b.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `entropy`: The steam entropy [J/Kg.K].
///
/// - Returns:
///   - The temperature [K].
/// 
pub fn temperature_ps_region3b(pressure: f64, entropy: f64) -> f64 {
    // Constant polynomial coefficients
    const II: [i32; 28] = [
        -12, -12, -12, -12, -8, -8, -8, -6, -6, -6, -5, -5, -5, -5, -5, -4, -3, -3, -2, 0, 2, 3, 4,
        5, 6, 8, 12, 14,
    ];

    const JJ: [i32; 28] = [
        1, 3, 4, 7, 0, 1, 3, 0, 2, 4, 0, 1, 2, 4, 6, 12, 1, 6, 2, 0, 1, 1, 0, 24, 0, 3, 1, 2,
    ];

    const NN: [f64; 28] = [
        0.527111701601660e0,
        -0.401317830052742e2,
        0.153020073134484e3,
        -0.224799398218827e4,
        -0.193993484669048e0,
        -0.140467557893768e1,
        0.426799878114024e2,
        0.752810643416743e0,
        0.226657238616417e2,
        -0.622873556909932e3,
        -0.660823667935396e0,
        0.841267087271658e0,
        -0.253717501764397e2,
        0.485708963532948e3,
        0.880531517490555e3,
        0.265015592794626e7,
        -0.359287150025783e0,
        -0.656991567673753e3,
        0.241768149185367e1,
        0.856873461222588e0,
        0.655143675313458e0,
        -0.213535213206406e0,
        0.562974957606348e-2,
        -0.316955725450471e15,
        -0.699997000152457e-3,
        0.119845803210767e-1,
        0.193848122022095e-4,
        -0.215095749182309e-4,
    ];
    //
    // Star pressure for the backward temperature(pressure, entropy)
    // equation in sub-region 3b in [MPa].
    //
    const REGION_3B_PS_PSTAR: f64 = 100.0;
    //
    // Star entropy for the backward temperature(pressure, entropy)
    // equation in sub-region 3b in [J/Kg].
    //
    const REGION_3B_PS_SSTAR: f64 = 5.3e3;
    //
    // Star temperature for the backward temperature(pressure, entropy)
    // equation in sub-region 3b in [J/Kg].
    //
    const REGION_3B_PS_TSTAR: f64 = 860.0;
    //
    // Constants.
    //
    const A: f64 = 0.760;
    const B: f64 = 0.818;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_3B_PS_PSTAR + A;
    let eta = entropy / REGION_3B_PS_SSTAR - B;

    let mut sum = 0.0;
    for i in 0..28 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3B_PS_TSTAR * sum
}

///
/// This function computes the specific volume with respect to pressure
/// and entropy (backward equation) for the sub-region 3a
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `entropy`: The steam entropy [J/Kg.K].
///
/// - Returns:
///   - The specific volume [m3/Kg].
/// 
pub fn specific_volume_ps_region3a(pressure: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 28] = [
        -12, -12, -12, -10, -10, -10, -10, -8, -8, -8, -8, -6, -5, -4, -3, -3, -2, -2, -1, -1, 0,
        0, 0, 1, 2, 4, 5, 6,
    ];

    const JJ: [i32; 28] = [
        10, 12, 14, 4, 8, 10, 20, 5, 6, 14, 16, 28, 1, 5, 2, 4, 3, 8, 1, 2, 0, 1, 3, 0, 0, 2, 2, 0,
    ];

    const NN: [f64; 28] = [
        0.795544074093975e2,
        -0.238261242984590e4,
        0.176813100617787e5,
        -0.110524727080379e-2,
        -0.153213833655326e2,
        0.297544599376982e3,
        -0.350315206871242e8,
        0.277513761062119e0,
        -0.523964271036888e0,
        -0.148011182995403e6,
        0.160014899374266e7,
        0.170802322663427e13,
        0.246866996006494e-3,
        0.165326084797980e1,
        -0.118008384666987e0,
        0.253798642355900e1,
        0.965127704669424e0,
        -0.282172420532826e2,
        0.203224612353823e0,
        0.110648186063513e1,
        0.526127948451280e0,
        0.277000018736321e0,
        0.108153340501132e1,
        -0.744127885357893e-1,
        0.164094443541384e-1,
        -0.680468275301065e-1,
        0.257988576101640e-1,
        -0.145749861944416e-3,
    ];
    //
    // Star pressure for the backward specific volume(pressure, entropy)
    // equation in sub-region 3a in [MPa].
    //
    const REGION_3A_PS_PSTAR: f64 = 100.0;
    //
    // Star entropy for the backward specific volume(pressure, entropy)
    // equation in sub-region 3a in [J/Kg].
    //
    const REGION_3A_PS_SSTAR: f64 = 4.4e3;
    //
    // Star specific volume for the backward specific volume(pressure, entropy)
    // equation in sub-region 3a in [J/Kg].
    //
    const REGION_3A_PS_VSTAR: f64 = 0.0028;
    //
    // Constants.
    //
    const A: f64 = 0.187;
    const B: f64 = 0.755;
    //
    // Compute the specific volume.
    //
    let pi = pressure / REGION_3A_PS_PSTAR + A;
    let eta = entropy / REGION_3A_PS_SSTAR - B;

    let mut sum = 0.0;
    for i in 0..28 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3A_PS_VSTAR * sum
}

///
/// This function computes the specific volume with respect to pressure
/// and entropy (backward equation) for the sub-region 3b.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `entropy`: The steam entropy [J/Kg.K].
///
/// - Returns:
///   - The specific volume [m3/Kg].
/// 
pub fn specific_volume_ps_region3b(pressure: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 31] = [
        -12, -12, -12, -12, -12, -12, -10, -10, -10, -10, -8, -5, -5, -5, -4, -4, -4, -4, -3, -2,
        -2, -2, -2, -2, -2, 0, 0, 0, 1, 1, 2,
    ];

    const JJ: [i32; 31] = [
        0, 1, 2, 3, 5, 6, 0, 1, 2, 4, 0, 1, 2, 3, 0, 1, 2, 3, 1, 0, 1, 2, 3, 4, 12, 0, 1, 2, 0, 2,
        2,
    ];

    const NN: [f64; 31] = [
        0.591599780322238e-4,
        -0.185465997137856e-2,
        0.104190510480013e-1,
        0.598647302038590e-2,
        -0.771391189901699e0,
        0.172549765557036e1,
        -0.467076079846526e-3,
        0.134533823384439e-1,
        -0.808094336805495e-1,
        0.508139374365767e0,
        0.128584643361683e-2,
        -0.163899353915435e1,
        0.586938199318063e1,
        -0.292466667918613e1,
        -0.614076301499537e-2,
        0.576199014049172e1,
        -0.121613320606788e2,
        0.167637540957944e1,
        -0.744135838773463e1,
        0.378168091437659e-1,
        0.401432203027688e1,
        0.160279837479185e2,
        0.317848779347728e1,
        -0.358362310304853e1,
        -0.115995260446827e7,
        0.199256573577909e0,
        -0.122270624794624e0,
        -0.191449143716586e2,
        -0.150448002905284e-1,
        0.146407900162154e2,
        -0.327477787188230e1,
    ];
    //
    // Star pressure for the backward specific volume(pressure, entropy)
    // equation in sub-region 3b in [MPa].
    //
    const REGION_3B_PS_PSTAR: f64 = 100.0;
    //
    // Star entropy for the backward specific volume(pressure, entropy)
    // equation in sub-region 3b in [J/Kg].
    //
    const REGION_3B_PS_SSTAR: f64 = 5.3e3;
    //
    // Star specific volume for the backward specific volume(pressure, entropy)
    // equation in sub-region 3b in [J/Kg].
    //
    const REGION_3B_PS_VSTAR: f64 = 0.0088;
    //
    // Constants.
    //
    const A: f64 = 0.298;
    const B: f64 = 0.816;
    //
    // Compute the specific volume.
    //
    let pi = pressure / REGION_3B_PS_PSTAR + A;
    let eta = entropy / REGION_3B_PS_SSTAR - B;

    let mut sum = 0.0;
    for i in 0..31 {
        sum += NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]);
    }

    REGION_3B_PS_VSTAR * sum
}
