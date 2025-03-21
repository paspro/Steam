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

//! Steam Region 5 Module
//!
//! This module implements the polynomials which compute the thermodynamic
//! properties of steam in region 5 according to the revised release on the
//! IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
//! Water and Steam, August 2007 (IAPWS-IF97)

use crate::steam_constants::*;
use crate::steam_units::*;

///
/// Constant coefficients "J" for the ideal gas part.
///
const J: [i32; 6] = [0, 1, -3, -2, -1, 2];

///
/// Constant coefficients "N" for the ideal gas part.
///
const N: [f64; 6] = [
    -0.13179983674201e2,
    0.68540841634434e1,
    -0.24805148933466e-1,
    0.36901534980333e0,
    -0.31161318213925e1,
    -0.32961626538917e0,
];

///
/// Constant coefficients "IR" for the ideal gas part.
///
const IR: [i32; 6] = [1, 1, 1, 2, 2, 3];

///
/// Constant coefficients "JR" for the ideal gas part.
///
const JR: [i32; 6] = [1, 2, 3, 3, 9, 7];

///
/// Constant coefficients "NR" for the ideal gas part.
///
const NR: [f64; 6] = [
    0.15736404855259e-2,
    0.90153761673944e-3,
    -0.50270077677648e-2,
    0.22440037409485e-5,
    -0.41163275453471e-5,
    0.37919454822955e-7,
];

///
/// Star pressure for region 5 in (MPa).
///
const REGION_5_PSTAR: f64 = 1.0;

///
/// Star temperature for region 5 in (K).
///
const REGION_5_TSTAR: f64 = 1000.0;

///
/// This function computes the dimensionless specific Gibbs free energy
/// for the ideal gas part.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The dimensionless specific Gibbs free energy.
///
fn gibbs_ideal(pi: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6).map(|i| N[i] * tau.powi(J[i])).sum();

    pi.ln() + sum
}

///
/// This function computes the derivative of the specific Gibbs free energy
/// with respect to the parameter tau for the ideal gas part.
///
/// - Arguments:
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d(Gibbs)/d(tau).
///
fn gibbs_ideal_grad_tau(tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .filter(|&i| J[i] != 0)
        .map(|i| N[i] * J[i] as f64 * tau.powi(J[i] - 1))
        .sum();

    sum
}

///
/// This function computes the second derivative of the specific Gibbs
/// free energy with respect to the parameter tau for the ideal gas part.
///
/// - Arguments:
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(Gibbs)/d(tau)2.
///
fn gibbs_ideal_grad2_tau(tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .filter(|&i| J[i] != 0 && J[i] != 1)
        .map(|i| N[i] * J[i] as f64 * (J[i] - 1) as f64 * tau.powi(J[i] - 2))
        .sum();

    sum
}

///
/// This function computes the dimensionless specific Gibbs free energy
/// for the residual part.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The dimensionless specific Gibbs free energy.
///
fn gibbs_residual(pi: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .map(|i| NR[i] * pi.powi(IR[i]) * tau.powi(JR[i]))
        .sum();

    sum
}

///
/// This function computes the derivative of the specific Gibbs free energy
/// with respect to the parameter pi for the residual part.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d(Gibbs)/d(pi).
///
fn gibbs_residual_grad_pi(pi: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .map(|i| NR[i] * IR[i] as f64 * pi.powi(IR[i] - 1) * tau.powi(JR[i]))
        .sum();

    sum
}

///
/// This function computes the second derivative of the specific Gibbs
/// free energy with respect to the parameter pi for the residual part.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(Gibbs)/d(pi)2.
///
fn gibbs_residual_grad2_pi(pi: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .filter(|&i| IR[i] > 1)
        .map(|i| NR[i] * IR[i] as f64 * (IR[i] - 1) as f64 * pi.powi(IR[i] - 2) * tau.powi(JR[i]))
        .sum();

    sum
}

///
/// This function computes the derivative of the specific Gibbs free energy
/// with respect to the parameter tau for the residual part.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d(Gibbs)/d(tau).
///
fn gibbs_residual_grad_tau(pi: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .map(|i| NR[i] * pi.powi(IR[i]) * JR[i] as f64 * tau.powi(JR[i] - 1))
        .sum();

    sum
}

///
/// This function computes the second derivative of the specific Gibbs
/// free energy with respect to the parameter tau for the residual part.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(Gibbs)/d(tau)2.
///
fn gibbs_residual_grad2_tau(pi: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .filter(|&i| JR[i] > 1)
        .map(|i| NR[i] * pi.powi(IR[i]) * JR[i] as f64 * (JR[i] - 1) as f64 * tau.powi(JR[i] - 2))
        .sum();

    sum
}

///
/// This function computes the second derivative of the specific Gibbs
/// free energy with respect to the parameters pi and tau for the
/// residual part.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(Gibbs)/d(pi)d(tau).
///
fn gibbs_residual_grad2_pi_tau(pi: f64, tau: f64) -> f64 {
    //
    // Polynomial expression.
    //
    let sum: f64 = (0..6)
        .map(|i| NR[i] * IR[i] as f64 * pi.powi(IR[i] - 1) * JR[i] as f64 * tau.powi(JR[i] - 1))
        .sum();

    sum
}

///
/// This function computes the specific internal energy with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The specific internal energy (J/Kg).
///
pub fn specific_internal_energy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the specific internal energy.
    //
    let gamtau = gibbs_ideal_grad_tau(tau) + gibbs_residual_grad_tau(pi, tau);
    IAPWS97_R * temperature * (tau * gamtau - ONE - pi * gibbs_residual_grad_pi(pi, tau))
}

///
/// This function computes the specific volume with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The specific volume (m3/Kg).
///
pub fn specific_volume(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the specific volume.
    //
    let res1 = ONE + gibbs_residual_grad_pi(pi, tau) * pi;
    IAPWS97_R * temperature * res1 / (pressure * MEGA)
}

///
/// This function computes the specific entropy with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The specific entropy (J/kg.K).
///
pub fn specific_entropy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the specific entropy.
    //
    let gam = gibbs_ideal(pi, tau) + gibbs_residual(pi, tau);
    let gamtau = gibbs_ideal_grad_tau(tau) + gibbs_residual_grad_tau(pi, tau);
    IAPWS97_R * (tau * gamtau - gam)
}

///
/// This function computes the specific enthalpy with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The specific enthalpy (J/Kg).
///
pub fn specific_enthalpy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the specific enthalpy.
    //
    let gamtau = gibbs_ideal_grad_tau(tau) + gibbs_residual_grad_tau(pi, tau);
    IAPWS97_R * temperature * tau * gamtau
}

///
/// This function computes the speed of sound with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The speed of sound (m/sec).
///
pub fn speed_of_sound(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the speed of sound.
    //
    let gp = gibbs_residual_grad_pi(pi, tau);
    let res1 = IAPWS97_R * temperature * (ONE + (TWO + pi * gp) * pi * gp);
    let res2 = ONE - pi * pi * gibbs_residual_grad2_pi(pi, tau);
    let res3 = (ONE + pi * gp - tau * pi * gibbs_residual_grad2_pi_tau(pi, tau)).powi(2);
    let res4 = tau * tau * (gibbs_ideal_grad2_tau(tau) + gibbs_residual_grad2_tau(pi, tau));

    (res1 / (res2 + res3 / res4)).sqrt()
}

///
/// This function computes the specific isobaric heat capacity with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The specific isobaric heat capacity (J/Kg.K).
///
pub fn specific_isobaric_heat_capacity(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the specific isobaric heat capacity.
    //
    -IAPWS97_R * tau * tau * (gibbs_ideal_grad2_tau(tau) + gibbs_residual_grad2_tau(pi, tau))
}

///
/// This function computes the specific isochoric heat capacity with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The specific isochoric heat capacity (J/Kg.K).
///
pub fn specific_isochoric_heat_capacity(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the specific isochoric heat capacity.
    //
    let res1 = -tau * tau * (gibbs_ideal_grad2_tau(tau) + gibbs_residual_grad2_tau(pi, tau));
    let res2 = (ONE + pi * gibbs_residual_grad_pi(pi, tau)
        - tau * pi * gibbs_residual_grad2_pi_tau(pi, tau))
    .powi(2);
    let res3 = ONE - pi * pi * gibbs_residual_grad2_pi(pi, tau);

    IAPWS97_R * (res1 - res2 / res3)
}

///
/// This function computes the ratio of specific heats with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The ratio of specific heats.
///
pub fn ratio_of_specific_heats(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the ratio of specific heats.
    //
    specific_isobaric_heat_capacity(pressure, temperature)
        / specific_isochoric_heat_capacity(pressure, temperature)
}

///
/// This function computes the specific Gibbs free energy with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The specific Gibbs free energy (J/kg).
///
pub fn specific_gibbs_free_energy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_5_PSTAR;
    let tau = REGION_5_TSTAR / temperature;
    //
    // Compute the specific Gibbs free energy.
    //
    IAPWS97_R * temperature * (gibbs_ideal(pi, tau) + gibbs_residual(pi, tau))
}
