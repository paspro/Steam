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

//! Steam Region 1 Module
//!
//! This module implements the polynomials which compute the thermodynamic
//! properties of steam in region 1 according to the revised release on the
//! IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
//! Water and Steam, August 2007 (IAPWS-IF97)

use crate::steam_constants::*;
use crate::steam_units::*;
use std::f64;

///
/// Constant coefficients "I".
///
const I: [i32; 34] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29,
    30, 31, 32,
];

///
/// Constant coefficients "J".
///
const J: [i32; 34] = [
    -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11,
    -6, -29, -31, -38, -39, -40, -41,
];

///
/// Constant coefficients "N".
///
const N: [f64; 34] = [
    0.14632971213167e0,
    -0.84548187169114e0,
    -0.37563603672040e1,
    0.33855169168385e1,
    -0.95791963387872e0,
    0.15772038513228e0,
    -0.16616417199501e-1,
    0.81214629983568e-3,
    0.28319080123804e-3,
    -0.60706301565874e-3,
    -0.18990068218419e-1,
    -0.32529748770505e-1,
    -0.21841717175414e-1,
    -0.52838357969930e-4,
    -0.47184321073267e-3,
    -0.30001780793026e-3,
    0.47661393906987e-4,
    -0.44141845330846e-5,
    -0.72694996297594e-15,
    -0.31679644845054e-4,
    -0.28270797985312e-5,
    -0.85205128120103e-9,
    -0.22425281908000e-5,
    -0.65171222895601e-6,
    -0.14341729937924e-12,
    -0.40516996860117e-6,
    -0.12734301741641e-8,
    -0.17424871230634e-9,
    -0.68762131295531e-18,
    0.14478307828521e-19,
    0.26335781662795e-22,
    -0.11947622640071e-22,
    0.18228094581404e-23,
    -0.93537087292458e-25,
];

///
/// Star pressure for region 1 in [MPa].
///
const REGION_1_PSTAR: f64 = 16.53;

///
/// Star temperature for region 1 in [K].
///
const REGION_1_TSTAR: f64 = 1386.0;

///
/// This function computes the dimensionless specific Gibbs free energy.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The dimensionless specific Gibbs free energy.
///
fn gibbs_nondim(pi: f64, tau: f64) -> f64 {
    //
    // Constants.
    //
    const A: f64 = 7.1;
    const B: f64 = 1.222;
    //
    // Polynomial expression.
    //
    let pia = A - pi;
    let taub = tau - B;

    let sum: f64 = (0..34)
        .map(|i| N[i] * pia.powi(I[i]) * taub.powi(J[i]))
        .sum();

    sum
}

///
/// This function computes the derivative of the specific Gibbs free energy
/// with respect to the parameter pi.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d(Gibbs)/d(pi).
///
fn gibbs_grad_pi(pi: f64, tau: f64) -> f64 {
    //
    // Constants.
    //
    const A: f64 = 7.1;
    const B: f64 = 1.222;
    //
    // Polynomial expression.
    //
    let pia = A - pi;
    let taub = tau - B;

    let sum: f64 = (0..34)
        .filter(|&i| I[i] > 0)
        .map(|i| -N[i] * I[i] as f64 * pia.powi(I[i] - 1) * taub.powi(J[i]))
        .sum();

    sum
}

///
/// This function computes the second derivative of the specific Gibbs
/// free energy with respect to the parameter pi.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(Gibbs)/d(pi)2.
///
fn gibbs_grad2_pi(pi: f64, tau: f64) -> f64 {
    //
    // Constants.
    //
    const A: f64 = 7.1;
    const B: f64 = 1.222;
    //
    // Polynomial expression.
    //
    let pia = A - pi;
    let taub = tau - B;

    let sum: f64 = (0..34)
        .filter(|&i| I[i] > 1)
        .map(|i| N[i] * I[i] as f64 * (I[i] - 1) as f64 * pia.powi(I[i] - 2) * taub.powi(J[i]))
        .sum();

    sum
}

///
/// This function computes the derivative of the specific Gibbs free energy
/// with respect to the parameter tau.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d(Gibbs)/d(tau).
///
fn gibbs_grad_tau(pi: f64, tau: f64) -> f64 {
    //
    // Constants.
    //
    const A: f64 = 7.1;
    const B: f64 = 1.222;
    //
    // Polynomial expression.
    //
    let pia = A - pi;
    let taub = tau - B;

    let sum: f64 = (0..34)
        .filter(|&i| J[i] != 0)
        .map(|i| N[i] * pia.powi(I[i]) * J[i] as f64 * taub.powi(J[i] - 1))
        .sum();

    sum
}

///
/// This function computes the second derivative of the specific Gibbs
/// free energy with respect to the parameter tau.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(Gibbs)/d(tau)2.
///
fn gibbs_grad2_tau(pi: f64, tau: f64) -> f64 {
    //
    // Constants.
    //
    const A: f64 = 7.1;
    const B: f64 = 1.222;
    //
    // Polynomial expression.
    //
    let pia = A - pi;
    let taub = tau - B;

    let sum: f64 = (0..34)
        .filter(|&i| J[i] > 1 || J[i] < 0)
        .map(|i| N[i] * pia.powi(I[i]) * J[i] as f64 * (J[i] - 1) as f64 * taub.powi(J[i] - 2))
        .sum();

    sum
}

///
/// This function computes the second derivative of the specific Gibbs
/// free energy with respect to the parameters pi and tau.
///
/// - Arguments:
///   - `pi`: Dimensionless pressure parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The derivative d2(Gibbs)/d(pi)d(tau).
///
fn gibbs_grad2_pi_tau(pi: f64, tau: f64) -> f64 {
    //
    // Constants.
    //
    const A: f64 = 7.1;
    const B: f64 = 1.222;
    //
    // Polynomial expression.
    //
    let pia = A - pi;
    let taub = tau - B;

    let sum: f64 = (0..34)
        .filter(|&i| I[i] > 0 && J[i] != 0)
        .map(|i| -N[i] * I[i] as f64 * pia.powi(I[i] - 1) * J[i] as f64 * taub.powi(J[i] - 1))
        .sum();

    sum
}

///
/// This function computes the specific internal energy with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific internal energy [J/Kg].
///
pub fn specific_internal_energy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the specific internal energy.
    //
    IAPWS97_R * temperature * (tau * gibbs_grad_tau(pi, tau) - pi * gibbs_grad_pi(pi, tau))
}

///
/// This function computes the specific volume with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific volume [m3/Kg].
///
pub fn specific_volume(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the specific volume.
    //
    IAPWS97_R * temperature * pi * gibbs_grad_pi(pi, tau) / (pressure * MEGA)
}

///
/// This function computes the specific entropy with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific entropy [J/kg.K].
///
pub fn specific_entropy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the specific entropy.
    //
    IAPWS97_R * (tau * gibbs_grad_tau(pi, tau) - gibbs_nondim(pi, tau))
}

///
/// This function computes the specific enthalpy with respect to
/// pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific enthalpy [J/Kg].
///
pub fn specific_enthalpy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the specific enthalpy.
    //
    IAPWS97_R * temperature * tau * gibbs_grad_tau(pi, tau)
}

///
/// This function computes the speed of sound with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The speed of sound [m/sec].
///
pub fn speed_of_sound(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the speed of sound.
    //
    let gp = gibbs_grad_pi(pi, tau);
    let res1 = IAPWS97_R * temperature;
    let res2 = (gp - tau * gibbs_grad2_pi_tau(pi, tau)).powi(2);
    let res3 = tau * tau * gibbs_grad2_tau(pi, tau);
    let res4 = res2 / res3;

    gp * (res1 / (res4 - gibbs_grad2_pi(pi, tau))).sqrt()
}

///
/// This function computes the specific isobaric heat capacity with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific isobaric heat capacity [J/Kg.K].
///
pub fn specific_isobaric_heat_capacity(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the specific isobaric heat capacity.
    //
    -IAPWS97_R * tau * tau * gibbs_grad2_tau(pi, tau)
}

///
/// This function computes the specific isochoric heat capacity with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific isochoric heat capacity [J/Kg.K].
///
pub fn specific_isochoric_heat_capacity(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the specific isochoric heat capacity.
    //
    let res1 = -tau * tau * gibbs_grad2_tau(pi, tau);
    let res2 = (gibbs_grad_pi(pi, tau) - tau * gibbs_grad2_pi_tau(pi, tau)).powi(2);

    IAPWS97_R * (res1 + res2 / gibbs_grad2_pi(pi, tau))
}

///
/// This function computes the ratio of specific heats with respect
/// to pressure and temperature.
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
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
///   - `pressure`: The steam pressure [MPa].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The specific Gibbs free energy [J/kg].
///
pub fn specific_gibbs_free_energy(pressure: f64, temperature: f64) -> f64 {
    //
    // Compute the dimensionless parameters pi and tau.
    //
    let pi = pressure / REGION_1_PSTAR;
    let tau = REGION_1_TSTAR / temperature;
    //
    // Compute the specific Gibbs free energy.
    //
    IAPWS97_R * temperature * gibbs_nondim(pi, tau)
}

///
/// This function computes the temperature with respect to pressure
/// and enthalpy (backward equation).
///
/// - Arguments:
///   - `pressure`: The steam pressure [MPa].
///   - `enthalpy`: The steam enthalpy [J/Kg].
///
/// - Returns:
///   - The temperature [K].
///
pub fn temperature_ph(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 20] = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 5, 6];

    const JJ: [i32; 20] = [
        0, 1, 2, 6, 22, 32, 0, 1, 2, 3, 4, 10, 32, 10, 32, 10, 32, 32, 32, 32,
    ];

    const NN: [f64; 20] = [
        -238.72489924521e0,
        404.21188637945e0,
        113.49746881718e0,
        -5.8457616048039e0,
        -1.528548241314e-4,
        -1.0866707695377e-6,
        -13.391744872602e0,
        43.211039183559e0,
        -54.010067170506e0,
        30.535892203916e0,
        -6.5964749423638e0,
        9.3965400878363e-3,
        1.157364750534e-7,
        -2.5858641282073e-5,
        -4.0644363084799e-9,
        6.6456186191635e-8,
        8.0670734103027e-11,
        -9.3477771213947e-13,
        5.8265442020601e-15,
        -1.5020185953503e-17,
    ];
    //
    // Star pressure for the backward temperature(pressure, enthalpy)
    // equation in region 1 in [MPa].
    //
    const REGION_1_TPH_PSTAR: f64 = 1.0;
    //
    // Star temperature for the backward temperature(pressure, enthalpy)
    // equation in region 1 in [K].
    //
    const REGION_1_TPH_TSTAR: f64 = 1.0;
    //
    // Star enthalpy for the backward temperature(pressure, enthalpy)
    // equation in region 1 in [J/kg].
    //
    const REGION_1_TPH_HSTAR: f64 = 2500.0e3;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_1_TPH_PSTAR;
    let e1 = ONE + (enthalpy / REGION_1_TPH_HSTAR);

    let sum: f64 = (0..20)
        .map(|i| NN[i] * pi.powi(II[i]) * e1.powi(JJ[i]))
        .sum();

    REGION_1_TPH_TSTAR * sum
}

///
/// This function computes the pressure with respect to enthalpy
/// and entropy (backward equation).
///
/// - Arguments:
///   - `enthalpy`: The steam enthalpy [J/Kg].
///   - `entropy`: The steam entropy [J/Kg.K].
///
/// - Returns:
///   - The pressure [MPa].
///
pub fn pressure_hs(enthalpy: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 19] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 5];

    const JJ: [i32; 19] = [0, 1, 2, 4, 5, 6, 8, 14, 0, 1, 4, 6, 0, 1, 10, 4, 1, 4, 0];

    const NN: [f64; 19] = [
        -0.691997014660582e0,
        -0.183612548787560e2,
        -0.928332409297335e1,
        0.659639569909906e2,
        -0.162060388912024e2,
        0.450620017338667e3,
        0.854680678224170e3,
        0.607523214001162e4,
        0.326487682621856e2,
        -0.269408844582931e2,
        -0.319947848334300e3,
        -0.928354307043320e3,
        0.303634537455249e2,
        -0.650540422444146e2,
        -0.430991316516130e4,
        -0.747512324096068e3,
        0.730000345529245e3,
        0.114284032569021e4,
        -0.436407041874559e3,
    ];
    //
    // Star pressure for the backward pressure(enthalpy, entropy)
    // equation in region 1 in [MPa].
    //
    const REGION_1_PHS_PSTAR: f64 = 100.0;
    //
    // Star enthalpy for the backward pressure(enthalpy, entropy)
    // equation in region 1 in [K].
    //
    const REGION_1_PHS_HSTAR: f64 = 3400.0e3;
    //
    // Star entropy for the backward pressure(enthalpy, entropy)
    // equation in region 1 in [J/kg].
    //
    const REGION_1_PHS_SSTAR: f64 = 7.6e3;
    //
    // Constants.
    //
    const A: f64 = 0.05;
    //
    // Compute the pressure.
    //
    let eta = enthalpy / REGION_1_PHS_HSTAR + A;
    let sigma = entropy / REGION_1_PHS_SSTAR + A;

    let sum: f64 = (0..19)
        .map(|i| NN[i] * eta.powi(II[i]) * sigma.powi(JJ[i]))
        .sum();

    REGION_1_PHS_PSTAR * sum
}
