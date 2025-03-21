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

//! Steam Region 2 Module
//!
//! This module implements the polynomials which compute the thermodynamic
//! properties of steam in region 2 according to the revised release on the
//! IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of
//! Water and Steam, August 2007 (IAPWS-IF97)

use crate::steam_constants::*;
use crate::steam_units::*;
use std::f64;

///
/// Constant coefficients "J" for the ideal gas part.
///
const J: [i32; 9] = [0, 1, -5, -4, -3, -2, -1, 2, 3];

///
/// Constant coefficients "n" for the ideal gas part.
///
const N: [f64; 9] = [
    -0.96927686500217e1,
    0.10086655968018e2,
    -0.56087911283020e-2,
    0.71452738081455e-1,
    -0.40710498223928e0,
    0.14240819171444e1,
    -0.43839511319450e1,
    -0.28408632460772e0,
    0.21268463753307e-1,
];

///
/// Constant coefficients "IR" for the residual part.
///
const IR: [i32; 43] = [
    1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10,
    16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24,
];

///
/// Constant coefficients "JR" for the residual part.
///
const JR: [i32; 43] = [
    0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4,
    10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58,
];

///
/// Constant coefficients "NR" for the residual part.
///
const NR: [f64; 43] = [
    -0.17731742473213e-2,
    -0.17834862292358e-1,
    -0.45996013696365e-1,
    -0.57581259083432e-1,
    -0.50325278727930e-1,
    -0.33032641670203e-4,
    -0.18948987516315e-3,
    -0.39392777243355e-2,
    -0.43797295650573e-1,
    -0.26674547914087e-4,
    0.20481737692309e-7,
    0.43870667284435e-6,
    -0.32277677238570e-4,
    -0.15033924542148e-2,
    -0.40668253562649e-1,
    -0.78847309559367e-9,
    0.12790717852285e-7,
    0.48225372718507e-6,
    0.22922076337661e-5,
    -0.16714766451061e-10,
    -0.21171472321355e-2,
    -0.23895741934104e2,
    -0.59059564324270e-17,
    -0.12621808899101e-5,
    -0.38946842435739e-1,
    0.11256211360459e-10,
    -0.82311340897998e1,
    0.19809712802088e-7,
    0.10406965210174e-18,
    -0.10234747095929e-12,
    -0.10018179379511e-8,
    -0.80882908646985e-10,
    0.10693031879409e0,
    -0.33662250574171e0,
    0.89185845355421e-24,
    0.30629316876232e-12,
    -0.42002467698208e-5,
    -0.59056029685639e-25,
    0.37826947613457e-5,
    -0.12768608934681e-14,
    0.73087610595061e-28,
    0.55414715350778e-16,
    -0.94369707241210e-6,
];

///
/// Star pressure for region 2 in (MPa).
///
const REGION_2_PSTAR: f64 = 1.0;

///
/// Star temperature for region 2 in (K).
///
const REGION_2_TSTAR: f64 = 540.0;

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
    let sum: f64 = (0..9).map(|i| N[i] * tau.powi(J[i])).sum();
    f64::ln(pi) + sum
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
    let sum: f64 = (0..9)
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
    let sum: f64 = (0..9)
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
    // Constants.
    //
    const B: f64 = 0.5;
    //
    // Polynomial expression.
    //
    let taub = tau - B;
    let sum: f64 = (0..43)
        .map(|i| NR[i] * pi.powi(IR[i]) * taub.powi(JR[i]))
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
    // Constants.
    //
    const B: f64 = 0.5;
    //
    // Polynomial expression.
    //
    let taub = tau - B;
    let sum: f64 = (0..43)
        .filter(|&i| IR[i] != 0)
        .map(|i| NR[i] * IR[i] as f64 * pi.powi(IR[i] - 1) * taub.powi(JR[i]))
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
    // Constants.
    //
    const B: f64 = 0.5;
    //
    // Polynomial expression.
    //
    let taub = tau - B;
    let sum: f64 = (0..43)
        .filter(|&i| IR[i] > 1)
        .map(|i| NR[i] * IR[i] as f64 * (IR[i] - 1) as f64 * pi.powi(IR[i] - 2) * taub.powi(JR[i]))
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
    // Constants.
    //
    const B: f64 = 0.5;
    //
    // Polynomial expression.
    //
    let taub = tau - B;
    let sum: f64 = (0..43)
        .filter(|&i| JR[i] != 0)
        .map(|i| NR[i] * pi.powi(IR[i]) * JR[i] as f64 * taub.powi(JR[i] - 1))
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
    // Constants.
    //
    const B: f64 = 0.5;
    //
    // Polynomial expression.
    //
    let taub = tau - B;
    let sum: f64 = (0..43)
        .filter(|&i| JR[i] > 1)
        .map(|i| NR[i] * pi.powi(IR[i]) * JR[i] as f64 * (JR[i] - 1) as f64 * taub.powi(JR[i] - 2))
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
    // Constants.
    //
    const B: f64 = 0.5;
    //
    // Polynomial expression.
    //
    let taub = tau - B;
    let sum: f64 = (0..43)
        .filter(|&i| IR[i] != 0 && JR[i] != 0)
        .map(|i| NR[i] * IR[i] as f64 * pi.powi(IR[i] - 1) * JR[i] as f64 * taub.powi(JR[i] - 1))
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
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
    //
    // Compute the specific internal energy.
    //
    let gamtau = gibbs_ideal_grad_tau(tau) + gibbs_residual_grad_tau(pi, tau);
    IAPWS97_R * temperature * (tau * gamtau - 1.0 - pi * gibbs_residual_grad_pi(pi, tau))
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
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
    //
    // Compute the specific volume.
    //
    let res1 = 1.0 + gibbs_residual_grad_pi(pi, tau) * pi;
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
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
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
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
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
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
    //
    // Compute the speed of sound.
    //
    let gp = gibbs_residual_grad_pi(pi, tau);
    let res1 = IAPWS97_R * temperature * (1.0 + (2.0 + pi * gp) * pi * gp);
    let res2 = 1.0 - pi * pi * gibbs_residual_grad2_pi(pi, tau);
    let res3 = (1.0 + pi * gp - tau * pi * gibbs_residual_grad2_pi_tau(pi, tau)).powi(2);
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
    // Compute the dimensionless parameters pi and tau
    //
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
    //
    // Compute the specific isobaric heat capacity
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
    // Compute the dimensionless parameters pi and tau
    //
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
    //
    // Compute the specific isochoric heat capacity
    //
    let res1 = -tau * tau * (gibbs_ideal_grad2_tau(tau) + gibbs_residual_grad2_tau(pi, tau));
    let res2 = (1.0 + pi * gibbs_residual_grad_pi(pi, tau)
        - tau * pi * gibbs_residual_grad2_pi_tau(pi, tau))
    .powi(2);
    let res3 = 1.0 - pi * pi * gibbs_residual_grad2_pi(pi, tau);
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
    // Compute the ratio of specific heats
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
    // Compute the dimensionless parameters pi and tau
    //
    let pi = pressure / REGION_2_PSTAR;
    let tau = REGION_2_TSTAR / temperature;
    //
    // Compute the specific Gibbs free energy
    //
    IAPWS97_R * temperature * (gibbs_ideal(pi, tau) + gibbs_residual(pi, tau))
}
///
/// Computes the temperature with respect to pressure and enthalpy (backward equation)
/// for the sub-region 2a of IAPWS97.
///
/// - Arguments:
///   - `pressure`: The steam pressure in MPa.
///   - `enthalpy`: The steam enthalpy in J/Kg.
///
/// - Returns:
///   - The temperature in K.
///
pub fn temperature_ph_region2a(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 34] = [
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5,
        5, 6, 6, 7,
    ];

    const JJ: [i32; 34] = [
        0, 1, 2, 3, 7, 20, 0, 1, 2, 3, 7, 9, 11, 18, 44, 0, 2, 7, 36, 38, 40, 42, 44, 24, 44, 12,
        32, 44, 32, 36, 42, 34, 44, 28,
    ];

    const NN: [f64; 34] = [
        1089.8952318288,
        849.51654495535,
        -107.81748091826,
        33.153654801263,
        -7.4232016790248,
        11.765048724356,
        1.844574935579,
        -4.1792700549624,
        6.2478196935812,
        -17.344563108114,
        -200.58176862096,
        271.96065473796,
        -455.11318285818,
        3091.9688604755,
        252266.40357872,
        -6.1707422868339e-3,
        -0.31078046629583,
        11.670873077107,
        128127984.04046,
        -985549096.23276,
        2822454697.3002,
        -3594897141.0703,
        1722734991.3197,
        -13551.334240775,
        12848734.66465,
        1.3865724283226,
        235988.32556514,
        -13105236.545054,
        7399.9835474766,
        -551966.9703006,
        3715408.5996233,
        19127.72923966,
        -415351.64835634,
        -62.459855192507,
    ];
    //
    // Constants.
    //
    const REGION_2_PH_PSTAR: f64 = 1.0; // Star pressure in MPa
    const REGION_2_PH_HSTAR: f64 = 2000.0e3; // Star enthalpy in J/Kg
    const B: f64 = 2.1;
    //
    // Calculate dimensionless parameters.
    //
    let pi = pressure / REGION_2_PH_PSTAR;
    let eta = enthalpy / REGION_2_PH_HSTAR - B;
    //
    // Calculate temperature using the polynomial equation.
    //
    let temperature: f64 = (0..II.len())
        .map(|i| NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]))
        .sum();

    temperature
}

///
/// Computes the temperature with respect to pressure and enthalpy (backward equation)
/// for the sub-region 2b of IAPWS97.
///
/// - Arguments:
///   - `pressure`: The steam pressure in MPa.
///   - `enthalpy`: The steam enthalpy in J/Kg.
///
/// - Returns:
///   - The temperature in K.
///
pub fn temperature_ph_region2b(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 38] = [
        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4,
        5, 5, 5, 6, 7, 7, 9, 9,
    ];

    const JJ: [i32; 38] = [
        0, 1, 2, 12, 18, 24, 28, 40, 0, 2, 6, 12, 18, 24, 28, 40, 2, 8, 18, 40, 1, 2, 12, 24, 2,
        12, 18, 24, 28, 40, 18, 24, 40, 28, 2, 28, 1, 40,
    ];

    const NN: [f64; 38] = [
        1489.5041079516,
        743.07798314034,
        -97.708318797837,
        2.4742464705674,
        -0.63281320016026,
        1.1385952129658,
        -0.47811863648625,
        8.5208123431544e-03,
        0.93747147377932,
        3.3593118604916,
        3.3809355601454,
        0.16844539671904,
        0.73875745236695,
        -0.47128737436186,
        0.15020273139707,
        -0.002176411421975,
        -0.021810755324761,
        -0.10829784403677,
        -0.046333324635812,
        7.1280351959551e-05,
        1.1032831789999e-04,
        1.8955248387902e-04,
        3.0891541160537e-03,
        1.3555504554949e-03,
        2.8640237477456e-07,
        -1.0779857357512e-05,
        -7.6462712454814e-05,
        1.4052392818316e-05,
        -3.1083814331434e-05,
        -1.0302738212103e-06,
        2.821728163504e-07,
        1.2704902271945e-06,
        7.3803353468292e-08,
        -1.1030139238909e-08,
        -8.1456365207833e-14,
        -2.5180545682962e-11,
        -1.7565233969407e-18,
        8.6934156344163e-15,
    ];
    //
    // Constants.
    //
    const REGION_2_PH_PSTAR: f64 = 1.0; // Star pressure in MPa
    const REGION_2_PH_HSTAR: f64 = 2000.0e3; // Star enthalpy in J/Kg
    const A: f64 = 2.0;
    const B: f64 = 2.6;
    //
    // Calculate dimensionless parameters.
    //
    let pi = pressure / REGION_2_PH_PSTAR - A;
    let eta = enthalpy / REGION_2_PH_HSTAR - B;
    //
    // Calculate temperature using the polynomial equation.
    //
    let temperature: f64 = (0..II.len())
        .map(|i| NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]))
        .sum();

    temperature
}

///
/// Computes the temperature with respect to pressure and enthalpy (backward equation)
/// for the sub-region 2c of IAPWS97.
///
/// - Arguments:
///   - `pressure`: The steam pressure in MPa.
///   - `enthalpy`: The steam enthalpy in J/Kg.
///
/// - Returns:
///   - The temperature in K.
///
pub fn temperature_ph_region2c(pressure: f64, enthalpy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 23] = [
        -7, -7, -6, -6, -5, -5, -2, -2, -1, -1, 0, 0, 1, 1, 2, 6, 6, 6, 6, 6, 6, 6, 6,
    ];

    const JJ: [i32; 23] = [
        0, 4, 0, 2, 0, 2, 0, 1, 0, 2, 0, 1, 4, 8, 4, 0, 1, 4, 10, 12, 16, 20, 22,
    ];

    const NN: [f64; 23] = [
        -3236839855524.2,
        7326335090218.1,
        358250899454.47,
        -583401318515.9,
        -10783068217.47,
        20825544563.171,
        610747.83564516,
        859777.2253558,
        -25745.72360417,
        31081.088422714,
        1208.2315865936,
        482.19755109255,
        3.7966001272486,
        -10.842984880077,
        -0.04536417267666,
        1.4559115658698e-13,
        1.126159740723e-12,
        -1.7804982240686e-11,
        1.2324579690832e-07,
        -1.1606921130984e-06,
        2.7846367088554e-05,
        -5.9270038474176e-04,
        1.2918582991878e-03,
    ];
    //
    // Constants.
    //
    const REGION_2_PH_PSTAR: f64 = 1.0; // Star pressure in MPa
    const REGION_2_PH_HSTAR: f64 = 2000.0e3; // Star enthalpy in J/Kg
    const A: f64 = 25.0;
    const B: f64 = 1.8;
    //
    // Calculate dimensionless parameters.
    //
    let pi = pressure / REGION_2_PH_PSTAR + A;
    let eta = enthalpy / REGION_2_PH_HSTAR - B;
    //
    // Calculate temperature using the polynomial equation.
    //
    let temperature: f64 = (0..II.len())
        .map(|i| NN[i] * pi.powi(II[i]) * eta.powi(JJ[i]))
        .sum();

    temperature
}

///
/// This function computes the temperature with respect to pressure
/// and entropy (backward equation) for the sub-region 2a.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `entropy`: The steam entropy (J/Kg.K).
///
/// - Returns:
///   - The temperature (K).
///
pub fn temperature_ps_region2a(pressure: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [f64; 46] = [
        -1.50, -1.50, -1.50, -1.50, -1.50, -1.50, -1.25, -1.25, -1.25, -1.00, -1.00, -1.00, -1.00,
        -1.00, -1.00, -0.75, -0.75, -0.50, -0.50, -0.50, -0.50, -0.25, -0.25, -0.25, -0.25, 0.25,
        0.25, 0.25, 0.25, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.75, 0.75, 0.75, 0.75, 1.00,
        1.00, 1.25, 1.25, 1.50, 1.50,
    ];

    const JJ: [i32; 46] = [
        -24, -23, -19, -13, -11, -10, -19, -15, -6, -26, -21, -17, -16, -9, -8, -15, -14, -26, -13,
        -9, -7, -27, -25, -11, -6, 1, 4, 8, 11, 0, 1, 5, 6, 10, 14, 16, 0, 4, 9, 17, 7, 18, 3, 15,
        5, 18,
    ];

    const NN: [f64; 46] = [
        -0.39235983861984e6,
        0.51526573827270e6,
        0.40482443161048e5,
        -0.32193790923902e3,
        0.96961424218694e2,
        -0.22867846371773e2,
        -0.44942914124357e6,
        -0.50118336020166e4,
        0.35684463560015e0,
        0.44235335848190e5,
        -0.13673388811708e5,
        0.42163260207864e6,
        0.22516925837475e5,
        0.47442144865646e3,
        -0.14931130797647e3,
        -0.19781126320452e6,
        -0.23554399470760e5,
        -0.19070616302076e5,
        0.55375669883164e5,
        0.38293691437363e4,
        -0.60391860580567e3,
        0.19363102620331e4,
        0.42660643698610e4,
        -0.59780638872718e4,
        -0.70401463926862e3,
        0.33836784107553e3,
        0.20862786635187e2,
        0.33834172656196e-1,
        -0.43124428414893e-4,
        0.16653791356412e3,
        -0.13986292055898e3,
        -0.78849547999872e0,
        0.72132411753872e-1,
        -0.59754839398283e-2,
        -0.12141358953904e-4,
        0.23227096733871e-6,
        -0.10538463566194e2,
        0.20718925496502e1,
        -0.72193155260427e-1,
        0.20749887081120e-6,
        -0.18340657911379e-1,
        0.29036272348696e-6,
        0.21037527893619e0,
        0.25681239729999e-3,
        -0.12799002933781e-1,
        -0.82198102652018e-5,
    ];
    //
    // Star pressure for the backward temperature(pressure, entropy)
    // equation in region 2 in (MPa).
    //
    const REGION_2_PS_PSTAR: f64 = 1.0;
    //
    // Star entropy for the backward temperature(pressure, entropy)
    // equation in region 2 in (J/Kg.K).
    //
    const REGION_2_PS_SSTAR: f64 = 2000.0;
    //
    // Constants.
    //
    const B: f64 = 2.0;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_2_PS_PSTAR;
    let sigma = entropy / REGION_2_PS_SSTAR - B;

    let sum: f64 = (0..46)
        .map(|i| NN[i] * pi.powf(II[i]) * sigma.powi(JJ[i]))
        .sum();

    sum
}

///
/// This function computes the temperature with respect to pressure
/// and entropy (backward equation) for the sub-region 2b.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `entropy`: The steam entropy (J/Kg.K).
///
/// - Returns:
///   - The temperature (K).
///
pub fn temperature_ps_region2b(pressure: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 44] = [
        -6, -6, -5, -5, -4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1, -1, 0, 0, 0, 0,
        0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5,
    ];

    const JJ: [i32; 44] = [
        0, 11, 0, 11, 0, 1, 11, 0, 1, 11, 12, 0, 1, 6, 10, 0, 1, 5, 8, 9, 0, 1, 2, 4, 5, 6, 9, 0,
        1, 2, 3, 7, 8, 0, 1, 5, 0, 1, 3, 0, 1, 0, 1, 2,
    ];

    const NN: [f64; 44] = [
        0.31687665083497e6,
        0.20864175881858e2,
        -0.39859399803599e6,
        -0.21816058518877e2,
        0.22369785194242e6,
        -0.27841703445817e4,
        0.99207436071480e1,
        -0.75197512299157e5,
        0.29708605951158e4,
        -0.34406878548526e1,
        0.38815564249115e0,
        0.17511295085750e5,
        -0.14237112854449e4,
        0.10943803364167e1,
        0.89971619308495e0,
        -0.33759740098958e4,
        0.47162885818355e3,
        -0.19188241993679e1,
        0.41078580492196e0,
        -0.33465378172097e0,
        0.13870034777505e4,
        -0.40663326195838e3,
        0.41727347159610e2,
        0.21932549434532e1,
        -0.10320050009077e1,
        0.35882943516703e0,
        0.52511453726066e-2,
        0.12838916450705e2,
        -0.28642437219381e1,
        0.56912683664855e0,
        -0.99962954584931e-1,
        -0.32632037778459e-2,
        0.23320922576723e-3,
        -0.15334809857450e0,
        0.29072288239902e-1,
        0.37534702741167e-3,
        0.17296691702411e-2,
        -0.38556050844504e-3,
        -0.35017712292608e-4,
        -0.14566393631492e-4,
        0.56420857267269e-5,
        0.41286150074605e-7,
        -0.20684671118824e-7,
        0.16409393674725e-8,
    ];
    //
    // Star pressure for the backward temperature(pressure, entropy)
    // equation in region 2 in (MPa).
    //
    const REGION_2_PS_PSTAR: f64 = 1.0;
    //
    // Star entropy for the backward temperature(pressure, entropy)
    // equation in region 2 in (J/Kg.K).
    //
    const REGION_2_PS_SSTAR: f64 = 0.7853e3;
    //
    // Constants.
    //
    const B: f64 = 10.0;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_2_PS_PSTAR;
    let sigma = B - entropy / REGION_2_PS_SSTAR;

    let sum: f64 = (0..44)
        .map(|i| NN[i] * pi.powi(II[i]) * sigma.powi(JJ[i]))
        .sum();

    sum
}

///
/// This function computes the temperature with respect to pressure
/// and entropy (backward equation) for the sub-region 2c.
///
/// - Arguments:
///   - `pressure`: The steam pressure (MPa).
///   - `entropy`: The steam entropy (J/Kg.K).
///
/// - Returns:
///   - The temperature (K).
///
pub fn temperature_ps_region2c(pressure: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 30] = [
        -2, -2, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7,
    ];

    const JJ: [i32; 30] = [
        0, 1, 0, 0, 1, 2, 3, 0, 1, 3, 4, 0, 1, 2, 0, 1, 5, 0, 1, 4, 0, 1, 2, 0, 1, 0, 1, 3, 4, 5,
    ];

    const NN: [f64; 30] = [
        0.90968501005365e3,
        0.24045667088420e4,
        -0.59162326387130e3,
        0.54145404128074e3,
        -0.27098308411192e3,
        0.97976525097926e3,
        -0.46966772959435e3,
        0.14399274604723e2,
        -0.19104204230429e2,
        0.53299167111971e1,
        -0.21252975375934e2,
        -0.31147334413760e0,
        0.60334840894623e0,
        -0.42764839702509e-1,
        0.58185597255259e-2,
        -0.14597008284753e-1,
        0.56631175631027e-2,
        -0.76155864584577e-4,
        0.22440342919332e-3,
        -0.12561095013413e-4,
        0.63323132660934e-6,
        -0.20541989675375e-5,
        0.36405370390082e-7,
        -0.29759897789215e-8,
        0.10136618529763e-7,
        0.59925719692351e-11,
        -0.20677870105164e-10,
        -0.20874278181886e-10,
        0.10162166825089e-9,
        -0.16429828281347e-9,
    ];
    //
    // Star pressure for the backward temperature(pressure, entropy)
    // equation in region 2 in (MPa).
    //
    const REGION_2_PS_PSTAR: f64 = 1.0;
    //
    // Star entropy for the backward temperature(pressure, entropy)
    // equation in region 2 in (J/Kg.K).
    //
    const REGION_2_PS_SSTAR: f64 = 2.9251e3;
    //
    // Constants.
    //
    const B: f64 = 2.0;
    //
    // Compute the temperature.
    //
    let pi = pressure / REGION_2_PS_PSTAR;
    let sigma = B - entropy / REGION_2_PS_SSTAR;

    let sum: f64 = (0..30)
        .map(|i| NN[i] * pi.powi(II[i]) * sigma.powi(JJ[i]))
        .sum();

    sum
}

///
/// This function computes the pressure with respect to enthalpy
/// and entropy (backward equation) for the sub-region 2a.
///
/// - Arguments:
///   - `enthalpy`: The steam enthalpy (J/Kg).
///   - `entropy`: The steam entropy (J/Kg.K).
///
/// - Returns:
///   - The pressure (MPa).
///
pub fn pressure_hs_region2a(enthalpy: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 29] = [
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 6, 7,
    ];

    const JJ: [i32; 29] = [
        1, 3, 6, 16, 20, 22, 0, 1, 2, 3, 5, 6, 10, 16, 20, 22, 3, 16, 20, 0, 2, 3, 6, 16, 16, 3,
        16, 3, 1,
    ];

    const NN: [f64; 29] = [
        -0.182575361923032e-1,
        -0.125229548799536e0,
        0.592290437320145e0,
        0.604769706185122e1,
        0.238624965444474e3,
        -0.298639090222922e3,
        0.512250813040750e-1,
        -0.437266515606486e0,
        0.413336902999504e0,
        -0.516468254574773e1,
        -0.557014838445711e1,
        0.128555037824478e2,
        0.114144108953290e2,
        -0.119504225652714e3,
        -0.284777985961560e4,
        0.431757846408006e4,
        0.112894040802650e1,
        0.197409186206319e4,
        0.151612444706319e4,
        0.141324451421235e-1,
        0.585501282219601e0,
        -0.297258075863012e1,
        0.594567314847319e1,
        -0.623656565798905e4,
        0.965986235133332e4,
        0.681500934948134e1,
        -0.633207286824489e4,
        -0.558919224465760e1,
        0.400645798472063e-1,
    ];
    //
    // Star pressure for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2a in (MPa).
    //
    const REGION_2A_PHS_PSTAR: f64 = 4.0;
    //
    // Star enthalpy for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2a in (K).
    //
    const REGION_2A_PHS_HSTAR: f64 = 4200.0e3;
    //
    // Star entropy for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2a in (J/kg).
    //
    const REGION_2A_PHS_SSTAR: f64 = 12.0e3;
    //
    // Constants.
    //
    const A: f64 = 0.5;
    const B: f64 = 1.2;
    //
    // Compute the pressure.
    //
    let eta = enthalpy / REGION_2A_PHS_HSTAR - A;
    let sigma = entropy / REGION_2A_PHS_SSTAR - B;

    let sum: f64 = (0..29)
        .map(|i| NN[i] * eta.powi(II[i]) * sigma.powi(JJ[i]))
        .sum();

    REGION_2A_PHS_PSTAR * sum.powi(4)
}

///
/// This function computes the pressure with respect to enthalpy
/// and entropy (backward equation) for the sub-region 2b.
///
/// - Arguments:
///   - `enthalpy`: The steam enthalpy (J/Kg).
///   - `entropy`: The steam entropy (J/Kg.K).
///
/// - Returns:
///   - The pressure (MPa).
///
pub fn pressure_hs_region2b(enthalpy: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 33] = [
        0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 8,
        8, 12, 14,
    ];

    const JJ: [i32; 33] = [
        0, 1, 2, 4, 8, 0, 1, 2, 3, 5, 12, 1, 6, 18, 0, 1, 7, 12, 1, 16, 1, 12, 1, 8, 18, 1, 16, 1,
        3, 14, 18, 10, 16,
    ];

    const NN: [f64; 33] = [
        0.801496989929495e-1,
        -0.543862807146111e0,
        0.337455597421283e0,
        0.890555451157450e1,
        0.313840736431485e3,
        0.797367065977789e0,
        -0.121616973556240e1,
        0.872803386937477e1,
        -0.169769781757602e2,
        -0.186552827328416e3,
        0.951159274344237e5,
        -0.189168510120494e2,
        -0.433407037194840e4,
        0.543212633012715e9,
        0.144793408386013e0,
        0.128024559637516e3,
        -0.672309534071268e5,
        0.336972380095287e8,
        -0.586634196762720e3,
        -0.221403224769889e11,
        0.171606668708389e4,
        -0.570817595806302e9,
        -0.312109693178482e4,
        0.207841384633010e7,
        0.305605946157786e13,
        0.322157004314333e4,
        0.326810259797295e12,
        -0.144104158934487e4,
        0.410694867802691e3,
        0.109077066873024e12,
        -0.247964654258893e14,
        0.188801906865134e10,
        -0.123651009018773e15,
    ];
    //
    // Star pressure for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2b in (MPa).
    //
    const REGION_2B_PHS_PSTAR: f64 = 100.0;
    //
    // Star enthalpy for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2b in (K).
    //
    const REGION_2B_PHS_HSTAR: f64 = 4100.0e3;
    //
    // Star entropy for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2b in (J/kg).
    //
    const REGION_2B_PHS_SSTAR: f64 = 7.9e3;
    //
    // Constants.
    //
    const A: f64 = 0.6;
    const B: f64 = 1.01;
    //
    // Compute the pressure.
    //
    let eta = enthalpy / REGION_2B_PHS_HSTAR - A;
    let sigma = entropy / REGION_2B_PHS_SSTAR - B;

    let sum: f64 = (0..33)
        .map(|i| NN[i] * eta.powi(II[i]) * sigma.powi(JJ[i]))
        .sum();

    REGION_2B_PHS_PSTAR * sum.powi(4)
}

///
/// This function computes the pressure with respect to enthalpy
/// and entropy (backward equation) for the sub-region 2c.
///
/// - Arguments:
///   - `enthalpy`: The steam enthalpy (J/Kg).
///   - `entropy`: The steam entropy (J/Kg.K).
///
/// - Returns:
///   - The pressure (MPa).
///
pub fn pressure_hs_region2c(enthalpy: f64, entropy: f64) -> f64 {
    //
    // Constant polynomial coefficients.
    //
    const II: [i32; 31] = [
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5, 6, 6, 10, 12,
        16,
    ];

    const JJ: [i32; 31] = [
        0, 1, 2, 3, 4, 8, 0, 2, 5, 8, 14, 2, 3, 7, 10, 18, 0, 5, 8, 16, 18, 18, 1, 4, 6, 14, 8, 18,
        7, 7, 10,
    ];

    const NN: [f64; 31] = [
        0.112225607199012e0,
        -0.339005953606712e1,
        -0.320503911730094e2,
        -0.197597305104900e3,
        -0.407693861553446e3,
        0.132943775222331e5,
        0.170846839774007e1,
        0.373694198142245e2,
        0.358144365815434e4,
        0.423014446424664e6,
        -0.751071025760063e9,
        0.523446127607898e2,
        -0.228351290812417e3,
        -0.960652417056937e6,
        -0.807059292526074e8,
        0.162698017225669e13,
        0.772465073604171e0,
        0.463929973837746e5,
        -0.137317885134128e8,
        0.170470392630512e13,
        -0.251104628187308e14,
        0.317748830835520e14,
        0.538685623675312e2,
        -0.553089094625169e5,
        -0.102861522421405e7,
        0.204249418756234e13,
        0.273918446626977e9,
        -0.263963146312685e16,
        -0.107890854108088e10,
        -0.296492620980124e11,
        -0.111754907323424e16,
    ];
    //
    // Star pressure for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2c in (MPa).
    //
    const REGION_2C_PHS_PSTAR: f64 = 100.0;
    //
    // Star enthalpy for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2c in (K).
    //
    const REGION_2C_PHS_HSTAR: f64 = 3500.0e3;
    //
    // Star entropy for the backward pressure(enthalpy, entropy)
    // equation in region 2, subregion 2c in (J/kg).
    //
    const REGION_2C_PHS_SSTAR: f64 = 5.9e3;
    //
    // Constants.
    //
    const A: f64 = 0.7;
    const B: f64 = 1.1;
    //
    // Compute the pressure.
    //
    let eta = enthalpy / REGION_2C_PHS_HSTAR - A;
    let sigma = entropy / REGION_2C_PHS_SSTAR - B;

    let sum: f64 = (0..31)
        .map(|i| NN[i] * eta.powi(II[i]) * sigma.powi(JJ[i]))
        .sum();

    REGION_2C_PHS_PSTAR * sum.powi(4)
}
