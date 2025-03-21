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

//! Steam Viscosity Module
//!
//! This module implements the polynomials which compute the viscosity
//! of steam according to the release on the IAPWS Formulation 2008 for
//! the Viscosity of Ordinary Water Substance (September 2008)
//!
//! @author Panos Asproulis
//! @date 2014
//! @copyright Panos Asproulis (2014-2016). All Rights Reserved.

use crate::steam_constants::*;

///
/// This function computes the viscosity in the dilute-gas limit.
///
/// - Arguments:
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The viscosity in the dilute-gas limit.
///
fn mu0(tau: f64) -> f64 {
    //
    // Constant coefficients "I".
    //
    const I: [i32; 4] = [0, 1, 2, 3];
    //
    // Constant coefficients "H".
    //
    const H: [f64; 4] = [1.67752, 2.20462, 0.6366564, -0.241605];
    //
    // Compute the viscosity in the dilute-gas limit.
    //
    let sum: f64 = (0..4).map(|i| H[i] * tau.powi(I[i])).sum();

    100.0 / (tau.sqrt() * sum)
}

///
/// This function computes the contribution to viscosity due to finite
/// density.
///
/// - Arguments:
///   - `del`: Dimensionless density parameter.
///   - `tau`: Dimensionless temperature parameter.
///
/// - Returns:
///   - The contribution to viscosity due to finite density.
///
fn mu1(del: f64, tau: f64) -> f64 {
    //
    // Constant coefficients "H".
    //
    const H: [[f64; 7]; 6] = [
        [
            5.20094e-1,
            2.22531e-1,
            -2.81378e-1,
            1.61913e-1,
            -3.25372e-2,
            0.0,
            0.0,
        ],
        [
            8.50895e-2,
            9.99115e-1,
            -9.06851e-1,
            2.57399e-1,
            0.0,
            0.0,
            0.0,
        ],
        [-1.08374e0, 1.88797e0, -7.72479e-1, 0.0, 0.0, 0.0, 0.0],
        [
            -2.89555e-1,
            1.26613e0,
            -4.89837e-1,
            0.0,
            6.98452e-2,
            0.0,
            -4.35673e-3,
        ],
        [0.0, 0.0, -2.57040e-1, 0.0, 0.0, 8.72102e-3, 0.0],
        [0.0, 1.20573e-1, 0.0, 0.0, 0.0, 0.0, -5.93264e-4],
    ];
    //
    // Compute the contribution to viscosity due to finite density.
    //
    let mut sum = 0.0;

    for i in 0..6 {
        let tau1 = (tau - 1.0).powi(i as i32);
        sum += (0..7)
            .filter(|&j| H[i][j] != 0.0)
            .map(|j| H[i][j] * tau1 * (del - 1.0).powi(j as i32))
            .sum::<f64>();
    }

    (del * sum).exp()
}

///
/// This function computes the viscosity as a function of density and
/// temperature.
///
/// - Arguments:
///   - `density`: The steam density (Kg/m3).
///   - `temperature`: The steam temperature (K).
///
/// - Returns:
///   - The viscosity in (Pa.sec).
///
pub fn viscosity(density: f64, temperature: f64) -> f64 {
    //
    // Star viscosity for the computation of viscosity(density,temperature)
    // in (Pa.sec).
    //
    const VISCOSITY_MUSTAR: f64 = 1.0e-6;
    //
    // Compute the viscosity.
    //
    let del = density / IAPWS97_RHOCRIT;
    let tau = IAPWS97_TCRIT / temperature;

    VISCOSITY_MUSTAR * mu0(tau) * mu1(del, tau)
}
