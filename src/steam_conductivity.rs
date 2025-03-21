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

//! Steam Conductivity Module
//!
//! This module implements the polynomials which compute the thermal
//! conductivity of steam according to the release on the IAPWS Formulation
//! 2011 for the Thermal Conductivity of Ordinary Water Substance
//! (September 2011)

use crate::steam_constants::*;

///
/// This function computes the thermal conductivity of steam as a function
/// of density and temperature.
///
/// - Arguments:
///   - `density`: The steam density [Kg/m3].
///   - `temperature`: The steam temperature [K].
///
/// - Returns:
///   - The thermal conductivity in [W/K.m].
///
pub fn thermal_conductivity(density: f64, temperature: f64) -> f64 {
    //
    // Star temperature for the computation of thermal
    // conductivity(density,temperature) in [K].
    //
    const THCOND_TSTAR: f64 = 647.26;
    //
    // Star density for the computation of thermal
    // conductivity(density,temperature) in [Kg/m3].
    //
    const THCOND_RHOSTAR: f64 = 317.7;
    //
    // Star thermal conductivity for the computation of thermal
    // conductivity(density,temperature) in [W/K.m].
    //
    const THCOND_KSTAR: f64 = 1.0;
    //
    // Constant coefficients "a".
    //
    const A: [f64; 4] = [0.0102811, 0.0299621, 0.0156146, -0.00422464];
    //
    // Constant coefficients "b".
    //
    const B: [f64; 3] = [-0.397070, 0.400302, 1.060000];
    //
    // Constant coefficients "bb".
    //
    const BB: [f64; 2] = [-0.171587, 2.392190];
    //
    // Constant coefficients "c".
    //
    const C: [f64; 6] = [
        0.64285700, -4.1171700, -6.17937, 0.00308976, 0.0822994, 10.0932,
    ];
    //
    // Constant coefficients "d".
    //
    const D: [f64; 4] = [0.0701309, 0.0118520, 0.00169937, -1.0200];
    //
    // Compute the heat conductivity.
    //
    let tbar = temperature / THCOND_TSTAR;
    let rhobar = density / THCOND_RHOSTAR;

    let troot = tbar.sqrt();
    let mut tpow = troot;
    let mut lam = ZERO;
    //
    // Calculate first term with polynomial of tbar.
    //
    for a in A.iter() {
        lam += a * tpow;
        tpow *= tbar;
    }
    //
    // Calculate second term.
    //
    let tmp = (rhobar + BB[1]) * (rhobar + BB[1]);
    lam += B[0] + B[1] * rhobar + B[2] * (BB[0] * tmp).exp();
    //
    // Calculate third term.
    //
    let dtbar = (tbar - ONE).abs() + C[3];
    let dtbarpow = dtbar.powf(THREE / FIVE);
    let q = TWO + C[4] / dtbarpow;

    let s = if tbar >= ONE {
        ONE / dtbar
    } else {
        C[5] / dtbarpow
    };

    let rhobar18 = rhobar.powf(1.8);
    let rhobarq = rhobar.powf(q);

    lam += (D[0] / tbar.powi(10) + D[1]) * rhobar18 * (C[0] * (ONE - rhobar * rhobar18)).exp()
        + D[2] * s * rhobarq * ((q / (ONE + q)) * (ONE - rhobar * rhobarq)).exp()
        + D[3] * (C[1] * troot.powi(3) + C[2] / rhobar.powi(5)).exp();

    THCOND_KSTAR * lam
}
