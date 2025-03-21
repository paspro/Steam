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

//! Integration tests for Steam Tables
//!
//! This file implements tests for the thermodynamic properties of steam
//! across all regions according to the IAPWS Industrial Formulation 1997
//! for the Thermodynamic Properties of Water and Steam, August 2007 (IAPWS-IF97)

use steam::steam_constants::*;
use steam::steam_numerics::*;
use steam::steam_regions::*;
use steam::steam_states::SteamState;

use approx::assert_relative_eq;
use rand::prelude::*;
use rayon::prelude::*;

#[test]
fn test_region1() {
    //
    // We are only working with region 1 in this test.
    //
    use steam::steam_region1::*;

    println!("\nRegion 1 Test - Computational Errors");
    println!("====================================");
    //
    // Input data.
    //
    let p = [3.0, 80.0, 3.0];
    let t = [300.0, 300.0, 500.0];
    //
    // Expected values.
    //
    let v_real = [0.100215168e-2, 0.971180894e-3, 0.120241800e-2];
    let h_real = [0.115331273e6, 0.184142828e6, 0.975542239e6];
    let u_real = [0.112324818e6, 0.106448356e6, 0.971934985e6];
    let s_real = [0.392294792e3, 0.368563852e3, 0.258041912e4];
    let cp_real = [0.417301218e4, 0.401008987e4, 0.465580682e4];
    let w_real = [0.150773921e4, 0.163469054e4, 0.124071337e4];

    for i in 0..3 {
        println!("============");
        println!("Test Case : {}", i + 1);
        println!("============");

        let mut ss = SteamState::new();
        let valid_state = ss.set_primary_properties(1, p[i], t[i]);

        println!("Steam state in region 1 ?                = {}", valid_state);

        let value = speed_of_sound(p[i], t[i]);
        let error = (value - w_real[i]).abs() / w_real[i].abs();
        println!("Error in speed of sound                  = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_volume();
        let error = (value - v_real[i]).abs() / v_real[i].abs();
        println!("Error in specific volume                 = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_enthalpy();
        let error = (value - h_real[i]).abs() / h_real[i].abs();
        println!("Error in specific enthalpy               = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_entropy();
        let error = (value - s_real[i]).abs() / s_real[i].abs();
        println!("Error in specific entropy                = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_internal_energy();
        let error = (value - u_real[i]).abs() / u_real[i].abs();
        println!("Error in specific internal energy        = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_isobaric_heat_capacity();
        let error = (value - cp_real[i]).abs() / cp_real[i].abs();
        println!("Error in specific isobaric heat capacity = {}", error);
        assert!(error <= 1e-2);

        let value = temperature_ph(p[i], h_real[i]);
        let error = (value - t[i]).abs() / t[i].abs();
        println!("Error in backward temperature(p,h)       = {}", error);
        assert!(error <= 1e-2);

        let value = pressure_hs(h_real[i], s_real[i]);
        let error = (value - p[i]).abs() / p[i].abs();
        println!("Error in backward pressure(h,s)          = {}", error);
        assert!(error <= 1e-2);

        println!("\nAttempting to compute temperature by inverting polynomial h = h(p,t)");
        let value = function_inverter_y(
            &specific_enthalpy,
            h_real[i],
            p[i],
            200.0,
            600.0,
            1e-5,
            100,
            3,
            true,
        );
        assert_relative_eq!(value, t[i], max_relative = 1e-4);

        println!("\nAttempting to compute enthalpy by inverting polynomial p = p(h,s)");
        let value = function_inverter_x(
            &pressure_hs,
            p[i],
            s_real[i],
            0.10e6,
            1.0e7,
            1e-5,
            100,
            3,
            true,
        );
        assert_relative_eq!(value, h_real[i], max_relative = 1e-4);
    }
}

#[test]
fn test_region2() {
    //
    // We are only working with region 2 in this test.
    //
    use steam::steam_region2::*;

    println!("\nRegion 2 Test - Computational Errors");
    println!("====================================");
    //
    // Input data.
    //
    let p = [0.0035, 0.0035, 30.0];
    let t = [300.0, 700.0, 700.0];
    //
    // Expected values.
    //
    let v_real = [0.394913866e2, 0.923015898e2, 0.542946619e-2];
    let h_real = [0.254991145e7, 0.333568375e7, 0.263149474e7];
    let u_real = [0.241169160e7, 0.301262819e7, 0.246861076e7];
    let s_real = [0.852238967e4, 0.101749996e5, 0.517540298e4];
    let cp_real = [0.191399162e4, 0.208141274e4, 0.103505092e5];
    let w_real = [0.427920172e3, 0.644289068e3, 0.480386523e3];

    for i in 0..3 {
        println!("============");
        println!("Test Case : {}", i + 1);
        println!("============");

        let mut ss = SteamState::new();
        let valid_state = ss.set_primary_properties(2, p[i], t[i]);

        let in2a = in_region2a(p[i], t[i]);
        let in2b = in_region2b(p[i], t[i]);
        let in2c = in_region2c(p[i], t[i]);

        println!("Steam state in region 2 ?                = {}", valid_state);
        println!("Steam state in region 2a ?               = {}", in2a);
        println!("Steam state in region 2b ?               = {}", in2b);
        println!("Steam state in region 2c ?               = {}", in2c);

        let value = speed_of_sound(p[i], t[i]);
        let error = (value - w_real[i]).abs() / w_real[i].abs();
        println!("Error in speed of sound                  = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_volume();
        let error = (value - v_real[i]).abs() / v_real[i].abs();
        println!("Error in specific volume                 = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_enthalpy();
        let error = (value - h_real[i]).abs() / h_real[i].abs();
        println!("Error in specific enthalpy               = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_entropy();
        let error = (value - s_real[i]).abs() / s_real[i].abs();
        println!("Error in specific entropy                = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_internal_energy();
        let error = (value - u_real[i]).abs() / u_real[i].abs();
        println!("Error in specific internal energy        = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_isobaric_heat_capacity();
        let error = (value - cp_real[i]).abs() / cp_real[i].abs();
        println!("Error in specific isobaric heat capacity = {}", error);
        assert!(error <= 1e-2);

        let (t_ph, t_ps, p_hs) = if in2a {
            (
                temperature_ph_region2a(p[i], h_real[i]),
                temperature_ps_region2a(p[i], s_real[i]),
                pressure_hs_region2a(h_real[i], s_real[i]),
            )
        } else if in2b {
            (
                temperature_ph_region2b(p[i], h_real[i]),
                temperature_ps_region2b(p[i], s_real[i]),
                pressure_hs_region2b(h_real[i], s_real[i]),
            )
        } else {
            (
                temperature_ph_region2c(p[i], h_real[i]),
                temperature_ps_region2c(p[i], s_real[i]),
                pressure_hs_region2c(h_real[i], s_real[i]),
            )
        };

        let error = (t[i] - t_ph).abs() / t_ph.abs();
        println!("Error in backward temperature(p,h)       = {}", error);
        assert!(error <= 1e-2);

        let error = (t[i] - t_ps).abs() / t_ps.abs();
        println!("Error in backward temperature(p,s)       = {}", error);
        assert!(error <= 1e-2);

        let error = (p[i] - p_hs).abs() / p_hs.abs();
        println!("Error in backward pressure(h,s)          = {}", error);
        assert!(error <= 1e-2);
    }
}

#[test]
fn test_region3() {
    //
    // We are only working with region 3 in this test.
    //
    use steam::steam_region3::*;

    println!("\nRegion 3 Test - Computational Errors");
    println!("====================================");
    //
    // Input data.
    //
    let rho = [500.0, 200.0, 500.0];
    let t = [650.0, 650.0, 750.0];
    //
    // Expected values.
    //
    let p_real = [0.255837018e2, 0.222930643e2, 0.783095639e2];
    let h_real = [0.186343019e7, 0.237512401e7, 0.225868845e7];
    let u_real = [0.181226279e7, 0.226365868e7, 0.210206932e7];
    let s_real = [0.405427273e4, 0.485438792e4, 0.446971906e4];
    let cp_real = [0.138935717e5, 0.446579342e5, 0.634165359e4];
    let w_real = [0.502005554e3, 0.383444594e3, 0.760696041e3];

    for i in 0..3 {
        println!("============");
        println!("Test Case : {}", i + 1);
        println!("============");

        let mut ss = SteamState::new();
        let valid_state = ss.set_primary_properties(3, rho[i], t[i]);

        println!("Steam state in region 3 ?                = {}", valid_state);

        let in3a = in_region3a(rho[i], t[i]);
        let in3b = in_region3b(rho[i], t[i]);

        println!("Steam state in region 3a ?               = {}", in3a);
        println!("Steam state in region 3b ?               = {}", in3b);

        let value = speed_of_sound(rho[i], t[i]);
        let error = (value - w_real[i]).abs() / w_real[i].abs();
        println!("Error in speed of sound                  = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_pressure();
        let error = (value - p_real[i]).abs() / p_real[i].abs();
        println!("Error in pressure                        = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_enthalpy();
        let error = (value - h_real[i]).abs() / h_real[i].abs();
        println!("Error in specific enthalpy               = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_entropy();
        let error = (value - s_real[i]).abs() / s_real[i].abs();
        println!("Error in specific entropy                = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_internal_energy();
        let error = (value - u_real[i]).abs() / u_real[i].abs();
        println!("Error in specific internal energy        = {}", error);
        assert!(error <= 1e-2);

        let value = ss.get_specific_isobaric_heat_capacity();
        let error = (value - cp_real[i]).abs() / cp_real[i].abs();
        println!("Error in specific isobaric heat capacity = {}", error);
        assert!(error <= 1e-2);

        let (t_ph, t_ps, v_ph, v_ps) = if in3a {
            (
                temperature_ph_region3a(p_real[i], h_real[i]),
                temperature_ps_region3a(p_real[i], s_real[i]),
                specific_volume_ph_region3a(p_real[i], h_real[i]),
                specific_volume_ps_region3a(p_real[i], s_real[i]),
            )
        } else {
            (
                temperature_ph_region3b(p_real[i], h_real[i]),
                temperature_ps_region3b(p_real[i], s_real[i]),
                specific_volume_ph_region3b(p_real[i], h_real[i]),
                specific_volume_ps_region3b(p_real[i], s_real[i]),
            )
        };

        let error = (t[i] - t_ph).abs() / t_ph.abs();
        println!("Error in backward temperature(p,h)       = {}", error);
        assert!(error <= 1e-2);

        let error = (t[i] - t_ps).abs() / t_ps.abs();
        println!("Error in backward temperature(p,s)       = {}", error);
        assert!(error <= 1e-2);

        let error = (1.0 / rho[i] - v_ph).abs() / v_ph.abs();
        println!("Error in backward specific volume(p,h)   = {}", error);
        assert!(error <= 1e-2);

        let error = (1.0 / rho[i] - v_ps).abs() / v_ps.abs();
        println!("Error in backward specific volume(p,s)   = {}", error);
        assert!(error <= 1e-2);

        println!("\nAttempting to compute temperature by inverting polynomial h = h(rho,t)");
        let value = function_inverter_y(
            &specific_enthalpy,
            h_real[i],
            rho[i],
            500.0,
            800.0,
            1e-5,
            100,
            3,
            true,
        );
        assert_relative_eq!(value, t[i], max_relative = 1e-4);

        println!("\nAttempting to compute density by inverting polynomial p = p(rho,t)");
        let value =
            function_inverter_x(&pressure, p_real[i], t[i], 100.0, 600.0, 1e-5, 100, 3, true);
        assert_relative_eq!(value, rho[i], max_relative = 1e-4);
    }
}

#[test]
fn test_region4() {
    //
    // We are only working with region 4 in this test.
    //
    use steam::steam_region4::*;

    println!("\nRegion 4 Test - Computational Errors");
    println!("====================================");
    //
    // Test case A - Test saturation pressure.
    //
    let p_a = [0.353658941e-2, 0.263889776e1, 0.123443146e2];
    let t_a = [300.0, 500.0, 600.0];

    for i in 0..3 {
        println!("=============");
        println!("Test Case A : {}", i + 1);
        println!("=============");

        let mut ss = SteamState::new();
        let valid_state = ss.set_primary_properties(4, t_a[i], 0.5);

        println!("Steam state in region 4 ?                = {}", valid_state);

        let psat = saturation_pressure(t_a[i]);
        let error = (psat - p_a[i]).abs() / p_a[i];
        println!("Error in saturation pressure             = {}", error);
        assert!(error <= 1e-2);
    }
    //
    // Test case B - Test saturation temperature.
    //
    let p_b = [0.1, 1.0, 10.0];
    let t_b = [0.372755919e3, 0.453035632e3, 0.584149488e3];

    for i in 0..3 {
        println!("=============");
        println!("Test Case B : {}", i + 1);
        println!("=============");

        let mut ss = SteamState::new();
        let valid_state = ss.set_primary_properties(4, t_b[i], 0.5);

        println!("Steam state in region 4 ?                = {}", valid_state);

        let tsat = saturation_temperature(p_b[i]);
        let error = (tsat - t_b[i]).abs() / t_b[i];
        println!("Error in saturation temperature          = {}", error);
        assert!(error <= 1e-2);
    }
}

#[test]
fn test_boundary23() {
    //
    // We are only working with steam boundaries in this test.
    //
    use steam::steam_boundaries::*;

    println!("\nBoundary 2-3 Test - Computational Errors");
    println!("========================================");

    let p_real = 0.165291643e2;
    let t_real = 0.623150000e3;

    let p = boundary23_pressure(t_real);
    let error = (p - p_real).abs() / p_real;
    println!("Error in pressure                          = {}", error);
    assert!(error <= 1e-2);

    let t = boundary23_temperature(p_real);
    let error = (t - t_real).abs() / t_real;
    println!("Error in temperature                       = {}", error);
    assert!(error <= 1e-2);
}

#[test]
fn test_boundary2bc() {
    //
    // We are only working with steam boundaries in this test.
    //
    use steam::steam_boundaries::*;

    println!("\nBoundary 2b-2c Test - Computational Errors");
    println!("==========================================");

    let p_real = 0.1000000000e3;
    let h_real = 0.3516004323e7;

    let p = boundary2bc_pressure(h_real);
    let error = (p - p_real).abs() / p_real;
    println!("Error in pressure                          = {}", error);
    assert!(error <= 1e-2);

    let h = boundary2bc_enthalpy(p_real);
    let error = (h - h_real).abs() / h_real;
    println!("Error in specific enthalpy                 = {}", error);
    assert!(error <= 1e-2);
}

#[test]
fn test_boundary2ab() {
    //
    // We are only working with steam boundaries in this test.
    //
    use steam::steam_boundaries::*;

    println!("\nBoundary 2a-2b Test - Computational Errors");
    println!("==========================================");

    let s_real = 7000.0;
    let h_real = 3376437.884;

    let h = boundary2ab_enthalpy(s_real);
    let error = (h - h_real).abs() / h_real;
    println!("Error in specific enthalpy                 = {}", error);
    assert!(error <= 1e-2);
}

#[test]
fn test_boundary3ab() {
    //
    // We are only working with steam boundaries in this test.
    //
    use steam::steam_boundaries::*;

    println!("\nBoundary 3a-3b Test - Computational Errors");
    println!("==========================================");

    let p_real = 25.0;
    let h_real = 2.095936454e6;

    let h = boundary3ab_enthalpy(p_real);
    let error = (h - h_real).abs() / h_real;
    println!("Error in specific enthalpy                 = {}", error);
    assert!(error <= 1e-2);
}

#[test]
fn test_cubic_interpolation() {
    //
    // We are only working with region 4 in this test.
    //
    use steam::steam_region4::*;

    println!("\nCubic Interpolation - Computational Errors");
    println!("==========================================");

    let imax = 1000;
    let tmin: f64 = 273.15;
    let tmax: f64 = 647.096;
    let dt: f64 = (tmax - tmin) / (imax - 1) as f64;
    //
    // Compute the table data.
    //
    let mut tsat = vec![0.0; imax];
    let mut psat = vec![0.0; imax];
    let mut pprime = vec![0.0; imax];

    for i in 0..imax {
        tsat[i] = tmin + dt * i as f64;
        psat[i] = saturation_pressure(tsat[i]);
        pprime[i] = saturation_pressure_gradient(tsat[i]);
    }
    //
    // Produce interpolated data at random locations.
    //
    let num = 10;
    //
    // Create a new random seed.
    //
    let mut rng = rand::rng();

    for _ in 0..num {
        let rf = rng.random::<f64>();
        let t = tmin + (tmax - tmin) * rf;

        let iloc = (((t - tsat[0]) / dt) as usize).min(imax - 2);
        let iloc_searched = binary_search_vector(&tsat, t, 0, imax as i32 - 1);

        let p_real = saturation_pressure(t);
        let p_int = interpolate_cubic(
            &[tsat[iloc], tsat[iloc + 1]],
            &[psat[iloc], psat[iloc + 1]],
            &[pprime[iloc], pprime[iloc + 1]],
            t,
        );

        let error: f64 = 100.0 * (p_int - p_real).abs() / p_real;
        assert!(error <= 1e-2);
        println!("Polynomial Pressure    = {} [MPa]", p_real);

        println!("Cubic with Exact Gradients:");
        println!("Interpolated Pressure  = {} [MPa]", p_int);
        println!("Interpolation Error    = {} [%]", error);

        let mut ptab: [f64; 4] = [0.0; 4];
        let mut ttab: [f64; 4] = [0.0; 4];

        ptab[1] = psat[iloc];
        ptab[2] = psat[iloc + 1];
        ttab[1] = tsat[iloc];
        ttab[2] = tsat[iloc + 1];

        if iloc == 0 {
            ptab[0] = psat[iloc];
            ttab[0] = tsat[iloc];
        } else {
            ptab[0] = psat[iloc - 1];
            ttab[0] = tsat[iloc - 1];
        }

        if iloc == imax - 2 {
            ptab[3] = psat[iloc + 1];
            ttab[3] = tsat[iloc + 1];
        } else {
            ptab[3] = psat[iloc + 2];
            ttab[3] = tsat[iloc + 2];
        }

        let p_int = interpolate_cubic_numerical(&ttab, &ptab, t);

        let error: f64 = 100.0 * (p_int - p_real).abs() / p_real;
        assert!(error <= 1e-2);

        println!("Cubic with Numerical Gradients:");
        println!("Interpolated Pressure  = {} [MPa]", p_int);
        println!("Interpolation Error    = {} [%]", error);

        println!("Interpolation index    = {} {}", iloc, iloc_searched);
        println!("-------------------------");
    }
}

#[test]
fn test_bilinear_interpolation() {
    //
    // We are only working with region 5 in this test.
    //
    use steam::steam_region5::*;

    println!("\nBilinear Interpolation - Computational Errors");
    println!("=============================================");

    let imax = 100;
    let jmax = 100;

    let pmin = 0.1e-2;
    let pmax = REGION_5_PMAX;
    let dp = (pmax - pmin) / (imax - 1) as f64;

    let tmin = REGION_2_TMAX;
    let tmax = IAPWS97_TMAX;
    let dt = (tmax - tmin) / (jmax - 1) as f64;
    //
    // Compute the table data.
    //
    let mut p = vec![0.0; imax];
    let mut t = vec![0.0; jmax];
    let mut h = vec![vec![0.0; jmax]; imax];

    for i in 0..imax {
        p[i] = pmin + dp * i as f64;
    }

    for j in 0..jmax {
        t[j] = tmin + dt * j as f64;
    }

    for i in 0..imax {
        for j in 0..jmax {
            h[i][j] = specific_enthalpy(p[i], t[j]);
        }
    }
    //
    // Produce interpolated data at random locations.
    //
    let num = 10;

    for _ in 0..num {
        let rf = rand::random::<f64>();
        let tt = tmin + (tmax - tmin) * rf;

        let rf = rand::random::<f64>();
        let pp = pmin + (pmax - pmin) * rf;

        let iloc = (((pp - p[0]) / dp) as usize).min(imax - 2);
        let jloc = (((tt - t[0]) / dt) as usize).min(jmax - 2);

        let iloc_searched = binary_search_vector(&p, pp, 0, imax as i32 - 1);
        let jloc_searched = binary_search_vector(&t, tt, 0, jmax as i32 - 1);

        let p_points = [p[iloc], p[iloc + 1]];
        let t_points = [t[jloc], t[jloc + 1]];
        let h_points = [
            [h[iloc][jloc], h[iloc][jloc + 1]],
            [h[iloc + 1][jloc], h[iloc + 1][jloc + 1]],
        ];
        let h_int = interpolate_bilinear(&p_points, &t_points, &h_points, pp, tt);

        let h_real = specific_enthalpy(pp, tt);
        let error = (h_int - h_real).abs() / h_real.abs();
        assert!(error <= 1e-2);

        println!("Polynomial Specific Enthalpy    = {} [J/Kg]", h_real);
        println!("Interpolation Specific Enthalpy = {} [J/Kg]", h_int);
        println!("Estimation Error                = {} [%]", error);
        println!(
            "Interpolation Indices           = ({},{}) - ({},{})",
            iloc, jloc, iloc_searched, jloc_searched
        );
        println!("----------------------------------");
    }
}

///
/// This function tests the generation of a table for the computation
/// of the specific enthalpy in region 5.
///
#[test]
fn test_table_generation() {
    use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
    //
    // We are only working with region 5 in this test.
    //
    use steam::steam_region5::*;
    use steam::steam_tables::*;

    println!();
    println!("Table Generation Test : Specific Enthalpy in R5");
    println!("===============================================");

    let error = 1.0e-05;

    let xmin = 0.1e-02;
    let xmax = REGION_5_PMAX;
    let ymin = REGION_2_TMAX;
    let ymax = IAPWS97_TMAX;
    //
    // Generate table with bilinear interpolation.
    //
    let (x, y, h) = generate_table(
        &specific_enthalpy,
        xmin,
        xmax,
        ymin,
        ymax,
        0.025,
        0.50,
        error,
        true,
    );
    //
    // Verify the accuracy of the table.
    //
    println!();
    println!("Table Verification - Interpolation Errors");
    println!("=========================================");
    println!("Testing 10E+06 random combinations of properties");

    let passed = AtomicBool::new(true);
    let num_errors = AtomicUsize::new(0);

    let imax = x.len();
    let jmax = y.len();
    //
    // Use Rayon for parallel execution.
    //
    (0..10_000_000).into_par_iter().for_each(|_| {
        //
        // Create a thread-local random number generator.
        //
        let mut thread_rng = rand::rng();
        let xp = xmin + (xmax - xmin) * thread_rng.random::<f64>();
        let yp = ymin + (ymax - ymin) * thread_rng.random::<f64>();

        let h_real = specific_enthalpy(xp, yp);

        let iloc = binary_search_vector(&x, xp, 0, (imax - 1) as i32);
        let jloc = binary_search_vector(&y, yp, 0, (jmax - 1) as i32);
        //
        // Create arrays for the interpolation.
        //
        let x_interp = [x[iloc as usize], x[(iloc + 1) as usize]];
        let y_interp = [y[jloc as usize], y[(jloc + 1) as usize]];
        let h_interp = [
            [
                h[iloc as usize][jloc as usize],
                h[iloc as usize][(jloc + 1) as usize],
            ],
            [
                h[(iloc + 1) as usize][jloc as usize],
                h[(iloc + 1) as usize][(jloc + 1) as usize],
            ],
        ];

        let h_int = interpolate_bilinear(&x_interp, &y_interp, &h_interp, xp, yp);

        let er = (h_int - h_real).abs() / h_real.abs();

        if er > error {
            passed.store(false, Ordering::Relaxed);
            num_errors.fetch_add(1, Ordering::Relaxed);
        }
    });

    if passed.load(Ordering::Relaxed) {
        println!("All computed errors are within the specified limit.");
        println!();
    } else {
        println!(
            "{} error(s) have been detected.",
            num_errors.load(Ordering::Relaxed)
        );
        println!();
    }

    assert!(passed.load(Ordering::Relaxed));
}

///
/// This function tests the computation of viscosity against known reference values.
///
#[test]
pub fn test_viscosity() {
    use steam::steam_viscosity::*;

    println!();
    println!("Viscosity Test - Computational Errors");
    println!("==========================================");
    //
    // Reference temperature values (K).
    //
    let t: [f64; 11] = [
        298.15, 298.15, 373.15, 433.15, 433.15, 873.15, 873.15, 873.15, 1173.15, 1173.15, 1173.15,
    ];
    //
    // Reference density values (kg/m³).
    //
    let rho: [f64; 11] = [
        998.0, 1200.0, 1000.0, 1.0, 1000.0, 1.0, 100.0, 600.0, 1.0, 100.0, 400.0,
    ];
    //
    // Reference viscosity values (Pa·s).
    //
    let mu: [f64; 11] = [
        889.735100e-06,
        1437.649467e-06,
        307.883622e-06,
        14.538324e-06,
        217.685358e-06,
        32.619287e-06,
        35.802262e-06,
        77.430195e-06,
        44.217245e-06,
        47.640433e-06,
        64.154608e-06,
    ];
    //
    // Test the computation of viscosity against reference values.
    //
    for i in 0..11 {
        let visc = viscosity(rho[i], t[i]);
        let error = (visc - mu[i]).abs() / mu[i];
        println!("Error in viscosity                         = {}", error);
        assert!(error <= 1e-2);
    }
}

#[test]
fn test_thermal_conductivity() {
    use steam::steam_conductivity::*;
    //
    // Reference temperature values (K).
    //
    const T: [f64; 4] = [298.15, 298.15, 298.15, 873.15];
    //
    // Reference density values (kg/m³).
    //
    const RHO: [f64; 4] = [0.0, 998.0, 1200.0, 0.0];
    //
    // Reference thermal conductivity values (W/m·K).
    //
    const K: [f64; 4] = [
        0.018310320893515347,
        0.6086385449226368,
        0.855468383940539,
        0.07984060251109305,
    ];

    println!();
    println!("Thermal Conductivity Test - Computational Errors");
    println!("================================================");
    //
    // Test each reference value.
    //
    for i in 0..4 {
        let hk = thermal_conductivity(RHO[i], T[i]);
        let error = (hk - K[i]).abs() / K[i];
        println!("Error in thermal conductivity              = {}", error);
        assert!(error <= 1e-2);
    }
}
