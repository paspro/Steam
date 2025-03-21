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

//! Steam Tables Module
//!
//! This module contains functions capable of generating steam tables
//! for thermodynamic properties at a specific level of accuracy.

use crate::steam_numerics::interpolate_bilinear;

///
/// This function creates a table suitable for using in order to obtain
/// interpolated properties that depend on two independent variables. The
/// table is constructed such as it results in interpolated values of a
/// specific maximum allowed error.
///
/// - Arguments:
///   - `poly`: A function pointer to the polynomial function that
///             generates the accurate, reference values to be used for
///             building the table. This should be of the form:
///             property = poly(x,y) where x, y are the independent variables
///   - `xmin`: The minimum value for the variable x to use
///   - `xmax`: The maximum value for the variable x to use
///   - `ymin`: The minimum value for the variable y to use
///   - `ymax`: The maximum value for the variable y to use
///   - `dx`: The initial fine step to use for variable x
///   - `dy`: The initial fine step to use for variable y
///   - `error`: The maximum allowed interpolation error (%).
///   - `logging`: If true then information will be printed during
///                the execution of the function
///
/// - Returns:
///   - A tuple containing:
///     - `x`: The vector with the values of the variable x forming the constructed table
///     - `y`: The vector with the values of the variable y forming the constructed table
///     - `val`: The 2D vector with the values of the property the supplied polynomial
///              function computes
///
#[allow(clippy::too_many_arguments)]
pub fn generate_table(
    poly: &dyn Fn(f64, f64) -> f64,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    dx: f64,
    dy: f64,
    error: f64,
    logging: bool,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    //
    // Compute the number of points to use for the table given the initial dx, dy values.
    //
    let imax = ((xmax - xmin) / dx).ceil() as usize + 1;
    let jmax = ((ymax - ymin) / dy).ceil() as usize + 1;
    let points = imax * jmax;

    if logging {
        println!();
        println!("====================================");
        println!("Dynamic Table Generation Algorithm");
        println!("====================================");
        println!("Generating the initial table.");
        println!("Number of columns       = {}", imax);
        println!("Number of rows          = {}", jmax);
        println!("Total table data points = {}", imax * jmax);
        println!("====================================");
        println!();
    }
    //
    // Generate the fine table.
    //
    let mut x_fine = vec![0.0; imax];
    let mut y_fine = vec![0.0; jmax];
    let mut val_fine = vec![vec![0.0; jmax]; imax];
    let mut i_used = vec![true; imax];
    let mut j_used = vec![true; jmax];
    //
    // Fill the fine table with values.
    //
    for i in 0..imax {
        x_fine[i] = xmin + (i as f64) * dx;
    }

    for j in 0..jmax {
        y_fine[j] = ymin + (j as f64) * dy;
    }

    for i in 0..imax {
        for j in 0..jmax {
            val_fine[i][j] = poly(x_fine[i], y_fine[j]);
        }
    }
    //
    // Check the initial table to verify that the interpolation
    // errors for this table are bounded by the maximum allowed
    // interpolation error. If this is not the case then the
    // optimisation process makes no sense.
    //
    println!("Checking the interpolation error for the initial table...");

    for i in 0..(imax - 1) {
        let xp = [x_fine[i], x_fine[i + 1]];

        for j in 0..(jmax - 1) {
            let yp = [y_fine[j], y_fine[j + 1]];

            let up = [
                [val_fine[i][j], val_fine[i][j + 1]],
                [val_fine[i + 1][j], val_fine[i + 1][j + 1]],
            ];

            for ik in 0..=2 {
                let xk = xp[0] + 0.5 * (ik as f64) * dx;
                for jk in 0..=2 {
                    let yk = yp[0] + 0.5 * (jk as f64) * dy;
                    //
                    // Compute the accurate value of the property.
                    //
                    let prop_real = poly(xk, yk);
                    //
                    // Compute the interpolated value of the property.
                    //
                    let prop_int = interpolate_bilinear(&xp, &yp, &up, xk, yk);
                    //
                    // Estimate the error and check it.
                    //
                    let error_estimate = 100.0 * (prop_int - prop_real).abs() / prop_real.abs();

                    if error_estimate > error {
                        println!("Encountered unacceptable error:");
                        println!("Error = {}", error_estimate);
                        println!("Location (i,j) = {}, {}", i, j);
                        println!("and (ik,jk)    = {}, {}", ik, jk);
                        println!("Properties have the values:");
                        println!("xp = {}, {}", xp[0], xp[1]);
                        println!("yp = {}, {}", yp[0], yp[1]);
                        println!(
                            "up = {}, {}, {}, {}",
                            up[0][0], up[0][1], up[1][0], up[1][1]
                        );
                        println!("xk, yk = {}, {}", xk, yk);
                        println!("Polynomial property   = {}", prop_real);
                        println!("Interpolated property = {}", prop_int);
                        println!("Program terminates");
                        panic!("Unacceptable interpolation error");
                    }
                }
            }
        }
    }

    println!("Interpolation errors are bounded.");
    //
    // Check the i and j lines alternatively to see
    // if they can be removed and still have an interpolation
    // error within the acceptable limit.
    //
    if logging {
        println!("Performing table optimisation...");
    }
    //
    // Start with the initial i and j lines that
    // can be checked for removal.
    //
    let mut icheck: usize = 2;
    let mut jcheck: usize = 2;
    let mut check_i = true;
    let mut check_j = true;
    //
    // Check the i and j lines for removal.
    //
    while check_i || check_j {
        //
        // Check the i-line for removal.
        //
        if check_i {
            //
            // Find the previous non-deleted i-line.
            //
            let mut icheck_m1 = icheck - 1;
            while !i_used[icheck_m1] {
                icheck_m1 -= 1;
            }
            //
            // Find the next non-deleted i-line.
            //
            let mut icheck_p1 = icheck + 1;
            while icheck_p1 < imax && !i_used[icheck_p1] {
                icheck_p1 += 1;
            }
            //
            // Set the xp values.
            //
            let xp = [x_fine[icheck_m1], x_fine[icheck_p1]];
            let inum = ((xp[1] - xp[0]) / dx) as usize;
            //
            // Iterate the j-lines and find pairs of
            // non-deleted successive lines.
            //
            let mut j_p1: usize = 1;
            let mut can_remove_line = true;

            'row_check: while can_remove_line && j_p1 < jmax {
                let j = j_p1;
                //
                // If this j-line has been deleted move to the next one.
                //
                if !j_used[j] {
                    continue;
                }
                //
                // Find the next non-deleted j-line.
                //
                j_p1 = j + 1;
                while j_p1 < jmax && !j_used[j_p1] {
                    j_p1 += 1;
                }

                if j_p1 >= jmax {
                    break;
                }
                //
                // Check the interpolation errors for the
                // interval (icheck-1, j) -> (icheck+1, j+1).
                //
                let yp = [y_fine[j], y_fine[j_p1]];
                let jnum = ((yp[1] - yp[0]) / dy) as usize;

                let up = [
                    [val_fine[icheck_m1][j], val_fine[icheck_m1][j_p1]],
                    [val_fine[icheck_p1][j], val_fine[icheck_p1][j_p1]],
                ];

                for i in 0..=inum {
                    let xk = xp[0] + (i as f64) * dx;
                    for j in 0..=jnum {
                        let yk = yp[0] + (j as f64) * dy;
                        //
                        // Compute the accurate value of the property.
                        //
                        let prop_real = poly(xk, yk);
                        //
                        // Compute the interpolated value of the property.
                        //
                        let prop_int = interpolate_bilinear(&xp, &yp, &up, xk, yk);
                        //
                        // Estimate the error and check it.
                        //
                        let error_estimate = 100.0 * (prop_int - prop_real).abs() / prop_real.abs();

                        if error_estimate > error {
                            can_remove_line = false;
                            break 'row_check;
                        }
                    }
                }
            }
            //
            // If the line can be removed then mark it as such.
            //
            if can_remove_line {
                i_used[icheck] = false;
            }
        }
        //
        // Check the j-line for removal.
        //
        if check_j {
            //
            // Find the previous non-deleted j-line
            //
            let mut jcheck_m1 = jcheck - 1;
            while !j_used[jcheck_m1] {
                jcheck_m1 -= 1;
            }
            //
            // Find the next non-deleted j-line.
            //
            let mut jcheck_p1 = jcheck + 1;
            while jcheck_p1 < jmax && !j_used[jcheck_p1] {
                jcheck_p1 += 1;
            }
            //
            // Set the yp values.
            //
            let yp = [y_fine[jcheck_m1], y_fine[jcheck_p1]];
            let jnum = ((yp[1] - yp[0]) / dy) as usize;
            //
            // Iterate the i-lines and find pairs of
            // non-deleted successive lines.
            //
            let mut i_p1: usize = 1;
            let mut can_remove_line = true;

            'column_check: while can_remove_line && i_p1 < imax {
                let i = i_p1;
                //
                // If this i-line has been deleted move to the next one.
                //
                if !i_used[i] {
                    continue;
                }
                //
                // Find the next non-deleted i-line.
                //
                i_p1 = i + 1;
                while i_p1 < imax && !i_used[i_p1] {
                    i_p1 += 1;
                }

                if i_p1 >= imax {
                    break;
                }
                //
                // Check the interpolation errors for the
                // interval (i, jcheck-1) -> (i+1, jcheck+1).
                //
                let xp = [x_fine[i], x_fine[i_p1]];
                let inum = ((xp[1] - xp[0]) / dx) as usize;

                let up = [
                    [val_fine[i][jcheck_m1], val_fine[i][jcheck_p1]],
                    [val_fine[i_p1][jcheck_m1], val_fine[i_p1][jcheck_p1]],
                ];

                for idx in 0..=inum {
                    let xk = xp[0] + (idx as f64) * dx;
                    for jdx in 0..=jnum {
                        let yk = yp[0] + (jdx as f64) * dy;
                        //
                        // Compute the accurate value of the property.
                        //
                        let prop_real = poly(xk, yk);
                        //
                        // Compute the interpolated value of the property.
                        //
                        let prop_int = interpolate_bilinear(&xp, &yp, &up, xk, yk);
                        //
                        // Estimate the error and check it.
                        //
                        let error_estimate = 100.0 * (prop_int - prop_real).abs() / prop_real.abs();

                        if error_estimate > error {
                            can_remove_line = false;
                            break 'column_check;
                        }
                    }
                }
            }
            //
            // If the line can be removed then mark it as such.
            //
            if can_remove_line {
                j_used[jcheck] = false;
            }
        }
        //
        // Move on to the next i and j lines to check for removal
        // that have not been already removed.
        //
        check_i = false;
        for i in (icheck + 1)..(imax - 1) {
            if i_used[i] {
                icheck = i;
                check_i = true;
                break;
            }
        }

        check_j = false;
        for j in (jcheck + 1)..(jmax - 1) {
            if j_used[j] {
                jcheck = j;
                check_j = true;
                break;
            }
        }
    }

    let imax_new = i_used.iter().filter(|&&used| used).count();
    let jmax_new = j_used.iter().filter(|&&used| used).count();
    let points_new = imax_new * jmax_new;

    if logging {
        println!();
        println!("====================================");
        println!("Table Optimisation Results");
        println!("====================================");
        println!("Eliminated {} columns", imax - imax_new);
        println!("Eliminated {} rows", jmax - jmax_new);
        println!("------------------------------------");
        println!("Optimised table has columns = {}", imax_new);
        println!("Optimised table has rows    = {}", jmax_new);
        println!("Total table data points     = {}", imax_new * jmax_new);
        println!(
            "Table size reduction        = {} %",
            100.0 * (points - points_new) as f64 / points as f64
        );
        println!("Maximum allowed error       = {} %", error);
        println!("====================================");
        println!();
    }
    //
    // Construct the table data.
    //
    let mut x = Vec::with_capacity(imax_new);
    let mut y = Vec::with_capacity(jmax_new);
    let mut val = vec![vec![0.0; jmax_new]; imax_new];

    let mut i_new = 0;
    for i in 0..imax {
        if i_used[i] {
            x.push(x_fine[i]);
            i_new += 1;
        }
    }

    let mut j_new = 0;
    for j in 0..jmax {
        if j_used[j] {
            y.push(y_fine[j]);
            j_new += 1;
        }
    }

    i_new = 0;
    for i in 0..imax {
        if i_used[i] {
            j_new = 0;
            for j in 0..jmax {
                if j_used[j] {
                    val[i_new][j_new] = val_fine[i][j];
                    j_new += 1;
                }
            }
            i_new += 1;
        }
    }

    (x, y, val)
}
