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

//! Steam Numerics Module
//!
//! This module provides numerical utilities for the steam calculations.

use crate::steam_constants::*;

///
/// This function reverses a generic equation f(x) by computing the value
/// of "x0" when the value of the equation f(x0) is known with the use of
/// an accelerated Secant root-finder algorithm operating on the function
/// g(x) = f(x) - f0, where f0 = f(x0) -> g(x0) = 0.
///
/// - Arguments:
///   - `f`: The function to invert.
///   - `f0`: The value of f0 = f(x0).
///   - `guess1`: A guess for the value of x0.
///   - `guess2`: A second guess for the value of x0.
///   - `tolerance`: The tolerance to utilise for the numerical method.
///   - `max_iterations`: The maximum number of allowed iterations.
///   - `n_order`: The interpolation order used by the method.
///   - `logging`: If true then write to the screen logging data.
///
/// - Returns:
///   - The value of "x0" when f(x0) is known.
///
#[allow(clippy::too_many_arguments)]
pub fn function_inverter(
    f: &dyn Fn(f64) -> f64,
    f0: f64,
    guess1: f64,
    guess2: f64,
    tolerance: f64,
    max_iterations: i32,
    n_order: i32,
    logging: bool,
) -> f64 {
    if n_order < 0 {
        eprintln!("\nError in Function: function_inverter");
        eprintln!("The order of interpolation must be a positive integer.");
        eprintln!();
    }
    if logging {
        println!();
        println!("=========================================");
        println!("Function f(x) Inverter Algorithm for x");
        println!("=========================================");
    }
    ///
    /// The solution context for the function inverter.
    ///
    struct SolutionContext<'a> {
        ///
        /// The function to invert.
        ///
        f: &'a dyn Fn(f64) -> f64,
        ///
        /// The value of f0 = f(x0).
        ///
        f0: f64,
        ///
        /// A guess for the value of x0.
        ///
        guess1: f64,
        ///
        /// A second guess for the value of x0.
        ///
        guess2: f64,
        ///
        /// The interpolation order used by the method.
        ///
        n_order: i32,
        ///
        /// The cache for the solution.
        ///
        cache: Vec<Vec<f64>>,
    }
    //
    // Implement the solution context.
    //
    impl<'a> SolutionContext<'a> {
        ///
        /// Create a new solution context.
        ///
        /// - Arguments:
        ///   - `f`: The function to invert.
        ///   - `f0`: The value of f0 = f(x0).
        ///   - `guess1`: A guess for the value of x0.
        ///   - `guess2`: A second guess for the value of x0.
        ///   - `n_order`: The interpolation order used by the method.
        ///   - `max_iterations`: The maximum number of allowed iterations.
        ///
        /// - Returns:
        ///   - The solution context.
        ///
        fn new(
            f: &'a dyn Fn(f64) -> f64,
            f0: f64,
            guess1: f64,
            guess2: f64,
            n_order: i32,
            max_iterations: i32,
        ) -> Self {
            //
            // Initialize cache with None values.
            //
            let mut cache = Vec::new();
            for _ in 0..(max_iterations + 1) {
                let row = vec![f64::NAN; (n_order + 1) as usize];
                cache.push(row);
            }
            Self {
                f,
                f0,
                guess1,
                guess2,
                n_order,
                cache,
            }
        }
        ///
        /// Compute the maximum value for the interpolation.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///
        /// - Returns:
        ///   - The maximum value for the interpolation.
        ///
        fn imax(&self, p: i32) -> i32 {
            if p == -1 || p == 0 {
                0
            } else if p <= self.n_order {
                p - 1
            } else {
                self.n_order
            }
        }
        ///
        /// Compute the solution for the function inverter.
        ///
        /// - Arguments:
        ///  - `p`: The value of p.
        ///  - `i`: The value of i.
        ///
        /// - Returns:
        ///   - The solution for the function inverter.
        ///
        fn solution(&mut self, p: i32, i: i32) -> f64 {
            //
            // Check cache first.
            //
            if p >= 0 && i >= 0 && !self.cache[p as usize][i as usize].is_nan() {
                return self.cache[p as usize][i as usize];
            }

            let result = if p == -1 {
                self.guess1
            } else if p == 0 {
                self.guess2
            } else if i == 0 {
                //
                // Secant step.
                //
                let x1 = self.solution(p - 1, self.imax(p - 1));
                let x2 = self.solution(p - 2, self.imax(p - 2));
                let f1 = (self.f)(x1) - self.f0;
                let f2 = (self.f)(x2) - self.f0;
                x1 - f1 * (x1 - x2) / (f1 - f2)
            } else {
                //
                // i-th order approximation.
                //
                let x1 = self.solution(p - 1, i - 1);
                let x2 = self.solution(p - 1, self.imax(p - 1));
                let x3 = self.solution(p, i - 1);
                let x4 = self.solution(p - i - 2, self.imax(p - i - 2));
                x3 + (x2 - x3) * (x1 - x3) / (x1 + x2 - x3 - x4)
            };
            //
            // Store in cache if possible.
            //
            if p >= 0 && i >= 0 {
                self.cache[p as usize][i as usize] = result;
            }

            result
        }
    }

    let mut context = SolutionContext::new(f, f0, guess1, guess2, n_order, max_iterations);
    let mut s_old = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        let s = context.solution(iter, n_order);
        let error = (s_old - s).abs() / s.abs();
        s_old = s;

        if logging {
            println!(" Iteration = {:<5} X = {:.5} Error = {:.5}", iter, s, error);
        }
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    if logging {
        println!("=========================================");
        println!();
    }

    s_old
}

///
/// This function reverses a generic equation f(x,y) by computing the value
/// of "x0" when the value of the equation f(x0,y0) and y0 are known using
/// an accelerated Secant root-finder algorithm operating on the function
/// g(x,y0) = f(x,y0) - f0, where f0 = f(x0,y0) -> g(x0,y0) = 0.
///
/// - Arguments:
///   - `f` - The function to invert.
///   - `f0` - The value of f0 = f(x0,y0).
///   - `y0` - The value y0.
///   - `guess1` - A guess for the value of x0.
///   - `guess2` - A second guess for the value of x0.
///   - `tolerance` - The tolerance to utilise for the numerical method.
///   - `max_iterations` - The maximum number of allowed iterations.
///   - `n_order` - The interpolation order used by the method.
///   - `logging` - If true then write to the screen logging data.
///
/// - Returns:
///   - The value of "x0" when f(x0,y0) is known.
///
#[allow(clippy::too_many_arguments)]
pub fn function_inverter_x(
    f: &dyn Fn(f64, f64) -> f64,
    f0: f64,
    y0: f64,
    guess1: f64,
    guess2: f64,
    tolerance: f64,
    max_iterations: i32,
    n_order: i32,
    logging: bool,
) -> f64 {
    if n_order < 0 {
        eprintln!("\nError in Function: function_inverter_x");
        eprintln!("The order of interpolation must be a positive integer.");
        eprintln!();
    }
    if logging {
        println!();
        println!("=========================================");
        println!("Function f(x,y0) Inverter Algorithm for x");
        println!("=========================================");
    }
    ///
    /// The solution context for the function inverter.
    ///
    struct SolutionContext<'a> {
        ///
        /// The function to invert.
        ///
        f: &'a dyn Fn(f64, f64) -> f64,
        ///
        /// The value of f0 = f(x0,y0).
        ///
        f0: f64,
        ///
        /// The value y0.
        ///
        y0: f64,
        ///
        /// A guess for the value of x0.
        ///     
        guess1: f64,
        ///
        /// A second guess for the value of x0.
        ///     
        guess2: f64,
        ///
        /// The interpolation order used by the method.
        ///
        n_order: i32,
        ///
        /// The cache for the solution.
        ///
        cache: Vec<Vec<f64>>,
    }
    //
    // Implement the solution context.
    //
    impl<'a> SolutionContext<'a> {
        ///
        /// Create a new solution context.
        ///
        /// - Arguments:
        ///   - `f`: The function to invert.
        ///   - `f0`: The value of f0 = f(x0,y0).
        ///   - `y0`: The value y0.
        ///   - `guess1`: A guess for the value of x0.
        ///   - `guess2`: A second guess for the value of x0.
        ///   - `n_order`: The interpolation order used by the method.
        ///   - `max_iterations`: The maximum number of allowed iterations.
        ///
        /// - Returns:
        ///   - The solution context.
        ///
        fn new(
            f: &'a dyn Fn(f64, f64) -> f64,
            f0: f64,
            y0: f64,
            guess1: f64,
            guess2: f64,
            n_order: i32,
            max_iterations: i32,
        ) -> Self {
            //
            // Initialize cache.
            //
            let mut cache = Vec::new();
            for _ in 0..(max_iterations + 1) {
                let row = vec![f64::NAN; (n_order + 1) as usize];
                cache.push(row);
            }
            Self {
                f,
                f0,
                y0,
                guess1,
                guess2,
                n_order,
                cache,
            }
        }
        ///
        /// Compute the maximum value for the interpolation.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///
        /// - Returns:
        ///   - The maximum value for the interpolation.
        ///
        fn imax(&self, p: i32) -> i32 {
            if p == -1 || p == 0 {
                0
            } else if p <= self.n_order {
                p - 1
            } else {
                self.n_order
            }
        }
        ///
        /// Compute the solution for the function inverter.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///   - `i`: The value of i.
        ///
        /// - Returns:
        ///   - The solution for the function inverter.
        ///
        fn solution(&mut self, p: i32, i: i32) -> f64 {
            //
            // Check cache first.
            //
            if p >= 0 && i >= 0 && !self.cache[p as usize][i as usize].is_nan() {
                return self.cache[p as usize][i as usize];
            }

            let result = if p == -1 {
                self.guess1
            } else if p == 0 {
                self.guess2
            } else if i == 0 {
                //
                // Secant step.
                //
                let x1 = self.solution(p - 1, self.imax(p - 1));
                let x2 = self.solution(p - 2, self.imax(p - 2));
                let f1 = (self.f)(x1, self.y0) - self.f0;
                let f2 = (self.f)(x2, self.y0) - self.f0;
                x1 - f1 * (x1 - x2) / (f1 - f2)
            } else {
                //
                // i-th order approximation.
                //
                let x1 = self.solution(p - 1, i - 1);
                let x2 = self.solution(p - 1, self.imax(p - 1));
                let x3 = self.solution(p, i - 1);
                let x4 = self.solution(p - i - 2, self.imax(p - i - 2));
                x3 + (x2 - x3) * (x1 - x3) / (x1 + x2 - x3 - x4)
            };

            //
            // Store in cache if possible.
            //
            if p >= 0 && i >= 0 {
                self.cache[p as usize][i as usize] = result;
            }

            result
        }
    }

    let mut context = SolutionContext::new(f, f0, y0, guess1, guess2, n_order, max_iterations);
    let mut s_old = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        let s = context.solution(iter, n_order);
        let error = (s_old - s).abs() / s.abs();
        s_old = s;

        if logging {
            println!(" Iteration = {:<5} X = {:.5} Error = {:.5}", iter, s, error);
        }
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    if logging {
        println!("=========================================");
        println!();
    }

    s_old
}

///
/// This function reverses a generic equation f(x,y) by computing the value
/// of "y0" when the value of the equation f(x0,y0) and x0 are known using
/// the numerical rootfinder algorithm of Newton-Raphson operating on the
/// function g(x0,y) = f(x0,y) - f0, where f0 = f(x0,y0) -> g(x0,y0) = 0.
///
/// - Arguments:
///   - `f`: The function to invert.
///   - `f0`: The value of f0 = f(x0,y0).
///   - `x0`: The value x0.
///   - `guess1`: A guess for the value of y0.
///   - `guess2`: A second guess for the value of y0.
///   - `tolerance`: The tolerance to utilise for the numerical method.
///   - `max_iterations`: The maximum number of allowed iterations.
///   - `n_order`: The interpolation order used by the method
///   - `logging`: If true then write to the screen logging data.
///
/// - Returns:
///   - The value of "y0" when f(x0,y0) is known.
///
#[allow(clippy::too_many_arguments)]
pub fn function_inverter_y(
    f: &dyn Fn(f64, f64) -> f64,
    f0: f64,
    x0: f64,
    guess1: f64,
    guess2: f64,
    tolerance: f64,
    max_iterations: i32,
    n_order: i32,
    logging: bool,
) -> f64 {
    if n_order < 0 {
        eprintln!("\nError in Function: function_inverter_y");
        eprintln!("The order of interpolation must be a positive integer.");
        eprintln!();
    }
    if logging {
        println!();
        println!("=========================================");
        println!("Function f(x0,y) Inverter Algorithm for y");
        println!("=========================================");
    }
    ///
    /// The solution context for the function inverter.
    ///
    struct SolutionContext<'a> {
        ///
        /// The function to invert.
        ///
        f: &'a dyn Fn(f64, f64) -> f64,
        ///
        /// The value of f0 = f(x0,y0).
        ///
        f0: f64,
        ///
        /// The value x0.
        ///
        x0: f64,
        ///
        /// A guess for the value of y0.
        ///     
        guess1: f64,
        ///
        /// A second guess for the value of y0.
        ///
        guess2: f64,
        ///
        /// The interpolation order used by the method.
        ///
        n_order: i32,
        ///
        /// The cache for the solution.
        ///
        cache: Vec<Vec<f64>>,
    }
    //
    // Implement the solution context.
    //
    impl<'a> SolutionContext<'a> {
        ///
        /// Create a new solution context.
        ///
        /// - Arguments:
        ///   - `f`: The function to invert.
        ///   - `f0`: The value of f0 = f(x0,y0).
        ///   - `x0`: The value x0.
        ///   - `guess1`: A guess for the value of y0.
        ///   - `guess2`: A second guess for the value of y0.
        ///   - `n_order`: The interpolation order used by the method.
        ///   - `max_iterations`: The maximum number of allowed iterations.
        ///
        fn new(
            f: &'a dyn Fn(f64, f64) -> f64,
            f0: f64,
            x0: f64,
            guess1: f64,
            guess2: f64,
            n_order: i32,
            max_iterations: i32,
        ) -> Self {
            //
            // Initialize cache.
            //
            let mut cache = Vec::new();
            for _ in 0..(max_iterations + 1) {
                let row = vec![f64::NAN; (n_order + 1) as usize];
                cache.push(row);
            }
            Self {
                f,
                f0,
                x0,
                guess1,
                guess2,
                n_order,
                cache,
            }
        }
        ///
        /// Compute the maximum value for the interpolation.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///
        /// - Returns:
        ///   - The maximum value for the interpolation.
        ///
        fn imax(&self, p: i32) -> i32 {
            if p == -1 || p == 0 {
                0
            } else if p <= self.n_order {
                p - 1
            } else {
                self.n_order
            }
        }
        ///
        /// Compute the solution for the function inverter.
        ///
        /// - Arguments:
        ///   - `p`: The value of p.
        ///   - `i`: The value of i.
        ///
        /// - Returns:
        ///   - The solution for the function inverter.
        ///
        fn solution(&mut self, p: i32, i: i32) -> f64 {
            //
            // Check cache first.
            //
            if p >= 0 && i >= 0 && !self.cache[p as usize][i as usize].is_nan() {
                return self.cache[p as usize][i as usize];
            }

            let result = if p == -1 {
                self.guess1
            } else if p == 0 {
                self.guess2
            } else if i == 0 {
                //
                // Secant step.
                //
                let y1 = self.solution(p - 1, self.imax(p - 1));
                let y2 = self.solution(p - 2, self.imax(p - 2));
                let f1 = (self.f)(self.x0, y1) - self.f0;
                let f2 = (self.f)(self.x0, y2) - self.f0;
                y1 - f1 * (y1 - y2) / (f1 - f2)
            } else {
                //
                // i-th order approximation.
                //
                let y1 = self.solution(p - 1, i - 1);
                let y2 = self.solution(p - 1, self.imax(p - 1));
                let y3 = self.solution(p, i - 1);
                let y4 = self.solution(p - i - 2, self.imax(p - i - 2));
                y3 + (y2 - y3) * (y1 - y3) / (y1 + y2 - y3 - y4)
            };
            //
            // Store in cache if possible.
            //
            if p >= 0 && i >= 0 {
                self.cache[p as usize][i as usize] = result;
            }

            result
        }
    }

    let mut context = SolutionContext::new(f, f0, x0, guess1, guess2, n_order, max_iterations);
    let mut s_old = guess2;

    for iter in (n_order + 1)..=max_iterations {
        //
        // Compute the new value.
        //
        let s = context.solution(iter, n_order);
        let error = (s_old - s).abs() / s.abs();
        s_old = s;

        if logging {
            println!(" Iteration = {:<5} Y = {:.5} Error = {:.5}", iter, s, error);
        }
        //
        // Check for convergence.
        //
        if error <= tolerance {
            break;
        }
    }

    if logging {
        println!("=========================================");
        println!();
    }

    s_old
}

///
/// This function performs a piecewise cubic interpolation in order to
/// compute the value of a property. The interpolation function has the
/// form:
///
/// s(x) = SUM(Cmj*(x-xj)^m), where the sum is over m = 0, 3
///        and Cmj are the interpolation coefficients computed based on
///        the supplied data in the form of vectors x, u(x), du(x)/dx.
///
/// - Arguments:
///   - `x`: The vector with the x-values (size=2).
///   - `u`: The vector with the property values (size=2).
///   - `uprime`: The vector with the gradient of the property.
///   - `xi`: The x-value to be used for computing the interpolated
///           property value which lies between x0 and x1.
///
/// - Returns:
///   - The interpolated value at location xi.
///
pub fn interpolate_cubic(x: &[f64; 2], u: &[f64; 2], uprime: &[f64; 2], xi: f64) -> f64 {
    //
    // Compute the interpolation coefficients.
    //
    let dx = x[1] - x[0];
    let du = u[1] - u[0];
    let dx2 = dx * dx;
    let dx3 = dx2 * dx;

    let coeff = [
        u[0],
        uprime[0],
        THREE * du / dx2 - (uprime[1] + TWO * uprime[0]) / dx,
        -TWO * du / dx3 + (uprime[1] + uprime[0]) / dx2,
    ];
    //
    // Perform the interpolation.
    //
    let x_diff = xi - x[0];
    coeff[0] + coeff[1] * x_diff + coeff[2] * x_diff.powi(2) + coeff[3] * x_diff.powi(3)
}

///
/// This function performs a cubic interpolation in order to compute
/// the value of a property.
///
/// - Arguments:
///   - `x`: The vector with the x-values (size=4).
///   - `u`: The vector with the property values (size=4).
///   - `xi`: The x-value to be used for computing the interpolated
///           property value which lies between x1 and x2.
///
/// - Returns:
///   - The interpolated value at location xi.
///
pub fn interpolate_cubic_numerical(x: &[f64; 4], u: &[f64; 4], xi: f64) -> f64 {
    //
    // Compute the interpolation coefficients.
    //
    let dx = [xi - x[0], xi - x[1], xi - x[2], xi - x[3]];

    let div = [
        (x[0] - x[1]) * (x[0] - x[2]) * (x[0] - x[3]),
        (x[1] - x[0]) * (x[1] - x[2]) * (x[1] - x[3]),
        (x[2] - x[0]) * (x[2] - x[1]) * (x[2] - x[3]),
        (x[3] - x[0]) * (x[3] - x[1]) * (x[3] - x[2]),
    ];

    let coeff = [
        dx[1] * dx[2] * dx[3] / div[0],
        dx[0] * dx[2] * dx[3] / div[1],
        dx[0] * dx[1] * dx[3] / div[2],
        dx[0] * dx[1] * dx[2] / div[3],
    ];
    //
    // Perform the interpolation.
    //
    u[0] * coeff[0] + u[1] * coeff[1] + u[2] * coeff[2] + u[3] * coeff[3]
}

///
/// This function performs a piecewise bilinear interpolation in order to
/// compute the value of a property that is a function of two parameters.
/// The interpolation function has the form:
///
/// s(x,y) = SUMn(SUMm(Cmnjk*(x-xj)^m*(y-yk)^n), where SUMm is over m = 0, 1,
///         SUMn is over n = 0, 1 and Cmnjk are the interpolation coefficients
///         computed based on the supplied data in the form of the table x, y,
///         u(x,y).
///
/// - Arguments:
///   - `x`: The vector with the x-values (dim=2).
///   - `y`: The vector with the y-values (dim=2).
///   - `u`: The array with the property values (dim=2x2).
///   - `xi`: The x-value to be used for computing the interpolated property value.
///   - `yi`: The y-value to be used for computing the interpolated property value.
///
/// - Returns:
///   - The interpolated value at location xi, yi.
///
pub fn interpolate_bilinear(
    x: &[f64; 2],
    y: &[f64; 2],
    u: &[[f64; 2]; 2],
    xi: f64,
    yi: f64,
) -> f64 {
    //
    // Compute the interpolation coefficients.
    //
    let dx21 = x[1] - x[0];
    let dy21 = y[1] - y[0];
    let dx2i = x[1] - xi;
    let dy2i = y[1] - yi;
    let dxi1 = xi - x[0];
    let dyi1 = yi - y[0];
    //
    // Perform the interpolation.
    //
    let u1d = (u[0][0] * dx2i + u[1][0] * dxi1) * dy2i;
    let u2d = (u[0][1] * dx2i + u[1][1] * dxi1) * dyi1;

    (u1d + u2d) / (dx21 * dy21)
}

///
/// This function performs a bicubic interpolation in order to
/// compute the value of a property that is a function of two parameters.
///
/// - Arguments:
///   - `x`: The vector with the x-values (dim=4).
///   - `y`: The vector with the y-values (dim=4).
///   - `u`: The array with the property values (dim=4x4).
///   - `xi`: The x-value to be used for computing the interpolated property value.
///   - `yi`: The y-value to be used for computing the interpolated property value.
///
/// - Returns:
///   - The interpolated value at location xi, yi.
///
pub fn interpolate_bicubic(x: &[f64; 4], y: &[f64; 4], u: &[[f64; 4]; 4], xi: f64, yi: f64) -> f64 {
    //
    // Perform a cubic interpolation for each column.
    //
    let mut uv = [0.0; 4];

    for i in 0..4 {
        let column = [u[i][0], u[i][1], u[i][2], u[i][3]];
        uv[i] = interpolate_cubic_numerical(y, &column, yi);
    }
    //
    // Perform a horizontal interpolation.
    //
    interpolate_cubic_numerical(x, &uv, xi)
}

///
/// This function implements a binary search algorithm in order to locate
/// a value in a vector. If the exact value does not exist it will return
/// the closest smaller value available in the table. Note that the
/// supplied vector must be sorted in ascending order.
///
/// - Arguments:
///   - `vec`: The vector to use for searching.
///   - `value`: The value to search for.
///   - `imin`: The minimum vector index to use for searching.
///   - `imax`: The maximum vector index to use for searching.
///
/// - Returns:
///   - The index with the location of the result.
///
pub fn binary_search_vector(vec: &[f64], value: f64, imin: i32, imax: i32) -> i32 {
    let mut high = imax;
    let mut low = imin;
    let mut med = imin + (imax - imin) / 2;

    while high != (low + 1) {
        let vmed = vec[med as usize];

        if vmed <= value {
            low = med;
        } else {
            high = med;
        }

        med = low + (high - low) / 2;
    }

    med
}

///
/// This function implements a binary search algorithm in order to locate
/// a value in a two dimensional array. If the exact value does not exist
/// it will return the closest smaller value available in the table.
///
/// - Arguments:
///   - `mat`: The matrix (2D) to use for searching.
///   - `value`: The value to search for.
///   - `imin`: The minimum array i-index to use for searching.
///   - `imax`: The maximum array i-index to use for searching.
///   - `jmin`: The minimum array j-index to use for searching.
///   - `jmax`: The maximum array j-index to use for searching.
///
/// - Returns:
///   - A tuple (i, j) with the indices for the location of the value.
///
pub fn binary_search_array(
    mat: &[Vec<f64>],
    value: f64,
    imin: i32,
    imax: i32,
    jmin: i32,
    jmax: i32,
) -> (i32, i32) {
    let mut j_index = binary_search_vector(&mat[imin as usize], value, jmin, jmax);
    let mut value_found = mat[imin as usize][j_index as usize];
    let mut error_previous = (value - value_found) / value_found.abs();

    let mut i = imin;
    let mut j = j_index;

    for i_index in (imin + 1)..=imax {
        j_index = binary_search_vector(&mat[i_index as usize], value, jmin, jmax);
        value_found = mat[i as usize][j as usize];
        let error_new = (value - value_found) / value_found.abs();

        if error_new >= 0.0 && error_new < error_previous {
            i = i_index;
            j = j_index;
        }

        error_previous = error_new;
    }

    (i, j)
}

///
/// This function implements the QuickSort algorithm for sorting the
/// elements of a vector to ascending order.
///
/// - Arguments:
///   - `vec`: The vector to sort.
///   - `indx`: A vector of indices that represents the modification
///             of the original vector. It can be used in order to
///             apply similar modifications to other associated vectors.
///
pub fn quicksort(vec: &mut [f64], indx: &mut Option<&mut [i32]>) {
    if vec.len() > 1 {
        let iq = partition(vec, indx);
        //
        // Sort the partitions.
        //
        if let Some(idx) = indx {
            let (left, right) = idx.split_at_mut(iq);
            quicksort(&mut vec[..iq], &mut Some(left));
            quicksort(&mut vec[iq..], &mut Some(right));
        } else {
            quicksort(&mut vec[..iq], &mut None);
            quicksort(&mut vec[iq..], &mut None);
        }
    }
}

///
/// This function partitions a vector for sorting purposes.
///
/// - Arguments:
///   - `vec`: The vector to partition.
///   - `indx`: An optional vector of indices to be partitioned in the same manner.
///
/// - Returns:
///   - The partition marker.
///
fn partition(vec: &mut [f64], indx: &mut Option<&mut [i32]>) -> usize {
    let x = vec[0];
    let mut i = 0;
    let mut j = vec.len();

    loop {
        j -= 1;
        while vec[j] > x {
            j -= 1;
        }

        i += 1;
        while i < j && vec[i] < x {
            i += 1;
        }

        if i < j {
            vec.swap(i, j);
            if let Some(idx) = indx {
                idx.swap(i, j);
            }
        } else {
            return if i == j { i + 1 } else { i };
        }
    }
}
