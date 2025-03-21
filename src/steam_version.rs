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

//! Steam Version Module
//!
//! This module defines the version of the Steam Polynomials library

///
/// Major version number.
///
pub const MAJOR_VERSION: i32 = 2;
///
/// Minor version number.
///
pub const MINOR_VERSION: i32 = 0;
///
/// Bug fix number
///
pub const BUGFIX_VERSION: i32 = 0;
///
/// Version status
///
pub const CODE_STATUS: &str = "release";
///
/// Returns the full version string
///
pub fn version_string() -> String {
    format!(
        "{}.{}.{}-{}",
        MAJOR_VERSION, MINOR_VERSION, BUGFIX_VERSION, CODE_STATUS
    )
}
