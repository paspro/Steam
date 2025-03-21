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

//! Steam Constants Module
//!
//! This module defines constants required for the computation of the steam
//! polynomials according to the revised release on the IAPWS Industrial
//! Formulation 1997 for the Thermodynamic Properties of Water and Steam,
//! August 2007 (IAPWS-IF97)

///
/// Zero constant (0.0).
///
pub const ZERO: f64 = 0.0;
///
/// One constant (1.0).
///
pub const ONE: f64 = 1.0;
///
/// Two constant (2.0).
/// 
pub const TWO: f64 = 2.0;
///
/// Three constant (3.0).
///
pub const THREE: f64 = 3.0;
///
/// Four constant (4.0).
/// 
pub const FOUR: f64 = 4.0;
///
/// Five constant (5.0).
/// 
pub const FIVE: f64 = 5.0;
///
/// One half constant (0.5).
///
pub const HALF: f64 = 0.5;
///
/// One quarter constant (0.25).
/// 
pub const QUARTER: f64 = 0.25;
///
/// One third constant (1.0 / 3.0).
/// 
pub const ONETHIRD: f64 = 1.0 / 3.0;
///
/// One sixth constant (1.0 / 6.0).
///
pub const ONESIXTH: f64 = 1.0 / 6.0;

//-----------------------------------------------------------------------
//
// Steam Constants
//
//-----------------------------------------------------------------------

///
/// Steam constant R in [J/Kg.K].
///
pub const IAPWS97_R: f64 = 461.526;
///
/// Steam maximum pressure in [MPa].
///
pub const IAPWS97_PMAX: f64 = 100.0;
///
/// Steam minimum temperature in [K].
///
pub const IAPWS97_TMIN: f64 = 273.15;
///
/// Steam maximum temperature in [K].
///
pub const IAPWS97_TMAX: f64 = 2273.15;
///
/// Steam critical temperature in [K].
///
pub const IAPWS97_TCRIT: f64 = 647.096;
///
/// Steam critical pressure in [MPa].
///
pub const IAPWS97_PCRIT: f64 = 22.064;
///
/// Steam critical density in [Kg/m3].
///
pub const IAPWS97_RHOCRIT: f64 = 322.0;

//-----------------------------------------------------------------------
//
// Constants for Region 1.
//
//-----------------------------------------------------------------------

///
/// Maximum temperature for region 1 in [K].
///
pub const REGION_1_TMAX: f64 = 623.15;
///
/// The value of pressure at the point of maximum temperature
/// that intersects the line that represents region 4 in [MPa].
///
pub const REGION_1_4_TMAX_P: f64 = 16.5292;

//-----------------------------------------------------------------------
//
// Constants for Region 2.
//
//-----------------------------------------------------------------------

///
/// Steam maximum temperature in region 2 in [K].
///
pub const REGION_2_TMAX: f64 = 1073.15;
///
/// Temperature at the point where the line that represents region 4
/// intersects region 2 in [K]
///
pub const REGION_2_4_T: f64 = 863.15;
///
/// Pressure at the boundary between the subregions 2a and 2b
/// in [MPa].
///
pub const REGION_2A_2B_P: f64 = 4.0;
///
/// Minimum pressure at the intersection of regions 2b and 2c in [MPa].
///
pub const REGION_2B_2C_PMIN: f64 = 6.5467;

//-----------------------------------------------------------------------
//
// Constants for Region 3.
//
//-----------------------------------------------------------------------

///
/// The entropy at the boundary of regions 3a and 3b in [J/Kg.K].
///
pub const REGION_3A_3B_S: f64 = 4.41202148223476e3;
///
/// The minimum enthalpy at the boundary of regions 3 and 4 in [J/Kg].
///
pub const REGION_3_4_HMIN: f64 = 1.670858218e6;
///
/// The maximum enthalpy at the boundary of regions 3 and 4 in [J/Kg].
///
pub const REGION_3_4_HMAX: f64 = 2.563592004e6;
///
/// The minimum entropy at the boundary of regions 3 and 4 in [J/Kg.K].
///
pub const REGION_3_4_SMIN: f64 = 3.778281340e6;
///
/// The maximum entropy at the boundary of regions 3 and 4 in [J/Kg.K].
///
pub const REGION_3_4_SMAX: f64 = 5.210887825e6;

//-----------------------------------------------------------------------
//
// Constants for Region 4.
//
//-----------------------------------------------------------------------

///
/// Maximum temperature for region 4 in [K].
///
pub const REGION_4_TMAX: f64 = 647.096;

//-----------------------------------------------------------------------
//
// Constants for Region 5.
//
//-----------------------------------------------------------------------

///
/// Maximum pressure in region 5 in [MPa].
///
pub const REGION_5_PMAX: f64 = 50.0;
