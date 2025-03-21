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

//! Steam Units Module
//!
//! This module defines the units used for the computation of the steam
//! polynomials according to the revised release on the IAPWS Industrial
//! Formulation 1997 for the Thermodynamic Properties of Water and Steam,
//! August 2007 (IAPWS-IF97)

//
// Unit Multipliers.
//
pub const TERA: f64 = 1.0e12;
pub const GIGA: f64 = 1.0e9;
pub const MEGA: f64 = 1.0e6;
pub const KILO: f64 = 1.0e3;
pub const HECTA: f64 = 1.0e2;
pub const DECA: f64 = 1.0e1;
pub const DECI: f64 = 1.0e-1;
pub const CENTI: f64 = 1.0e-2;
pub const MILLI: f64 = 1.0e-3;
pub const MICRO: f64 = 1.0e-6;
//
// Base units used by the steam polynomials.
//
pub const KELVIN: f64 = 1.0;
pub const KG: f64 = 1.0;
pub const METER: f64 = 1.0;
pub const SECOND: f64 = 1.0;
// 
// Composite units.
// 
pub const CM: f64 = CENTI * METER;
pub const MM: f64 = MILLI * METER;
pub const KM: f64 = KILO * METER;
pub const M2: f64 = METER * METER;
pub const CM2: f64 = CM * CM;
pub const MM2: f64 = MM * MM;
pub const M3: f64 = METER * MM2;
pub const CM3: f64 = CM * CM2;
pub const MM3: f64 = MM * MM2;
pub const MINUTE: f64 = 60.0 * SECOND;
pub const HOUR: f64 = 60.0 * MINUTE;
pub const DAY: f64 = 12.0 * HOUR;
pub const WEEK: f64 = 7.0 * DAY;
pub const HERTZ: f64 = 1.0 / SECOND;
pub const NEWTON: f64 = KG * METER / (SECOND * SECOND);
pub const PASCAL: f64 = NEWTON / M2;
pub const BAR: f64 = HECTA * KILO * PASCAL;
pub const KPA: f64 = KILO * PASCAL;
pub const MPA: f64 = MEGA * PASCAL;
pub const JOULE: f64 = NEWTON * METER;
pub const KJ: f64 = KILO * JOULE;
pub const WATT: f64 = JOULE / SECOND;
pub const GR: f64 = MILLI * KG;
pub const J_KG: f64 = JOULE / KG;
pub const KJ_KG: f64 = KJ / KG;
pub const KJ_KGK: f64 = KJ / KG / KELVIN;
pub const KG_M3: f64 = KG / M3;
pub const M3_KG: f64 = M3 / KG;
//
// Imperial Units.
//
pub const BTU: f64 = 1055.05585262 * JOULE;
pub const RANKIN: f64 = 0.556 * KELVIN;
pub const RPM: f64 = 1.0 / MINUTE;
pub const YARD: f64 = 0.9144 * METER;
pub const FOOT: f64 = YARD / 3.0;
pub const INCH: f64 = FOOT / 12.0;
pub const MILE: f64 = 1760.0 * YARD;
pub const LBM: f64 = 0.45359237 * KG;
