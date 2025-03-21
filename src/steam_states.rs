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

//! Steam States Module
//!
//! This module implements a steam state representation
//!
//! @author Panos Asproulis
//! @date 2014
//! @copyright Panos Asproulis (2014-2016). All Rights Reserved.

use crate::steam_conductivity;
use crate::steam_constants::{ONE, ZERO};
use crate::steam_region1;
use crate::steam_region2;
use crate::steam_region3;
use crate::steam_region4;
use crate::steam_region5;
use crate::steam_regions;
use crate::steam_viscosity;

///
/// This struct implements a steam state.
///
#[derive(Debug, Clone)]
pub struct SteamState {
    ///
    /// The region the steam state is located.
    ///
    region: usize,
    ///
    /// The steam pressure in [MPa].
    ///
    pressure: f64,
    ///
    /// The steam temperature in [K].
    ///
    temperature: f64,
    ///
    /// The steam density in [kg/m3].
    ///
    density: f64,
    ///
    /// The steam quality in [%].
    ///
    quality: f64,
    ///
    /// The specific volume of the steam state in [m3/kg].
    ///
    v: f64,
    ///
    /// The specific internal energy of the steam state in [J/kg].
    ///
    u: f64,
    ///
    /// The specific enthalpy of the steam state in [J/kg].
    ///
    h: f64,
    ///
    /// The specific entropy of the steam state in [J/Kg.k].
    ///
    s: f64,
    ///
    /// The specific isobaric heat capacity of the steam state in [J/Kg.K].
    ///
    cp: f64,
    ///
    /// The specific isochoric heat capacity of the steam state in [J/Kg.K].
    ///
    cv: f64,
    ///
    /// The ratio of specific heats.
    ///
    ratio_cp_cv: f64,
    ///
    /// Viscosity in [Pa.sec].
    ///
    mu: f64,
    ///
    /// Thermal Conductivity in [W/K.m].
    ///
    tk: f64,
}

//
// Implement the default trait for the SteamState struct.
//
impl Default for SteamState {
    //
    // Create a default SteamState.
    //
    fn default() -> Self {
        SteamState {
            region: 0,
            pressure: 0.0,
            temperature: 0.0,
            density: 0.0,
            quality: 0.0,
            v: -1.0,
            u: -1.0,
            h: -1.0,
            s: -1.0,
            cp: -1.0,
            cv: -1.0,
            ratio_cp_cv: -1.0,
            mu: -1.0,
            tk: -1.0,
        }
    }
}

//
// Implement the SteamState struct.
//
impl SteamState {
    ///
    /// Creates a new SteamState with default values.
    ///
    pub fn new() -> Self {
        Self::default()
    }

    ///
    /// Sets the primary properties of the steam state.
    /// These properties depend on the region the steam state belongs to.
    /// For regions 1, 2 and 5 these are pressure and temperature. For
    /// region 3 they are density and temperature and for region 4 they are
    /// temperature and steam quality.
    ///
    /// - Arguments:
    ///   - `region`: The region of this steam state.
    ///   - `prop1`: The first primary property (as described above).
    ///   - `prop2`: The second primary property (as described above).
    ///
    /// - Returns:
    ///   - `bool`: True if the specified thermodynamic properties are valid
    ///             for the region of this state.
    ///
    pub fn set_primary_properties(&mut self, region: usize, prop1: f64, prop2: f64) -> bool {
        self.region = region;
        let valid;

        match region {
            1 => {
                self.pressure = prop1;
                self.temperature = prop2;
                valid = steam_regions::in_region1(prop1, prop2);
            }
            2 => {
                self.pressure = prop1;
                self.temperature = prop2;
                valid = steam_regions::in_region2(prop1, prop2);
            }
            3 => {
                self.density = prop1;
                self.temperature = prop2;
                let p = steam_region3::pressure(prop1, prop2);
                valid = steam_regions::in_region3(p, prop2);
            }
            4 => {
                self.temperature = prop1;
                self.quality = prop2;
                let p = steam_region4::saturation_pressure(prop1);
                valid = steam_regions::in_region4(p, prop1, 1.0e-4);
            }
            5 => {
                self.pressure = prop1;
                self.temperature = prop2;
                valid = steam_regions::in_region5(prop1, prop2);
            }
            _ => {
                //
                // Invalid region.
                //
                valid = false;
            }
        }
        //
        // Set the other thermodynamic properties to invalid values.
        //
        self.v = -1.0;
        self.u = -1.0;
        self.h = -1.0;
        self.s = -1.0;
        self.cp = -1.0;
        self.cv = -1.0;
        self.ratio_cp_cv = -1.0;
        self.mu = -1.0;
        self.tk = -1.0;

        valid
    }

    ///
    /// Returns the specific volume of the steam state.
    ///
    /// - Returns:
    ///   - The specific volume [m3/kg].
    ///
    pub fn get_specific_volume(&mut self) -> f64 {
        if self.v > ZERO {
            return self.v;
        }

        self.v = match self.region {
            1 => steam_region1::specific_volume(self.pressure, self.temperature),
            2 => steam_region2::specific_volume(self.pressure, self.temperature),
            3 => ONE / self.density,
            4 => steam_region4::specific_volume(self.temperature, self.quality),
            5 => steam_region5::specific_volume(self.pressure, self.temperature),
            _ => ZERO,
        };

        self.v
    }

    ///
    /// Returns the pressure of the steam state.
    ///
    /// - Returns
    ///   - The pressure [MPa].
    ///
    pub fn get_pressure(&mut self) -> f64 {
        match self.region {
            1 | 2 | 5 => self.pressure,
            3 => {
                self.pressure = steam_region3::pressure(self.density, self.temperature);
                self.pressure
            }
            4 => {
                self.pressure = steam_region4::saturation_pressure(self.temperature);
                self.pressure
            }
            _ => ZERO,
        }
    }

    ///
    /// Returns the density of the steam state.
    ///
    /// - Returns:
    ///   - The density [Kg/m3].
    ///
    pub fn get_density(&mut self) -> f64 {
        ONE / self.get_specific_volume()
    }

    ///
    /// Returns the temperature of the steam state.
    ///
    /// - Returns:
    ///   - The temperature [K].
    ///
    pub fn get_temperature(&self) -> f64 {
        self.temperature
    }

    ///
    /// Returns the specific internal energy of the steam state.
    ///
    /// - Returns:
    ///   - The specific internal energy [J/kg].
    ///
    pub fn get_specific_internal_energy(&mut self) -> f64 {
        if self.u > ZERO {
            return self.u;
        }

        self.u = match self.region {
            1 => steam_region1::specific_internal_energy(self.pressure, self.temperature),
            2 => steam_region2::specific_internal_energy(self.pressure, self.temperature),
            3 => steam_region3::specific_internal_energy(self.density, self.temperature),
            4 => steam_region4::specific_internal_energy(self.temperature, self.quality),
            5 => steam_region5::specific_internal_energy(self.pressure, self.temperature),
            _ => ZERO,
        };

        self.u
    }

    ///
    /// Returns the specific enthalpy of the steam state.
    ///
    /// - Returns:
    ///   - The specific enthalpy [J/kg].
    ///
    pub fn get_specific_enthalpy(&mut self) -> f64 {
        if self.h > ZERO {
            return self.h;
        }

        self.h = match self.region {
            1 => steam_region1::specific_enthalpy(self.pressure, self.temperature),
            2 => steam_region2::specific_enthalpy(self.pressure, self.temperature),
            3 => steam_region3::specific_enthalpy(self.density, self.temperature),
            4 => steam_region4::specific_enthalpy(self.temperature, self.quality),
            5 => steam_region5::specific_enthalpy(self.pressure, self.temperature),
            _ => ZERO,
        };

        self.h
    }

    ///
    /// Returns the specific entropy of the steam state.
    ///
    /// - Returns:
    ///   - The specific entropy [J/kg.K].
    ///
    pub fn get_specific_entropy(&mut self) -> f64 {
        if self.s > ZERO {
            return self.s;
        }

        self.s = match self.region {
            1 => steam_region1::specific_entropy(self.pressure, self.temperature),
            2 => steam_region2::specific_entropy(self.pressure, self.temperature),
            3 => steam_region3::specific_entropy(self.density, self.temperature),
            4 => steam_region4::specific_entropy(self.temperature, self.quality),
            5 => steam_region5::specific_entropy(self.pressure, self.temperature),
            _ => ZERO,
        };

        self.s
    }

    ///
    /// Returns the specific isobaric heat capacity of the steam state.
    ///
    /// - Returns:
    ///   - The specific isobaric heat capacity [J/kg.K].
    ///
    pub fn get_specific_isobaric_heat_capacity(&mut self) -> f64 {
        if self.cp > ZERO {
            return self.cp;
        }

        self.cp = match self.region {
            1 => steam_region1::specific_isobaric_heat_capacity(self.pressure, self.temperature),
            2 => steam_region2::specific_isobaric_heat_capacity(self.pressure, self.temperature),
            3 => steam_region3::specific_isobaric_heat_capacity(self.density, self.temperature),
            4 => steam_region4::specific_isobaric_heat_capacity(self.temperature, self.quality),
            5 => steam_region5::specific_isobaric_heat_capacity(self.pressure, self.temperature),
            _ => ZERO,
        };

        self.cp
    }

    ///
    /// Returns the specific isochoric heat capacity of the steam state.
    ///
    /// - Returns:
    ///   - The specific isochoric heat capacity [J/kg.K].
    ///
    pub fn get_specific_isochoric_heat_capacity(&mut self) -> f64 {
        if self.cv > ZERO {
            return self.cv;
        }

        self.cv = match self.region {
            1 => steam_region1::specific_isochoric_heat_capacity(self.pressure, self.temperature),
            2 => steam_region2::specific_isochoric_heat_capacity(self.pressure, self.temperature),
            3 => steam_region3::specific_isochoric_heat_capacity(self.density, self.temperature),
            4 => steam_region4::specific_isochoric_heat_capacity(self.temperature, self.quality),
            5 => steam_region5::specific_isochoric_heat_capacity(self.pressure, self.temperature),
            _ => ZERO,
        };

        self.cv
    }

    ///
    /// Returns the ratio of specific heats of the steam state.
    ///
    /// - Returns:
    ///   - The ratio of specific heats.
    ///
    pub fn get_ratio_of_specific_heats(&mut self) -> f64 {
        self.ratio_cp_cv = self.get_specific_isobaric_heat_capacity()
            / self.get_specific_isochoric_heat_capacity();

        self.ratio_cp_cv
    }

    ///
    /// Returns the viscosity of the steam state.
    ///
    /// - Returns:
    ///   - The viscosity [Pa.sec].
    ///
    pub fn get_viscosity(&mut self) -> f64 {
        if self.mu > ZERO {
            return self.mu;
        }
        self.mu = steam_viscosity::viscosity(self.get_density(), self.temperature);
        self.mu
    }

    ///
    /// Returns the thermal conductivity of the steam state.
    ///
    /// - Returns:
    ///   - The thermal conductivity [W/K.m].
    ///
    pub fn get_thermal_conductivity(&mut self) -> f64 {
        if self.tk > ZERO {
            return self.tk;
        }
        self.tk = steam_conductivity::thermal_conductivity(self.get_density(), self.temperature);
        self.tk
    }

    ///
    /// Prints the state information to standard output.
    ///
    pub fn print_state(&mut self) {
        println!();
        println!("Steam State - Thermodynamic Properties");
        println!("======================================");
        println!();

        println!("Steam Region                     = {}", self.region);
        println!(
            "Pressure                         = {} [MPa]",
            self.get_pressure()
        );
        println!(
            "Density                          = {} [Kg/m3]",
            self.get_density()
        );
        println!(
            "Temperature                      = {} [K]",
            self.get_temperature()
        );
        println!(
            "Specific Volume                  = {} [m3/Kg]",
            self.get_specific_volume()
        );
        println!(
            "Specific Internal Energy         = {} [J/Kg]",
            self.get_specific_internal_energy()
        );
        println!(
            "Specific Enthalpy                = {} [J/Kg]",
            self.get_specific_enthalpy()
        );
        println!(
            "Specific Entropy                 = {} [J/Kg.K]",
            self.get_specific_entropy()
        );
        println!(
            "Specific Isobaric Heat Capacity  = {} [J/Kg.K]",
            self.get_specific_isobaric_heat_capacity()
        );
        println!(
            "Specific Isochoric Heat Capacity = {} [J/Kg.K]",
            self.get_specific_isochoric_heat_capacity()
        );
        println!(
            "Ratio of Specific Heats          = {}",
            self.get_ratio_of_specific_heats()
        );
        println!(
            "Viscosity                        = {} [Pa.sec]",
            self.get_viscosity()
        );
        println!(
            "Thermal Conductivity             = {} [W/K.m]",
            self.get_thermal_conductivity()
        );

        println!();
    }
}
