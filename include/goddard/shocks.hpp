#pragma once
#include "cantera/core.h"
#include <utility>
namespace Goddard {

struct ShockResult {
    bool valid;
    double mach_in;
    double mach_out;
    double static_pressure_ratio;
    double static_temperature_ratio;
    double total_pressure_ratio;
};

struct ObliqueShockResult {
    bool valid;
    double mach_in;
    double mach_out;
    ShockResult shock;
    double beta;
    double theta;
};

/**
 * Solve for post-shock conditions for a calorically perfect gas.
 */
ShockResult normal_shock(double mach, double gamma);
/**
 * Solve for post-shock conditions for a general gas, assuming
 * frozen composition (no chemical reactions). 
 */
ShockResult normal_shock(Cantera::ThermoPhase& gas, double mach);

/**
 * Solve for perfect gas properties across a reflected shock.
 */
ShockResult reflected_shock(double mach, double gamma);

/** 
 * Solve for frozen gas properties across a reflected shock. 
 */
ShockResult reflected_shock(Cantera::ThermoPhase& gas, double mach);

/**
 * Solve for oblique shock angles for a perfect gas oblique shock with known deflection angle.
 * @returns First entry is the weak shock angle, second entry is the strong shock angle.
 */
std::pair<double,double> oblique_shock_wave_angle(double mach, double deflection_angle, double gamma);

/**
 * Solve for deflection angle of a perfect gas given the angle of the oblique shock.
 */
double oblique_shock_deflection_angle(double mach, double wave_angle, double gamma);

double oblique_shock_max_deflection(double mach, double gamma);


/**
 * Get the full set of shock relations for a perfect gas and a specified shock angle relative
 * to the flow.
 */
ObliqueShockResult oblique_shock_from_wave_angle(
    double mach, double wave_angle, double gamma);
    
ObliqueShockResult oblique_shock_from_wave_angle(
    Cantera::ThermoPhase& gas, double mach, double wave_angle);

/**
 * Get the full set of shock relations for a perfect gas with a specified deflection angle.
 */
ObliqueShockResult oblique_shock_from_deflection(
    double mach, double deflection_angle, double gamma, bool weak = true);

/**
 * Get the full set of shock relations for a general gas with frozen composition with a 
 * specified deflection angle.
 */
ObliqueShockResult oblique_shock_from_deflection(
    Cantera::ThermoPhase& gas, double mach, double deflection_angle, bool weak = true);



}