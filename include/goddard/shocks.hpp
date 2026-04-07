#pragma once
#include <vector>
#include <utility>
#include "goddard/gas.hpp"

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

// --- Perfect gas free functions (pure, no state) ---

ShockResult normal_shock(double mach, double gamma);
ShockResult reflected_shock(double mach, double gamma);

std::pair<double,double> oblique_shock_wave_angle(double mach, double deflection_angle, double gamma);
double oblique_shock_deflection_angle(double mach, double wave_angle, double gamma);
double oblique_shock_max_deflection(double mach, double gamma);

ObliqueShockResult oblique_shock_from_wave_angle(
    double mach, double wave_angle, double gamma);
ObliqueShockResult oblique_shock_from_deflection(
    double mach, double deflection_angle, double gamma, bool weak = true);

// --- ShockSolver: chemistry-dispatched solver with state management ---

class ShockSolver {
public:
    ShockSolver(Gas gas, SolverOptions options = {});

    ShockResult normal_shock(double mach);
    ShockResult reflected_shock(double mach);
    ObliqueShockResult oblique_shock_from_wave_angle(double mach, double wave_angle);
    ObliqueShockResult oblique_shock_from_deflection(double mach, double deflection, bool weak = true);

    const Gas& pre_shock_state() const;
    const Gas& post_shock_state() const;

private:
    mutable Gas m_gas;
    SolverOptions m_options;
    std::vector<double> m_pre_shock_state;
    std::vector<double> m_post_shock_state;
};

}
