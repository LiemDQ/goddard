#pragma once

namespace Goddard {

enum class GasChemistry {
    PERFECT_GAS,    // constant gamma, no Cantera dependency
    FROZEN,         // composition fixed, variable thermodynamic properties via Cantera
    EQUILIBRIUM,    // full chemical equilibrium at each point via Cantera
    KINETIC         // finite-rate chemistry via Cantera kinetics
};

struct SolverOptions {
    double abstol = 1e-6;
    double reltol = 1e-5;
    int max_iterations = 100;
};

} // namespace Goddard
