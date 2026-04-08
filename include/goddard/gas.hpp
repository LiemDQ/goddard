#pragma once
#include <memory>
#include <vector>
#include <string>
#include "cantera/core.h"
#include "goddard/chemistry.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/thermo.hpp"

namespace Goddard {

/**
 * Main Goddard class for querying thermodynamic information. `Gas` wraps Cantera's `Solution`
 * object and adds several convenience methods for processes involving reactive flow.
 */
class Gas {
public:
    explicit Gas(std::shared_ptr<Cantera::Solution> gas,
                 GasChemistry chemistry = GasChemistry::FROZEN);
    Gas(Cantera::Solution& gas,
        GasChemistry chemistry = GasChemistry::FROZEN);
    
    Gas(Cantera::Solution&& gas,
        GasChemistry chemistry = GasChemistry::FROZEN);

    Gas(const Gas& gas);
    Gas operator=(const Gas& gas);

    // State setters
    void set_state_TP(double T, double P);
    void set_state_TPX(double T, double P, const std::string& composition);
    void set_state_HP(double H, double P);
    void set_state_SP(double S, double P);
    void save_state(std::vector<double>& state) const;
    void restore_state(const std::vector<double>& state);

    // Basic thermodynamic properties
    double temperature() const;
    double pressure() const;
    double density() const;
    double enthalpy_mass() const;
    double entropy_mass() const;
    double cp_mass() const;
    double cv_mass() const;
    double molecular_weight() const;

    // Chemistry-aware derived properties
    double gamma_s() const;
    double speed_of_sound() const;
    double stagnation_enthalpy(double velocity) const;
    double stagnation_pressure(double velocity) const;
    double isenthalpic_velocity(double H_stagnation) const;
    double isenthalpic_velocity() const;
    double mach(double velocity) const;

    // Expansion properties 
    ExpansionProperties expansion_properties() const;
    void equilibrate(const std::string& XY);

    // Snapshot
    ThermoStateInfo snapshot() const;

    // Reference state
    void set_stagnation_enthalpy(double H);
    double get_stagnation_enthalpy() const;
    void set_reference_entropy(double S);
    double get_reference_entropy() const;

    // Access to underlying Cantera objects
    std::shared_ptr<Cantera::Solution> solution() const;
    std::shared_ptr<Cantera::ThermoPhase> thermo() const;
    std::shared_ptr<Cantera::Kinetics> kinetics() const;
    std::shared_ptr<Cantera::Transport> transport() const;

    GasChemistry chemistry;

private:
    std::shared_ptr<Cantera::Solution> m_sol;
    double m_H_stagnation;
    double m_S0;
};

} // namespace Goddard
