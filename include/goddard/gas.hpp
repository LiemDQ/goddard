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
    Gas(const Cantera::Solution& gas,
        GasChemistry chemistry = GasChemistry::FROZEN);
    
    Gas(Cantera::Solution&& gas,
        GasChemistry chemistry = GasChemistry::FROZEN);

    Gas(const Gas& gas);
    Gas operator=(const Gas& gas);

    static Gas create(const std::string& filename, 
        GasChemistry chemistry = GasChemistry::FROZEN,
        const std::string& phase_name = "");


    // State setters
    using Composition = Cantera::Composition;
    void set_state_TD(double T, double D);
    void set_state_TP(double T, double P);
    void set_state_TPX(double T, double P, const std::string& composition);
    void set_state_TPX(double T, double P, const Composition& composition);
    void set_state_TPY(double T, double P, const std::string& composition);
    void set_state_TPY(double T, double P, const Composition& composition);
    void set_state_HP(double H, double P);
    void set_state_SP(double S, double P);
    void set_state_UV(double U, double V);

    // Saving/loading state

    std::vector<double> save_state() const;
    void copy_state(std::vector<double>& state) const;
    void restore_state(const std::vector<double>& state);
    ThermodynamicState snapshot() const;

    // Basic thermodynamic properties
    double temperature() const;
    double pressure() const;
    double density() const;
    double enthalpy_mass() const;
    double entropy_mass() const;
    double cp_mass() const;
    double cv_mass() const;
    double molecular_weight() const;
    size_t num_species() const;
    std::vector<std::string> species_names() const;

    // Chemistry-aware derived properties
    double gamma_s() const;
    double speed_of_sound() const;
    double stagnation_enthalpy(double velocity) const;
    double stagnation_pressure(double velocity) const;
    double isenthalpic_velocity(double H_stagnation) const;
    double isenthalpic_velocity() const;
    double mach(double velocity) const;
    double cstar() const;
    double isp() const;
    double ivac() const;

    // Mixture properties
    double fuel_fraction(
        const std::string& fuel_comp, 
        const std::string& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar) const;
    double fuel_fraction(
        const Composition& fuel_comp, 
        const Composition& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar) const;
    
    double equivalence_ratio() const;
    double equivalence_ratio(
        const std::string& fuel_comp, 
        const std::string& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar) const;
    double equivalence_ratio(
        const Composition& fuel_comp, 
        const Composition& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar) const;
    
    double stoich_OF_ratio(
        const std::string& fuel_comp, 
        const std::string& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar) const;
    double stoich_OF_ratio(
        const Composition& fuel_comp, 
        const Composition& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar) const;
    
    void set_fuel_fraction(
        double fuel_frac,
        const std::string& fuel_comp, 
        const std::string& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar);
    void set_fuel_fraction(
        double fuel_frac,
        const Composition& fuel_comp, 
        const Composition& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar);

    void set_equivalence_ratio(
        double phi,
        const std::string& fuel_comp, 
        const std::string& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar);
    void set_equivalence_ratio(
        double phi,
        const Composition& fuel_comp, 
        const Composition& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar);

    void set_OF_ratio(
        double OF,
        const std::string& fuel_comp, 
        const std::string& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar);
    void set_OF_ratio(
        double OF,
        const Composition& fuel_comp, 
        const Composition& ox_comp,
        Cantera::ThermoBasis basis = Cantera::ThermoBasis::molar);


    // Expansion properties 
    ExpansionProperties expansion_properties() const;
    void equilibrate(const std::string& XY);
    

    // Reference state
    /**
     * Set stagnation enthalpy.The stagnation enthalpy is used to calculate 
     * stagnation values and must be set explicitly. 
     * It does not change when the thermodynamic state of `Gas` changes.
     */
    void set_stagnation_enthalpy(double H);
    /**
     * Get stagnation enthalpy. The stagnation enthalpy is used to calculate 
     * stagnation values and must be set explicitly. 
     * It does not change when the thermodynamic state of `Gas` changes.
     */
    double get_stagnation_enthalpy() const;
    void set_reference_entropy(double S);
    double get_reference_entropy() const;
    

    void set_current_state_as_reference();

    // Access to underlying Cantera objects
    std::shared_ptr<Cantera::Solution> solution() const;
    std::shared_ptr<Cantera::ThermoPhase> thermo() const;
    std::shared_ptr<Cantera::Kinetics> kinetics() const;
    std::shared_ptr<Cantera::Transport> transport() const;

    /**
     * Generate report string of underlying Cantera `ThermoPhase` object.
     */
    std::string report(bool show_thermo = true, double threshold = -1e-14) const;

    GasChemistry chemistry;

private:
    std::shared_ptr<Cantera::Solution> m_sol;
    double m_H_stagnation;
    double m_S0;
};

} // namespace Goddard
