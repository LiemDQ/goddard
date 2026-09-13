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
 * 
 * `Gas` is itself a wrapper class that carries very little data, and is cheap to copy. 
 */
class Gas {
public:
    /**
     * @note The created `Gas` will reference the same `Solution` object pointed to
     * by the `shared_ptr`.  
     */
    explicit Gas(std::shared_ptr<Cantera::Solution> gas,
                 GasChemistry chemistry = GasChemistry::FROZEN);
    
    /**
     * @note This constructor copies the `Solution` object. 
     */
    Gas(const Cantera::Solution& gas,
        GasChemistry chemistry = GasChemistry::FROZEN);
    
    /**
     * @note This constructor will reference the same `Solution` object. 
     */
    Gas(Cantera::Solution&& gas,
        GasChemistry chemistry = GasChemistry::FROZEN);

    /**
     * Create `Gas` from input data file. 
     */
    Gas(const std::string& infile, 
        const std::string& phase_name,
        GasChemistry chemistry = GasChemistry::FROZEN);
    
    /**
     * Create a perfect gas object from the adiabatic index.
     */
    Gas(double gamma);

    Gas(const Gas& gas);
    Gas operator=(const Gas& gas);

    static Gas create(const std::string& infile, 
        const std::string& phase_name,
        GasChemistry chemistry = GasChemistry::FROZEN);
    
    static Gas create_from_elements(const std::string& infile,
        const std::string& name, 
        const std::vector<std::string>& elements,
        GasChemistry chemistry = GasChemistry::FROZEN);
    
    static Gas create_from_species(const std::string& infile, 
        const std::string& name,
        const std::vector<std::string>& species,
        GasChemistry chemistry = GasChemistry::FROZEN);

    // State setters
    /** Set state from temperature [K] and density [kg/m^3]. */
    void set_state_TD(double T, double D);
    /** Set state from temperature [K] and pressure [Pa]. */
    void set_state_TP(double T, double P);
    /** Set state from temperature [K], pressure [Pa], and mole-fraction composition string. */
    void set_state_TPX(double T, double P, const std::string& composition);
    /** Set state from temperature [K], pressure [Pa], and mole-fraction map. */
    void set_state_TPX(double T, double P, const Composition& composition);
    void set_state_TPX(double T, double P, const double* composition);
    /** Set state from temperature [K], pressure [Pa], and mass-fraction composition string. */
    void set_state_TPY(double T, double P, const std::string& composition);
    /** Set state from temperature [K], pressure [Pa], and mass-fraction map. */
    void set_state_TPY(double T, double P, const Composition& composition);
    void set_state_TPY(double T, double P, const double* composition);
    /** Set state from specific enthalpy [J/kg] and pressure [Pa]. */
    void set_state_HP(double H, double P);
    /** Set state from specific entropy [J/(kg.K)] and pressure [Pa]. */
    void set_state_SP(double S, double P);
    /** Set state from specific internal energy [J/kg] and specific volume [m^3/kg]. */
    void set_state_UV(double U, double V);

    // Saving/loading state

    /** Snapshot the current thermodynamic state into an opaque vector that can later be passed to restore_state(). */
    std::vector<double> save_state() const;
    void copy_state(std::vector<double>& state) const;
    /** Restore a thermodynamic state previously captured by save_state(). */
    void restore_state(const std::vector<double>& state);
    /** Return a ThermodynamicState struct containing all current properties. */
    ThermodynamicState snapshot() const;
    std::string name() const;
    void set_name(const std::string& name);

    // Basic thermodynamic properties
    /** Current temperature [K]. */
    double temperature() const;
    /** Current pressure [Pa]. */
    double pressure() const;
    /** Current mass density [kg/m^3]. */
    double density() const;
    /** Specific enthalpy [J/kg] at the current state. */
    double enthalpy_mass() const;
    /** Specific entropy [J/(kg.K)] at the current state. */
    double entropy_mass() const;
    /** Constant-pressure specific heat [J/(kg.K)] at the current state. */
    double cp_mass() const;
    /** Constant-volume specific heat [J/(kg.K)] at the current state. */
    double cv_mass() const;
    /** Mean molecular weight of the mixture [kg/kmol]. */
    double molecular_weight() const;
    std::vector<double> mole_fractions() const;
    std::vector<double> mass_fractions() const;

    size_t num_species() const;
    std::vector<std::string> species_names() const;

    // Chemistry-aware derived properties
    /**
     * Effective ratio of specific heats used for compressible-flow calculations.
     *
     * The result depends on the active GasChemistry mode: PERFECT_GAS uses cp/cv,
     * FROZEN holds composition fixed, EQUILIBRIUM accounts for shifting equilibrium.
     */
    double gamma_s() const;
    /** Speed of sound [m/s] at the current state, consistent with gamma_s(). */
    double speed_of_sound() const;
    /**
     * Stagnation (total) enthalpy h0 = h + v^2/2.
     * @param velocity flow velocity [m/s]
     * @return stagnation enthalpy [J/kg]
     */
    double stagnation_enthalpy(double velocity) const;
    /**
     * Stagnation pressure obtained by isentropic deceleration from the current state.
     * @param velocity flow velocity [m/s]
     * @return stagnation pressure [Pa]
     */
    double stagnation_pressure(double velocity) const;
    /**
     * Velocity that produces a given stagnation enthalpy from the current static state.
     * @param H_stagnation target stagnation enthalpy [J/kg]
     * @return velocity [m/s]
     */
    double isenthalpic_velocity(double H_stagnation) const;
    /**
     * Velocity computed from the stored reference stagnation enthalpy
     * (see set_stagnation_enthalpy()).
     * @return velocity [m/s]
     */
    double isenthalpic_velocity() const;
    double area_per_mdot(double velocity) const;
    /**
     * Mach number at the current state for a given velocity.
     * @param velocity flow velocity [m/s]
     */
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
    /**
     * Compute expansion-related thermodynamic derivatives at the current state:
     * gamma_s, dlnV/dlnT|_P, dlnV/dlnP|_T, and cp.
     */
    ExpansionProperties expansion_properties() const;
    /**
     * Equilibrate the gas mixture, holding two thermodynamic properties constant.
     * @param XY two-letter property pair to hold constant, e.g. "HP", "TP", "SP".
     * @param solver Cantera equilibrium solver name (default "gibbs").
     */
    void equilibrate(const std::string& XY, const std::string& solver = "gibbs");
    

    // Reference state
    /**
     * Set the stored reference stagnation enthalpy [J/kg].
     *
     * The reference value is used by isenthalpic_velocity() and other stagnation
     * helpers and must be set explicitly. It does not change when the
     * thermodynamic state of `Gas` changes.
     */
    void set_stagnation_enthalpy(double H);
    /**
     * Get the stored reference stagnation enthalpy [J/kg].
     *
     * See set_stagnation_enthalpy() for the semantics of the reference value.
     */
    double get_stagnation_enthalpy() const;
    /** Set the stored reference specific entropy [J/(kg.K)]. */
    void set_reference_entropy(double S);
    /** Get the stored reference specific entropy [J/(kg.K)]. */
    double get_reference_entropy() const;
    

    void set_current_state_as_reference();

    // Access to underlying Cantera objects
    /** Underlying Cantera Solution handle. Pass to solver constructors that accept one. */
    std::shared_ptr<Cantera::Solution> solution() const;
    std::shared_ptr<Cantera::ThermoPhase> thermo() const;
    std::shared_ptr<Cantera::Kinetics> kinetics() const;
    std::shared_ptr<Cantera::Transport> transport() const;

    /**
     * Generate report string of underlying Cantera `ThermoPhase` object.
     */
    std::string report(bool show_thermo = true, double threshold = -1e-14) const;
    
    /**
     * True if underlying Cantera Solution object has been initialized.
     */
    inline bool has_cantera_sln() const {
        return m_sol != nullptr;
    }

    GasChemistry chemistry;

private:
    std::shared_ptr<Cantera::Solution> m_sol;
    double m_H_stagnation;
    double m_S0;
    double m_gamma;

    inline void check_for_valid_cantera() const {
        if (!has_cantera_sln()) throw std::runtime_error("No Cantera Solution owned by this Gas object.");
    }
};

} // namespace Goddard
