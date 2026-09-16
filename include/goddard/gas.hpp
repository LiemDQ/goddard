#pragma once
#include <memory>
#include <utility>
#include <vector>
#include <string>
#include "cantera/core.h"
#include "eigen3/Eigen/Dense"
#include "goddard/chemistry.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/thermo.hpp"

namespace Goddard {

/** Candidate condensed species of a `Gas`. Implementation detail; see `condensed.hpp`. */
class CondensedPhaseSet;

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
     * Create a `Gas` that carries a set of candidate condensed species.
     *
     * @note Both the `Solution` and the condensed set are referenced, not copied. Copies of the
     * resulting `Gas` share both; `clone()` deep-copies both.
     */
    Gas(std::shared_ptr<Cantera::Solution> gas,
        std::shared_ptr<CondensedPhaseSet> condensed,
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

    /**
     * Create a new Gas with a deep copy of the underlying Solution object and of the condensed
     * species set, if any. The chemistry mode and stored reference stagnation enthalpy and
     * entropy are preserved.
     *
     * @note Copying a `Gas` normally is shallow: the copy references the same `Solution` and the
     * same condensed species set.
     */
    Gas clone() const;

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

    /**
     * Snapshot the current thermodynamic state into an opaque vector that can later be passed to
     * restore_state().
     *
     * The vector is the Cantera state of the gas phase followed by one entry per candidate
     * condensed species (its amount in kmol per kg of mixture). With no candidates attached it is
     * the Cantera state alone, exactly as before.
     */
    std::vector<double> save_state() const;
    void copy_state(std::vector<double>& state) const;
    /**
     * Restore a thermodynamic state previously captured by save_state().
     *
     * Accepts either the extended length described in `save_state()` or the bare Cantera state
     * length, in which case all condensed amounts are set to zero.
     *
     * @throws std::invalid_argument if the vector has neither length.
     */
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

    /**
     * Equilibrate at fixed temperature [K] and pressure [Pa].
     * @note With condensed candidates attached this throws `NotImplementedError` until the
     * multiphase solver lands (work package B).
     */
    void equilibrate_TP(double T, double P);
    /**
     * Equilibrate at fixed specific enthalpy [J/kg] and pressure [Pa].
     * @note See `equilibrate_TP()` for the condensed-phase restriction.
     */
    void equilibrate_HP(double H, double P);
    /**
     * Equilibrate at fixed specific entropy [J/(kg.K)] and pressure [Pa].
     * @note See `equilibrate_TP()` for the condensed-phase restriction.
     */
    void equilibrate_SP(double S, double P);

    /**
     * Set the composition from elemental amounts, without equilibrating.
     *
     * The state becomes a "basis" composition holding the requested elements: each element is
     * assigned to a single single-element species (the homonuclear diatomic if the phase has one,
     * otherwise the monatomic, otherwise any species made of that element alone). This is the
     * starting point a subsequent `equilibrate_*` call refines.
     *
     * @param element_moles Amount of each element [kmol per kg of mixture], in the element order
     *                      of `element_names()`.
     * @param T Temperature [K] of the resulting state.
     * @param P Pressure [Pa] of the resulting state.
     *
     * @throws FmtError if the phase contains no single-element species for one of the requested
     * elements.
     * @note See `equilibrate_TP()` for the condensed-phase restriction.
     */
    void set_element_moles(const Eigen::ArrayXd& element_moles, double T, double P);

    // Condensed species
    /**
     * @name Condensed species
     *
     * A `Gas` optionally carries a set of candidate condensed (solid or liquid) product species
     * and the current amount of each. The set is shared by copies of the `Gas` exactly as the
     * underlying `Solution` is, and deep-copied by `clone()`.
     *
     * @warning Until the multiphase solver lands (work package B) the property getters above
     * (`density()`, `enthalpy_mass()`, `entropy_mass()`, `cp_mass()`, `molecular_weight()`,
     * `mole_fractions()`, `mass_fractions()`, `gamma_s()`, ...) report **gas-phase** values even
     * when condensed phases are present. Use `gas_mass_fraction()`,
     * `mixture_molecular_weight()` and `mixture_mass_fractions()` for mixture quantities.
     * @{
     */
    /**
     * Add named condensed species from a data file as candidates.
     *
     * @param infile Cantera YAML data file containing the species.
     * @param names Species names, matching the data file exactly.
     * @throws FmtError if a name is absent from the file, or if a species contains an element
     * the gas phase does not have.
     */
    void add_condensed_species(const std::string& infile, const std::vector<std::string>& names);
    /** Add every species of `infile` whose elements are a subset of this phase's elements. */
    void add_all_condensed_species(const std::string& infile);
    /** True if any candidate condensed species is attached. */
    bool has_condensed_candidates() const;
    /** True if any candidate condensed species is present in a nonzero amount. */
    bool has_condensed_phases() const;
    /** Names of all candidate condensed species, in candidate order. */
    std::vector<std::string> condensed_species_names() const;

    /** Amount of each candidate condensed species [kmol per kg of mixture], in candidate order. */
    std::vector<double> condensed_moles() const;
    /**
     * Set the amount of each candidate condensed species [kmol per kg of mixture].
     * @throws std::invalid_argument if `moles` does not have one entry per candidate.
     */
    void set_condensed_moles(const std::vector<double>& moles);

    /** Gas mass fraction w_g = 1 - sum_k n_k M_k [-]. 1 for a gas-only mixture. */
    double gas_mass_fraction() const;
    /**
     * Mixture molecular weight [kg/kmol], CEA's "MW": 1 / (n + sum_k n_k) with n the moles of
     * gas per kg of mixture. `molecular_weight()` remains CEA's "M" = 1/n, the gas-phase value.
     */
    double mixture_molecular_weight() const;

    /** True if the mixture sits exactly at a condensed phase transition (a pinned state). */
    bool at_phase_transition() const;
    /**
     * Indices of the two coexisting polymorphs at a pinned phase transition, lower-temperature
     * polymorph first. The indices refer to the condensed species that are currently *present*
     * (nonzero moles), in candidate order. `{-1, -1}` when no group is pinned.
     */
    std::pair<long, long> pinned_polymorphs() const;
    /**
     * Mass fractions of the whole mixture [-]: the gas-phase mass fractions scaled by
     * `gas_mass_fraction()`, followed by n_k M_k for every candidate condensed species in
     * candidate order. Sums to 1.
     */
    std::vector<double> mixture_mass_fractions() const;

    /** Reference-state molar enthalpy divided by R*T [-] of each *present* condensed species. */
    Eigen::ArrayXd condensed_enthalpy_RT() const;
    /** Reference-state molar heat capacity divided by R [-] of each *present* condensed species. */
    Eigen::ArrayXd condensed_cp_R() const;
    /** Molar mass [kg/kmol] of each *present* condensed species. */
    Eigen::ArrayXd condensed_molar_masses() const;
    /**
     * Stoichiometric coefficients of the *present* condensed species: a |C| x l array whose
     * rows are species and columns are elements, in the element order of `element_names()`.
     */
    Eigen::ArrayXXd condensed_stoich_coeffs() const;
    /** @} */

    /** Element names of the gas phase, in the phase's element order. */
    std::vector<std::string> element_names() const;
    /**
     * Amount of each element [kmol per kg of mixture], gas plus condensed, in the element order
     * of `element_names()`.
     */
    Eigen::ArrayXd element_moles() const;

    /** Settings for the equilibrium solvers used by `equilibrate_TP/HP/SP`. */
    EquilibriumOptions equilibrium_options;


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
    // `ThermoArray` stores a clone of the condensed species set of the `Gas` it is built from.
    // `CondensedPhaseSet` is an implementation detail, so it is reached through friendship rather
    // than through a public accessor.
    friend class ThermoArray;

    std::shared_ptr<Cantera::Solution> m_sol;
    /** Candidate condensed species and their amounts. Null means the `Gas` is gas-only. */
    std::shared_ptr<CondensedPhaseSet> m_condensed;
    double m_H_stagnation = 0.0;
    double m_S0 = 0.0;
    double m_gamma = 0.0;

    inline void check_for_valid_cantera() const {
        if (!has_cantera_sln()) throw std::runtime_error("No Cantera Solution owned by this Gas object.");
    }
};

} // namespace Goddard
