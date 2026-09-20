#pragma once

#include "cantera/base/AnyMap.h"
#include "cantera/core.h"
#include "eigen3/Eigen/Dense"

#include <memory>
#include <string>
#include <vector>

namespace Goddard {

/**
 * @brief Bookkeeping for one condensed (solid or liquid) product species.
 *
 * Each condensed species is carried as its own single-species Cantera phase, following the
 * Gordon & McBride convention that a condensed product is a pure phase of fixed composition.
 * The condensed volume is neglected everywhere in Goddard, so only the reference-state
 * enthalpy, entropy and heat capacity of the phase are ever used.
 */
struct CondensedSpecies {
    /** Species name exactly as written in the source data file, e.g. "AL2O3(a)". */
    std::string name;
    /** Single-species `fixed-stoichiometry` phase holding this species' thermo data. */
    std::shared_ptr<Cantera::ThermoPhase> phase;
    /** Molar mass [kg/kmol]. */
    double molar_mass = 0.0;
    /** Lower bound of the species' thermodynamic data range [K]. */
    double T_min = 0.0;
    /** Upper bound of the species' thermodynamic data range [K]. */
    double T_max = 0.0;
    /** Index of the polymorph group this species belongs to (see PolymorphGroup). */
    int group = -1;
    /**
     * Number of atoms of each element in one mole of this species [-], indexed in the
     * element order of the gas phase the set was built against.
     */
    Eigen::ArrayXd element_atoms;
};

/**
 * @brief A set of condensed species with identical elemental composition (polymorphs).
 *
 * Members are ordered by increasing `T_min`. Where two consecutive members have touching
 * data ranges (`T_max[i] == T_min[i+1]`) the shared temperature is a phase-transition
 * temperature at which both polymorphs may coexist.
 */
struct PolymorphGroup {
    /** Indices into `CondensedPhaseSet::species`, sorted by increasing `T_min`. */
    std::vector<size_t> members;
    /** Transition temperatures [K] between consecutive members with contiguous data ranges. */
    std::vector<double> transition_temperatures;
};

/**
 * @brief Candidate condensed species of a mixture, together with their current amounts.
 *
 * `CondensedPhaseSet` is an implementation detail held by `Gas` through a `std::shared_ptr`:
 * a null pointer means the `Gas` is gas-only. There is no equilibrium calculation; it only
 * owns the candidate phases, their temperature ranges, polymorph grouping and the current
 * number of moles of each candidate.
 */
class CondensedPhaseSet {
public:
    /**
     * Build a set from named species in an already-loaded data file root node.
     *
     * @param root_node Root node of a Cantera YAML data file (see `load_root_node`); the
     *                  species are looked up in its `species` list.
     * @param names     Species names to add. Must match the data file exactly.
     * @param gas       Gas phase whose element order the stoichiometry is expressed in. Every
     *                  element of every requested species must exist in this phase.
     */
    CondensedPhaseSet(const Cantera::AnyMap& root_node,
                      const std::vector<std::string>& names,
                      const Cantera::ThermoPhase& gas);

    /**
     * Names of all species in `root_node` whose elements are a subset of the elements of `gas`.
     */
    static std::vector<std::string> compatible_species(const Cantera::AnyMap& root_node,
                                                       const Cantera::ThermoPhase& gas);

    /**
     * Append further candidate species, re-deriving the polymorph groups. Moles of species
     * already present are preserved; new species start with zero moles.
     */
    void add_species(const Cantera::AnyMap& root_node,
                     const std::vector<std::string>& names,
                     const Cantera::ThermoPhase& gas);

    /** Deep copy, including freshly built phase objects, moles and pinned state. */
    std::shared_ptr<CondensedPhaseSet> clone() const;

    /** Index of a candidate by name, or `Cantera::npos` if it is not a candidate. */
    size_t species_index(const std::string& name) const;

    /** Number of candidate species. */
    inline size_t size() const { return species.size(); }

    /** Names of all candidate species, in candidate order. */
    std::vector<std::string> names() const;

    /** True if the thermodynamic data range of candidate `k` contains `T` (closed interval). */
    bool in_range(size_t k, double T) const;

    /**
     * Candidates that may be offered to the equilibrium solver at temperature `T`.
     *
     * A candidate is offered when `in_range(k, T)` holds. Exactly at a phase-transition
     * temperature both polymorphs of a group are in range; only the lower one is offered,
     * because offering both makes Cantera's `MultiPhaseEquil` oscillate (the two have equal
     * chemical potential there).
     */
    std::vector<size_t> offered_at(double T) const;

    /**
     * The member of group `group_index` that `offered_at()` would offer at temperature `T`: the
     * first member whose data range contains `T`, or `Cantera::npos` if none does.
     */
    size_t offered_member(size_t group_index, double T) const;

    /**
     * Re-derive `pinned_group` from the current amounts.
     *
     * A group is pinned exactly when two of its polymorphs are simultaneously present, which can
     * only happen at a phase-transition temperature. This is how the pinned flag is recovered by
     * `Gas::restore_state()`, so it is never stored in the state vector.
     */
    void update_pinned_group();

    /** Reference-state molar enthalpy of candidate `k` divided by R*T [-]. */
    double enthalpy_RT(size_t k, double T) const;
    /** Reference-state molar entropy of candidate `k` divided by R [-]. */
    double entropy_R(size_t k, double T) const;
    /** Reference-state molar heat capacity of candidate `k` divided by R [-]. */
    double cp_R(size_t k, double T) const;

    /** Candidate species, in the order they were added. */
    std::vector<CondensedSpecies> species;
    /** Polymorph groups. Every candidate belongs to exactly one group, possibly alone. */
    std::vector<PolymorphGroup> groups;
    /** Current amount of each candidate [kmol per kg of mixture]. One entry per candidate. */
    std::vector<double> moles;
    /**
     * Index into `groups` of the group currently pinned at a phase transition, or -1 when the
     * mixture is not at a pinned transition. Set by the multiphase solver.
     */
    int pinned_group = -1;

private:
    /** Rebuild `groups` and the `group` index of every species from the current species list. */
    void build_groups();

    /** Species nodes as read from the data file, aligned with `species`; used by `clone()`. */
    std::vector<Cantera::AnyMap> m_species_nodes;
    /** Element names of the gas phase the set was built against, in gas element order. */
    std::vector<std::string> m_element_names;
};

/**
 * @brief Build a single-species `fixed-stoichiometry` (Cantera `StoichSubstance`) phase.
 *
 * Condensed species in the NASA databases have no equation of state. When `species_node` has
 * no `equation-of-state` entry, a `constant-volume` one with the given density is injected.
 * The density only enters through the neglected `(P - P_ref) V / (R T)` term (~1e-8 at rocket
 * pressures with the default), so its exact value is immaterial.
 *
 * @param species_node Species node from a Cantera YAML data file. Taken by value; the copy is
 *                     modified before use.
 * @param density      Density [kg/m^3] used when the species node has no equation of state.
 */
std::shared_ptr<Cantera::ThermoPhase> build_condensed_phase(Cantera::AnyMap species_node,
                                                            double density = 1.0e6);

/** Phase node for a single-species `fixed-stoichiometry` phase containing `species_name`. */
Cantera::AnyMap make_condensed_phase_node(const std::string& species_name);

} // namespace Goddard
