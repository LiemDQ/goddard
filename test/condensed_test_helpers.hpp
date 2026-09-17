#pragma once
/**
 * @file condensed_test_helpers.hpp
 * @brief Multiphase equilibrium helpers for condensed-species tests.
 *
 * These solve a multiphase TP equilibrium with Cantera's `MultiPhaseEquil` directly and load the
 * answer into a `Gas`, so tests can build converged condensed states without depending on
 * Goddard's own multiphase solver (work package B).
 */

#include "goddard/condensed.hpp"
#include "goddard/config.h"
#include "goddard/gas.hpp"
#include "goddard/utils.hpp"

#include "cantera/base/AnyMap.h"
#include "cantera/core.h"
#include "cantera/equil/MultiPhase.h"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Goddard {
namespace test_helpers {

/** CEA-derived condensed species data used by the condensed-phase tests. */
inline const char* condensed_data_file() { return DATA_DIR "/nasa9_condensed.yaml"; }

/** CEA-derived gas species data used by the condensed-phase tests. */
inline const char* gas_data_file() { return DATA_DIR "/nasa9_gas.yaml"; }

/**
 * Parsed data file root node, cached per file name.
 *
 * `nasa9_condensed.yaml` holds 750 species; re-parsing it for every finite-difference step would
 * dominate the runtime of these tests.
 */
inline const Cantera::AnyMap& cached_root_node(const std::string& infile) {
    static std::map<std::string, Cantera::AnyMap> roots;
    auto found = roots.find(infile);
    if (found == roots.end()) {
        found = roots.emplace(infile, load_root_node(infile)).first;
    }
    return found->second;
}

/** A phase set mirroring the candidate condensed species of `gas`, in candidate order. */
inline CondensedPhaseSet candidate_set(const Gas& gas,
                                       const std::string& infile = condensed_data_file())
{
    return CondensedPhaseSet(cached_root_node(infile), gas.condensed_species_names(),
                             *gas.thermo());
}

/**
 * Solve a multiphase TP equilibrium for `gas` and write the result back into it.
 *
 * The elemental composition is taken from the current state of `gas` (gas phase plus condensed
 * amounts), normalised to 1 kg of mixture. Only candidates whose data range contains `T` are
 * offered: `MultiPhaseEquil` throws when an out-of-range phase carries moles, and it oscillates
 * when both polymorphs of a group are offered at their shared transition temperature. On return
 * the gas phase holds the equilibrium composition and the condensed amounts are in kmol per kg
 * of mixture.
 *
 * @param gas Mixture to equilibrate; must carry condensed candidates.
 * @param T Temperature [K].
 * @param P Pressure [Pa].
 * @param restrict_to Names of the only candidates that may be offered. Empty means every
 *                    in-range candidate is offered.
 */
inline void solve_multiphase_TP(Gas& gas, double T, double P,
                                const std::vector<std::string>& restrict_to = {})
{
    CondensedPhaseSet candidates = candidate_set(gas);
    const std::vector<double> current_moles = gas.condensed_moles();
    const double gas_moles = gas.gas_mass_fraction() / gas.thermo()->meanMolecularWeight();

    std::shared_ptr<Cantera::ThermoPhase> gas_phase = gas.thermo();
    gas_phase->setState_TP(T, P);

    std::vector<size_t> offered = candidates.offered_at(T);
    if (!restrict_to.empty()) {
        std::vector<size_t> filtered;
        for (size_t k : offered) {
            if (std::find(restrict_to.begin(), restrict_to.end(), candidates.species[k].name)
                != restrict_to.end()) {
                filtered.push_back(k);
            }
        }
        offered = filtered;
    }

    std::vector<double> initial_mole_fractions(gas_phase->nSpecies());
    gas_phase->getMoleFractions(initial_mole_fractions.data());

    // `MultiPhaseEquil` starts from the mole fractions of the phases (estimate_equil = 0), which
    // fails when the initial state is far from equilibrium (an unburnt reactant mixture, say).
    // Then let it build its own estimate (-1), which is the same ladder the solver of work
    // package B will use.
    Cantera::MultiPhase mix;
    double mixture_mass = 0.0;
    for (int estimate_equil : {0, -1}) {
        gas_phase->setState_TPX(T, P, initial_mole_fractions.data());

        mix = Cantera::MultiPhase();
        mix.addPhase(gas_phase, gas_moles);
        for (size_t k : offered) {
            candidates.species[k].phase->setState_TP(T, P);
            mix.addPhase(candidates.species[k].phase, current_moles[k]);
        }
        mix.init();
        mix.setState_TP(T, P);

        try {
            mix.equilibrate("TP", "gibbs", 1e-9, 20000, 100, estimate_equil, 0);
        } catch (const Cantera::CanteraError&) {
            if (estimate_equil == -1) throw;
            continue;
        }

        // Mass of the whole mixture, so the amounts can be normalised to 1 kg.
        mixture_mass = 0.0;
        for (size_t i = 0; i < mix.nPhases(); i++) {
            mixture_mass += mix.phaseMoles(i) * mix.phase(i).meanMolecularWeight();
        }
        break;
    }

    std::vector<double> mole_fractions(gas_phase->nSpecies());
    gas_phase->getMoleFractions(mole_fractions.data());
    gas.set_state_TPX(T, P, mole_fractions.data());

    std::vector<double> new_moles(current_moles.size(), 0.0);
    for (size_t i = 0; i < offered.size(); i++) {
        new_moles[offered[i]] = mix.phaseMoles(i + 1) / mixture_mass;
    }
    gas.set_condensed_moles(new_moles);
}

/** Specific enthalpy of the whole mixture [J/kg]: w_g h_gas + sum_C n_k H_k. */
inline double mixture_enthalpy(const Gas& gas) {
    const double T = gas.temperature();
    double enthalpy = gas.gas_mass_fraction() * gas.thermo()->enthalpy_mass();

    const CondensedPhaseSet candidates = candidate_set(gas);
    const std::vector<double> moles = gas.condensed_moles();
    for (size_t k = 0; k < moles.size(); k++) {
        if (moles[k] <= 0.0) continue;
        enthalpy += moles[k] * candidates.enthalpy_RT(k, T) * Cantera::GasConstant * T;
    }
    return enthalpy;
}

/** Specific entropy of the whole mixture [J/(kg.K)]: w_g s_gas + sum_C n_k S_k. */
inline double mixture_entropy(const Gas& gas) {
    const double T = gas.temperature();
    double entropy = gas.gas_mass_fraction() * gas.thermo()->entropy_mass();

    const CondensedPhaseSet candidates = candidate_set(gas);
    const std::vector<double> moles = gas.condensed_moles();
    for (size_t k = 0; k < moles.size(); k++) {
        if (moles[k] <= 0.0) continue;
        entropy += moles[k] * candidates.entropy_R(k, T) * Cantera::GasConstant;
    }
    return entropy;
}

/** Specific volume of the whole mixture [m^3/kg], neglecting condensed volume: n R T / P. */
inline double mixture_specific_volume(const Gas& gas) {
    const double gas_moles = gas.gas_mass_fraction() / gas.thermo()->meanMolecularWeight();
    return gas_moles * Cantera::GasConstant * gas.temperature() / gas.pressure();
}

} // namespace test_helpers
} // namespace Goddard
