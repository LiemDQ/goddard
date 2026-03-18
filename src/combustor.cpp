#include "goddard/combustor.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/error.hpp"
#include "cantera/core.h"
#include <utility>
#include <iostream>
#include <cassert>
namespace Goddard {

// ---- BaseCombustor ----

BaseCombustor::BaseCombustor(
    std::shared_ptr<Cantera::Solution> thermo
) : m_thermo(std::move(thermo))
{}

ThermoArray BaseCombustor::combust(ThermoArray& states, const CombustorOptions& options) {
    //combustion is adiabatic and isobaric for subsonic flames
    switch (options.type) {
        case CombustorType::INFINITE_AREA: {
            states.equilibrate("HP", "gibbs");
            break;
        }
        default: throw NotImplementedError("Finite area combustors are not implemented.");
    }
    return states;
}

void BaseCombustor::assign_mole_frac_row_entries(
    Eigen::ArrayXXd& matrix,
    long row_idx,
    const Cantera::Composition& composition,
    double coeff) const {

    auto row = matrix.row(row_idx);
    auto product_thermo = m_thermo->thermo();

    for (auto&& entry:composition){
        long entry_idx = product_thermo->speciesIndex(entry.first);
        row(entry_idx) += entry.second * coeff;
    }
}

// ---- Combustor ----

Combustor::Combustor(
    std::shared_ptr<Cantera::Solution> thermo,
    std::vector<double>& fuel,
    std::vector<double>& oxidizer
) : BaseCombustor(std::move(thermo)), fuel_state(fuel), oxidizer_state(oxidizer)
{}

ThermoArray Combustor::solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options) {
    Eigen::ArrayXXd mole_fracs = generate_mole_fraction_matrix(mr);

    //NOTE: this assumes the species layout in "products" is the same as in the feed object
    ThermoArray combustion_states = ThermoArray(m_thermo, {temperatures.size(), pressures.size(), mr.molar_ratio().size() });
    combustion_states.TPX(temperatures, pressures, mole_fracs);

    return combust(combustion_states, options);
}

ThermoArray Combustor::solve(const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options) {
    auto thermo = m_thermo->thermo();

    Eigen::ArrayXd mass_ox = mr.oxidizer_mass_frac();
    Eigen::ArrayXd mass_f = mr.fuel_mass_frac();
    Eigen::ArrayXXd mass_fracs = generate_mass_fraction_matrix(mr);

    thermo->restoreState(oxidizer_state);
    double oxidizer_enthalpy = thermo->enthalpy_mass();

    thermo->restoreState(fuel_state);
    double fuel_enthalpy = thermo->enthalpy_mass();

    Eigen::ArrayXd enthalpies = oxidizer_enthalpy * mass_ox + fuel_enthalpy * mass_f;

    ThermoArray combustion_states = ThermoArray(m_thermo, {enthalpies.size(), pressures.size(), mr.molar_ratio().size()});
    combustion_states.HPY(enthalpies, pressures, mass_fracs);

    return combust(combustion_states, options);
}

/**
 * Create a matrix with the mole fractions of each component participating in the reaction.
 */
Eigen::ArrayXXd Combustor::generate_mole_fraction_matrix(const MixtureRatios& mr) const {
    auto thermo = m_thermo->thermo();

    Eigen::ArrayXd mole_ratios = mr.molar_ratio();
    Eigen::ArrayXd moles_ox = mole_ratios / (1 + mole_ratios);
    Eigen::ArrayXd moles_f = 1 - moles_ox;

    //create the matrix.
    //this needs to include all species from both oxidizer and fuel
    //this may be problematic if there are unused species in the reactants
    Eigen::ArrayXXd mole_frac_matrix = Eigen::ArrayXXd::Zero(mole_ratios.size(), thermo->nSpecies());
    thermo->restoreState(oxidizer_state);
    Cantera::Composition ox_species = thermo->getMoleFractionsByName();

    thermo->restoreState(fuel_state);
    Cantera::Composition fuel_species = thermo->getMoleFractionsByName();

    for (int i = 0; i < mole_frac_matrix.rows(); i++) {
        assign_mole_frac_row_entries(mole_frac_matrix, i, ox_species, moles_ox[i]);
        assign_mole_frac_row_entries(mole_frac_matrix, i, fuel_species, moles_f[i]);
    }

    assert((mole_frac_matrix >= 0).all() && "Mole fractions must be nonnegative");

    return mole_frac_matrix;
}


/**
 * Create a matrix with the mass fractions of each component participating in the reaction.
 */
Eigen::ArrayXXd Combustor::generate_mass_fraction_matrix(const MixtureRatios& mr) const {
    auto thermo = m_thermo->thermo();

    Eigen::ArrayXd moles_ox = mr.oxidizer_mass_frac();
    Eigen::ArrayXd moles_f = mr.fuel_mass_frac();

    //create the matrix.
    //this needs to include all species from both oxidizer and fuel
    //this may be problematic if there are unused species in the reactants
    Eigen::ArrayXXd mass_frac_matrix = Eigen::ArrayXXd::Zero(moles_ox.size(), thermo->nSpecies());
    thermo->restoreState(oxidizer_state);
    Cantera::Composition ox_species = thermo->getMassFractionsByName();

    thermo->restoreState(fuel_state);
    Cantera::Composition fuel_species = thermo->getMassFractionsByName();

    for (int i = 0; i < mass_frac_matrix.rows(); i++) {
        assign_mole_frac_row_entries(mass_frac_matrix, i, ox_species, moles_ox[i]);
        assign_mole_frac_row_entries(mass_frac_matrix, i, fuel_species, moles_f[i]);
    }

    assert((mass_frac_matrix >= 0).all() && "Mole fractions must be nonnegative");

    return mass_frac_matrix;
}

// ---- DilutedCombustor ----

DilutedCombustor::DilutedCombustor(
    std::shared_ptr<Cantera::Solution> thermo,
    std::vector<double>& fuel,
    std::vector<double>& oxidizer,
    std::vector<double>& flue
) : BaseCombustor(std::move(thermo)),
    fuel_state(fuel), oxidizer_state(oxidizer), flue_state(flue)
{}

ThermoArray DilutedCombustor::solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, double recirculation_ratio, const CombustorOptions& options) {
    Eigen::ArrayXXd mole_fracs = generate_mole_fraction_matrix(mr, recirculation_ratio);

    ThermoArray combustion_states = ThermoArray(m_thermo, {temperatures.size(), pressures.size(), mr.molar_ratio().size()});
    combustion_states.TPX(temperatures, pressures, mole_fracs);

    return combust(combustion_states, options);
}

ThermoArray DilutedCombustor::solve(const Eigen::ArrayXd& pressures, const MixtureRatios& mr, double dilution_ratio, const CombustorOptions& options) {
    auto thermo = m_thermo->thermo();

    double r = dilution_ratio;
    double wReactant = 1.0 / (1.0 + r);
    double wDilution = r / (1.0 + r);

    // Mass fractions of fuel and oxidizer within the fresh feed, scaled by wFresh
    Eigen::ArrayXd massFuelTotal = mr.fuel_mass_frac() * wReactant;
    Eigen::ArrayXd massOxTotal = mr.oxidizer_mass_frac() * wReactant;

    Eigen::ArrayXXd mass_fracs = generate_mass_fraction_matrix(mr, dilution_ratio);

    thermo->restoreState(fuel_state);
    double fuelEnthalpy = thermo->enthalpy_mass();

    thermo->restoreState(oxidizer_state);
    double oxidizerEnthalpy = thermo->enthalpy_mass();

    thermo->restoreState(flue_state);
    double flueEnthalpy = thermo->enthalpy_mass();

    Eigen::ArrayXd enthalpies = fuelEnthalpy * massFuelTotal
                              + oxidizerEnthalpy * massOxTotal
                              + flueEnthalpy * wDilution;

    ThermoArray combustion_states = ThermoArray(m_thermo, {enthalpies.size(), pressures.size(), mr.molar_ratio().size()});
    combustion_states.HPY(enthalpies, pressures, mass_fracs);

    return combust(combustion_states, options);
}

/**
 * Create a mole fraction matrix blending fuel, oxidizer, and flue gas streams.
 * The fuel/oxidizer split is determined by MixtureRatios (molar basis).
 * The fresh vs flue split is determined by dilution_ratio (mass basis),
 * converted to a molar basis using the mean molecular weights of each stream.
 */
Eigen::ArrayXXd DilutedCombustor::generate_mole_fraction_matrix(const MixtureRatios& mr, double recirculation_ratio) const {
    auto thermo = m_thermo->thermo();

    double r = recirculation_ratio;

    // Compute molar weights of each stream to convert mass-based recirculation ratio to moles
    thermo->restoreState(fuel_state);
    double mwFuel = thermo->meanMolecularWeight();
    Cantera::Composition fuel_species = thermo->getMoleFractionsByName();

    thermo->restoreState(oxidizer_state);
    double mwOx = thermo->meanMolecularWeight();
    Cantera::Composition ox_species = thermo->getMoleFractionsByName();

    thermo->restoreState(flue_state);
    double mwFlue = thermo->meanMolecularWeight();
    Cantera::Composition flue_species = thermo->getMoleFractionsByName();

    // Molar ratio of oxidizer to fuel (from MixtureRatios)
    Eigen::ArrayXd moleRatios = mr.molar_ratio();

    // Mass of each stream per unit total mass (fuel + ox + flue):
    //   m_fuel = 1/(OF+1) per unit fresh feed, scaled by 1/(1+r) for total
    //   m_ox   = OF/(OF+1) per unit fresh feed, scaled by 1/(1+r) for total
    //   m_flue = r/(1+r)
    // Convert to moles: n_i = m_i / MW_i
    // Then normalize to get mole fractions of each stream.
    Eigen::ArrayXd OF = mr.OF_ratio();
    Eigen::ArrayXd nFuel = (1.0 / (OF + 1.0)) / ((1.0 + r) * mwFuel);
    Eigen::ArrayXd nOx = (OF / (OF + 1.0)) / ((1.0 + r) * mwOx);
    double nFlueBase = r / ((1.0 + r) * mwFlue);
    Eigen::ArrayXd nFlue = Eigen::ArrayXd::Constant(nFuel.size(), nFlueBase);

    Eigen::ArrayXd nTotal = nFuel + nOx + nFlue;
    Eigen::ArrayXd moleFracFuel = nFuel / nTotal;
    Eigen::ArrayXd moleFracOx = nOx / nTotal;
    Eigen::ArrayXd moleFracFlue = nFlue / nTotal;

    Eigen::ArrayXXd mole_frac_matrix = Eigen::ArrayXXd::Zero(moleRatios.size(), thermo->nSpecies());

    for (int i = 0; i < mole_frac_matrix.rows(); i++) {
        assign_mole_frac_row_entries(mole_frac_matrix, i, fuel_species, moleFracFuel[i]);
        assign_mole_frac_row_entries(mole_frac_matrix, i, ox_species, moleFracOx[i]);
        assign_mole_frac_row_entries(mole_frac_matrix, i, flue_species, moleFracFlue[i]);
    }

    assert((mole_frac_matrix >= 0).all() && "Mole fractions must be nonnegative");

    return mole_frac_matrix;
}

/**
 * Create a mass fraction matrix blending fuel, oxidizer, and flue gas streams.
 * Fresh feed mass fractions are scaled by 1/(1+r), flue gas by r/(1+r).
 */
Eigen::ArrayXXd DilutedCombustor::generate_mass_fraction_matrix(const MixtureRatios& mr, double recirculation_ratio) const {
    auto thermo = m_thermo->thermo();

    double r = recirculation_ratio;
    double wFresh = 1.0 / (1.0 + r);
    double wFlue = r / (1.0 + r);

    Eigen::ArrayXd massFuel = mr.fuel_mass_frac() * wFresh;
    Eigen::ArrayXd massOx = mr.oxidizer_mass_frac() * wFresh;

    Eigen::ArrayXXd mass_frac_matrix = Eigen::ArrayXXd::Zero(massFuel.size(), thermo->nSpecies());

    thermo->restoreState(oxidizer_state);
    Cantera::Composition ox_species = thermo->getMassFractionsByName();

    thermo->restoreState(fuel_state);
    Cantera::Composition fuel_species = thermo->getMassFractionsByName();

    thermo->restoreState(flue_state);
    Cantera::Composition flue_species = thermo->getMassFractionsByName();

    for (int i = 0; i < mass_frac_matrix.rows(); i++) {
        assign_mole_frac_row_entries(mass_frac_matrix, i, fuel_species, massFuel[i]);
        assign_mole_frac_row_entries(mass_frac_matrix, i, ox_species, massOx[i]);
        assign_mole_frac_row_entries(mass_frac_matrix, i, flue_species, wFlue);
    }

    if (!((mass_frac_matrix >= 0).all())) {
        throw ConvergenceError("Mass fractions must be nonnegative.");
    }

    return mass_frac_matrix;
}

} //namespace Goddard
