#include "goddard/combustor.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/error.hpp"
#include "cantera/core.h"
#include <utility>
#include <iostream>
#include <cassert>
namespace Goddard {

Combustor::Combustor(
    std::shared_ptr<Cantera::Solution> thermo,
    std::vector<double>& fuel,
    std::vector<double>& oxidizer
) : fuel_state(fuel), oxidizer_state(oxidizer), m_thermo(std::move(thermo))
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

ThermoArray Combustor::combust(ThermoArray& states, const CombustorOptions& options) {
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

void Combustor::assign_mole_frac_row_entries(
    Eigen::ArrayXXd& matrix, 
    long row_idx, 
    const Cantera::Composition& composition, 
    double coeff) const {

    auto row = matrix.row(row_idx);
    auto product_thermo = m_thermo->thermo();

    for (auto&& entry:composition){
        long entry_idx = product_thermo->speciesIndex(entry.first);
        row(entry_idx) = entry.second * coeff;
    }
}

} //nameplace Goddard