#include "goddard/combustor.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/error.hpp"
#include "cantera/core.h"
#include <cassert>
namespace Goddard {

Combustor::Combustor(
    std::shared_ptr<Cantera::Solution>& fuel, 
    std::shared_ptr<Cantera::Solution>& oxidizer, 
    std::shared_ptr<Cantera::Solution>& products
) : m_fuel_sln(fuel), m_oxidizer_sln(oxidizer), m_product_sln(products)
{}

ThermoArray Combustor::solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options) {
    Eigen::ArrayXXd mole_fracs = generate_mole_fraction_matrix(mr);
    

    //NOTE: this assumes the species layout in "products" is the same as in the feed object
    ThermoArray combustion_states = ThermoArray(m_product_sln, {temperatures.size(), pressures.size(), mr.molar_ratio().size() });
    combustion_states.TPX(temperatures, pressures, mole_fracs);

    //combustion is adiabatic and isobaric for subsonic flames
    switch (options.type) {
        case CombustorType::INFINITE_AREA: {
            combustion_states.equilibrate("HP", "gibbs");
            break;
        }
        default: throw NotImplementedError("Finite area combustors");
    }
    
    return combustion_states;
}

/**
 * Create a matrix with the mole fractions of each component participating in the reaction.
 */
Eigen::ArrayXXd Combustor::generate_mole_fraction_matrix(const MixtureRatios& mr) const {
    
    auto oxidizer_thermo = m_oxidizer_sln->thermo();
    auto fuel_thermo = m_fuel_sln->thermo();
    Eigen::ArrayXd mole_ratios = mr.molar_ratio();
    Eigen::ArrayXd moles_ox = mole_ratios / (1 + mole_ratios);
    Eigen::ArrayXd moles_f = 1 - moles_ox;

    //create the matrix.
    //this needs to include all species from both oxidizer and fuel
    //this may be problematic if there are unused species in the reactants
    Eigen::ArrayXXd mole_frac_matrix = Eigen::ArrayXXd::Zero(mole_ratios.size(), m_product_sln->thermo()->nSpecies()); 
    Cantera::Composition ox_species = oxidizer_thermo->getMoleFractionsByName();
    Cantera::Composition fuel_species = fuel_thermo->getMoleFractionsByName();
    
    for (int i = 0; i < mole_frac_matrix.rows(); i++) {
        assign_mole_frac_row_entries(mole_frac_matrix, i, ox_species, moles_ox[i]);
        assign_mole_frac_row_entries(mole_frac_matrix, i, fuel_species, moles_f[i]);
    }

    assert((mole_frac_matrix >= 0).all() && "Mole fractions must be nonnegative");

    return mole_frac_matrix;
}

void Combustor::assign_mole_frac_row_entries(
    Eigen::ArrayXXd& matrix, 
    long row_idx, 
    const Cantera::Composition& composition, 
    double coeff) const {

    auto row = matrix.row(row_idx);
    auto product_thermo = m_product_sln->thermo();

    for (auto&& entry:composition){
        long entry_idx = product_thermo->speciesIndex(entry.first);
        row(entry_idx) = entry.second * coeff;
    }
}

} //nameplace Goddard