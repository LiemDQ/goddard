#include "goddard/combustor.hpp"
#include "goddard/thermoarray.hpp"

namespace Goddard {

Combustor::Combustor(std::shared_ptr<Cantera::Solution> fuel, std::shared_ptr<Cantera::Solution> oxidizer, Eigen::ArrayXd&& temperatures, Eigen::ArrayXd&& pressures, MixtureRatios&& mr)
: fuel(fuel), oxidizer(oxidizer), temperatures(temperatures), pressures(pressures), mr(mr)
{
    feed.addPhase(fuel->thermo().get(), 1.0);
    feed.addPhase(oxidizer->thermo().get(), 1.0);

}

ThermoArray Combustor::solve(CombustionOptions options) {
    auto mole_fracs = generate_mole_fraction_matrix();
    //todo: make "reactant" solution instead of using fuel solution 
    ThermoArray combustion_states = ThermoArray(products, {temperatures.size(), pressures.size(), mr.molar_ratio().size() });
    combustion_states.TPX(temperatures, pressures, mole_fracs);
    combustion_states.equilibrate("HP", "gibbs");

    return combustion_states;
}

/**
 * Create a matrix with the mole fractions of each component participating in the reaction.
 */
Eigen::ArrayXXd Combustor::generate_mole_fraction_matrix() {
    
    auto oxidizer_thermo = oxidizer->thermo();
    auto fuel_thermo = fuel->thermo();
    Eigen::ArrayXd mole_ratios = mr.molar_ratio() / (oxidizer_thermo->meanMolecularWeight() / fuel_thermo->meanMolecularWeight());
    Eigen::ArrayXd moles_ox = mole_ratios / (1 + mole_ratios);
    Eigen::ArrayXd moles_f = 1 - moles_ox;

    //create the matrix.
    //TODO: should this be an ArrayXXd instead of a vector<ArrayXd>?
    //this needs to include all species from both oxidizer and fuel
    // auto mole_frac_matrix = std::vector<Eigen::ArrayXd>(mole_ratios.size(), Eigen::ArrayXd::Zero(feed.nSpecies())); 
    auto mole_frac_matrix = Eigen::ArrayXXd(mole_ratios.size(), feed.nSpecies());
    for (size_t i = 0; i < mole_frac_matrix.cols(); i++) {
        feed.setPhaseMoles(0, moles_f[i]);
        feed.setPhaseMoles(1, moles_ox[i]);
        feed.getMoleFractions(mole_frac_matrix.row(i).data());
    }

    return mole_frac_matrix;
}

} //nameplace Goddard