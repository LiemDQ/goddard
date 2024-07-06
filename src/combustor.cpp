#include "goddard/combustor.hpp"
#include "goddard/thermo_states.hpp"

namespace Goddard {

Combustor::Combustor(std::shared_ptr<Cantera::Solution> fuel, std::shared_ptr<Cantera::Solution> oxidizer, Eigen::ArrayXd&& pressures, MixtureRatios&& mr)
: fuel(std::move(fuel)), oxidizer(std::move(oxidizer)), pressures(pressures), mr(mr)
{
    feed.addPhase(fuel->thermo().get(), 1.0);
    feed.addPhase(oxidizer->thermo().get(), 1.0);

}

ThermoArray Combustor::solve() {
    auto mole_fracs = this->generate_mole_fraction_tensor();
    //todo: make "reactant" solution instead of using fuel solution 
    ThermoArray combustion_states = ThermoArray(fuel, {feed_temperatures.size(), pressures.size(), mr.molar_ratio().size() });
    combustion_states.TPX(feed_temperatures, pressures, mole_fracs);
    combustion_states.equilibrate("HP", "gibbs");

    return combustion_states;
}

/**
 * Create a n-D tensor with the mole fractions of each component participating in the reaction.
 */
std::vector<Eigen::ArrayXd> Combustor::generate_mole_fraction_tensor() {
    
    auto oxidizer_thermo = oxidizer->thermo();
    auto fuel_thermo = fuel->thermo();
    Eigen::ArrayXd mole_ratios = mr.molar_ratio() / (oxidizer_thermo->meanMolecularWeight() / fuel_thermo->meanMolecularWeight());
    Eigen::ArrayXd moles_ox = mole_ratios / (1 + mole_ratios);
    Eigen::ArrayXd moles_f = 1 - moles_ox;

    //this needs to include all species from both oxidizer and fuel
    auto mole_frac_tensor = std::vector<Eigen::ArrayXd>(mole_ratios.size(), Eigen::ArrayXd::Zero(feed.nSpecies())); 
    
    for (size_t i = 0; i < mole_ratios.size(); i++){
        double fuel_moles = moles_f[i];
        double ox_moles = moles_ox[i];
        feed.setPhaseMoles(0, fuel_moles);
        feed.setPhaseMoles(1, ox_moles);
        feed.getMoleFractions(mole_frac_tensor[i].data());
        
    }

    return mole_frac_tensor;
}

} //nameplace Goddard