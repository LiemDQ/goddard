#include "goddard/mixture_ratio.hpp"

namespace Goddard {

MixtureRatio::MixtureRatio(double OF, 
        const std::shared_ptr<Cantera::Solution>& fuel, 
        const std::shared_ptr<Cantera::Solution>& oxidizer)
        : OF_ratio(OF) 
        {
            M_fuel = fuel->thermo()->meanMolecularWeight();
            M_ox = oxidizer->thermo()->meanMolecularWeight();
        }

MixtureRatios::MixtureRatios(double OF, 
        const std::shared_ptr<Cantera::Solution>& fuel, 
        const std::shared_ptr<Cantera::Solution>& oxidizer)
        : OF_ratio(1) 
        {
            OF_ratio << OF;
            M_fuel = fuel->thermo()->meanMolecularWeight();
            M_ox = oxidizer->thermo()->meanMolecularWeight();
        }


MixtureRatios::MixtureRatios(const Eigen::ArrayXd& OF,
    const std::shared_ptr<Cantera::Solution>& fuel, 
    const std::shared_ptr<Cantera::Solution>& oxidizer)
    : OF_ratio(OF) 
    {
        M_fuel = fuel->thermo()->meanMolecularWeight();
        M_ox = oxidizer->thermo()->meanMolecularWeight();
    }

} //namespace Goddard