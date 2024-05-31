#include "goddard/mixture_ratio.hpp"

namespace Goddard {

MixtureRatio::MixtureRatio(double OF, 
        const std::shared_ptr<Cantera::Solution>& fuel, 
        const std::shared_ptr<Cantera::Solution>& oxidizer)
        : OF_ratio(OF) 
        {
            this->M_fuel = fuel->thermo()->meanMolecularWeight();
            this->M_ox = oxidizer->thermo()->meanMolecularWeight();
        }

}