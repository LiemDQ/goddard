#pragma once
#include "cantera/core.h"
#include <memory>
namespace Goddard {

class MixtureRatio {
    public: 
    MixtureRatio(double OF, double M_fuel, double M_oxidizer): OF_ratio(OF), M_fuel(M_fuel), M_ox(M_oxidizer) {}
    MixtureRatio(double OF, 
        const std::shared_ptr<Cantera::Solution>& fuel, 
        const std::shared_ptr<Cantera::Solution>& oxidizer);
    inline double fuel_mole_frac() {return 1 - molar_ratio() / (1 + molar_ratio());}
    inline double oxidizer_mole_frac() { return molar_ratio() / (1 + molar_ratio());}
    inline double OF_to_molar_ratio(double OF) {return OF/(M_ox / M_fuel);}
    inline double molar_ratio() {return OF_to_molar_ratio(OF_ratio);}

    private:
    
    double OF_ratio;
    double M_fuel;
    double M_ox;
    
    // double phi() {}
    // double equiv_ratio() {}
    // double fuel_percent() {}

};

}