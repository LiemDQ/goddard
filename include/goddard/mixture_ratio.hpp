#pragma once
#include "cantera/core.h"
#include <memory>
#include <vector>
#include <eigen3/Eigen/Core>
namespace Goddard {

class MixtureRatio {
    public: 
    MixtureRatio(double OF, double fuel_molar_mass, double oxidizer_molar_mass): OF_ratio(OF), M_fuel(fuel_molar_mass), M_ox(oxidizer_molar_mass) {}
    MixtureRatio(double OF, 
        const std::shared_ptr<Cantera::Solution>& fuel, 
        const std::shared_ptr<Cantera::Solution>& oxidizer);
        
    inline double fuel_mole_frac() {return 1 - molar_ratio() / (1 + molar_ratio());}
    inline double oxidizer_mole_frac() { return molar_ratio() / (1 + molar_ratio());}
    inline double OF_to_molar_ratio(double OF) const {return OF/(M_ox / M_fuel);}
    inline double molar_ratio() const {return OF_to_molar_ratio(OF_ratio);}

    private:
    double OF_ratio;
    double M_fuel;
    double M_ox;
    // double phi() {}
    // double equiv_ratio() {}
    // double fuel_percent() {}

};

class MixtureRatios {
    public: 
    MixtureRatios(double OF, double fuel_molar_mass, double oxidizer_molar_mass): OF_ratio(1), M_fuel(fuel_molar_mass), M_ox(oxidizer_molar_mass) { OF_ratio << OF;}
    MixtureRatios(const Eigen::ArrayXd& OF, double fuel_molar_mass, double oxidizer_molar_mass): OF_ratio(OF), M_fuel(fuel_molar_mass), M_ox(oxidizer_molar_mass) {}

    MixtureRatios(double OF, 
        const std::shared_ptr<Cantera::Solution>& fuel, 
        const std::shared_ptr<Cantera::Solution>& oxidizer);
    MixtureRatios(const Eigen::ArrayXd& OF, 
        const std::shared_ptr<Cantera::Solution>& fuel, 
        const std::shared_ptr<Cantera::Solution>& oxidizer);
        
    inline Eigen::ArrayXd fuel_mole_frac() {return 1.0 - molar_ratio() / (1.0 + molar_ratio());}
    inline Eigen::ArrayXd oxidizer_mole_frac() { return molar_ratio() / (1 + molar_ratio());}
    inline Eigen::ArrayXd OF_to_molar_ratio(const Eigen::ArrayXd& OF) const {return OF/(M_ox / M_fuel);}
    inline Eigen::ArrayXd molar_ratio() const {return OF_to_molar_ratio(OF_ratio);}

    private:
    Eigen::ArrayXd OF_ratio;
    double M_fuel;
    double M_ox;
    // double phi() {}
    // double equiv_ratio() {}
    // double fuel_percent() {}

};

}