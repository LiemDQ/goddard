#include "cantera/core.h"
#include "cantera/base/SolutionArray.h"
#include "goddard/equilibrium.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/global.hpp"
#include <iostream>

int main() {
    Goddard::setup_defaults();
    
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    double temp = 2400.0; //K
    double pressure = 50.0*Cantera::OneAtm;
    auto gas = sln->thermo();
    
    gas->setState_TPX(temp, pressure, "H2O:1, N2:1, O2:1. AR:0.1"); //completely random composition lol
    auto slnarr = Cantera::SolutionArray::create(sln, 10);


    auto coeffs = Goddard::get_stoichiometric_coeffs(*sln->thermo());

    auto manager = Goddard::ThermoArray(sln, std::vector<long>{1,1,1});

    return 0;
}