#pragma once
#include "cantera/core.h"

#include "eigen3/Eigen/Dense"
#include "goddard/mixture_ratio.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"
#include <memory>
#include <vector>
#include <utility>
#include <string>

namespace Goddard {

enum class CombustorType {
    INFINITE_AREA,
    FINITE_MASS_FLUX,
    FINITE_CONTRACTION_RATIO
};

struct CombustionOptions {
    CombustorType combustor_type = CombustorType::INFINITE_AREA;
    bool include_transport = false;
    bool include_ionized = false;
    double reltol = DEFAULT_RELTOL;
};

/**
 * @brief Handles isobaric combustion reactions.
 */
class Combustor {
    public:
    //this may lead to lots of unnecessary copies
    Combustor(
        std::shared_ptr<Cantera::Solution>& fuel, 
        std::shared_ptr<Cantera::Solution>& oxidizer, 
        std::shared_ptr<Cantera::Solution>& products
    );


    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustionOptions& options = {});
    Eigen::ArrayXXd generate_mole_fraction_matrix(const MixtureRatios& mr);
    

    inline std::vector<std::string> get_combustion_species() {return product_sln->thermo()->speciesNames();}
    inline std::shared_ptr<Cantera::Solution> get_fuel() {return fuel_sln;}
    inline std::shared_ptr<Cantera::Solution> get_oxidizer() {return oxidizer_sln;}
    inline std::shared_ptr<Cantera::Solution> get_products() { return product_sln; }

    protected:
    
    std::shared_ptr<Cantera::Solution> fuel_sln;
    std::shared_ptr<Cantera::Solution> oxidizer_sln;
    std::shared_ptr<Cantera::Solution> product_sln;

    
    void assign_mole_frac_row_entries(
        Eigen::ArrayXXd& matrix, 
        long row_idx, 
        const Cantera::Composition& composition, 
        double coeff = 1.0);


};

} //namespace Goddard