#pragma once
#include "cantera/core.h"

#include "eigen3/Eigen/Dense"
#include "goddard/mixture_ratio.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"
#include "goddard/case_options.hpp"

#include <memory>
#include <vector>
#include <utility>
#include <string>

namespace Goddard {

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


    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options = {});
    Eigen::ArrayXXd generate_mole_fraction_matrix(const MixtureRatios& mr) const;
    

    inline std::vector<std::string> get_combustion_species() {return m_product_sln->thermo()->speciesNames();}
    inline std::shared_ptr<Cantera::Solution> get_fuel() {return m_fuel_sln;}
    inline std::shared_ptr<Cantera::Solution> get_oxidizer() {return m_oxidizer_sln;}
    inline std::shared_ptr<Cantera::Solution> get_products() { return m_product_sln; }

    protected:
    
    std::shared_ptr<Cantera::Solution> m_fuel_sln;
    std::shared_ptr<Cantera::Solution> m_oxidizer_sln;
    std::shared_ptr<Cantera::Solution> m_product_sln;

    
    void assign_mole_frac_row_entries(
        Eigen::ArrayXXd& matrix, 
        long row_idx, 
        const Cantera::Composition& composition, 
        double coeff = 1.0) const;


};

} //namespace Goddard