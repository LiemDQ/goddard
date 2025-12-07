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
    FINITE_CONTRACTION_RATIO,
    NONE
};

struct CombustorOptions {
    CombustorType type;
    std::vector<double> pressures;
    double mass_flux;
    double contraction_ratio;
};

/**
 * @brief Handles isobaric combustion reactions.
 */
class Combustor {
    public:
    //this may lead to lots of unnecessary copies
    Combustor(
        std::shared_ptr<Cantera::Solution> thermo,
        std::vector<double>& fuel_state, 
        std::vector<double>& oxidizer_state
    );


    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options = {});
    ThermoArray solve(const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options = {});

    Eigen::ArrayXXd generate_mole_fraction_matrix(const MixtureRatios& mr) const;
    Eigen::ArrayXXd generate_mass_fraction_matrix(const MixtureRatios& mr) const;

    inline std::vector<std::string> get_combustion_species() {return m_thermo->thermo()->speciesNames();}
    std::vector<double> fuel_state, oxidizer_state;

    protected:
    std::shared_ptr<Cantera::Solution> m_thermo;
    
    void assign_mole_frac_row_entries(
        Eigen::ArrayXXd& matrix, 
        long row_idx, 
        const Cantera::Composition& composition, 
        double coeff = 1.0) const;
    
    ThermoArray combust(ThermoArray& states, const CombustorOptions& options);

};

} //namespace Goddard