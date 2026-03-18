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
 * @brief Base class for isobaric combustion reactions.
 * Provides shared utilities for stream mixing and equilibration.
 */
class BaseCombustor {
    public:
    explicit BaseCombustor(std::shared_ptr<Cantera::Solution> thermo);
    virtual ~BaseCombustor() = default;

    inline std::vector<std::string> get_combustion_species() {return m_thermo->thermo()->speciesNames();}

    protected:
    std::shared_ptr<Cantera::Solution> m_thermo;

    void assign_mole_frac_row_entries(
        Eigen::ArrayXXd& matrix,
        long row_idx,
        const Cantera::Composition& composition,
        double coeff = 1.0) const;

    ThermoArray combust(ThermoArray& states, const CombustorOptions& options);
};

/**
 * @brief Handles isobaric combustion reactions with fuel and oxidizer streams.
 */
class Combustor : public BaseCombustor {
    public:
    Combustor(
        std::shared_ptr<Cantera::Solution> thermo,
        std::vector<double>& fuel_state,
        std::vector<double>& oxidizer_state
    );

    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options = {});
    ThermoArray solve(const Eigen::ArrayXd& pressures, const MixtureRatios& mr, const CombustorOptions& options = {});

    Eigen::ArrayXXd generate_mole_fraction_matrix(const MixtureRatios& mr) const;
    Eigen::ArrayXXd generate_mass_fraction_matrix(const MixtureRatios& mr) const;

    std::vector<double> fuel_state, oxidizer_state;
};

/**
 * @brief Handles isobaric combustion with fuel, oxidizer, and recirculated flue gas streams.
 * The recirculation ratio is defined as r = m_flue / (m_fuel + m_ox).
 * The flue gas state is user-supplied (fixed), not iteratively solved.
 */
class DilutedCombustor : public BaseCombustor {
    public:
    DilutedCombustor(
        std::shared_ptr<Cantera::Solution> thermo,
        std::vector<double>& fuel_state,
        std::vector<double>& oxidizer_state,
        std::vector<double>& flue_state
    );

    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures, const MixtureRatios& mr, double recirculation_ratio, const CombustorOptions& options = {});
    ThermoArray solve(const Eigen::ArrayXd& pressures, const MixtureRatios& mr, double recirculation_ratio, const CombustorOptions& options = {});

    Eigen::ArrayXXd generate_mole_fraction_matrix(const MixtureRatios& mr, double recirculation_ratio) const;
    Eigen::ArrayXXd generate_mass_fraction_matrix(const MixtureRatios& mr, double recirculation_ratio) const;

    std::vector<double> fuel_state, oxidizer_state, flue_state;
};

} //namespace Goddard
