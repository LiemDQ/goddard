#pragma once
#include "cantera/core.h"

#include "eigen3/Eigen/Dense"
#include "goddard/gas.hpp"
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

enum class MixtureRatioType {
    FUEL_FRAC, // fuel fraction
    OF_RATIO,  // oxidizer-to-fuel
    PHI_RATIO, // equivalence ratio
};

struct CombustorOptions {
    CombustorType type = CombustorType::INFINITE_AREA;
    MixtureRatioType mixture_type = MixtureRatioType::OF_RATIO;
    std::vector<double> pressures;
    double mass_flux = 0.0;
    double contraction_ratio = 0.0;
};

/**
 * @brief Base class for isobaric combustion reactions.
 * Provides shared utilities for stream mixing and equilibration.
 */
class BaseCombustor {
    public:
    explicit BaseCombustor(Gas gas);
    virtual ~BaseCombustor() = default;

    inline std::vector<std::string> get_combustion_species() {return m_gas.thermo()->speciesNames();}

    protected:
    mutable Gas m_gas;

    ThermoArray combust(ThermoArray& states, const CombustorOptions& options);

    void set_mixture_composition(double value, MixtureRatioType type,
        const Composition& fuel, const Composition& oxidizer) const;
};

/**
 * @brief Handles isobaric combustion reactions with fuel and oxidizer streams.
 */
class Combustor : public BaseCombustor {
    public:
    Combustor(Gas gas, const std::string& fuel, const std::string& oxidizer);
    Combustor(Gas gas, const Composition& fuel, const Composition& oxidizer);

    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures,
        const Eigen::ArrayXd& mixture_ratios, const CombustorOptions& options = {});
    ThermoArray solve(double fuel_temperature, double oxidizer_temperature,
        const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
        const CombustorOptions& options = {});

    Eigen::ArrayXXd generate_mole_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const;
    Eigen::ArrayXXd generate_mass_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const;

    private:
    Composition m_fuel_composition;
    Composition m_oxidizer_composition;
};

/**
 * @brief Handles isobaric combustion with fuel, oxidizer, and recirculated flue gas streams.
 * The recirculation ratio is defined as r = m_flue / (m_fuel + m_ox).
 * The flue gas state is user-supplied (fixed), not iteratively solved.
 */
class DilutedCombustor : public BaseCombustor {
    public:
    DilutedCombustor(Gas gas, const std::string& fuel, const std::string& oxidizer, const std::string& flue);
    DilutedCombustor(Gas gas, const Composition& fuel, const Composition& oxidizer, const Composition& flue);

    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures,
        const Eigen::ArrayXd& mixture_ratios, double recirculation_ratio,
        const CombustorOptions& options = {});
    ThermoArray solve(double fuel_temperature, double oxidizer_temperature, double flue_temperature,
        const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
        double recirculation_ratio, const CombustorOptions& options = {});

    Eigen::ArrayXXd generate_mole_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type,
        double recirculation_ratio) const;
    Eigen::ArrayXXd generate_mass_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type,
        double recirculation_ratio) const;

    private:
    Composition m_fuel_composition;
    Composition m_oxidizer_composition;
    Composition m_flue_composition;
};

} //namespace Goddard
