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

/**
 * Thermodynamic constraint held fixed while the reactants burn to equilibrium.
 */
enum class CombustionProcess {
    ISOBARIC,  ///< Constant enthalpy and pressure (HP).
    ISOCHORIC, ///< Constant internal energy and specific volume (UV).
};

struct CombustorOptions {
    CombustorType type = CombustorType::INFINITE_AREA;
    MixtureRatioType mixture_type = MixtureRatioType::OF_RATIO;
    /**
     * Pressures [Pa]. For `CombustionProcess::ISOBARIC` these are the chamber pressures. For
     * `CombustionProcess::ISOCHORIC` they are the initial pressures of the unburnt reactants;
     * the chamber pressure is the result of the constant-volume combustion.
     */
    std::vector<double> pressures;
    double mass_flux = 0.0;
    double contraction_ratio = 0.0;
    /**
     * Combustion constraint. With `ISOCHORIC` in a `RocketProblem`, the constant-volume
     * equilibrium state is used as the stagnation state of a steady isentropic nozzle
     * expansion, so the reported performance is that idealization. It is not a
     * Chapman-Jouguet detonation model.
     */
    CombustionProcess process = CombustionProcess::ISOBARIC;
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

    /**
     * Construct from reactant streams given as `Gas` objects.
     *
     * The fuel and oxidizer streams carry their own temperature, pressure and composition, set
     * beforehand with `Gas::set_state_TPX`/`set_state_TPY`. Their species need not be product
     * species: element amounts and enthalpies are transferred to `products` by element name.
     *
     * @param products Product gas, which also defines the element set of the problem.
     * @param fuel Fuel stream at its own state.
     * @param oxidizer Oxidizer stream at its own state.
     *
     * @note Not implemented yet (work package C); `solve(pressures, mixture_ratios, options)`
     * throws `NotImplementedError`.
     */
    Combustor(Gas products, Gas fuel, Gas oxidizer);

    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures,
        const Eigen::ArrayXd& mixture_ratios, const CombustorOptions& options = {});
    ThermoArray solve(double fuel_temperature, double oxidizer_temperature,
        const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
        const CombustorOptions& options = {});

    /**
     * Burn the reactant streams given to the `Gas`-stream constructor over a grid of pressures
     * and mixture ratios.
     *
     * @param pressures Chamber pressures [Pa].
     * @param mixture_ratios Mixture ratios, interpreted according to `options.mixture_type`.
     * @param options Combustor settings.
     * @return Array of combustion states with shape (mixture ratios, pressures).
     *
     * @note Not implemented yet (work package C).
     */
    ThermoArray solve(const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
        const CombustorOptions& options = {});

    Eigen::ArrayXXd generate_mole_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const;
    Eigen::ArrayXXd generate_mass_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const;

    private:
    Composition m_fuel_composition;
    Composition m_oxidizer_composition;
    // Reactant streams of the `Gas`-stream constructor. `Gas` has no default constructor, so the
    // placeholder perfect-gas value is what the other constructors leave here; `m_use_reactant_gases`
    // says whether the streams are meaningful.
    Gas m_fuel_gas{1.4};
    Gas m_oxidizer_gas{1.4};
    bool m_use_reactant_gases = false;
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
