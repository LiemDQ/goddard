#pragma once
#include "cantera/core.h"

#include "eigen3/Eigen/Dense"
#include "goddard/gas.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"

#include <memory>
#include <optional>
#include <vector>
#include <utility>
#include <string>

namespace Goddard {

/** Chamber model of a rocket combustor. */
enum class CombustorType {
    /// Infinite-area chamber: the combustion state is the stagnation state of the nozzle.
    INFINITE_AREA,
    /// Finite-area chamber given by its mass flux mdot/A_c (`CombustorOptions::mass_flux`).
    FINITE_MASS_FLUX,
    /// Finite-area chamber given by its contraction ratio A_c/A_t (`CombustorOptions::contraction_ratio`).
    FINITE_CONTRACTION_RATIO
};

/** How a mixture ratio value is interpreted. */
enum class MixtureRatioType {
    /// Fuel mass fraction m_fuel/(m_fuel + m_ox) [-], in [0, 1].
    FUEL_FRAC,
    /// Oxidizer-to-fuel mass ratio m_ox/m_fuel [-].
    OF_RATIO,
    /// Equivalence ratio [-]: the fuel-to-oxidizer ratio over its stoichiometric value.
    PHI_RATIO,
};

/**
 * Thermodynamic constraint held fixed while the reactants burn to equilibrium.
 */
enum class CombustionProcess {
    ISOBARIC,  ///< Constant enthalpy and pressure (HP).
    ISOCHORIC, ///< Constant internal energy and specific volume (UV).
};

/** Settings of a combustor solve. */
struct CombustorOptions {
    /// Chamber model.
    CombustorType type = CombustorType::INFINITE_AREA;
    /// How the mixture ratios passed to the solver are interpreted.
    MixtureRatioType mixture_type = MixtureRatioType::OF_RATIO;
    /**
     * Pressures [Pa]. For `CombustionProcess::ISOBARIC` these are the chamber pressures. For
     * `CombustionProcess::ISOCHORIC` they are the initial pressures of the unburnt reactants;
     * the chamber pressure is the result of the constant-volume combustion. Read by
     * `RocketProblem`; the `Combustor::solve` overloads take their pressures as an argument.
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
    std::vector<double> pressures;
    /** Mass flux through a `FINITE_MASS_FLUX` chamber, mdot/A_c [kg/(m^2 s)]. */
    double mass_flux = 0.0;
    /** Contraction ratio A_c/A_t [-] of a `FINITE_CONTRACTION_RATIO` chamber. */
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
 * @brief Base class for isobaric (HP) and isochoric (UV) combustion.
 * Provides shared utilities for stream mixing and equilibration.
 */
class BaseCombustor {
    public:
    /**
     * @param gas product gas. The solver works on its own copy of `gas`; the caller's `Gas` is
     * not modified.
     */
    explicit BaseCombustor(const Gas& gas);
    virtual ~BaseCombustor() = default;

    /** Names of the gas-phase product species, in phase order. */
    inline std::vector<std::string> get_combustion_species() {return m_gas.thermo()->speciesNames();}

    protected:
    mutable Gas m_gas;

    /**
     * Equilibrate reactant states in place and return them.
     *
     * @throws NotImplementedError for `CombustionProcess::ISOCHORIC` when `states` carries
     * candidate condensed species.
     */
    ThermoArray combust(ThermoArray& states, const CombustorOptions& options);

    /**
     * Array of reactant states with the given shape, carrying the candidate condensed species of
     * `m_gas` with every amount set to zero (the reactants are all gas phase). Sets the condensed
     * amounts of `m_gas` to zero as well.
     */
    ThermoArray reactant_states(const std::vector<long>& shape) const;

    /**
     * Check that the combustor options describe a supported combustor.
     *
     * The combustor itself is the same for every `CombustorType`: for `FINITE_MASS_FLUX` and
     * `FINITE_CONTRACTION_RATIO` it produces the injector-face state, and the finite-area chamber
     * is solved afterwards by `Nozzle::solve_finite_area_chamber`.
     *
     * @throws std::invalid_argument for a finite-area type with `contraction_ratio` <= 1 or
     * `mass_flux` <= 0, or combined with `CombustionProcess::ISOCHORIC`.
     */
    static void validate_options(const CombustorOptions& options);

    /**
     * Reject an empty input grid before any entry of it is read.
     *
     * @param values Grid of input values, e.g. pressures [Pa] or mixture ratios [-].
     * @param name Name of the grid, for the error message.
     * @throws std::invalid_argument if `values` is empty.
     */
    static void require_nonempty(const Eigen::ArrayXd& values, const std::string& name);

    void set_mixture_composition(double value, MixtureRatioType type,
        const Composition& fuel, const Composition& oxidizer) const;
};

/**
 * @brief Combustion of a fuel and an oxidizer.
 *
 * The reactants are given in one of two ways, and each has its own `solve`:
 * - as compositions of product species (the first two constructors), burnt with the
 *   `solve` overloads that take reactant temperatures;
 * - as reactant `Gas` streams carrying their own state (the third constructor), burnt with the
 *   `solve` overload that takes only pressures and mixture ratios.
 */
class Combustor : public BaseCombustor {
    public:
    /**
     * Construct from fuel and oxidizer compositions, given as product species of `gas`.
     * Solve with the `solve` overloads that take reactant temperatures.
     */
    Combustor(Gas gas, const std::string& fuel, const std::string& oxidizer);
    /**
     * Construct from fuel and oxidizer mole-fraction maps of product species of `gas`.
     * Solve with the `solve` overloads that take reactant temperatures.
     *
     * @param gas Product gas. The solver works on its own copy.
     * @param fuel Fuel composition, mole fractions by species name.
     * @param oxidizer Oxidizer composition, mole fractions by species name.
     */
    Combustor(Gas gas, const Composition& fuel, const Composition& oxidizer);

    /**
     * Construct from reactant streams given as `Gas` objects.
     *
     * The fuel and oxidizer streams carry their own temperature, pressure and composition, set
     * beforehand with `Gas::set_state_TPX`/`set_state_TPY`. Their species need not be product
     * species: element amounts [kmol per kg of stream] and specific enthalpies [J/kg] are
     * transferred to `products` by element name. A stream built from a condensed reactant such as
     * `H2(L)` or `RP-1` is an ideal-gas phase whose NASA9 polynomials give the correct enthalpy
     * and element amounts; its density and entropy are meaningless and are only read on the
     * isochoric path, which is therefore restricted to gaseous reactants.
     *
     * @param products Product gas, which also defines the element set of the problem.
     * @param fuel Fuel stream at its own state.
     * @param oxidizer Oxidizer stream at its own state.
     */
    Combustor(Gas products, Gas fuel, Gas oxidizer);

    /**
     * Burn the fuel and oxidizer compositions given to the composition constructor over a grid
     * of reactant temperatures, pressures and mixture ratios.
     *
     * @param temperatures Reactant temperatures [K], applied to both streams.
     * @param pressures Pressures [Pa]: chamber pressures for `ISOBARIC`, initial reactant
     *                  pressures for `ISOCHORIC`.
     * @param mixture_ratios Mixture ratios, interpreted according to `options.mixture_type`.
     * @param options Combustor settings.
     * @return Array of combustion states with shape (temperatures, pressures, mixture ratios).
     */
    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures,
        const Eigen::ArrayXd& mixture_ratios, const CombustorOptions& options = {});
    /**
     * Burn the fuel and oxidizer compositions given to the composition constructor, with the
     * fuel and the oxidizer at their own temperatures, over a grid of pressures and mixture
     * ratios.
     *
     * @param fuel_temperature Fuel temperature [K].
     * @param oxidizer_temperature Oxidizer temperature [K].
     * @param pressures Pressures [Pa]: chamber pressures for `ISOBARIC`, initial reactant
     *                  pressures for `ISOCHORIC`.
     * @param mixture_ratios Mixture ratios, interpreted according to `options.mixture_type`.
     * @param options Combustor settings.
     * @return Array of combustion states with shape (mixture ratios, pressures).
     */
    ThermoArray solve(double fuel_temperature, double oxidizer_temperature,
        const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
        const CombustorOptions& options = {});

    /**
     * Burn the reactant streams given to the `Gas`-stream constructor over a grid of pressures
     * and mixture ratios.
     *
     * For every mixture ratio the fuel mass fraction f follows from `options.mixture_type`
     * (`OF_RATIO`: f = 1/(1+O/F); `FUEL_FRAC`: f is the value itself). The element amounts
     * b [kmol/kg of mixture] and the specific enthalpy h [J/kg] of the mixture are the mass
     * blends f*x_fuel + (1-f)*x_oxidizer of the two stream states, with the element amounts
     * mapped by element name onto the element order of the product gas. Each combination is then
     * equilibrated at constant enthalpy and pressure (`ISOBARIC`) or at constant internal energy
     * and specific volume (`ISOCHORIC`), warm-started from the previous pressure of the same
     * mixture ratio.
     *
     * @param pressures Pressures [Pa]: chamber pressures for `ISOBARIC`, initial reactant
     *                  pressures for `ISOCHORIC`.
     * @param mixture_ratios Mixture ratios, interpreted according to `options.mixture_type`.
     * @param options Combustor settings.
     * @return Array of combustion states with shape (mixture ratios, pressures).
     *
     * @throws NotImplementedError if the combustor was not built from reactant `Gas` streams, if
     * `options.mixture_type` is `PHI_RATIO` (CEA's valence rule is not implemented on this path),
     * or if the process is `ISOCHORIC` while the
     * product gas carries candidate condensed species.
     * @throws std::invalid_argument if a reactant contains an element the product gas does not have.
     * @throws std::invalid_argument if `pressures` or `mixture_ratios` is empty.
     */
    ThermoArray solve(const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
        const CombustorOptions& options = {});

    /**
     * Unburnt reactant mole fractions [-] for each mixture ratio, from the compositions given to
     * the composition constructor.
     *
     * @param mixture_ratios Mixture ratios, interpreted according to `type`.
     * @param type How `mixture_ratios` are interpreted.
     * @return Array with one row per mixture ratio and one column per product species.
     * @throws std::invalid_argument if a mixture ratio is out of range or gives negative fractions.
     */
    Eigen::ArrayXXd generate_mole_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const;
    /**
     * Unburnt reactant mass fractions [-] for each mixture ratio; see
     * `generate_mole_fraction_matrix`.
     */
    Eigen::ArrayXXd generate_mass_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const;

    private:
    Composition m_fuel_composition;
    Composition m_oxidizer_composition;
    // Reactant streams of the `Gas`-stream constructor; empty for the product-species constructors.
    std::optional<Gas> m_fuel_gas;
    std::optional<Gas> m_oxidizer_gas;

    /**
     * Element amounts of a reactant stream [kmol per kg of stream] re-indexed onto the element
     * order of the product gas. Product elements the stream does not have are zero.
     */
    Eigen::ArrayXd stream_element_moles(const Gas& stream) const;
};

/**
 * @brief Combustion of fuel, oxidizer and recirculated flue gas streams.
 * The recirculation ratio is defined as r = m_flue / (m_fuel + m_ox).
 * The flue gas state is user-supplied (fixed), not iteratively solved.
 */
class DilutedCombustor : public BaseCombustor {
    public:
    /**
     * Construct from fuel, oxidizer and flue gas compositions given as mole-fraction strings of
     * product species of `gas`, e.g. "H2:1".
     */
    DilutedCombustor(Gas gas, const std::string& fuel, const std::string& oxidizer, const std::string& flue);
    /**
     * Construct from fuel, oxidizer and flue gas mole-fraction maps of product species of `gas`.
     *
     * @param gas Product gas. The solver works on its own copy.
     * @param fuel Fuel composition, mole fractions by species name.
     * @param oxidizer Oxidizer composition, mole fractions by species name.
     * @param flue Recirculated flue gas composition, mole fractions by species name.
     */
    DilutedCombustor(Gas gas, const Composition& fuel, const Composition& oxidizer, const Composition& flue);

    /**
     * Burn the diluted mixture with every stream at the same temperature, over a grid of
     * temperatures, pressures and mixture ratios.
     *
     * @param temperatures Reactant temperatures [K], applied to all three streams.
     * @param pressures Pressures [Pa]: chamber pressures for `ISOBARIC`, initial reactant
     *                  pressures for `ISOCHORIC`.
     * @param mixture_ratios Fuel/oxidizer mixture ratios of the fresh feed, interpreted according
     *                       to `options.mixture_type`.
     * @param recirculation_ratio Flue gas mass over fresh feed mass, m_flue/(m_fuel + m_ox) [-].
     * @param options Combustor settings.
     * @return Array of combustion states with shape (temperatures, pressures, mixture ratios).
     */
    ThermoArray solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures,
        const Eigen::ArrayXd& mixture_ratios, double recirculation_ratio,
        const CombustorOptions& options = {});
    /**
     * Burn the diluted mixture with each stream at its own temperature, over a grid of
     * pressures and mixture ratios. The mixture enthalpy is the mass-weighted blend of the
     * stream enthalpies.
     *
     * @param fuel_temperature Fuel temperature [K].
     * @param oxidizer_temperature Oxidizer temperature [K].
     * @param flue_temperature Flue gas temperature [K].
     * @param pressures Pressures [Pa]: chamber pressures for `ISOBARIC`, initial reactant
     *                  pressures for `ISOCHORIC`.
     * @param mixture_ratios Fuel/oxidizer mixture ratios of the fresh feed, interpreted according
     *                       to `options.mixture_type`.
     * @param recirculation_ratio Flue gas mass over fresh feed mass, m_flue/(m_fuel + m_ox) [-].
     * @param options Combustor settings.
     * @return Array of combustion states with shape (mixture ratios, pressures).
     */
    ThermoArray solve(double fuel_temperature, double oxidizer_temperature, double flue_temperature,
        const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
        double recirculation_ratio, const CombustorOptions& options = {});

    /**
     * Unburnt mole fractions [-] of the fresh feed blended with flue gas, for each mixture ratio.
     *
     * @param mixture_ratios Fuel/oxidizer mixture ratios of the fresh feed, interpreted according
     *                       to `type`.
     * @param type How `mixture_ratios` are interpreted.
     * @param recirculation_ratio Flue gas mass over fresh feed mass [-].
     * @return Array with one row per mixture ratio and one column per product species.
     * @throws std::invalid_argument if a mixture ratio is out of range or gives negative fractions.
     */
    Eigen::ArrayXXd generate_mole_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type,
        double recirculation_ratio) const;
    /**
     * Unburnt mass fractions [-] of the fresh feed blended with flue gas; see
     * `generate_mole_fraction_matrix`.
     */
    Eigen::ArrayXXd generate_mass_fraction_matrix(
        const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type,
        double recirculation_ratio) const;

    private:
    Composition m_fuel_composition;
    Composition m_oxidizer_composition;
    Composition m_flue_composition;
};

} //namespace Goddard
