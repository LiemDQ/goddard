#include "goddard/combustor.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/error.hpp"
#include "cantera/core.h"
#include "cantera/base/stringUtils.h"
#include <algorithm>
#include <iterator>
#include <utility>
#include <iostream>
#include <cassert>
namespace Goddard {

// ---- BaseCombustor ----

BaseCombustor::BaseCombustor(Gas gas)
    : m_gas(std::move(gas))
{}

void BaseCombustor::validate_options(const CombustorOptions& options) {
    switch (options.type) {
        case CombustorType::INFINITE_AREA:
            break;
        case CombustorType::FINITE_MASS_FLUX:
            if (options.mass_flux <= 0.0) {
                throw std::invalid_argument(
                    "BaseCombustor: FINITE_MASS_FLUX requires mass_flux > 0.");
            }
            if (options.process == CombustionProcess::ISOCHORIC) {
                throw std::invalid_argument(
                    "BaseCombustor: finite-area combustors do not support ISOCHORIC combustion.");
            }
            break;
        case CombustorType::FINITE_CONTRACTION_RATIO:
            if (options.contraction_ratio <= 1.0) {
                throw std::invalid_argument(
                    "BaseCombustor: FINITE_CONTRACTION_RATIO requires contraction_ratio > 1.");
            }
            if (options.process == CombustionProcess::ISOCHORIC) {
                throw std::invalid_argument(
                    "BaseCombustor: finite-area combustors do not support ISOCHORIC combustion.");
            }
            break;
        case CombustorType::NONE:
            throw NotImplementedError("CombustorType::NONE is not implemented.");
    }
}

ThermoArray BaseCombustor::combust(ThermoArray& states, const CombustorOptions& options) {
    validate_options(options);
    switch (options.process) {
        case CombustionProcess::ISOBARIC:
            try {
                states.equilibrate("HP", "gibbs");
            } catch (const Cantera::CanteraError&) {
                if (states.num_condensed() > 0) {
                    throw;
                }
                // "gibbs" tests its enthalpy residual relative to H and fails when the target
                // enthalpy is close to zero (e.g. gaseous H2/O2 at 298.15 K); see
                // `Gas::equilibrate_HP`. Locations already solved are HP equilibria, so
                // re-solving them with "vcs" leaves them unchanged.
                states.equilibrate("HP", "vcs");
            }
            break;
        case CombustionProcess::ISOCHORIC:
            // Cantera's "gibbs" (MultiPhaseEquil) solver does not support UV, and
            // "element_potential" fails to converge for H2/O2 at initial pressures of
            // ~10 bar and above. "vcs" converges across O/F 1-32 and 1e3-1e7 Pa.
            states.equilibrate("UV", "vcs");
            break;
    }
    return states;
}

void BaseCombustor::set_mixture_composition(double value, MixtureRatioType type,
    const Composition& fuel, const Composition& oxidizer) const {
    switch (type) {
        case MixtureRatioType::OF_RATIO:
            m_gas.set_OF_ratio(value, fuel, oxidizer);
            break;
        case MixtureRatioType::PHI_RATIO:
            m_gas.set_equivalence_ratio(value, fuel, oxidizer);
            break;
        case MixtureRatioType::FUEL_FRAC:
            m_gas.set_fuel_fraction(value, fuel, oxidizer);
            break;
    }
}

// ---- Combustor ----

Combustor::Combustor(Gas gas, const std::string& fuel_comp, const std::string& ox_comp)
    : BaseCombustor(std::move(gas))
{
    m_fuel_composition = Cantera::parseCompString(fuel_comp, m_gas.species_names());
    m_oxidizer_composition = Cantera::parseCompString(ox_comp, m_gas.species_names());
}

Combustor::Combustor(Gas gas, const Composition& fuel, const Composition& oxidizer)
    : BaseCombustor(std::move(gas)),
      m_fuel_composition(fuel), m_oxidizer_composition(oxidizer)
{}

Combustor::Combustor(Gas products, Gas fuel, Gas oxidizer)
    : BaseCombustor(std::move(products)),
      m_fuel_gas(std::move(fuel)), m_oxidizer_gas(std::move(oxidizer))
{}

Eigen::ArrayXd Combustor::stream_element_moles(const Gas& stream) const {
    const std::vector<std::string> product_elements = m_gas.element_names();
    const std::vector<std::string> stream_elements = stream.element_names();
    const Eigen::ArrayXd stream_moles = stream.element_moles();

    Eigen::ArrayXd mapped = Eigen::ArrayXd::Zero(static_cast<long>(product_elements.size()));
    for (size_t m = 0; m < stream_elements.size(); m++) {
        const double amount = stream_moles(static_cast<long>(m));
        const auto match = std::find(product_elements.begin(), product_elements.end(),
            stream_elements[m]);
        if (match == product_elements.end()) {
            // Elements the reactant phase declares but does not actually contain are ignored:
            // a stream built from a shared species list carries the elements of every listed
            // species, not only of the ones it is made of.
            if (amount > 0.0) {
                throw FmtError(
                    "Combustor: reactant stream contains element '{}', which is absent from the "
                    "product gas.", stream_elements[m]);
            }
            continue;
        }
        mapped(static_cast<long>(std::distance(product_elements.begin(), match))) = amount;
    }
    return mapped;
}

namespace {

/** Mass fraction of fuel in the mixture from a mixture ratio and its interpretation. */
double fuel_mass_fraction(double mixture_ratio, MixtureRatioType type) {
    switch (type) {
        case MixtureRatioType::OF_RATIO:
            return 1.0 / (1.0 + mixture_ratio);
        case MixtureRatioType::FUEL_FRAC:
            return mixture_ratio;
        case MixtureRatioType::PHI_RATIO:
            break;
    }
    throw NotImplementedError(
        "Combustor: equivalence ratios are not implemented for reactant Gas streams; CEA's "
        "valence rule is needed to define the stoichiometric ratio of an arbitrary reactant.");
}

} // namespace

ThermoArray Combustor::solve(const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
    const CombustorOptions& options) {

    if (!m_fuel_gas || !m_oxidizer_gas) {
        throw NotImplementedError(
            "Combustor::solve(pressures, mixture_ratios) requires the combustor to be built from "
            "reactant Gas streams. Use solve(fuel_T, oxidizer_T, pressures, mixture_ratios) for "
            "reactants given as product-species compositions.");
    }
    validate_options(options);
    if (options.process == CombustionProcess::ISOCHORIC && m_gas.has_condensed_candidates()) {
        throw NotImplementedError(
            "Constant-volume combustion with candidate condensed species is not implemented.");
    }

    const Eigen::ArrayXd fuel_elements = stream_element_moles(*m_fuel_gas);
    const Eigen::ArrayXd oxidizer_elements = stream_element_moles(*m_oxidizer_gas);
    const double fuel_enthalpy = m_fuel_gas->enthalpy_mass();
    const double oxidizer_enthalpy = m_oxidizer_gas->enthalpy_mass();
    // Specific volume and internal energy of each stream at its own state. Only meaningful for
    // gaseous reactants, and only read on the isochoric path.
    const double fuel_volume = 1.0 / m_fuel_gas->density();
    const double oxidizer_volume = 1.0 / m_oxidizer_gas->density();
    const double fuel_energy = fuel_enthalpy - m_fuel_gas->pressure() * fuel_volume;
    const double oxidizer_energy =
        oxidizer_enthalpy - m_oxidizer_gas->pressure() * oxidizer_volume;

    const long n_ratios = mixture_ratios.size();
    const long n_pressures = pressures.size();
    ThermoArray combustion_states(m_gas, {n_ratios, n_pressures});

    for (long i = 0; i < n_ratios; i++) {
        const double fuel_fraction = fuel_mass_fraction(mixture_ratios(i), options.mixture_type);
        const Eigen::ArrayXd element_moles =
            fuel_fraction * fuel_elements + (1.0 - fuel_fraction) * oxidizer_elements;
        const double enthalpy =
            fuel_fraction * fuel_enthalpy + (1.0 - fuel_fraction) * oxidizer_enthalpy;

        double temperature_guess = m_gas.equilibrium_options.T_default;
        for (long j = 0; j < n_pressures; j++) {
            m_gas.set_element_moles(element_moles, temperature_guess, pressures(j));
            switch (options.process) {
                case CombustionProcess::ISOBARIC:
                    m_gas.equilibrate_HP(enthalpy, pressures(j));
                    break;
                case CombustionProcess::ISOCHORIC: {
                    const double energy =
                        fuel_fraction * fuel_energy + (1.0 - fuel_fraction) * oxidizer_energy;
                    const double volume =
                        fuel_fraction * fuel_volume + (1.0 - fuel_fraction) * oxidizer_volume;
                    m_gas.thermo()->setState_UV(energy, volume);
                    // See `BaseCombustor::combust` for why the UV problem uses the "vcs" solver.
                    m_gas.equilibrate("UV", "vcs");
                    break;
                }
            }
            // Warm start of the next pressure at the same mixture ratio.
            temperature_guess = m_gas.temperature();
            combustion_states.set_state(combustion_states.flat_index(i, j), m_gas.save_state());
        }
    }

    return combustion_states;
}

ThermoArray Combustor::solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures,
    const Eigen::ArrayXd& mixture_ratios, const CombustorOptions& options) {

    Eigen::ArrayXXd mole_fracs = generate_mole_fraction_matrix(mixture_ratios, options.mixture_type);

    ThermoArray combustion_states(m_gas.solution(), {temperatures.size(), pressures.size(), mixture_ratios.size()});
    combustion_states.TPX(temperatures, pressures, mole_fracs);

    return combust(combustion_states, options);
}

ThermoArray Combustor::solve(double fuel_temperature, double oxidizer_temperature,
    const Eigen::ArrayXd& pressures, const Eigen::ArrayXd& mixture_ratios,
    const CombustorOptions& options) {

    auto thermo = m_gas.thermo();
    MixtureRatioType type = options.mixture_type;

    // Compute reactant enthalpies at their respective temperatures
    // Use first pressure as reference (enthalpy is P-independent for ideal gas)
    double ref_pressure = pressures[0];

    thermo->setState_TPX(fuel_temperature, ref_pressure, m_fuel_composition);
    double fuel_enthalpy = thermo->enthalpy_mass();

    thermo->setState_TPX(oxidizer_temperature, ref_pressure, m_oxidizer_composition);
    double oxidizer_enthalpy = thermo->enthalpy_mass();

    Eigen::ArrayXXd mass_fracs = generate_mass_fraction_matrix(mixture_ratios, type);

    // Compute fuel mass fraction for each mixture ratio to blend enthalpies
    long n_compositions = mixture_ratios.size();
    long n_pressures = pressures.size();
    Eigen::ArrayXd enthalpies(n_compositions);

    for (long i = 0; i < n_compositions; i++) {
        set_mixture_composition(mixture_ratios[i], type, m_fuel_composition, m_oxidizer_composition);
        double fuel_mass_frac = thermo->mixtureFraction(
            m_fuel_composition, m_oxidizer_composition, Cantera::ThermoBasis::mass);
        enthalpies[i] = fuel_enthalpy * fuel_mass_frac + oxidizer_enthalpy * (1.0 - fuel_mass_frac);
    }

    ThermoArray combustion_states(m_gas.solution(), {n_compositions, n_pressures});

    for (long i = 0; i < n_compositions; i++) {
        Eigen::ArrayXd row = mass_fracs.row(i);
        for (long j = 0; j < n_pressures; j++) {
            thermo->setMassFractions(row.data());
            thermo->setState_HP(enthalpies[i], pressures[j]);
            combustion_states.set_state(combustion_states.flat_index(i, j), m_gas.save_state());
        }
    }

    return combust(combustion_states, options);
}

Eigen::ArrayXXd Combustor::generate_mole_fraction_matrix(
    const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const {

    auto thermo = m_gas.thermo();
    long n_species = static_cast<long>(thermo->nSpecies());
    long n_ratios = mixture_ratios.size();

    Eigen::ArrayXXd mole_frac_matrix = Eigen::ArrayXXd::Zero(n_ratios, n_species);
    std::vector<double> mole_fracs(n_species);

    for (long i = 0; i < n_ratios; i++) {
        set_mixture_composition(mixture_ratios[i], type, m_fuel_composition, m_oxidizer_composition);
        thermo->getMoleFractions(mole_fracs.data());
        for (long j = 0; j < n_species; j++) {
            mole_frac_matrix(i, j) = mole_fracs[j];
        }
    }

    assert((mole_frac_matrix >= 0).all() && "Mole fractions must be nonnegative");
    return mole_frac_matrix;
}

Eigen::ArrayXXd Combustor::generate_mass_fraction_matrix(
    const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type) const {

    auto thermo = m_gas.thermo();
    long n_species = static_cast<long>(thermo->nSpecies());
    long n_ratios = mixture_ratios.size();

    Eigen::ArrayXXd mass_frac_matrix = Eigen::ArrayXXd::Zero(n_ratios, n_species);
    std::vector<double> mass_fracs(n_species);

    for (long i = 0; i < n_ratios; i++) {
        set_mixture_composition(mixture_ratios[i], type, m_fuel_composition, m_oxidizer_composition);
        thermo->getMassFractions(mass_fracs.data());
        for (long j = 0; j < n_species; j++) {
            mass_frac_matrix(i, j) = mass_fracs[j];
        }
    }

    assert((mass_frac_matrix >= 0).all() && "Mass fractions must be nonnegative");
    return mass_frac_matrix;
}

// ---- DilutedCombustor ----

DilutedCombustor::DilutedCombustor(Gas gas, const std::string& fuel_comp,
    const std::string& ox_comp, const std::string& dilute_comp)
    : BaseCombustor(std::move(gas))
{
    m_fuel_composition = Cantera::parseCompString(fuel_comp, m_gas.species_names());
    m_oxidizer_composition = Cantera::parseCompString(ox_comp, m_gas.species_names());
    m_flue_composition = Cantera::parseCompString(dilute_comp, m_gas.species_names());
}

DilutedCombustor::DilutedCombustor(Gas gas, const Composition& fuel,
    const Composition& oxidizer, const Composition& flue)
    : BaseCombustor(std::move(gas)),
      m_fuel_composition(fuel), m_oxidizer_composition(oxidizer), m_flue_composition(flue)
{}

ThermoArray DilutedCombustor::solve(const Eigen::ArrayXd& temperatures, const Eigen::ArrayXd& pressures,
    const Eigen::ArrayXd& mixture_ratios, double recirculation_ratio, const CombustorOptions& options) {

    Eigen::ArrayXXd mole_fracs = generate_mole_fraction_matrix(mixture_ratios, options.mixture_type, recirculation_ratio);

    ThermoArray combustion_states(m_gas.solution(), {temperatures.size(), pressures.size(), mixture_ratios.size()});
    combustion_states.TPX(temperatures, pressures, mole_fracs);

    return combust(combustion_states, options);
}

ThermoArray DilutedCombustor::solve(double fuel_temperature, double oxidizer_temperature,
    double flue_temperature, const Eigen::ArrayXd& pressures,
    const Eigen::ArrayXd& mixture_ratios, double recirculation_ratio,
    const CombustorOptions& options) {

    auto thermo = m_gas.thermo();
    MixtureRatioType type = options.mixture_type;
    double r = recirculation_ratio;
    double w_reactant = 1.0 / (1.0 + r);
    double w_dilution = r / (1.0 + r);

    double ref_pressure = pressures[0];

    thermo->setState_TPX(fuel_temperature, ref_pressure, m_fuel_composition);
    double fuel_enthalpy = thermo->enthalpy_mass();

    thermo->setState_TPX(oxidizer_temperature, ref_pressure, m_oxidizer_composition);
    double oxidizer_enthalpy = thermo->enthalpy_mass();

    thermo->setState_TPX(flue_temperature, ref_pressure, m_flue_composition);
    double flue_enthalpy = thermo->enthalpy_mass();

    Eigen::ArrayXXd mass_fracs = generate_mass_fraction_matrix(mixture_ratios, type, recirculation_ratio);

    long n_compositions = mixture_ratios.size();
    long n_pressures = pressures.size();
    Eigen::ArrayXd enthalpies(n_compositions);

    for (long i = 0; i < n_compositions; i++) {
        // Get fuel mass fraction within the fresh feed
        set_mixture_composition(mixture_ratios[i], type, m_fuel_composition, m_oxidizer_composition);
        double fuel_mass_frac_fresh = thermo->mixtureFraction(
            m_fuel_composition, m_oxidizer_composition, Cantera::ThermoBasis::mass);
        double ox_mass_frac_fresh = 1.0 - fuel_mass_frac_fresh;

        enthalpies[i] = fuel_enthalpy * fuel_mass_frac_fresh * w_reactant
                      + oxidizer_enthalpy * ox_mass_frac_fresh * w_reactant
                      + flue_enthalpy * w_dilution;
    }

    ThermoArray combustion_states(m_gas.solution(), {n_compositions, n_pressures});

    for (long i = 0; i < n_compositions; i++) {
        Eigen::ArrayXd row = mass_fracs.row(i);
        for (long j = 0; j < n_pressures; j++) {
            thermo->setMassFractions(row.data());
            thermo->setState_HP(enthalpies[i], pressures[j]);
            combustion_states.set_state(combustion_states.flat_index(i, j), m_gas.save_state());
        }
    }

    return combust(combustion_states, options);
}

Eigen::ArrayXXd DilutedCombustor::generate_mole_fraction_matrix(
    const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type,
    double recirculation_ratio) const {

    auto thermo = m_gas.thermo();
    long n_species = static_cast<long>(thermo->nSpecies());
    long n_ratios = mixture_ratios.size();
    double r = recirculation_ratio;

    // Get molecular weight and mole fractions for flue gas
    thermo->setState_TPX(300.0, 101325.0, m_flue_composition);
    double mw_flue = thermo->meanMolecularWeight();
    std::vector<double> flue_mole_fracs(n_species);
    thermo->getMoleFractions(flue_mole_fracs.data());

    Eigen::ArrayXXd mole_frac_matrix = Eigen::ArrayXXd::Zero(n_ratios, n_species);
    std::vector<double> fresh_mole_fracs(n_species);

    for (long i = 0; i < n_ratios; i++) {
        // Set fresh feed composition
        set_mixture_composition(mixture_ratios[i], type, m_fuel_composition, m_oxidizer_composition);
        double mw_fresh = thermo->meanMolecularWeight();
        thermo->getMoleFractions(fresh_mole_fracs.data());

        // Convert mass-based recirculation ratio to molar basis
        // n_fresh = m_fresh / MW_fresh, n_flue = m_flue / MW_flue
        // where m_fresh = 1/(1+r), m_flue = r/(1+r)
        double n_fresh = 1.0 / ((1.0 + r) * mw_fresh);
        double n_flue = r / ((1.0 + r) * mw_flue);
        double n_total = n_fresh + n_flue;

        double mole_frac_fresh = n_fresh / n_total;
        double mole_frac_flue = n_flue / n_total;

        for (long j = 0; j < n_species; j++) {
            mole_frac_matrix(i, j) = fresh_mole_fracs[j] * mole_frac_fresh
                                   + flue_mole_fracs[j] * mole_frac_flue;
        }
    }

    assert((mole_frac_matrix >= 0).all() && "Mole fractions must be nonnegative");
    return mole_frac_matrix;
}

Eigen::ArrayXXd DilutedCombustor::generate_mass_fraction_matrix(
    const Eigen::ArrayXd& mixture_ratios, MixtureRatioType type,
    double recirculation_ratio) const {

    auto thermo = m_gas.thermo();
    long n_species = static_cast<long>(thermo->nSpecies());
    long n_ratios = mixture_ratios.size();
    double r = recirculation_ratio;
    double w_fresh = 1.0 / (1.0 + r);
    double w_flue = r / (1.0 + r);

    // Get mass fractions for flue gas
    thermo->setState_TPX(300.0, 101325.0, m_flue_composition);
    std::vector<double> flue_mass_fracs(n_species);
    thermo->getMassFractions(flue_mass_fracs.data());

    Eigen::ArrayXXd mass_frac_matrix = Eigen::ArrayXXd::Zero(n_ratios, n_species);
    std::vector<double> fresh_mass_fracs(n_species);

    for (long i = 0; i < n_ratios; i++) {
        set_mixture_composition(mixture_ratios[i], type, m_fuel_composition, m_oxidizer_composition);
        thermo->getMassFractions(fresh_mass_fracs.data());

        for (long j = 0; j < n_species; j++) {
            mass_frac_matrix(i, j) = fresh_mass_fracs[j] * w_fresh
                                   + flue_mass_fracs[j] * w_flue;
        }
    }

    if (!((mass_frac_matrix >= 0).all())) {
        throw ConvergenceError("Mass fractions must be nonnegative.");
    }

    return mass_frac_matrix;
}

} //namespace Goddard
