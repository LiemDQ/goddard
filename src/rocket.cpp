#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "goddard/error.hpp"
#include "goddard/rocket.hpp"

namespace Goddard {

RocketPerformance rocket_performance(const ThermodynamicState& stagnation,
                                     const ThermodynamicState& throat,
                                     const ThermodynamicState& exit)
{
    RocketPerformance performance;
    performance.pressure_ratio = stagnation.pressure / exit.pressure;
    // c* = P_0 A_t / mdot, with mdot/A_t = rho_t u_t. The mass flux is stationary at the throat,
    // so a throat solved to Mach 1 within a tolerance gives c* to second order in that tolerance;
    // rho_t a_t would carry the Mach error to first order.
    const double throat_velocity = std::sqrt(2.0 * (stagnation.enthalpy - throat.enthalpy));
    const double throat_mass_flux = throat.density * throat_velocity;
    performance.cstar = stagnation.pressure / throat_mass_flux;
    // Energy: u_e = sqrt(2 (h_0 - h_e))
    const double exit_velocity = std::sqrt(2.0 * (stagnation.enthalpy - exit.enthalpy));
    // Mass: rho u is constant times area, so Ae/At = rho_t u_t / (rho_e u_e)
    performance.area_ratio = throat_mass_flux / (exit.density * exit_velocity);
    // Thrust coefficient at optimum expansion
    performance.CF = exit_velocity / performance.cstar;
    performance.isp = exit_velocity;
    performance.ivac = exit_velocity
        + (exit.pressure / stagnation.pressure) * performance.area_ratio * performance.cstar;
    performance.mach_number = exit_velocity / exit.speed_of_sound;
    return performance;
}

namespace {

/** Reject a reactant stream that was left unset or given a nonphysical state. */
void validate_stream(const PhaseSpecification& stream, const std::string& name) {
    if (stream.composition.empty()) {
        throw std::invalid_argument(std::format(
            "RocketSolver: RocketPoint.{} has no composition; set its temperature, pressure and "
            "composition.", name));
    }
    if (!(stream.T > 0.0) || !std::isfinite(stream.T)) {
        throw std::invalid_argument(std::format(
            "RocketSolver: RocketPoint.{} must have a positive temperature, got {:.6g} K.",
            name, stream.T));
    }
    if (!(stream.P > 0.0) || !std::isfinite(stream.P)) {
        throw std::invalid_argument(std::format(
            "RocketSolver: RocketPoint.{} must have a positive pressure, got {:.6g} Pa.",
            name, stream.P));
    }
}

/** Reject a mixture ratio outside its domain, as `Combustor` does. */
void validate_mixture_ratio(double value, MixtureRatioType type) {
    bool valid = false;
    switch (type) {
        case MixtureRatioType::OF_RATIO:
        case MixtureRatioType::PHI_RATIO:
            valid = value >= 0.0 && std::isfinite(value);
            break;
        case MixtureRatioType::FUEL_FRAC:
            valid = value >= 0.0 && value <= 1.0;
            break;
    }
    if (!valid) {
        throw std::invalid_argument(std::format(
            "RocketSolver: mixture ratio {:.6g} is out of range; O/F and equivalence ratios must "
            "be nonnegative and finite, and fuel fractions must lie in [0, 1].", value));
    }
}

/** Status of the first failed station of a result, or ok when every station was solved. */
SolveStatus first_failure(const RocketResult& result) {
    for (const NozzleStation* station : {&result.chamber, &result.stagnation, &result.combustion_end}) {
        if (!station->status.ok()) {
            return station->status;
        }
    }
    if (!result.throat.status.ok()) {
        return result.throat.status;
    }
    for (const NozzleStation& exit : result.exits) {
        if (!exit.status.ok()) {
            return exit.status;
        }
    }
    return {};
}

} // namespace

RocketSolver::RocketSolver(const Gas& products, std::string reactant_file)
    : m_products(products.clone()), m_reactant_file(std::move(reactant_file)) {}

void RocketSolver::check_implemented() const {
    if (include_transport) {
        throw NotImplementedError("RocketSolver: transport properties are not implemented.");
    }
    if (include_ionized_species) {
        throw NotImplementedError("RocketSolver: ionized species are not implemented.");
    }
    if (trace_cutoff != 0.0) {
        throw NotImplementedError("RocketSolver: a trace cutoff is not implemented.");
    }
}

const Gas& RocketSolver::reactant_phase(const Composition& composition) const {
    std::vector<std::string> species;
    species.reserve(composition.size());
    for (const auto& entry : composition) {
        species.push_back(entry.first);
    }
    // `Composition` is an ordered map, so the species list is a canonical key.
    auto found = m_reactant_phases.find(species);
    if (found == m_reactant_phases.end()) {
        found = m_reactant_phases.emplace(species,
            Gas::create_from_species(m_reactant_file, "reactants", species)).first;
    }
    return found->second;
}

void RocketSolver::validate(const RocketPoint& point) const {
    check_implemented();
    validate_stream(point.fuel, "fuel");
    validate_stream(point.oxidizer, "oxidizer");
    if (!(point.chamber_pressure > 0.0) || !std::isfinite(point.chamber_pressure)) {
        throw std::invalid_argument(std::format(
            "RocketSolver: RocketPoint.chamber_pressure must be positive and finite, got {:.6g} Pa.",
            point.chamber_pressure));
    }
    validate_mixture_ratio(point.mixture_ratio, point.combustor.mixture_type);
    BaseCombustor::validate_options(point.combustor);
    Nozzle::validate(point.nozzle, point.expansion_type, point.expansion_ratios);

    if (point.combustor.type != CombustorType::INFINITE_AREA
        && point.nozzle.chemistry == GasChemistry::FROZEN && point.nozzle.frozen_NFZ == 0) {
        throw NotImplementedError(
            "RocketSolver: a finite-area chamber frozen at the combustion end (frozen_NFZ = 0) is "
            "not implemented.");
    }
    if (point.combustor.process == CombustionProcess::ISOCHORIC
        && m_products.has_condensed_candidates()) {
        throw NotImplementedError(
            "RocketSolver: ISOCHORIC combustion with candidate condensed species is not "
            "implemented.");
    }

    if (m_reactant_file.empty()) {
        // The compositions are product species.
        const auto& thermo = *m_products.thermo();
        for (const PhaseSpecification* stream : {&point.fuel, &point.oxidizer}) {
            for (const auto& entry : stream->composition) {
                if (thermo.speciesIndex(entry.first, false) == Cantera::npos) {
                    throw std::invalid_argument(std::format(
                        "RocketSolver: reactant species '{}' is not a product species, and the "
                        "solver has no reactant file.", entry.first));
                }
            }
        }
    } else {
        if (point.combustor.mixture_type == MixtureRatioType::PHI_RATIO) {
            throw NotImplementedError(
                "RocketSolver: equivalence ratios are not implemented with a reactant file; use "
                "OF_RATIO or FUEL_FRAC.");
        }
        // Builds (and keeps) the reactant phases, which rejects unknown species.
        reactant_phase(point.fuel.composition);
        reactant_phase(point.oxidizer.composition);
    }
}

Gas RocketSolver::reactant_gas(const PhaseSpecification& stream) const {
    Gas gas = reactant_phase(stream.composition).clone();
    gas.set_state_TPX(stream.T, stream.P, stream.composition);
    return gas;
}

ThermoArray RocketSolver::combust(const RocketPoint& point) {
    const Eigen::ArrayXd pressures = Eigen::ArrayXd::Constant(1, point.chamber_pressure);
    const Eigen::ArrayXd ratios = Eigen::ArrayXd::Constant(1, point.mixture_ratio);

    // Consecutive points usually share their reactants, so the last combustor is kept.
    if (!m_combustor || point.fuel != m_combustor_fuel || point.oxidizer != m_combustor_oxidizer) {
        m_combustor.reset();
        if (m_reactant_file.empty()) {
            m_combustor.emplace(m_products, point.fuel.composition, point.oxidizer.composition);
        } else {
            m_combustor.emplace(m_products, reactant_gas(point.fuel), reactant_gas(point.oxidizer));
        }
        m_combustor_fuel = point.fuel;
        m_combustor_oxidizer = point.oxidizer;
    }
    if (m_reactant_file.empty()) {
        return m_combustor->solve(point.fuel.T, point.oxidizer.T, pressures, ratios, point.combustor);
    }
    return m_combustor->solve(pressures, ratios, point.combustor);
}

Nozzle& RocketSolver::nozzle_for(const NozzleOptions& options) {
    if (!m_nozzle || m_nozzle_options != options) {
        m_nozzle.reset();
        m_nozzle.emplace(m_products, options);
        m_nozzle_options = options;
    }
    return *m_nozzle;
}

RocketResult RocketSolver::solve_validated(const RocketPoint& point) {
    RocketResult result;
    result.point = point;
    const std::size_t num_exits = point.expansion_ratios.size();

    const ThermoArray combustion = combust(point);
    const SolveStatus& combustion_status = combustion.status_at(0);
    if (!combustion_status.ok()) {
        // Nothing to expand: every station is downstream of the failed combustion.
        result.chamber = NozzleStation::failed(combustion_status);
        const SolveStatus downstream = SolveStatus::not_computed(combustion_status);
        result.stagnation = NozzleStation::failed(downstream);
        result.combustion_end = NozzleStation::failed(downstream);
        result.throat = ThroatCondition::failed(downstream);
        result.exits.assign(num_exits, NozzleStation::failed(downstream));
    } else {
        const std::vector<double> state = combustion.get_state(0);
        Nozzle& nozzle = nozzle_for(point.nozzle);
        switch (point.combustor.type) {
            case CombustorType::INFINITE_AREA: {
                nozzle.set_inlet_state(state);
                NozzleResults expansion = nozzle.solve(point.expansion_type, point.expansion_ratios);
                result.chamber = expansion.inlet;
                result.stagnation = expansion.inlet;
                result.combustion_end = std::move(expansion.inlet);
                result.throat = std::move(expansion.throat);
                result.exits = std::move(expansion.expansions);
                if (result.chamber.status.ok()) {
                    result.contraction_ratio = std::numeric_limits<double>::infinity();
                    result.mass_flux = 0.0;
                }
                break;
            }
            case CombustorType::FINITE_MASS_FLUX:
            case CombustorType::FINITE_CONTRACTION_RATIO: {
                const double value = point.combustor.type == CombustorType::FINITE_CONTRACTION_RATIO
                    ? point.combustor.contraction_ratio : point.combustor.mass_flux;
                FiniteAreaChamber chamber = nozzle.solve_finite_area_chamber(
                    state, point.combustor.type, value);
                result.exits = nozzle.solve_stations(chamber, point.expansion_type,
                                                     point.expansion_ratios);
                result.chamber = std::move(chamber.injector);
                result.stagnation = std::move(chamber.stagnation);
                result.combustion_end = std::move(chamber.combustion_end);
                result.throat = std::move(chamber.throat);
                result.contraction_ratio = chamber.contraction_ratio;
                result.mass_flux = chamber.mass_flux;
                break;
            }
        }
    }

    // Every station is measured from the stagnation enthalpy of the point, which the adiabatic
    // expansion conserves.
    const double stagnation_enthalpy = result.stagnation.status.ok()
        ? result.stagnation.thermo.enthalpy : NaN;
    for (NozzleStation* station : {&result.chamber, &result.stagnation, &result.combustion_end}) {
        station->thermo.stagnation_enthalpy = stagnation_enthalpy;
    }
    result.throat.thermo.stagnation_enthalpy = stagnation_enthalpy;
    for (NozzleStation& exit : result.exits) {
        exit.thermo.stagnation_enthalpy = stagnation_enthalpy;
    }

    result.performance.reserve(num_exits);
    for (const NozzleStation& exit : result.exits) {
        RocketPerformance performance;
        if (!result.stagnation.status.ok()) {
            performance.status = result.stagnation.status;
        } else if (!result.throat.status.ok()) {
            performance.status = result.throat.status;
        } else if (!exit.status.ok()) {
            performance.status = exit.status;
        } else {
            performance = rocket_performance(result.stagnation.thermo, result.throat.thermo,
                                             exit.thermo);
        }
        result.performance.push_back(std::move(performance));
    }

    result.status = first_failure(result);
    return result;
}

RocketResult RocketSolver::solve(const RocketPoint& point) {
    validate(point);
    return solve_validated(point);
}

std::vector<RocketResult> RocketSolver::solve(const std::vector<RocketPoint>& points) {
    // Every point is checked before any is solved.
    for (std::size_t i = 0; i < points.size(); i++) {
        try {
            validate(points[i]);
        } catch (const std::invalid_argument& error) {
            throw std::invalid_argument(std::format("point {}: {}", i, error.what()));
        }
    }
    std::vector<RocketResult> results;
    results.reserve(points.size());
    for (const RocketPoint& point : points) {
        results.push_back(solve_validated(point));
    }
    return results;
}

} // namespace Goddard
