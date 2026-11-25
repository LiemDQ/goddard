#include "cantera/core.h"

#include "goddard/thermoarray.hpp"
#include "goddard/error.hpp"
#include "goddard/solve_guard.hpp"
#include "goddard/utils.hpp"

#include <algorithm>
#include <format>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Goddard {

using Eigen::ArrayXd;
using Eigen::ArrayXXd;

namespace {

/** Number of entries of `shape`, after checking that it is a valid array shape. */
size_t checked_size(const std::vector<long>& shape) {
	if (shape.empty() || shape.size() > 3) {
		throw std::invalid_argument(std::format(
			"ThermoArray: the shape must have one to three dimensions, got {}.", shape.size()));
	}
	size_t array_size = 1;
	for (long extent : shape) {
		if (extent < 0) {
			throw std::invalid_argument("ThermoArray: shape dimensions must be non-negative.");
		}
		array_size *= static_cast<size_t>(extent);
	}
	return array_size;
}

} // namespace

ThermoArray::ThermoArray(const Gas& gas, const std::vector<long>& shape) :
	m_shape(shape),
	m_gas(gas.clone()),
	m_states(checked_size(shape), m_gas.save_state()),
	m_statuses(m_states.size()) {}

void ThermoArray::reshape(const std::vector<long>& shape) {
	if (checked_size(shape) != m_states.size()) {
		throw std::invalid_argument(std::format(
			"ThermoArray::reshape: the new shape has {} entries, the array {}.",
			checked_size(shape), m_states.size()));
	}
	m_shape = shape;
}

int ThermoArray::flat_index(long i, long j, long k) const {
	const long indices[3] = {i, j, k};
	long loc = 0;
	for (size_t dim = 0; dim < 3; dim++) {
		const long extent = dim < m_shape.size() ? m_shape[dim] : 1;
		if (indices[dim] < 0 || indices[dim] >= extent) {
			throw std::out_of_range(std::format(
				"ThermoArray index {} out of range for dimension {} of extent {}",
				indices[dim], dim, extent));
		}
		loc = loc * extent + indices[dim];
	}
	return static_cast<int>(loc);
}

int ThermoArray::flat_index(const std::vector<long>& indices) const {
	if (indices.size() != m_shape.size()) {
		throw std::invalid_argument(std::format(
			"ThermoArray: {} indices given for an array of {} dimensions.",
			indices.size(), m_shape.size()));
	}
	return flat_index(indices[0], indices.size() > 1 ? indices[1] : 0,
		indices.size() > 2 ? indices[2] : 0);
}

std::vector<std::string> ThermoArray::species_names() const {
	return m_gas.species_names();
}

size_t ThermoArray::num_condensed() const {
	return m_gas.condensed_species_names().size();
}

std::vector<std::string> ThermoArray::condensed_species_names() const {
	return m_gas.condensed_species_names();
}

void ThermoArray::check_location(int loc) const {
	if (loc < 0 || loc >= size()) {
		throw std::out_of_range(std::format(
			"ThermoArray location {} out of range for an array of {} entries.", loc, size()));
	}
}

void ThermoArray::load(int loc) const {
	m_gas.restore_state(m_states[static_cast<size_t>(loc)]);
}

std::string ThermoArray::entry_label(int loc) const {
	// Row-major: peel the indices off from the last dimension.
	std::vector<long> indices(m_shape.size());
	long remainder = loc;
	for (size_t dim = m_shape.size(); dim-- > 0;) {
		const long extent = std::max(m_shape[dim], 1L);
		indices[dim] = remainder % extent;
		remainder /= extent;
	}
	std::string label = "entry (";
	for (size_t dim = 0; dim < indices.size(); dim++) {
		label += (dim == 0 ? "" : ", ") + std::to_string(indices[dim]);
	}
	return label + ")";
}

std::vector<double> ThermoArray::get_state(int loc) const {
	// A failed entry holds the state it had before the failed solve, not a solution.
	status_at(loc).raise_if_failed();
	return m_states[static_cast<size_t>(loc)];
}

void ThermoArray::set_state(int loc, const std::vector<double>& state) {
	check_location(loc);
	// Restoring checks the length and normalizes a bare gas state to the full layout.
	m_gas.restore_state(state);
	m_states[static_cast<size_t>(loc)] = m_gas.save_state();
	set_status(loc, {});
}

std::vector<double> ThermoArray::get_condensed_moles(int loc) const {
	status_at(loc).raise_if_failed();
	load(loc);
	return m_gas.condensed_moles();
}

void ThermoArray::set_condensed_moles(int loc, const std::vector<double>& moles) {
	check_location(loc);
	if (moles.size() != num_condensed()) {
		throw std::invalid_argument(std::format(
			"ThermoArray::set_condensed_moles: expected {} values, got {}.",
			num_condensed(), moles.size()));
	}
	load(loc);
	m_gas.set_condensed_moles(moles);
	m_states[static_cast<size_t>(loc)] = m_gas.save_state();
}

Gas ThermoArray::gas_at(int loc) const {
	status_at(loc).raise_if_failed();
	load(loc);
	return m_gas.clone();
}

ThermodynamicState ThermoArray::snapshot(int loc, CompositionModel model) const {
	status_at(loc).raise_if_failed();
	load(loc);
	return m_gas.snapshot(model);
}

ArrayXd ThermoArray::evaluate(const std::function<double(const Gas&)>& property) const {
	ArrayXd values(size());
	for (int loc = 0; loc < size(); loc++) {
		if (!m_statuses[static_cast<size_t>(loc)].ok()) {
			values(loc) = NaN;
			continue;
		}
		load(loc);
		values(loc) = property(m_gas);
	}
	return values;
}

ArrayXd ThermoArray::temperature() const {
	return evaluate([](const Gas& gas) { return gas.temperature(); });
}

ArrayXd ThermoArray::pressure() const {
	return evaluate([](const Gas& gas) { return gas.pressure(); });
}

ArrayXd ThermoArray::density() const {
	return evaluate([](const Gas& gas) { return gas.density(); });
}

ArrayXd ThermoArray::enthalpy_mass() const {
	return evaluate([](const Gas& gas) { return gas.enthalpy_mass(); });
}

ArrayXd ThermoArray::entropy_mass() const {
	return evaluate([](const Gas& gas) { return gas.entropy_mass(); });
}

ArrayXd ThermoArray::internal_energy_mass() const {
	return evaluate([](const Gas& gas) {
		return gas.enthalpy_mass() - gas.pressure() / gas.density();
	});
}

ArrayXd ThermoArray::cp_mass() const {
	return evaluate([](const Gas& gas) { return gas.cp_mass(); });
}

ArrayXd ThermoArray::cv_mass() const {
	return evaluate([](const Gas& gas) { return gas.cv_mass(); });
}

ArrayXd ThermoArray::molecular_weight() const {
	return evaluate([](const Gas& gas) { return gas.molecular_weight(); });
}

ArrayXd ThermoArray::mixture_molecular_weight() const {
	return evaluate([](const Gas& gas) { return gas.mixture_molecular_weight(); });
}

ArrayXd ThermoArray::gas_mass_fraction() const {
	return evaluate([](const Gas& gas) { return gas.gas_mass_fraction(); });
}

ArrayXd ThermoArray::gamma_s(CompositionModel model) const {
	return evaluate([model](const Gas& gas) { return gas.gamma_s(model); });
}

ArrayXd ThermoArray::speed_of_sound(CompositionModel model) const {
	return evaluate([model](const Gas& gas) { return gas.speed_of_sound(model); });
}

ArrayXXd ThermoArray::composition_rows(size_t columns,
	const std::function<std::vector<double>(const Gas&)>& fractions) const {
	ArrayXXd rows = ArrayXXd::Constant(size(), static_cast<long>(columns), NaN);
	for (int loc = 0; loc < size(); loc++) {
		if (!m_statuses[static_cast<size_t>(loc)].ok()) {
			continue;
		}
		load(loc);
		const std::vector<double> values = fractions(m_gas);
		for (size_t k = 0; k < columns; k++) {
			rows(loc, static_cast<long>(k)) = values[k];
		}
	}
	return rows;
}

ArrayXXd ThermoArray::mole_fractions() const {
	return composition_rows(m_gas.num_species(), [](const Gas& gas) { return gas.mole_fractions(); });
}

ArrayXXd ThermoArray::mass_fractions() const {
	return composition_rows(m_gas.num_species(), [](const Gas& gas) { return gas.mass_fractions(); });
}

ArrayXXd ThermoArray::mixture_mass_fractions() const {
	return composition_rows(m_gas.num_species() + num_condensed(),
		[](const Gas& gas) { return gas.mixture_mass_fractions(); });
}

const SolveStatus& ThermoArray::status_at(int loc) const {
	check_location(loc);
	return m_statuses[static_cast<size_t>(loc)];
}

void ThermoArray::set_status(int loc, SolveStatus status) {
	check_location(loc);
	m_statuses[static_cast<size_t>(loc)] = std::move(status);
}

int ThermoArray::num_failed() const {
	return static_cast<int>(std::count_if(m_statuses.begin(), m_statuses.end(),
		[](const SolveStatus& status) { return !status.ok(); }));
}

SolveStatus ThermoArray::status() const {
	const int failed = num_failed();
	if (failed == 0) {
		return {};
	}
	const auto first = std::find_if(m_statuses.begin(), m_statuses.end(),
		[](const SolveStatus& status) { return !status.ok(); });
	SolveStatus overall = *first;
	if (failed > 1) {
		overall.message += std::format(" ({} of {} entries failed; this is the first.)", failed, size());
	}
	return overall;
}

Eigen::Array<bool, Eigen::Dynamic, 1> ThermoArray::ok() const {
	Eigen::Array<bool, Eigen::Dynamic, 1> flags(size());
	for (int loc = 0; loc < size(); loc++) {
		flags(loc) = m_statuses[static_cast<size_t>(loc)].ok();
	}
	return flags;
}

void ThermoArray::equilibrate_location(int loc, const std::string& XY, const std::string& solver) {
	load(loc);
	SolveStatus status = capture_failure(entry_label(loc), [&] { m_gas.equilibrate(XY, solver); });
	if (status.ok()) {
		m_states[static_cast<size_t>(loc)] = m_gas.save_state();
	}
	// A failed entry keeps the state it had before the solve.
	set_status(loc, std::move(status));
}

void ThermoArray::equilibrate(const std::string& XY, const std::string& solver) {
	validate_equilibrium_request(XY, solver);
	for (int loc = 0; loc < size(); loc++) {
		equilibrate_location(loc, XY, solver);
	}
}

int ThermoArray::retry_failed(const std::string& XY, const std::string& solver) {
	validate_equilibrium_request(XY, solver);
	for (int loc = 0; loc < size(); loc++) {
		if (!m_statuses[static_cast<size_t>(loc)].ok()) {
			equilibrate_location(loc, XY, solver);
		}
	}
	return num_failed();
}

void ThermoArray::check_shape(const std::vector<long>& extents, const std::string& setter) const {
	if (extents == m_shape) {
		return;
	}
	auto format_shape = [](const std::vector<long>& shape) {
		std::string text;
		for (long extent : shape) {
			text += (text.empty() ? "" : ", ") + std::to_string(extent);
		}
		return "(" + text + ")";
	};
	throw std::length_error(std::format("ThermoArray::{}: the arguments need an array of shape {}, "
		"but the array has shape {}.", setter, format_shape(extents), format_shape(m_shape)));
}

void ThermoArray::broadcast(const std::string& setter, const ArrayXd& var1, const ArrayXd& var2,
	const ArrayXXd* compositions, CompositionInput composition,
	const std::function<void(double, double)>& set) {

	const long rows = compositions ? compositions->rows() : 1;
	if (compositions) {
		check_shape({var1.size(), var2.size(), rows}, setter);
		if (compositions->cols() != static_cast<long>(m_gas.num_species())) {
			throw std::invalid_argument(std::format(
				"ThermoArray::{}: the composition matrix has {} columns, the phase {} species.",
				setter, compositions->cols(), m_gas.num_species()));
		}
	} else {
		check_shape({var1.size(), var2.size()}, setter);
	}

	auto thermo = m_gas.thermo();
	for (long i = 0; i < var1.size(); i++) {
		for (long j = 0; j < var2.size(); j++) {
			for (long k = 0; k < rows; k++) {
				const int loc = compositions ? flat_index(i, j, k) : flat_index(i, j);
				load(loc);
				translate_backend_errors(
					[&] { return std::format("ThermoArray::{}: setting {} from the values ({:.6g}, {:.6g})",
						setter, entry_label(loc), var1(i), var2(j)); },
					[&] {
						if (compositions) {
							// ArrayXXd is column-major: copy the row into contiguous storage.
							const ArrayXd row = compositions->row(k);
							if (composition == CompositionInput::MASS) {
								thermo->setMassFractions(row.data());
							} else {
								thermo->setMoleFractions(row.data());
							}
						}
						set(var1(i), var2(j));
					});
				m_states[static_cast<size_t>(loc)] = m_gas.save_state();
				set_status(loc, {});
			}
		}
	}
}

void ThermoArray::TD(const ArrayXd& Ts, const ArrayXd& Ds) {
	auto thermo = m_gas.thermo();
	broadcast("TD", Ts, Ds, nullptr, CompositionInput::NONE,
		[&](double T, double D) { thermo->setState_TD(T, D); });
}

void ThermoArray::TV(const ArrayXd& Ts, const ArrayXd& Vs) {
	auto thermo = m_gas.thermo();
	broadcast("TV", Ts, Vs, nullptr, CompositionInput::NONE,
		[&](double T, double V) { thermo->setState_TV(T, V); });
}

void ThermoArray::TP(const ArrayXd& Ts, const ArrayXd& Ps) {
	auto thermo = m_gas.thermo();
	broadcast("TP", Ts, Ps, nullptr, CompositionInput::NONE,
		[&](double T, double P) { thermo->setState_TP(T, P); });
}

void ThermoArray::TPX(const ArrayXd& Ts, const ArrayXd& Ps, const ArrayXXd& xs) {
	auto thermo = m_gas.thermo();
	broadcast("TPX", Ts, Ps, &xs, CompositionInput::MOLE,
		[&](double T, double P) { thermo->setState_TP(T, P); });
}

void ThermoArray::TPY(const ArrayXd& Ts, const ArrayXd& Ps, const ArrayXXd& ys) {
	auto thermo = m_gas.thermo();
	broadcast("TPY", Ts, Ps, &ys, CompositionInput::MASS,
		[&](double T, double P) { thermo->setState_TP(T, P); });
}

void ThermoArray::HP(const ArrayXd& Hs, const ArrayXd& Ps) {
	auto thermo = m_gas.thermo();
	broadcast("HP", Hs, Ps, nullptr, CompositionInput::NONE,
		[&](double H, double P) { thermo->setState_HP(H, P); });
}

void ThermoArray::HPX(const ArrayXd& Hs, const ArrayXd& Ps, const ArrayXXd& xs) {
	auto thermo = m_gas.thermo();
	broadcast("HPX", Hs, Ps, &xs, CompositionInput::MOLE,
		[&](double H, double P) { thermo->setState_HP(H, P); });
}

void ThermoArray::HPY(const ArrayXd& Hs, const ArrayXd& Ps, const ArrayXXd& ys) {
	auto thermo = m_gas.thermo();
	broadcast("HPY", Hs, Ps, &ys, CompositionInput::MASS,
		[&](double H, double P) { thermo->setState_HP(H, P); });
}

void ThermoArray::SP(const ArrayXd& Ss, const ArrayXd& Ps) {
	auto thermo = m_gas.thermo();
	broadcast("SP", Ss, Ps, nullptr, CompositionInput::NONE,
		[&](double S, double P) { thermo->setState_SP(S, P); });
}

void ThermoArray::SPX(const ArrayXd& Ss, const ArrayXd& Ps, const ArrayXXd& xs) {
	auto thermo = m_gas.thermo();
	broadcast("SPX", Ss, Ps, &xs, CompositionInput::MOLE,
		[&](double S, double P) { thermo->setState_SP(S, P); });
}

void ThermoArray::SPY(const ArrayXd& Ss, const ArrayXd& Ps, const ArrayXXd& ys) {
	auto thermo = m_gas.thermo();
	broadcast("SPY", Ss, Ps, &ys, CompositionInput::MASS,
		[&](double S, double P) { thermo->setState_SP(S, P); });
}

void ThermoArray::SH(const ArrayXd& Ss, const ArrayXd& Hs) {
	auto thermo = m_gas.thermo();
	broadcast("SH", Ss, Hs, nullptr, CompositionInput::NONE,
		[&](double S, double H) { thermo->setState_SH(S, H); });
}

void ThermoArray::UV(const ArrayXd& Us, const ArrayXd& Vs) {
	auto thermo = m_gas.thermo();
	broadcast("UV", Us, Vs, nullptr, CompositionInput::NONE,
		[&](double U, double V) { thermo->setState_UV(U, V); });
}

} // namespace Goddard
