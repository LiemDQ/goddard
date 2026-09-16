#include "cantera/core.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/SolutionArray.h"

#include "goddard/thermoarray.hpp"
#include "goddard/condensed.hpp"
#include "goddard/error.hpp"
#include "goddard/gas.hpp"
#include "goddard/utils.hpp"

#include <vector>
#include <exception>
#include <utility>
#include <memory>
#include <string>
#include <cassert>
#include <iostream>
#include <stdexcept>
namespace Goddard {

using Cantera::SolutionArray;
using Cantera::Solution;
using Cantera::ThermoPhase;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;


namespace {

// Copy of `sol` holding only its thermodynamic model, in the same state as `sol`.
std::shared_ptr<Solution> private_copy(Solution& sol) {
	std::shared_ptr<Solution> copy = sol.clone({}, false, false);
	std::vector<double> state(sol.thermo()->stateSize());
	sol.thermo()->saveState(state);
	copy->thermo()->restoreState(state);
	return copy;
}

long shape_size(const std::vector<long>& shape) {
	long array_size = 1;
	for (long dim : shape) {
		if (dim < 0) {
			throw std::invalid_argument("ThermoArray shape dimensions must be non-negative.");
		}
		array_size *= dim;
	}
	return array_size;
}

}

// SolutionArray::create fills every entry with the current Solution state, so each entry starts valid.
ThermoArray::ThermoArray(const std::shared_ptr<Solution>& sol, int len, const Cantera::AnyMap& meta) : 
	m_solution(private_copy(*sol)), 
	m_states(SolutionArray::create(m_solution, len, meta)) {}

ThermoArray::ThermoArray(const std::shared_ptr<Solution>& sol, const std::vector<long>& shape) : 
	m_solution(private_copy(*sol)), 
	m_states(SolutionArray::create(m_solution, static_cast<int>(shape_size(shape)), {})) {
	
	if (!shape.empty()) {
		m_states->setApiShape(shape);
		m_shape_is_set = true;
	}
}

ThermoArray::ThermoArray(const Gas& gas, const std::vector<long>& shape) :
	ThermoArray(gas.solution(), shape) {

	if (gas.m_condensed) {
		m_condensed = gas.m_condensed->clone();
	}
}

size_t ThermoArray::num_condensed() const {
	return m_condensed ? m_condensed->size() : 0;
}

std::vector<std::string> ThermoArray::condensed_species_names() const {
	if (!m_condensed) return {};
	return m_condensed->names();
}

void ThermoArray::reshape(const std::vector<long>& shape) {
	// setApiShape resizes the storage. Entries that remain keep their data, so the buffered
	// location still matches the Solution state if it is in range; new entries are zero-filled.
	shape_size(shape);
	m_states->setApiShape(shape);
	m_shape_is_set = !shape.empty();
}

int ThermoArray::flat_index(long i, long j, long k) const {
	const std::vector<long> data_shape = shape();
	const std::vector<long> indices = {i, j, k};
	long loc = 0;
	long stride = 1;
	for (size_t dim = 0; dim < indices.size(); dim++) {
		long extent = dim < data_shape.size() ? data_shape[dim] : 1;
		if (indices[dim] < 0 || indices[dim] >= extent) {
			throw std::out_of_range("ThermoArray index " + std::to_string(indices[dim])
				+ " out of range for dimension " + std::to_string(dim)
				+ " of extent " + std::to_string(extent));
		}
		loc += indices[dim] * stride;
		stride *= extent;
	}
	return static_cast<int>(loc);
}

std::vector<double> ThermoArray::get_state(int loc) const {
	// TODO(WP-B): return the extended state vector (Cantera state ++ condensed moles).
	if (num_condensed() > 0) {
		throw NotImplementedError(
			"ThermoArray::get_state: arrays carrying condensed species are not implemented yet (WP-B).");
	}
	return m_states->getState(loc);
}

void ThermoArray::set_state(int loc, const std::vector<double>& state) {
	// TODO(WP-B): accept the extended state vector (Cantera state ++ condensed moles).
	if (num_condensed() > 0) {
		throw NotImplementedError(
			"ThermoArray::set_state: arrays carrying condensed species are not implemented yet (WP-B).");
	}
	m_states->setState(loc, state);
}

ArrayXXd ThermoArray::temperature(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::temperature, slice);
}

ArrayXXd ThermoArray::pressure(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::pressure, slice);
}

ArrayXXd ThermoArray::internal_energy_mass(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::intEnergy_mass, slice);
}

ArrayXXd ThermoArray::internal_energy_mole(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::intEnergy_mole, slice);
}

ArrayXXd ThermoArray::enthalpy_mass(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::enthalpy_mass, slice);
}

ArrayXXd ThermoArray::enthalpy_mole(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::enthalpy_mole, slice);
}

ArrayXXd ThermoArray::entropy_mass(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::entropy_mass, slice);
}

ArrayXXd ThermoArray::entropy_mole(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::entropy_mole, slice);
}

ArrayXXd ThermoArray::mean_molecular_weight(int slice) const {
	return retrieve_thermo_data(&ThermoPhase::meanMolecularWeight, slice);
}


void ThermoArray::equilibrate(const std::string& XY, const std::string& solver, double rtol, int max_steps, int max_iter, int estimate_equil, int log_level){

	for (int loc = 0; loc < size(); loc++){
		m_states->setLoc(loc);
		m_solution->thermo()->equilibrate(XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level);
		m_states->updateState(loc);
	}
}

void ThermoArray::TD(const ArrayXd& Ts, const ArrayXd& Ds) {
	update_states(&ThermoPhase::setState_TD, Ts, Ds);
}

void ThermoArray::TV(const ArrayXd& Ts, const ArrayXd& Vs) {
	update_states(&ThermoPhase::setState_TV, Ts, Vs);
}

void ThermoArray::TP(const ArrayXd& Ts, const ArrayXd& Ps) {
	update_states(&ThermoPhase::setState_TP, Ts, Ps);
}

void ThermoArray::TPX(const ArrayXd& Ts, const ArrayXd& Ps, const ArrayXXd& xs){
	update_states_with_composition(&ThermoPhase::setState_TPX, Ts, Ps, xs);
}

void ThermoArray::TPY(const ArrayXd& Ts, const ArrayXd& Ps, const ArrayXXd& ys){
	update_states_with_composition(&ThermoPhase::setState_TPY, Ts, Ps, ys);
}

void ThermoArray::HP(const ArrayXd& Hs, const ArrayXd& Ps) {
	update_states(&ThermoPhase::setState_HP, Hs, Ps);
}

void ThermoArray::HPX(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs) {
	update_states_with_mole_composition(&ThermoPhase::setState_HP, Hs, Ps, xs);
}
void ThermoArray::HPY(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys) {
	update_states_with_mass_composition(&ThermoPhase::setState_HP, Hs, Ps, ys);
}

void ThermoArray::SP(const ArrayXd& Ss, const ArrayXd& Ps) {
	update_states(&ThermoPhase::setState_SP, Ss, Ps);
}

void ThermoArray::SPX(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs) {
	update_states_with_mole_composition(&ThermoPhase::setState_SP, Ss, Ps, xs);
}

void ThermoArray::SPY(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys) {
	update_states_with_mass_composition(&ThermoPhase::setState_SP, Ss, Ps, ys);
}


void ThermoArray::SH(const ArrayXd& Ss, const ArrayXd& Hs) {
	update_states(&ThermoPhase::setState_SH, Ss, Hs);
}


void ThermoArray::UV(const ArrayXd& Us, const ArrayXd& Vs) {
	update_states(&ThermoPhase::setState_UV, Us, Vs);
}
	
void ThermoArray::check_ndim(int expected_ndim) {
	if (ndim() != expected_ndim) {
		throw std::length_error("Operation requires a " + std::to_string(expected_ndim)
			+ "-D ThermoArray but the array has " + std::to_string(ndim()) + " dimensions.");
	}
}

void ThermoArray::check_dimensionality(size_t len, size_t dim){
	const auto& shape = m_states->apiShape();
	if (len != static_cast<size_t>(shape[dim])){
		throw std::length_error("Provided vector for dimension " + std::to_string(dim) 
			+ " was of length " + std::to_string(len) + " but expected length " + std::to_string(shape[dim]));
	}
}

ArrayXXd ThermoArray::retrieve_thermo_data(double (Cantera::ThermoPhase::*f)(void) const, int slice) const {
	if (!m_shape_is_set) {
		throw std::runtime_error("Attempted to retrieve data from ThermoArray before setting its shape.");
	}
	const std::vector<long> data_shape = shape();
	const long rows = data_shape[0];
	const long cols = ndim() >= 2 ? data_shape[1] : 1;
	const long slices = ndim() >= 3 ? data_shape[2] : 1;
	if (slice < 0 || slice >= slices) {
		throw std::out_of_range("ThermoArray slice " + std::to_string(slice)
			+ " out of range for " + std::to_string(slices) + " slices.");
	}

	auto fn = std::mem_fn(f);
	ArrayXXd data(rows, cols);
	for (long j = 0; j < cols; j++) {
		for (long i = 0; i < rows; i++) {
			// Loading an entry leaves the private Solution equal to the buffered entry.
			m_states->setLoc(flat_index(i, j, slice));
			data(i, j) = fn(m_solution->thermo());
		}
	}
	return data;
}

void ThermoArray::update_states(void (ThermoPhase::*f)(double, double), const ArrayXd& var1, const ArrayXd& var2){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	_update_states([&](double v1, double v2){fn(m_states->thermo(),v1, v2);}, var1, var2);
}

void ThermoArray::update_states(void (ThermoPhase::*f)(double, double, double), const ArrayXd& var1, const ArrayXd& var2, double tol){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	_update_states([&](double v1, double v2){fn(m_states->thermo(),v1, v2, tol);}, var1, var2);
}

void ThermoArray::update_states_with_mole_composition(void (ThermoPhase::*f)(double, double, double), const ArrayXd& var1, const ArrayXd& var2, const ArrayXXd& var3, double tol){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	auto update_f = [&](double v1, double v2, const double* v3){
		m_states->thermo()->setMoleFractions(v3);
		fn(m_states->thermo(), v1, v2, tol);
	};
	_update_states_with_composition(update_f, var1, var2, var3);
}

void ThermoArray::update_states_with_mass_composition(void (ThermoPhase::*f)(double, double, double), const ArrayXd& var1, const ArrayXd& var2, const ArrayXXd& var3, double tol) {
	auto fn = std::mem_fn(f);
	auto update_f = [&](double v1, double v2, const double* v3){
		m_states->thermo()->setMassFractions(v3);
		fn(m_states->thermo(), v1, v2, tol);
	};
	_update_states_with_composition(update_f, var1, var2, var3);
}

void ThermoArray::update_states_with_composition(void (ThermoPhase::*f)(double, double, const double*),const ArrayXd& var1, const ArrayXd& var2, const ArrayXXd& var3){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	_update_states_with_composition([&](double v1, double v2, const double* v3){fn(m_states->thermo(), v1, v2, v3);}, var1, var2, var3);
}


template <typename Func>
void ThermoArray::_update_states(Func&& f, const ArrayXd& arr1, const ArrayXd& arr2){

	size_t len1 = arr1.size();
	size_t len2 = arr2.size();
	if (!m_shape_is_set) {
		m_states->setApiShape({static_cast<long>(len1), static_cast<long>(len2)});
		m_shape_is_set = true;
	} else{
		check_ndim(2);
		check_dimensionality(len1, 0);
		check_dimensionality(len2, 1);
	}

	int loc = 0;
	for (double v2: arr2){
		for (double v1: arr1){
			f(v1, v2);
			m_states->updateState(loc);
			loc++;			
		}
	}
}

/**
 * @brief Updates thermodynamic state, including compositions, 
 * by broadcasting across the specified arrays and matrix. 
 */
template <typename Func>
void ThermoArray::_update_states_with_composition(
	Func&& f, 
	const ArrayXd& arr1, 
	const ArrayXd& arr2, 
	const ArrayXXd& arr3
){

	size_t len1 = arr1.size();
	size_t len2 = arr2.size();
	size_t len3 = arr3.rows();
	if (!m_shape_is_set){
		m_states->setApiShape({static_cast<long>(len1), static_cast<long>(len2), static_cast<long>(len3)});
		m_shape_is_set = true;
	} else {
		check_ndim(3);
		check_dimensionality(len1, 0);
		check_dimensionality(len2, 1);
		check_dimensionality(len3, 2);
	}

	int loc = 0;
	for (int i = 0; i < arr3.rows(); i++){
		// Evaluate row into contiguous storage. ArrayXXd is column-major,
		// so a rowwise view's data() pointer is strided and cannot be
		// passed directly to Cantera functions that expect contiguous arrays.
		ArrayXd row = arr3.row(i);
		for (double v2: arr2){
			for (double v1: arr1){
				f(v1, v2, row.data());
				m_states->updateState(loc);
				loc++;
			}
		}
	}
}
 
std::shared_ptr<Solution> ThermoArray::copy_original_solution(){

	return copy_solution(*m_solution);
}

} //namespace Goddard