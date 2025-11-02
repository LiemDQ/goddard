#include "cantera/core.h"
#include "cantera/thermo.h"
#include "cantera/base/SolutionArray.h"

#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"

#include <vector>
#include <exception>
#include <utility>
#include <memory>
#include <string>
#include <cassert>
namespace Goddard {

using Cantera::SolutionArray;
using Cantera::Solution;
using Cantera::ThermoPhase;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;


ThermoArray::ThermoArray(std::shared_ptr<Solution> sol, int len, const Cantera::AnyMap& meta) : 
	m_solution(sol), 
	m_states(SolutionArray::create(sol, len, meta)),
	m_orig_solution_state(m_solution->thermo()->stateSize()) {

	m_solution->thermo()->saveState(m_orig_solution_state);
}

ThermoArray::ThermoArray(std::shared_ptr<Solution> sol, const std::vector<long>& shape) : 
	m_solution(sol), 
	m_states(SolutionArray::create(std::move(sol), static_cast<int>(shape.size()), {})),
	m_orig_solution_state(m_solution->thermo()->stateSize()) {
	
	m_solution->thermo()->saveState(m_orig_solution_state);
	m_states->setApiShape(shape);
	m_shape_is_set = true;
}

void ThermoArray::reshape(const std::vector<long>& shape) {
	long array_size = 1;
	for (long val: shape) {
		array_size *= val;
	}
	m_states->resize(static_cast<int>(array_size));
	m_states->setApiShape(shape);
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
		m_states->thermo()->restoreState(m_states->getState(loc));
		m_states->thermo()->equilibrate(XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level);
		m_states->updateState(loc);
	}

	m_states->thermo()->restoreState(m_orig_solution_state);
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
	
void ThermoArray::check_dimensionality(size_t len, size_t dim){
	const auto& shape = m_states->apiShape();
	if (len != static_cast<size_t>(shape[dim])){
		throw std::length_error("Provided vector for dimension " + std::to_string(dim) 
			+ " was of length " + std::to_string(len) + " but expected length " + std::to_string(shape[dim]));
	}
}

ArrayXXd ThermoArray::retrieve_thermo_data(double (Cantera::ThermoPhase::*f)(void) const, int slice) const {
	std::vector<double> old_state(m_solution->thermo()->stateSize());
	m_solution->thermo()->saveState(old_state);

	auto fn = std::mem_fn(f);
	
	long data_size = 0;
	int initial_index = 0;

	//determine the slice of data to extract, if the array is 3D.
	if (ndim() > 2) {
		const auto& data_shape = shape();
		data_size = data_shape[0]*data_shape[1];
		initial_index = static_cast<int>(data_size * slice);
	} else {
		data_size = size();
	}

	std::vector<double> retrieved_data(data_size); //TODO: use an XXd array directly?
	for (int loc = initial_index; loc <= data_size+initial_index; loc++) {
		m_solution->thermo()->restoreState(m_states->getState(loc));
		retrieved_data.push_back(fn(m_solution->thermo()));
	}
	return reshape_thermo_data(retrieved_data);
}

ArrayXXd ThermoArray::reshape_thermo_data(const std::vector<double>& vec) const {
	if (!m_shape_is_set) {
		throw std::runtime_error("Attempted to retrieve data from ThermoArray before setting its shape.");
	}
	const auto& data_shape = shape();
	
	long rows = data_shape[0];
	long cols = 0;

	if (ndim() == 1) {
		cols = 1;
	} else if (ndim() >= 2) {
		cols = data_shape[1];
	}

	assert(cols*rows == vec.size() && "Data vector and Eigen matrix sizes do not match.");

	ArrayXXd data_array = Eigen::Map<const ArrayXXd>(vec.data(), rows, cols);
	
	return data_array;
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
		check_dimensionality(len1, 0);
		check_dimensionality(len2, 1);
		check_dimensionality(len3, 2);
	}

	int loc = 0;
	for (auto xs: arr3.rowwise()){
		for (double v2: arr2){
			for (double v1: arr1){
				f(v1, v2, xs.data());
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