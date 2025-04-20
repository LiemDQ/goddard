#include "cantera/core.h"
#include "cantera/thermo.h"
#include "cantera/base/SolutionArray.h"

#include "goddard/thermoarray.hpp"
#include <vector>
#include <exception>
#include <utility>
#include <memory>
#include <string>

namespace Goddard {

using Cantera::SolutionArray;
using Cantera::Solution;
using Cantera::ThermoPhase;


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

void ThermoArray::equilibrate(const std::string& XY, const std::string& solver, double rtol, int max_steps, int max_iter, int estimate_equil, int log_level){

	for (int loc = 0; loc < size(); loc++){
		m_states->thermo()->restoreState(m_states->getState(loc));
		m_states->thermo()->equilibrate(XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level);
		m_states->updateState(loc);
	}

	m_states->thermo()->restoreState(m_orig_solution_state);
}

void ThermoArray::TD(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ds) {
	update_states(&ThermoPhase::setState_TD, Ts, Ds);
}

void ThermoArray::TV(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Vs) {
	update_states(&ThermoPhase::setState_TV, Ts, Vs);
}

void ThermoArray::TP(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps) {
	update_states(&ThermoPhase::setState_TP, Ts, Ps);
}

void ThermoArray::TPX(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs){
	update_states_with_composition(&ThermoPhase::setState_TPX, Ts, Ps, xs);
}

void ThermoArray::TPY(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys){
	update_states_with_composition(&ThermoPhase::setState_TPY, Ts, Ps, ys);
}

void ThermoArray::HP(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps) {
	update_states(&ThermoPhase::setState_HP, Hs, Ps);
}

void ThermoArray::SP(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps) {
	update_states(&ThermoPhase::setState_SP, Ss, Ps);
}

void ThermoArray::SH(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Hs) {
	update_states(&ThermoPhase::setState_SH, Ss, Hs);
}


void ThermoArray::UV(const Eigen::ArrayXd& Us, const Eigen::ArrayXd& Vs) {
	update_states(&ThermoPhase::setState_UV, Us, Vs);
}
	
void ThermoArray::check_dimensionality(size_t len, size_t dim){
	const auto& shape = m_states->apiShape();
	if (len != static_cast<size_t>(shape[dim])){
		throw std::length_error("Provided vector for dimension " + std::to_string(dim) 
			+ " was of length " + std::to_string(len) + " but expected length " + std::to_string(shape[dim]));
	}
}

void ThermoArray::update_states(void (ThermoPhase::*f)(double, double), const Eigen::ArrayXd& var1, const Eigen::ArrayXd& var2){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	_update_states([&](double v1, double v2){fn(m_states->thermo(),v1, v2);}, var1, var2);
}

void ThermoArray::update_states(void (ThermoPhase::*f)(double, double, double), const Eigen::ArrayXd& var1, const Eigen::ArrayXd& var2, double tol){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	_update_states([&](double v1, double v2){fn(m_states->thermo(),v1, v2, tol);}, var1, var2);
}

void ThermoArray::update_states_with_composition(void (ThermoPhase::*f)(double, double, double), const Eigen::ArrayXd& var1, const Eigen::ArrayXd& var2, const Eigen::ArrayXXd& var3, double tol){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	auto update_f = [&](double v1, double v2, const double* v3){
		fn(m_states->thermo(), v1, v2, tol);
		m_states->thermo()->setMoleFractions(v3);
	};
	_update_states_with_composition(update_f, var1, var2, var3);
}

void ThermoArray::update_states_with_composition(void (ThermoPhase::*f)(double, double, const double*),const Eigen::ArrayXd& var1, const Eigen::ArrayXd& var2, const Eigen::ArrayXXd& var3){
	//function pointer signature is needed so compiler can resolve which overloaded function to use
	auto fn = std::mem_fn(f);
	_update_states_with_composition([&](double v1, double v2, const double* v3){fn(m_states->thermo(), v1, v2, v3);}, var1, var2, var3);
}

template <typename Func>
void ThermoArray::_update_states(Func&& f, const Eigen::ArrayXd& arr1, const Eigen::ArrayXd& arr2){

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
	const Eigen::ArrayXd& arr1, 
	const Eigen::ArrayXd& arr2, 
	const Eigen::ArrayXXd& arr3
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
	std::shared_ptr<Solution> new_sln = Cantera::newSolution(m_solution->source(), m_solution->name());
	auto thermo = m_solution->thermo();
	std::vector<double> Xs(thermo->nSpecies());
	thermo->getMoleFractions(Xs.data());
	new_sln->thermo()->setState_TPX(thermo->temperature(), thermo->pressure(), Xs.data());
	//kinetics and transport are probably not needed

	return new_sln;
}

} //namespace Goddard