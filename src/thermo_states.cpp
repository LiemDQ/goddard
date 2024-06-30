#include "cantera/core.h"
#include "cantera/thermo.h"
#include "goddard/thermo_states.hpp"
#include <vector>
#include <exception>
#include <utility>
#include <memory>
#include <string>

namespace Goddard {

using Cantera::SolutionArray;
using Cantera::Solution;
using Cantera::ThermoPhase;


ThermoStateManager::ThermoStateManager(std::shared_ptr<Solution> sol, int len) : 
	solution(sol), 
	states(SolutionArray::create(std::move(sol), len, {})),
	orig_solution_state(solution->thermo()->stateSize()) {

	solution->thermo()->saveState(orig_solution_state);
}

ThermoStateManager::ThermoStateManager(std::shared_ptr<Solution> sol, const std::vector<long>& shape) : 
	solution(sol), 
	states(SolutionArray::create(std::move(sol), static_cast<int>(shape.size()), {})),
	orig_solution_state(solution->thermo()->stateSize()) {

	solution->thermo()->saveState(orig_solution_state);
	states->setApiShape(shape);
}

void ThermoStateManager::equilibrate(const std::string& XY, const std::string& solver, double rtol, int max_steps, int max_iter, int estimate_equil, int log_level){
	//not sure if this will work; this is based on what they did in python implementation
	for (size_t loc = 0; loc < this->size(); loc++){
		states->setLoc(loc);
		states->thermo()->equilibrate(XY, solver, rtol, max_steps, max_iter, estimate_equil, log_level);
		states->updateState(loc);
	}
}

void ThermoStateManager::TP(const std::vector<double>& Ts, const std::vector<double>& Ps) {
	this->update_states(&ThermoPhase::setState_TP, Ts, Ps);
}

void ThermoStateManager::TPX(const std::vector<double>& Ts, const std::vector<double>& Ps, const std::vector<std::vector<double>>& xs){
	this->update_states_with_composition(&ThermoPhase::setState_TPX, Ts, Ps, xs);
}

void ThermoStateManager::HP(const std::vector<double>& Hs, const std::vector<double>& Ps) {
	this->update_states(&ThermoPhase::setState_HP, Hs, Ps);
}

void ThermoStateManager::SP(const std::vector<double>& Ss, const std::vector<double>& Ps) {
	this->update_states(&ThermoPhase::setState_SP, Ss, Ps);
}

void ThermoStateManager::SPX(const std::vector<double>& Ss, const std::vector<double>& Ps, const std::vector<std::vector<double>>& xs) {
	this->update_states_with_composition(&ThermoPhase::setState_SP, Ss, Ps, xs);
}
	
void ThermoStateManager::check_dimensionality(size_t len, size_t dim){
	const auto& shape = this->states->apiShape();
	if (len != static_cast<size_t>(shape[dim])){
		throw std::length_error("Provided vector for dimension " + std::to_string(dim) 
			+ " was of length " + std::to_string(len) + " but expected length " + std::to_string(shape[dim]));
	}
}

void ThermoStateManager::update_states(void (ThermoPhase::*f)(double, double), const std::vector<double>& var1, const std::vector<double>& var2){
	auto fn = std::mem_fn(f);
	this->_update_states([&](double v1, double v2){fn(states->thermo(),v1, v2);}, var1, var2);
}

void ThermoStateManager::update_states(void (ThermoPhase::*f)(double, double, double), const std::vector<double>& var1, const std::vector<double>& var2, double tol){
	auto fn = std::mem_fn(f);
	this->_update_states([&](double v1, double v2){fn(states->thermo(),v1, v2, tol);}, var1, var2);
}

void ThermoStateManager::update_states_with_composition(void (ThermoPhase::*f)(double, double, double), const std::vector<double>& var1, const std::vector<double>& var2, const std::vector<std::vector<double>>& var3, double tol){
	auto fn = std::mem_fn(f);
	auto update_f = [&](double v1, double v2, const double* v3){
		fn(states->thermo(), v1, v2, tol);
		states->thermo()->setMoleFractions(v3);
	};
	this->_update_states(update_f, var1, var2, var3);
}

void ThermoStateManager::update_states_with_composition(void (ThermoPhase::*f)(double, double, const double*), const std::vector<double>& var1, const std::vector<double>& var2, const std::vector<std::vector<double>>& var3){
	auto fn = std::mem_fn(f);
	this->_update_states([&](double v1, double v2, const double* v3){fn(states->thermo(), v1, v2, v3);}, var1, var2, var3);
}

template <typename Func>
void ThermoStateManager::_update_states(Func&& f, const std::vector<double>& var1, const std::vector<double>& var2){

	size_t len1 = var1.size();
	size_t len2 = var2.size();
	if (!shape_is_set) {
		states->setApiShape({static_cast<long>(len1), static_cast<long>(len2)});
		shape_is_set = true;
	} else{
		this->check_dimensionality(len1, 0);
		this->check_dimensionality(len2, 1);
	}

	int loc = 0;
	for (size_t idx2 = 0; idx2 < len2; idx2++){
		for (size_t idx1 = 0; idx1 < len1; idx1++){
			states->setLoc(loc);
			f(var1[idx1], var2[idx2]);
			states->updateState(loc);
			loc++;			
		}
	}
}

template <typename Func>
void ThermoStateManager::_update_states(Func&& f, const std::vector<double>& var1, const std::vector<double>& var2, const std::vector<std::vector<double>>& var3) {

	size_t len1 = var1.size();
	size_t len2 = var2.size();
	size_t len3 = var3.size();
	if (!shape_is_set){
		states->setApiShape({static_cast<long>(len1), static_cast<long>(len2), static_cast<long>(len3)});
		shape_is_set = true;
	} else {
		this->check_dimensionality(len1, 0);
		this->check_dimensionality(len2, 1);
		this->check_dimensionality(len3, 2);
	}

	int loc = 0;
	for (size_t idx3 = 0; idx3 < len3;  idx3++){
		for (size_t idx2 = 0; idx2 < len2; idx2++){
			for (size_t idx1 = 0; idx1 < len1; idx1++){
				states->setLoc(loc);
				f(var1[idx1], var2[idx2], var3[idx3].data());
				states->updateState(loc);
				loc++;			
			}
		}
	}
}
 
std::shared_ptr<Solution> ThermoStateManager::copy_original_state(){
	std::shared_ptr<Solution> new_sln = Cantera::newSolution(solution->source(), solution->name());
	auto thermo = solution->thermo();
	std::vector<double> Xs(thermo->nSpecies());
	thermo->getMoleFractions(Xs.data());
	new_sln->thermo()->setState_TPX(thermo->temperature(), thermo->pressure(), Xs.data());
	//kinetics and transport are probably not needed

	return new_sln;
}

} //namespace Goddard