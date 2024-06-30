#pragma once 

#include "cantera/thermo.h"
#include "cantera/base/SolutionArray.h"

#include <memory>
#include <string>
#include <vector>

namespace Goddard {

/**
 * Wrapper around `Cantera::SolutionArray` with a higher-level API. 
 * Supports ND arrays.
*/
class ThermoStateManager {
	public:
	ThermoStateManager(const std::shared_ptr<Cantera::Solution> sol, int len);
	ThermoStateManager(const std::shared_ptr<Cantera::Solution> sol, const std::vector<long>& shape);

	void equilibrate(const std::string& XY, 
		const std::string& solver="auto", 
		double rtol=1e-9,
		int max_steps=50000,
		int max_iter=100,
		int estimate_equil=0,
		int log_level=0);

	inline int size() const {return states->size();}
	inline bool is_shape_set() const {return shape_is_set;}

	void TP(const std::vector<double>& Ts, const std::vector<double>& Ps);	
	void TPX(const std::vector<double>& Ts, const std::vector<double>& Ps, const std::vector<std::vector<double>>& xs);
	void HP(const std::vector<double>& Hs, const std::vector<double>& Ps);
	void SP(const std::vector<double>& Ss, const std::vector<double>& Ps);
	void SPX(const std::vector<double>& Ss, const std::vector<double>& Ps, const std::vector<std::vector<double>>& xs);
	

	private:
	void check_dimensionality(size_t len, size_t dim);
	
	void update_states(void (Cantera::ThermoPhase::*f)(double, double), 
		const std::vector<double>& var1,
		const std::vector<double>& var2);

	void update_states(void (Cantera::ThermoPhase::*f)(double, double, double), 
		const std::vector<double>& var1,
		const std::vector<double>& var2,
		double tol = 1e-9);

	void update_states_with_composition(void (Cantera::ThermoPhase::*f)(double, double, double), 
		const std::vector<double>& var1,
		const std::vector<double>& var2, 
		const std::vector<std::vector<double>>& var3,
		double tol = 1e-9);

	void update_states_with_composition(void (Cantera::ThermoPhase::*f)(double, double, const double*),
		const std::vector<double>& var1,
		const std::vector<double>& var2,
		const std::vector<std::vector<double>>& var3);

	template <typename Func>
	void _update_states(Func&& f, const std::vector<double>& var1, const std::vector<double>& var2);

	template <typename Func>
	void _update_states(Func&& f, const std::vector<double>& var1, const std::vector<double>& var2, const std::vector<std::vector<double>>& var3);
	
	std::shared_ptr<Cantera::Solution> copy_original_state();
	
	std::shared_ptr<Cantera::Solution> solution;
	std::shared_ptr<Cantera::SolutionArray> states; //TODO: maybe make these public
	std::vector<double> orig_solution_state;
	bool shape_is_set = false;
};

}