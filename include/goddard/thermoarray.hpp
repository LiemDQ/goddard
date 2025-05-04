#pragma once 

#include "cantera/thermo.h"
#include "cantera/base/SolutionArray.h"
#include "eigen3/Eigen/Dense"
#include <memory>
#include <string>
#include <vector>

namespace Goddard {

/**
 * @brief Wrapper around `Cantera::SolutionArray` with a higher-level API for 
 * broadcasting thermodynamic operations.
 * Supports ND arrays.
 * 
*/
class ThermoArray {
	public:
	/**
	 * NOTE: `ThermoArray` is a straightforward extension of Cantera' `SolutionArray`. The most efficient
	 * solution, in terms of performance and repeated code, would be to make it a derived class of `SolutionArray`.
	 * Unfortunately, as of 2025-04-19 all of `SolutionArray`'s constructors are private, which makes it impossible to
	 * inherit from. The only way to construct a `SolutionArray` class is to call the `create` method.
	 * 
	 * The alternative for now is to implement `ThermoArray` as a wrapper class that contains a pointer to a `SolutionArray`.
	 */
	ThermoArray(std::shared_ptr<Cantera::Solution> sol, int len, const Cantera::AnyMap& meta={});
	ThermoArray(std::shared_ptr<Cantera::Solution> sol, const std::vector<long>& shape);
	
	static std::shared_ptr<ThermoArray> create(const std::shared_ptr<Cantera::Solution>& sol, int size=0, const Cantera::AnyMap& meta={}) {
		return std::shared_ptr<ThermoArray>(new ThermoArray(sol, size, meta));
	}

	void reshape(const std::vector<long>& shape);
	inline std::vector<long> shape() const {return m_states->apiShape();}
	inline int size() const {return m_states->size();}
	inline int ndim() const {return m_states->apiNdim();}
	inline bool is_shape_set() const {return m_shape_is_set;}

	/**
	 * @brief Get a pointer to the underlying `SolutionArray` object. 
	 */
	inline std::shared_ptr<Cantera::SolutionArray> solutionarray() {return m_states;}
	
	/**
	 * @brief Get a pointer to the underlying `Solution` object.
	 */
	inline std::shared_ptr<Cantera::Solution> solution() {return m_solution;}
	
	
	void equilibrate(const std::string& XY, 
		const std::string& solver="auto", 
		double rtol=1e-6,
		int max_steps=50000,
		int max_iter=100,
		int estimate_equil=0,
		int log_level=0);

	/** Thermodynamic state methods
		TODO: Possibly enable dynamic resizing of the array based on the arguments
		provided to the method
	*/

	/**
	 * @brief Set the temperature and density of the array.
	 */
	void TD(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ds);

	/**
	 * @brief Set the enthalpy and pressure of the array.
	 */
	void TV(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Vs);

	/**
	 * @brief Set the temperature and pressure of the array.
	 */
	void TP(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps);
	/**
	 * @brief Set the temperature, pressure, and mole fractions of the array. 
	 * The mole fraction matrix columns should represent species and 
	 * the rows should represent distinct compositions. 
	 * */	
	void TPX(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs);
	void TPY(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys);

	/**
	 * @brief Set the enthalpy and pressure of the array.
	 */
	void HP(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps);

	/**
	 * @brief Set the entropy and pressure of the array.
	 */
	void SP(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps);
	
	/**
	 * @brief Set the entropy and the enthalpy of the array.
	 */
	void SH(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Hs);

	/**
	 * Set the internal energy and the volume of the array.
	 */
	void UV(const Eigen::ArrayXd& Us, const Eigen::ArrayXd& Vs);
	

	private:
	void check_dimensionality(size_t len, size_t dim);
	
	void update_states(void (Cantera::ThermoPhase::*f)(double, double), 
		const Eigen::ArrayXd& var1,
		const Eigen::ArrayXd& var2);

	void update_states(void (Cantera::ThermoPhase::*f)(double, double, double), 
		const Eigen::ArrayXd& var1,
		const Eigen::ArrayXd& var2,
		double tol = 1e-9);

	void update_states_with_composition(void (Cantera::ThermoPhase::*f)(double, double, double), 
		const Eigen::ArrayXd& var1,
		const Eigen::ArrayXd& var2, 
		const Eigen::ArrayXXd& var3,
		double tol = 1e-9);

	void update_states_with_composition(void (Cantera::ThermoPhase::*f)(double, double, const double*),
		const Eigen::ArrayXd& var1,
		const Eigen::ArrayXd& var2,
		const Eigen::ArrayXXd& var3);
	
		
	/**
	 * @brief Updates thermodynamic state by broadcasting a function `f` with values in `var1` and `var2`. 
	 * 
	 * This template function is defined in the source file because it is a private method that is only used 
	 * within the same source file.
	 */
	template <typename Func>
	void _update_states(Func&& f, const Eigen::ArrayXd& arr1, const Eigen::ArrayXd& arr2);

	template <typename Func>
	void _update_states_with_composition(Func&& f, const Eigen::ArrayXd& arr1, const Eigen::ArrayXd& arr2, const Eigen::ArrayXXd& arr3);
	
	/**
	 * @brief Create a copy of the underlying solution that the array is derived from.
	 * 
	 */
	std::shared_ptr<Cantera::Solution> copy_original_solution();
	

	std::shared_ptr<Cantera::Solution> m_solution; //Underlying solution that SolutionArray is derived from.
	std::shared_ptr<Cantera::SolutionArray> m_states; 
	std::vector<double> m_orig_solution_state;
	bool m_shape_is_set = false;
};

}