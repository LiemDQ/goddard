#pragma once 

#include "cantera/core.h"
#include "cantera/base/SolutionArray.h"
#include "cantera/base/AnyMap.h"
#include "eigen3/Eigen/Dense"
#include <memory>
#include <string>
#include <vector>

// NOLINTBEGIN(readability-identifier-naming)
namespace Goddard {

class Gas;
class CondensedPhaseSet;

/**
 * @brief Wrapper around `Cantera::SolutionArray` with a higher-level API for 
 * broadcasting thermodynamic operations.
 * Supports up to 3D arrays.
 * 
 * If the third dimension is used, it is assumed to be used to set a composition. 
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
	 *
	 * Entries are stored with the first dimension varying fastest: the entry at indices (i, j, k) of an array with
	 * shape (n0, n1, n2) is at flat location `i + j*n0 + k*n0*n1` (see `flat_index`).
	 *
	 * The array holds its own copy of the `Solution` passed to the constructor, so changes to that `Solution` made
	 * elsewhere (e.g. through a `Gas` or `Nozzle` sharing it) do not affect the stored states. Copies of a
	 * `ThermoArray` share the same storage.
	 */
	ThermoArray(const std::shared_ptr<Cantera::Solution>& sol, int len, const Cantera::AnyMap& meta={});
	/**
	 * Create an array with the given shape. Every entry is initialized to the current state of `sol`.
	 * An empty `shape` leaves the shape unset; the first setter call then sets it.
	 */
	ThermoArray(const std::shared_ptr<Cantera::Solution>& sol, const std::vector<long>& shape);

	/**
	 * Create an array from a `Gas`, carrying over its candidate condensed species.
	 *
	 * The array stores its own clone of the gas's condensed species set, so the candidate list and
	 * the phase objects are independent of the `Gas` it was built from. Every entry additionally
	 * carries its own amount of each candidate, initialized to the amounts of `gas`.
	 */
	ThermoArray(const Gas& gas, const std::vector<long>& shape);

	static std::shared_ptr<ThermoArray> create(const std::shared_ptr<Cantera::Solution>& sol, int size=0, const Cantera::AnyMap& meta={}) {
		return std::shared_ptr<ThermoArray>(new ThermoArray(sol, size, meta));
	}

	void reshape(const std::vector<long>& shape);
	inline std::vector<long> shape() const {return m_states->apiShape();}

	//Size of ThermoArray (number of entries)
	inline int size() const {return m_states->size();}
	inline int ndim() const {return m_states->apiNdim();}
	inline bool is_shape_set() const {return m_shape_is_set;}

	/**
	 * Flat storage location of the entry at indices (i, j, k). Indices beyond the array's number of
	 * dimensions must be zero.
	 */
	int flat_index(long i, long j = 0, long k = 0) const;

	/** Number of candidate condensed species carried by the array. Zero for a gas-only array. */
	size_t num_condensed() const;

	/** Names of the candidate condensed species, in candidate order. */
	std::vector<std::string> condensed_species_names() const;

	/**
	 * State vector of the entry at flat location `loc`: the Cantera state of the gas phase followed
	 * by one entry per candidate condensed species (kmol per kg of mixture), exactly as
	 * `Gas::save_state()` lays it out.
	 */
	std::vector<double> get_state(int loc) const;

	/**
	 * Set the entry at flat location `loc` from a state vector of the same phase.
	 *
	 * Accepts either the extended length returned by `get_state()` or the bare Cantera state
	 * length, in which case the condensed amounts of that entry are set to zero.
	 *
	 * @throws std::invalid_argument if the vector has neither length.
	 */
	void set_state(int loc, const std::vector<double>& state);

	/**
	 * Amounts of the candidate condensed species at flat location `loc` [kmol per kg of mixture],
	 * in candidate order. Empty for a gas-only array.
	 */
	std::vector<double> get_condensed_moles(int loc) const;

	/**
	 * Set the amounts of the candidate condensed species at flat location `loc`
	 * [kmol per kg of mixture].
	 *
	 * @throws std::invalid_argument if `moles` does not have one entry per candidate.
	 */
	void set_condensed_moles(int loc, const std::vector<double>& moles);

	/**
	 * Property getters. The result has shape (n0, 1) for a 1-D array and (n0, n1) for 2-D and 3-D arrays,
	 * where element (i, j) is the entry at `flat_index(i, j, slice)`. `slice` selects the index along the third
	 * dimension and must be 0 for arrays with fewer than 3 dimensions.
	 */
	Eigen::ArrayXXd temperature(int slice = 0) const;
	Eigen::ArrayXXd pressure(int slice = 0) const;
	Eigen::ArrayXXd internal_energy_mass(int slice = 0) const;
	Eigen::ArrayXXd internal_energy_mole(int slice = 0) const;
	Eigen::ArrayXXd enthalpy_mass(int slice = 0) const;
	Eigen::ArrayXXd enthalpy_mole(int slice = 0) const;
	Eigen::ArrayXXd entropy_mass(int slice = 0) const;
	Eigen::ArrayXXd entropy_mole(int slice = 0) const;
	Eigen::ArrayXXd mean_molecular_weight(int slice = 0) const;


	/**
	 * Equilibrate every entry, holding the two properties named by `XY` constant.
	 *
	 * When the array carries condensed species each location is equilibrated through a `Gas` built
	 * on the array's phase and candidate set, so the condensed amounts take part and are written
	 * back; `rtol` and `max_steps` are forwarded to the Gibbs solver and the remaining arguments
	 * are ignored.
	 *
	 * @throws std::invalid_argument if `solver` is "vcs" and the array carries condensed species.
	 */
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
	 * @brief Set the enthalpy, pressure, and mole fractions of the array. 
	 * The mole fraction matrix columns should represent species and 
	 * the rows should represent distinct compositions. 
	 * */
	void HPX(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs);
	void HPY(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys);

	/**
	 * @brief Set the entropy and pressure of the array.
	 */
	void SP(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps);

	/**
	 * @brief Set the enthalpy, pressure, and mole fractions of the array. 
	 * The mole fraction matrix columns should represent species and 
	 * the rows should represent distinct compositions. 
	 * */
	void SPX(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs);
	void SPY(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys);

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
	void check_ndim(int expected_ndim);

	/**
	 * @brief retrieve thermodynamic state values from SolutionArray.
	 */
	Eigen::ArrayXXd retrieve_thermo_data(double (Cantera::ThermoPhase::*f)(void) const, int slice = 0) const;
	
	
	void update_states(void (Cantera::ThermoPhase::*f)(double, double), 
		const Eigen::ArrayXd& var1,
		const Eigen::ArrayXd& var2);

	void update_states(void (Cantera::ThermoPhase::*f)(double, double, double), 
		const Eigen::ArrayXd& var1,
		const Eigen::ArrayXd& var2,
		double tol = 1e-9);

	void update_states_with_mole_composition(void (Cantera::ThermoPhase::*f)(double, double, double), 
		const Eigen::ArrayXd& var1,
		const Eigen::ArrayXd& var2, 
		const Eigen::ArrayXXd& var3,
		double tol = 1e-9);
		
	void update_states_with_mass_composition(void (Cantera::ThermoPhase::*f)(double, double, double), 
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
	

	// Private copy of the constructor's Solution, used by `m_states` as its working state.
	//
	// Invariant: the state of `m_solution` always equals the stored state at the SolutionArray's buffered
	// location. Cantera 3.2's `SolutionArray::getState`/`setLoc` skip restoring the stored state when asked for
	// the buffered location (Cantera issue #2067, fixed after 3.2.0), so they return whatever state the Solution
	// holds. Keeping the Solution private and leaving it equal to the buffered entry after every operation makes
	// those calls correct. Every method that changes the Solution's state must write it back with
	// `updateState` at the location it loaded.
	std::shared_ptr<Cantera::Solution> m_solution;
	std::shared_ptr<Cantera::SolutionArray> m_states;
	// Candidate condensed species carried over from the `Gas` the array was built from; null for a
	// gas-only array. `Cantera::SolutionArray` has no room for them, so the amounts live in a side
	// table of `size() * num_condensed()` entries, location major: the amount of candidate `k` at
	// flat location `loc` is `m_condensed_moles[loc * num_condensed() + k]`.
	std::shared_ptr<CondensedPhaseSet> m_condensed;
	std::vector<double> m_condensed_moles;
	bool m_shape_is_set = false;
};

}
// NOLINTEND(readability-identifier-naming)