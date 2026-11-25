#pragma once

#include "eigen3/Eigen/Dense"
#include <functional>
#include <string>
#include <vector>

#include "goddard/chemistry.hpp"
#include "goddard/gas.hpp"
#include "goddard/status.hpp"
#include "goddard/thermo.hpp"

// NOLINTBEGIN(readability-identifier-naming)
namespace Goddard {

/**
 * @brief Array of thermodynamic states of one `Gas`, with one to three dimensions.
 *
 * Each entry is a `Gas` state vector (`Gas::save_state`), condensed amounts included, and its
 * properties are evaluated through a private clone of the `Gas` the array was built from. Every
 * property is therefore the mixture property `Gas` reports, condensed phases included.
 *
 * Entries are stored in row-major order, last index fastest, as numpy stores an array: the entry
 * at indices (i, j, k) of an array with shape (n0, n1, n2) is at flat location
 * `(i * n1 + j) * n2 + k` (see `flat_index`). The property getters return one value per entry in
 * this order.
 *
 * Every entry carries the status of the last solve that wrote it. A failed entry has no state:
 * its properties are NaN and its state accessors throw its status as `SolveError`.
 *
 * Copying a `ThermoArray` copies its states.
 */
class ThermoArray {
	public:
	/**
	 * Create an array with the given shape, every entry in the current state of `gas`.
	 *
	 * The array evaluates its entries through its own clone of `gas`, so later changes to `gas`
	 * do not affect it.
	 *
	 * @param gas Phase, species set and candidate condensed species of the entries.
	 * @param shape Extent of each dimension: one to three non-negative values.
	 * @throws std::invalid_argument if `shape` has no or more than three dimensions, or a
	 * negative extent.
	 */
	ThermoArray(const Gas& gas, const std::vector<long>& shape);

	/** Extent of each dimension [-]. */
	const std::vector<long>& shape() const {return m_shape;}
	/** Number of entries [-]. */
	int size() const {return static_cast<int>(m_statuses.size());}
	/** Number of dimensions, 1 to 3 [-]. */
	int ndim() const {return static_cast<int>(m_shape.size());}

	/**
	 * Change the shape, keeping the entries in storage order.
	 *
	 * @throws std::invalid_argument if `shape` is not a valid shape (see the constructor) or
	 * does not have `size()` entries.
	 */
	void reshape(const std::vector<long>& shape);

	/**
	 * Flat storage location of the entry at indices (i, j, k). Indices beyond the array's number
	 * of dimensions must be zero.
	 *
	 * @throws std::out_of_range if an index is outside its dimension.
	 */
	int flat_index(long i, long j = 0, long k = 0) const;

	/**
	 * Flat storage location of the entry at `indices`, one index per dimension.
	 *
	 * @throws std::invalid_argument if `indices` does not have `ndim()` entries.
	 * @throws std::out_of_range if an index is outside its dimension.
	 */
	int flat_index(const std::vector<long>& indices) const;

	/** Names of the gas-phase species, in the order of the composition getters. */
	std::vector<std::string> species_names() const;
	/** Number of candidate condensed species carried by the array. Zero for a gas-only array. */
	size_t num_condensed() const;
	/** Names of the candidate condensed species, in candidate order. */
	std::vector<std::string> condensed_species_names() const;

	/**
	 * State vector of the entry at flat location `loc`, laid out as `Gas::save_state()`: the
	 * Cantera state of the gas phase followed by one entry per candidate condensed species
	 * (kmol per kg of mixture).
	 *
	 * @throws SolveError with the entry's status if the entry failed (see `status_at`).
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	std::vector<double> get_state(int loc) const;

	/**
	 * Set the entry at flat location `loc` from a state vector of the array's `Gas`, and mark it
	 * ok. Accepts the lengths `Gas::restore_state` accepts; the bare Cantera state length sets
	 * the condensed amounts to zero.
	 *
	 * @throws std::invalid_argument if the vector has neither length.
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	void set_state(int loc, const std::vector<double>& state);

	/**
	 * Amounts of the candidate condensed species at flat location `loc` [kmol per kg of mixture],
	 * in candidate order. Empty for a gas-only array.
	 *
	 * @throws SolveError with the entry's status if the entry failed.
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	std::vector<double> get_condensed_moles(int loc) const;

	/**
	 * Set the amounts of the candidate condensed species at flat location `loc`
	 * [kmol per kg of mixture].
	 *
	 * @throws std::invalid_argument if `moles` does not have one entry per candidate.
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	void set_condensed_moles(int loc, const std::vector<double>& moles);

	/**
	 * Independent `Gas` in the state of the entry at flat location `loc`.
	 *
	 * @throws SolveError with the entry's status if the entry failed.
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	Gas gas_at(int loc) const;

	/**
	 * `Gas::snapshot` of the entry at flat location `loc`.
	 *
	 * @param loc Flat location.
	 * @param model Composition model of the derived properties.
	 * @throws SolveError with the entry's status if the entry failed.
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	ThermodynamicState snapshot(int loc, CompositionModel model) const;

	/**
	 * Evaluate `property` on the `Gas` of every entry, in storage order. Failed entries are NaN.
	 * The `Gas` passed to `property` is the array's private evaluator in the entry's state; it
	 * must not be kept.
	 */
	Eigen::ArrayXd evaluate(const std::function<double(const Gas&)>& property) const;

	/**
	 * @name Property getters
	 * One value per entry in storage order (see `flat_index`), NaN for a failed entry. Each is
	 * the `Gas` property of the same name: a mixture property, per kg of gas plus condensed
	 * phases.
	 * @{
	 */
	/** Temperature [K]. */
	Eigen::ArrayXd temperature() const;
	/** Pressure [Pa]. */
	Eigen::ArrayXd pressure() const;
	/** Density [kg/m^3]. */
	Eigen::ArrayXd density() const;
	/** Specific enthalpy [J/kg]. */
	Eigen::ArrayXd enthalpy_mass() const;
	/** Specific entropy [J/(kg.K)]. */
	Eigen::ArrayXd entropy_mass() const;
	/** Specific internal energy [J/kg], h - P/rho. */
	Eigen::ArrayXd internal_energy_mass() const;
	/** Frozen constant-pressure specific heat [J/(kg.K)]. */
	Eigen::ArrayXd cp_mass() const;
	/** Frozen constant-volume specific heat [J/(kg.K)]. */
	Eigen::ArrayXd cv_mass() const;
	/** Molecular weight [kg/kmol], CEA's "M" (see `Gas::molecular_weight`). */
	Eigen::ArrayXd molecular_weight() const;
	/** Mixture molecular weight [kg/kmol], CEA's "MW" (see `Gas::mixture_molecular_weight`). */
	Eigen::ArrayXd mixture_molecular_weight() const;
	/** Mass fraction of the gas phase [-]. */
	Eigen::ArrayXd gas_mass_fraction() const;
	/** Isentropic exponent [-] in `model` (see `Gas::gamma_s`). */
	Eigen::ArrayXd gamma_s(CompositionModel model) const;
	/** Speed of sound [m/s] in `model` (see `Gas::speed_of_sound`). */
	Eigen::ArrayXd speed_of_sound(CompositionModel model) const;
	/** @} */

	/**
	 * Mole fractions of the gas-phase species [-]: one row per entry in storage order, one column
	 * per species of `species_names`. A failed entry's row is NaN.
	 */
	Eigen::ArrayXXd mole_fractions() const;
	/** Mass fractions of the gas-phase species [-]; laid out as `mole_fractions`. */
	Eigen::ArrayXXd mass_fractions() const;
	/**
	 * Mass fractions of the whole mixture [-] (see `Gas::mixture_mass_fractions`): one row per
	 * entry, one column per gas species followed by one per candidate condensed species.
	 */
	Eigen::ArrayXXd mixture_mass_fractions() const;

	/**
	 * Outcome of the most recent solve at flat location `loc`: ok for an entry that was set
	 * directly or equilibrated successfully.
	 *
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	const SolveStatus& status_at(int loc) const;

	/**
	 * Outcome over the whole array: ok when every entry is, otherwise the status of the first
	 * failed entry (in storage order), whose message also gives the number of failed entries.
	 */
	SolveStatus status() const;

	/** Number of entries whose solve failed [-]. */
	int num_failed() const;

	/** Whether each entry is ok, one value per entry in storage order. */
	Eigen::Array<bool, Eigen::Dynamic, 1> ok() const;

	/**
	 * Record the outcome of a solve at flat location `loc`. Used by the solvers that fill the
	 * array entry by entry; `set_state` and the bulk setters reset an entry to ok.
	 *
	 * @throws std::out_of_range if `loc` is not a location of the array.
	 */
	void set_status(int loc, SolveStatus status);

	/**
	 * Equilibrate every entry with `Gas::equilibrate`, holding the two properties named by `XY`
	 * constant. The solver settings are the `equilibrium_options` of the `Gas` the array was
	 * built from.
	 *
	 * An entry whose equilibrium solve fails keeps the state it had before the call and records
	 * the failure (see `status_at`); the other entries are unaffected. Nothing is thrown for a
	 * failed solve.
	 *
	 * @throws std::invalid_argument for an `XY` or `solver` that `Gas::equilibrate` rejects.
	 * @throws NotImplementedError for a property pair `Gas::equilibrate` does not support with
	 * the array's condensed species.
	 */
	void equilibrate(const std::string& XY, const std::string& solver = "auto");

	/**
	 * Equilibrate again the entries whose last solve failed, e.g. with another solver. Entries
	 * that are ok are left alone. See `equilibrate` for the arguments.
	 *
	 * @return Number of entries still failed [-].
	 */
	int retry_failed(const std::string& XY, const std::string& solver = "auto");

	/**
	 * @name Bulk setters
	 * The two-variable setters need a 2-D array of shape (len(var1), len(var2)) and set the
	 * entry (i, j) from (var1[i], var2[j]). The setters with a composition need a 3-D array of
	 * shape (len(var1), len(var2), rows of the composition matrix) and set the entry (i, j, k)
	 * from (var1[i], var2[j]) and row k, whose columns are the species of `species_names`.
	 * A setter that takes no composition keeps each entry's composition and condensed amounts.
	 *
	 * Each setter resets the status of every entry it writes to ok. A state the thermodynamic
	 * backend rejects (e.g. a nonpositive temperature) throws `SolveError`: the entries set before
	 * it are written and ok, and the entry that was rejected and those after it keep their
	 * previous state and status.
	 *
	 * @throws std::length_error if the array's shape does not match the arguments.
	 * @{
	 */
	/** Temperature [K] and density [kg/m^3]. */
	void TD(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ds);
	/** Temperature [K] and specific volume [m^3/kg]. */
	void TV(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Vs);
	/** Temperature [K] and pressure [Pa]. */
	void TP(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps);
	/** Temperature [K], pressure [Pa] and mole fractions [-]. */
	void TPX(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs);
	/** Temperature [K], pressure [Pa] and mass fractions [-]. */
	void TPY(const Eigen::ArrayXd& Ts, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys);
	/** Specific enthalpy [J/kg] and pressure [Pa]. */
	void HP(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps);
	/** Specific enthalpy [J/kg], pressure [Pa] and mole fractions [-]. */
	void HPX(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs);
	/** Specific enthalpy [J/kg], pressure [Pa] and mass fractions [-]. */
	void HPY(const Eigen::ArrayXd& Hs, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys);
	/** Specific entropy [J/(kg.K)] and pressure [Pa]. */
	void SP(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps);
	/** Specific entropy [J/(kg.K)], pressure [Pa] and mole fractions [-]. */
	void SPX(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& xs);
	/** Specific entropy [J/(kg.K)], pressure [Pa] and mass fractions [-]. */
	void SPY(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Ps, const Eigen::ArrayXXd& ys);
	/** Specific entropy [J/(kg.K)] and specific enthalpy [J/kg]. */
	void SH(const Eigen::ArrayXd& Ss, const Eigen::ArrayXd& Hs);
	/** Specific internal energy [J/kg] and specific volume [m^3/kg]. */
	void UV(const Eigen::ArrayXd& Us, const Eigen::ArrayXd& Vs);
	/** @} */

	private:
	/** How a bulk setter writes the evaluator's state from one entry's values. */
	enum class CompositionInput { NONE, MOLE, MASS };

	void check_location(int loc) const;
	/** Put the evaluator in the state of the entry at `loc`, whatever its status. */
	void load(int loc) const;
	/** "entry (i, j, k)" for flat location `loc`, with as many indices as dimensions. */
	std::string entry_label(int loc) const;
	/** Equilibrate the entry at `loc`; on failure keep its state and record the failure. */
	void equilibrate_location(int loc, const std::string& XY, const std::string& solver);
	/** Check that the array has `extents` as its shape. */
	void check_shape(const std::vector<long>& extents, const std::string& setter) const;
	/**
	 * Broadcast `set(v1, v2)` over a 2-D array, or over a 3-D array with the rows of
	 * `compositions` set first as `composition` says.
	 */
	void broadcast(const std::string& setter, const Eigen::ArrayXd& var1, const Eigen::ArrayXd& var2,
		const Eigen::ArrayXXd* compositions, CompositionInput composition,
		const std::function<void(double, double)>& set);
	Eigen::ArrayXXd composition_rows(size_t columns,
		const std::function<std::vector<double>(const Gas&)>& fractions) const;

	std::vector<long> m_shape;
	// Evaluator: a private clone of the constructor's `Gas`. Every operation loads the state of
	// the entry it reads or writes first, so the evaluator's own state carries no meaning.
	mutable Gas m_gas;
	// State vector of each entry (`Gas::save_state` layout), in storage order.
	std::vector<std::vector<double>> m_states;
	// Outcome of the last solve at each entry.
	std::vector<SolveStatus> m_statuses;
};

}
// NOLINTEND(readability-identifier-naming)
