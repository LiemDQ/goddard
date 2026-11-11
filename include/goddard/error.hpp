#pragma once
#include <stdexcept>
#include <string>

#include "goddard/status.hpp"

namespace Goddard {

class NotImplementedError: public std::logic_error {
    using std::logic_error::logic_error;

    public:
    const char * what() const noexcept override;

};

/**
 * A numerical solve failed.
 *
 * Thrown by functions that have no result to carry a `SolveStatus`: the state-changing methods
 * of `Gas` (state setters, `equilibrate_*`) and scalar-valued solves. Solvers that return a
 * result catch it at their boundary and report it as the result's `status`;
 * `SolveStatus::raise_if_failed` turns a failed status back into this exception.
 *
 * Errors of the thermodynamic backend (Cantera) surface as a `SolveError` with reason
 * `THERMO_FAILURE`.
 */
class SolveError: public std::runtime_error {
public:
    /** @param status failed status; `what()` is its `describe()`. */
    explicit SolveError(SolveStatus status);
    SolveError(FailureReason reason, const std::string& message);

    /** The failure: reason, message, stage and iteration data. */
    const SolveStatus& status() const { return m_status; }
    /** Kind of failure. */
    FailureReason reason() const { return m_status.reason; }

private:
    SolveStatus m_status;
};

/**
 * An iterative solve did not converge: a `SolveError` with reason `NOT_CONVERGED`. Carries the
 * iteration count, the tolerance and the last residual (or step) when the solver reports them.
 */
class ConvergenceError: public SolveError {

public:
    /**
     * @param what description of the iteration that failed
     * @param num_iters iterations performed; -1 if not reported
     * @param tol convergence tolerance applied; NaN if not reported
     * @param residual last residual or step; NaN if not reported
     */
    ConvergenceError(const std::string& what, int num_iters = -1,
        double tol = NaN, double residual = NaN);
    /** @param status failed status with reason `NOT_CONVERGED`. */
    explicit ConvergenceError(SolveStatus status);

    /** Iterations performed before giving up; -1 if not reported. */
    int iterations() const { return status().iterations; }
    /** Convergence tolerance the solver applied; NaN if not reported. */
    double tolerance() const { return status().tolerance; }
    /** Last residual or step size; NaN if not reported. */
    double residual() const { return status().residual; }
};

} // namespace Goddard
