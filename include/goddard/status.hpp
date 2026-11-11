#pragma once
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace Goddard {

/** Value of a result field that was not computed. */
inline constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

/**
 * Why a solve did not produce a result.
 *
 * The reason says what kind of failure it was, so that a script can decide what to do with the
 * point; `SolveStatus::message` says what happened in detail.
 */
enum class FailureReason {
    /// The solve succeeded.
    NONE,
    /**
     * An iteration ran out of steps, stalled or produced a non-finite step. A better starting
     * point, a looser tolerance or a higher iteration limit may help.
     */
    NOT_CONVERGED,
    /**
     * The inputs admit no solution in this model: a detached oblique shock, a wave slower than
     * the Chapman-Jouguet speed, condensed products in a gas-only solver. Retrying does not help.
     */
    NO_SOLUTION,
    /**
     * A state left the range of the thermodynamic data or of a lookup table, e.g. a temperature
     * outside the fit range of a condensed species.
     */
    OUT_OF_DATA_RANGE,
    /**
     * The thermodynamic backend could not set a state or solve a chemical equilibrium, e.g. no
     * temperature matches the requested enthalpy.
     */
    THERMO_FAILURE,
    /**
     * The solve stopped early. The part of the result computed before the stop is valid; the
     * rest is missing.
     */
    INCOMPLETE,
    /// The solve was not attempted because a stage it depends on failed.
    NOT_COMPUTED
};

/** Name of a `FailureReason`, e.g. "NOT_CONVERGED". */
std::string_view to_string(FailureReason reason);

/**
 * Outcome of a solve. Every result a solver returns carries one, as its `status`.
 *
 * A solver reports a numerical failure here and does not throw, so a sweep keeps its other
 * points. When `ok()` is false, the outputs that depend on the failed solve are NaN (or empty).
 * Invalid arguments and unimplemented options are not solve failures: they throw.
 */
struct SolveStatus {
    /// Kind of failure; `NONE` when the solve succeeded.
    FailureReason reason = FailureReason::NONE;
    /// What went wrong, with the values involved. Empty when the solve succeeded.
    std::string message;
    /// Part of the solve that failed, e.g. "throat" or "exit station 2". Empty when it succeeded.
    std::string stage;
    /// Iterations performed by the failed iteration [-]; -1 if not reported.
    int iterations = -1;
    /// Convergence tolerance the failed iteration applied; NaN if not reported.
    double tolerance = NaN;
    /// Last residual or step of the failed iteration; NaN if not reported.
    double residual = NaN;
    /// Conditions that did not stop the solve but may affect the result.
    std::vector<std::string> warnings;

    /** True when the solve succeeded. */
    bool ok() const { return reason == FailureReason::NONE; }

    /**
     * One-line description: stage, message, reason and iteration data. "ok" when the solve
     * succeeded.
     */
    std::string describe() const;

    /**
     * Throw the failure as an exception, for code that should stop at the first failed solve.
     * Does nothing when the solve succeeded.
     *
     * @throws ConvergenceError if `reason` is `NOT_CONVERGED`.
     * @throws SolveError for any other failure.
     */
    void raise_if_failed() const;

    /** Failed status with the given reason and message. */
    static SolveStatus failure(FailureReason reason, std::string message);

    /**
     * Status of a solve that was skipped because the stage `upstream` depends on failed. The
     * message names the upstream stage and repeats its message. If `upstream` is itself
     * `NOT_COMPUTED`, it already names the failure that caused it and is returned unchanged.
     */
    static SolveStatus not_computed(const SolveStatus& upstream);
};

} // namespace Goddard
