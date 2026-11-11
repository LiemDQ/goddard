#pragma once
#include <string>
#include <string_view>
#include <utility>

#include "goddard/error.hpp"
#include "goddard/status.hpp"

/**
 * Conversion of numerical failures into `SolveStatus` values at solver boundaries. Internal: not
 * part of the Python API.
 *
 * Inside a solver a failure travels as an exception (`SolveError`, or a Cantera error from a call
 * that does not go through `Gas`). The public solve function runs its body under
 * `capture_failure`, which turns those into the status of the result. Everything else
 * (`std::invalid_argument`, `NotImplementedError`, logic errors) is a usage or programming error
 * and passes through.
 */

namespace Goddard {

/**
 * Status for the exception being handled. Call it only inside a `catch (...)` block.
 *
 * A `SolveError` gives its own status and a Cantera error a `THERMO_FAILURE`; `stage` is recorded
 * unless the status already names one. Any other exception is rethrown.
 */
SolveStatus status_from_current_exception(std::string_view stage);

/**
 * Run `body` and report a numerical failure as a status.
 *
 * @param stage part of the solve `body` computes, recorded in a failed status
 * @param body callable taking no arguments; its return value is ignored
 * @return an ok status if `body` returns, the failure if it throws a `SolveError` or a Cantera
 * error. Other exceptions propagate.
 */
template <typename Body>
SolveStatus capture_failure(std::string_view stage, Body&& body) {
    try {
        std::forward<Body>(body)();
    } catch (...) {
        return status_from_current_exception(stage);
    }
    return {};
}

/** Record `stage` on a failed status that does not name one yet. */
inline void set_stage_if_unset(SolveStatus& status, std::string_view stage) {
    if (!status.ok() && status.stage.empty()) {
        status.stage = std::string(stage);
    }
}

/**
 * Rethrow the exception being handled, turning a Cantera error into a `SolveError` with reason
 * `THERMO_FAILURE` whose message starts with `operation`. Call it only inside a `catch (...)`
 * block. Any other exception is rethrown unchanged.
 */
[[noreturn]] void rethrow_backend_error(std::string_view operation);

/**
 * As `rethrow_backend_error`, for calls whose failure means a bad argument (an unknown species,
 * phase or file): a Cantera error becomes `std::invalid_argument`.
 */
[[noreturn]] void rethrow_backend_usage_error(std::string_view operation);

/**
 * Run a call into the thermodynamic backend, rethrowing a Cantera error as a `SolveError` with
 * reason `THERMO_FAILURE`.
 *
 * @param operation what was being done, with the values involved, e.g.
 * "Gas::set_state_TP(T = 300 K, P = 1e5 Pa)"
 * @return the value `call` returns
 */
template <typename Call>
decltype(auto) translate_backend_errors(std::string_view operation, Call&& call) {
    try {
        return std::forward<Call>(call)();
    } catch (...) {
        rethrow_backend_error(operation);
    }
}

/**
 * Run a call into the thermodynamic backend whose failure means a bad argument (an unknown
 * species, phase or file), rethrowing a Cantera error as `std::invalid_argument`.
 *
 * @param operation what was being done, e.g. "Gas::create"
 * @return the value `call` returns
 */
template <typename Call>
decltype(auto) translate_backend_usage_errors(std::string_view operation, Call&& call) {
    try {
        return std::forward<Call>(call)();
    } catch (...) {
        rethrow_backend_usage_error(operation);
    }
}

} // namespace Goddard
