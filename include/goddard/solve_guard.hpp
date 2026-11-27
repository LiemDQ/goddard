#pragma once
#include <format>
#include <string>
#include <string_view>
#include <type_traits>
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

class Gas;

/**
 * Range check of a state a solver reports: `Gas::check_temperature_range()` of `gas`, with its
 * warning added to `status`. Call it at every state a result reports, inside the body of the
 * solve's `capture_failure`: a state far outside the data throws `SolveError`
 * (`OUT_OF_DATA_RANGE`), which then becomes the status of the result like any other failure.
 */
void check_reported_state(const Gas& gas, SolveStatus& status);

/**
 * Add the warnings of `part`, a part of a composite result, to the status of the whole, each
 * prefixed with `label` (e.g. "incident shock: ...").
 */
inline void add_part_warnings(SolveStatus& status, std::string_view label, const SolveStatus& part) {
    for (const std::string& warning : part.warnings) {
        status.warnings.push_back(std::format("{}: {}", label, warning));
    }
}

/**
 * `capture_failure` for a solve whose result is the state `gas` holds when `body` returns: the
 * status carries the range check of that state (`check_reported_state`).
 */
template <typename Body>
SolveStatus capture_failure(std::string_view stage, const Gas& gas, Body&& body) {
    SolveStatus status;
    try {
        std::forward<Body>(body)();
        check_reported_state(gas, status);
    } catch (...) {
        return status_from_current_exception(stage);
    }
    return status;
}

/**
 * Reject an equilibrium request the backend does not know, before any solve is attempted: a
 * misspelt property pair or solver name is a bad argument, not a failed solve.
 *
 * @param XY pair of properties held constant: "TP", "TV", "HP", "SP", "SV" or "UV"
 * @param solver "auto", "element_potential", "gibbs" or "vcs"
 * @throws std::invalid_argument
 */
void validate_equilibrium_request(const std::string& XY, const std::string& solver);

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

namespace detail {

/**
 * Text of an operation description: a string-like value as it is, or the result of a callable
 * returning `std::string`, which is called only here.
 */
template <typename Operation>
std::string operation_text(Operation&& operation) {
    if constexpr (std::is_invocable_v<Operation>) {
        return std::string(std::forward<Operation>(operation)());
    } else {
        return std::string(std::string_view(operation));
    }
}

} // namespace detail

/**
 * Run a call into the thermodynamic backend, rethrowing a Cantera error as a `SolveError` with
 * reason `THERMO_FAILURE`.
 *
 * @param operation what was being done, with the values involved, e.g.
 * "Gas::set_state_TP(T = 300 K, P = 1e5 Pa)". Either a string, or a callable returning one that
 * is invoked only when `call` fails, so that a message built with `std::format` costs nothing on
 * the success path.
 * @return the value `call` returns
 */
template <typename Operation, typename Call>
decltype(auto) translate_backend_errors(Operation&& operation, Call&& call) {
    try {
        return std::forward<Call>(call)();
    } catch (...) {
        rethrow_backend_error(detail::operation_text(std::forward<Operation>(operation)));
    }
}

/**
 * Run a call into the thermodynamic backend whose failure means a bad argument (an unknown
 * species, phase or file), rethrowing a Cantera error as `std::invalid_argument`.
 *
 * @param operation what was being done, e.g. "Gas::create"; a string or a callable returning
 * one, as for `translate_backend_errors`
 * @return the value `call` returns
 */
template <typename Operation, typename Call>
decltype(auto) translate_backend_usage_errors(Operation&& operation, Call&& call) {
    try {
        return std::forward<Call>(call)();
    } catch (...) {
        rethrow_backend_usage_error(detail::operation_text(std::forward<Operation>(operation)));
    }
}

} // namespace Goddard
