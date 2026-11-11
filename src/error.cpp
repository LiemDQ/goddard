#include "goddard/error.hpp"
#include "goddard/solve_guard.hpp"
#include "goddard/status.hpp"

#include "cantera/base/ctexceptions.h"

#include <cmath>
#include <format>
#include <string>
#include <utility>

namespace Goddard {

std::string_view to_string(FailureReason reason) {
    switch (reason) {
        case FailureReason::NONE:              return "NONE";
        case FailureReason::NOT_CONVERGED:     return "NOT_CONVERGED";
        case FailureReason::NO_SOLUTION:       return "NO_SOLUTION";
        case FailureReason::OUT_OF_DATA_RANGE: return "OUT_OF_DATA_RANGE";
        case FailureReason::THERMO_FAILURE:    return "THERMO_FAILURE";
        case FailureReason::INCOMPLETE:        return "INCOMPLETE";
        case FailureReason::NOT_COMPUTED:      return "NOT_COMPUTED";
    }
    return "UNKNOWN";
}

std::string SolveStatus::describe() const {
    if (ok()) {
        return "ok";
    }
    std::string text;
    if (!stage.empty()) {
        text += stage + ": ";
    }
    text += message;
    text += std::format(" ({}", to_string(reason));
    if (iterations >= 0) {
        text += std::format(", {} iterations", iterations);
    }
    if (!std::isnan(residual)) {
        text += std::format(", residual {:.3g}", residual);
    }
    if (!std::isnan(tolerance)) {
        text += std::format(", tolerance {:.3g}", tolerance);
    }
    text += ")";
    return text;
}

void SolveStatus::raise_if_failed() const {
    switch (reason) {
        case FailureReason::NONE:
            return;
        case FailureReason::NOT_CONVERGED:
            throw ConvergenceError(*this);
        case FailureReason::NO_SOLUTION:
        case FailureReason::OUT_OF_DATA_RANGE:
        case FailureReason::THERMO_FAILURE:
        case FailureReason::INCOMPLETE:
        case FailureReason::NOT_COMPUTED:
            throw SolveError(*this);
    }
}

SolveStatus SolveStatus::failure(FailureReason reason, std::string message) {
    SolveStatus status;
    status.reason = reason;
    status.message = std::move(message);
    return status;
}

SolveStatus SolveStatus::not_computed(const SolveStatus& upstream) {
    SolveStatus status;
    status.reason = FailureReason::NOT_COMPUTED;
    status.message = upstream.stage.empty()
        ? std::format("not computed because an earlier stage failed: {}", upstream.message)
        : std::format("not computed because the {} solve failed: {}", upstream.stage,
                      upstream.message);
    return status;
}

const char * NotImplementedError::what() const noexcept {
    return std::logic_error::what();
}

SolveError::SolveError(SolveStatus status)
    : std::runtime_error(status.describe()), m_status(std::move(status)) {}

SolveError::SolveError(FailureReason reason, const std::string& message)
    : SolveError(SolveStatus::failure(reason, message)) {}

namespace {

SolveStatus convergence_status(const std::string& what, int num_iters, double tol, double residual) {
    SolveStatus status = SolveStatus::failure(FailureReason::NOT_CONVERGED, what);
    status.iterations = num_iters;
    status.tolerance = tol;
    status.residual = residual;
    return status;
}

/** Text of a Cantera error without its banner: "<message> (<procedure>)". */
std::string backend_message(const Cantera::CanteraError& error) {
    std::string message = error.getMessage();
    // An input-file error is followed by an excerpt of the file, which does not belong in a
    // one-line message.
    const size_t excerpt = message.find("\n|  Line |");
    if (excerpt != std::string::npos) {
        message.erase(excerpt);
    }
    // Cantera messages may end in a newline or span several lines.
    while (!message.empty() && (message.back() == '\n' || message.back() == ' ')) {
        message.pop_back();
    }
    for (char& c : message) {
        if (c == '\n') {
            c = ' ';
        }
    }
    const std::string procedure = error.getMethod();
    if (procedure.empty()) {
        return message;
    }
    return std::format("{} (in {})", message, procedure);
}

} // namespace

ConvergenceError::ConvergenceError(const std::string& what, int num_iters, double tol, double residual)
    : SolveError(convergence_status(what, num_iters, tol, residual)) {}

ConvergenceError::ConvergenceError(SolveStatus status)
    : SolveError(std::move(status)) {}

SolveStatus status_from_current_exception(std::string_view stage) {
    SolveStatus status;
    try {
        throw;
    } catch (const SolveError& error) {
        status = error.status();
    } catch (const Cantera::CanteraError& error) {
        status = SolveStatus::failure(FailureReason::THERMO_FAILURE, backend_message(error));
    }
    if (status.stage.empty()) {
        status.stage = std::string(stage);
    }
    return status;
}

void rethrow_backend_error(std::string_view operation) {
    try {
        throw;
    } catch (const Cantera::CanteraError& error) {
        throw SolveError(FailureReason::THERMO_FAILURE,
            std::format("{}: {}", operation, backend_message(error)));
    }
}

void rethrow_backend_usage_error(std::string_view operation) {
    try {
        throw;
    } catch (const Cantera::CanteraError& error) {
        throw std::invalid_argument(std::format("{}: {}", operation, backend_message(error)));
    }
}

} // namespace Goddard
