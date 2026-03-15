#include "goddard/error.hpp"
#include <string>

namespace Goddard {

const char * NotImplementedError::what() const noexcept {
    return std::logic_error::what();
}

ConvergenceError::ConvergenceError(const std::string& what, int num_iters, double tol, double residual)
    : std::runtime_error(what), m_num_iters(num_iters), m_tol(tol), m_residual(residual) {}

const char * ConvergenceError::what() const noexcept {
    return std::runtime_error::what();
}

}