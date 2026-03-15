#pragma once
#include <stdexcept>
#include <string>

namespace Goddard {

class NotImplementedError: public std::logic_error {
    using std::logic_error::logic_error;

    public:
    virtual const char * what() const noexcept override;

};


class ConvergenceError: public std::runtime_error {

    public:
    ConvergenceError(const std::string& what, int num_iters = -1, double tol = -1.0, double residual = -1.0);
    virtual const char * what() const noexcept override;

    private:
    int m_num_iters;
    double m_tol;
    double m_residual;
};

}