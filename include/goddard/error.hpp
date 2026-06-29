#pragma once
#include <stdexcept>
#include <string>
#include <string_view>
#include <format>

namespace Goddard {

class NotImplementedError: public std::logic_error {
    using std::logic_error::logic_error;

    public:
    const char * what() const noexcept override;

};


class ConvergenceError: public std::runtime_error {

public:
    ConvergenceError(const std::string& what, int num_iters = -1, double tol = -1.0, double residual = -1.0);
    const char * what() const noexcept override;

private:
    int m_num_iters;
    double m_tol;
    double m_residual;
};

class FmtError: public std::runtime_error {
public:

    template <typename... Args>
    FmtError(std::string_view fmt, Args&&... args)
        : std::runtime_error(std::vformat(fmt, std::make_format_args(args...))) {}

    ~FmtError() noexcept override = default;

};

} // namespace Goddard
