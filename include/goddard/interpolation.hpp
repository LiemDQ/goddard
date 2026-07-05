#pragma once
#include <Eigen/Dense>
#include <vector>
#include <utility>

namespace Goddard {
/** Perform a cubic spline interpolation. 
 * t is the fractional distance of the interpolated point between l1 and r1 (0.0 <= t <= 1.0) 
 * */ 
constexpr double cubic_interpolate(double l2, double l1, double r1, double r2, double t) {
    return l1 + 0.5 * t * (r1 - l2 + t * (2.0 * l2 - 5.0 * l1 + 4.0 * r1 - r2 + t * (3.0 * (l1 - r1) + r2 - l2)));
}

enum class BoundsHandling {
    ERROR,
    CLAMP,
};

struct BoundsOptions {
    BoundsHandling x_low, x_high;
    BoundsHandling y_low, y_high;
};

/**
 * Performs bicubic spline interpolation over a 2D grid.
 * @warning Assumes grid points are evenly spaced. 
 */
class BicubicInterpolator {
public:
    const double x_min, x_step;
    const double y_min, y_step;
    const size_t x_size, y_size;

private:
    // Declared after x_size/y_size: members initialize in declaration order
    // (not initializer-list order), and x_size/y_size are computed from the
    // `values` constructor parameter before it is moved into m_values below.
    Eigen::MatrixXd m_values;
    BoundsOptions m_bound_opts;

    auto handle_boundary(double val, double low, double high, BoundsHandling option) const -> double;

public:
    BicubicInterpolator(double x_start, double x_spacing,
                        double y_start, double y_spacing,
                        Eigen::MatrixXd values,
                        BoundsOptions bounds = {})
    : x_min(x_start), x_step(x_spacing),
      y_min(y_start), y_step(y_spacing),
      x_size(values.rows()), y_size(values.cols()),
      m_values(std::move(values)), m_bound_opts(bounds) 
    {}
    
    auto interpolate(double x, double y) const -> double;

    auto x_max() const -> double;
    auto y_max() const -> double;


};


} // namespace Goddard