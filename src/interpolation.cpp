#include "goddard/interpolation.hpp"
#include "goddard/error.hpp"
#include <algorithm>
#include <format>

namespace Goddard {

auto BicubicInterpolator::handle_boundary(double val, double low, double high, BoundsHandling option) const -> double {
    switch (option) {
        case BoundsHandling::ERROR: {
            throw std::invalid_argument(std::format("Value {} is out of bounds {}, {}", val, low, high));
        }
        case BoundsHandling::CLAMP: {
            return std::clamp(val, low, high);
        }
        default: {
            //unreachable
            std::runtime_error("Invalid value for BoundsHandling.");
        }
    }
}

auto BicubicInterpolator::interpolate(double x, double y) const -> double {
    if (x > x_max()) {
        x = handle_boundary(x, x_min, x_max(), m_bound_opts.x_high);
    }
    else if (x < x_min) {
        x = handle_boundary(x, x_min, x_max(), m_bound_opts.x_low);        
    }
    if (y > y_max()) {
        y = handle_boundary(y, y_min, y_max(), m_bound_opts.y_high);
    }
    else if (y < y_min) {
        y = handle_boundary(y, y_min, y_max(), m_bound_opts.y_low);
    }

    double x_grid = (x - x_min) / x_step;
    double y_grid = (y - y_min) / y_step;

    size_t idx_x1 = static_cast<size_t>(std::floor(x_grid));
    size_t idx_y1 = static_cast<size_t>(std::floor(y_grid));

    // Edge case on the exact boundary
    if (idx_x1 >= x_size - 1) idx_x1 = x_size - 2;
    if (idx_y1 >= y_size - 1) idx_y1 = y_size - 2;

    double t = x_grid - idx_x1;
    double u = y_grid - idx_y1;
    
    size_t ix[4], iy[4];
    for (int i = 0; i < 4; ++i) {
        long long rx = static_cast<long long>(idx_x1) - 1 + i;
        long long ry = static_cast<long long>(idx_y1) - 1 + i;
        ix[i] = std::clamp(rx, 0LL, static_cast<long long>(x_size - 1));
        iy[i] = std::clamp(ry, 0LL, static_cast<long long>(y_size - 1));
    }

    // x-axis interpolatation
    double row_arr[4];
    for (int j = 0; j < 4; ++j) {
        row_arr[j] = cubic_interpolate(
            m_values(ix[0], iy[j]),
            m_values(ix[1], iy[j]),
            m_values(ix[2], iy[j]),
            m_values(ix[3], iy[j]),
            t
        );
    }

    return cubic_interpolate(row_arr[0], row_arr[1], row_arr[2], row_arr[3], u);
}

auto BicubicInterpolator::x_max() const -> double {
    return x_min + x_step * (x_size - 1);
}


auto BicubicInterpolator::y_max() const -> double {
    return y_min + y_step * (y_size - 1);
}

} // namespace Goddard