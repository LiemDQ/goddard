#include "goddard/interpolation.hpp"
#include <Eigen/Dense>
#include <stdexcept>
#include "gtest/gtest.h"

using namespace Goddard;

namespace {

// Synthetic grid dimensions chosen so that every axis has at least 3 "interior"
// segments whose 4-point stencil never needs boundary clamping (needed to test
// exact off-grid recovery without boundary effects contaminating the result).
constexpr size_t GRID_ROWS = 6; // bound to the interpolator's x axis
constexpr size_t GRID_COLS = 8; // bound to the interpolator's y axis

constexpr double GRID_A = 1.0;
constexpr double GRID_B = 2.0;
constexpr double GRID_C = -1.5;

// f(x, y) = a + b*x + c*y on integer grid coordinates (x_min=0, x_step=1, etc.)
Eigen::MatrixXd make_bilinear_grid() {
    Eigen::MatrixXd m(GRID_ROWS, GRID_COLS);
    for (size_t i = 0; i < GRID_ROWS; i++) {
        for (size_t j = 0; j < GRID_COLS; j++) {
            m(i, j) = GRID_A + GRID_B * static_cast<double>(i) + GRID_C * static_cast<double>(j);
        }
    }
    return m;
}

// f(x, y) = x^2, independent of y.
Eigen::MatrixXd make_quadratic_in_x_grid() {
    Eigen::MatrixXd m(GRID_ROWS, GRID_COLS);
    for (size_t i = 0; i < GRID_ROWS; i++) {
        for (size_t j = 0; j < GRID_COLS; j++) {
            m(i, j) = static_cast<double>(i) * static_cast<double>(i);
        }
    }
    return m;
}

constexpr BoundsOptions ALL_ERROR{
    .x_low = BoundsHandling::ERROR, .x_high = BoundsHandling::ERROR,
    .y_low = BoundsHandling::ERROR, .y_high = BoundsHandling::ERROR,
};

BicubicInterpolator make_interpolator(BoundsOptions opts, const Eigen::MatrixXd& values) {
    return BicubicInterpolator(0.0, 1.0, 0.0, 1.0, values, opts);
}

} // namespace

// ============================================================
// cubic_interpolate (free function)
// ============================================================

TEST(CubicInterpolateTest, ReturnsL1AtT0) {
    EXPECT_DOUBLE_EQ(cubic_interpolate(3.0, 5.0, 8.0, 1.0, 0.0), 5.0);
}

TEST(CubicInterpolateTest, ReturnsR1AtT1) {
    EXPECT_DOUBLE_EQ(cubic_interpolate(3.0, 5.0, 8.0, 1.0, 1.0), 8.0);
}

TEST(CubicInterpolateTest, ExactOnLinearData) {
    // l2, l1, r1, r2 in arithmetic progression with common difference d.
    double d = 3.0;
    double l1 = 5.0;
    double l2 = l1 - d;
    double r1 = l1 + d;
    double r2 = l1 + 2 * d;
    for (double t : {0.0, 0.2, 0.5, 0.7, 1.0}) {
        EXPECT_NEAR(cubic_interpolate(l2, l1, r1, r2, t), l1 + t * d, 1e-12)
            << "Catmull-Rom cubic must reproduce linear data exactly at t=" << t;
    }
}

// ============================================================
// BicubicInterpolator construction and grid geometry
// ============================================================

TEST(BicubicInterpolatorTest, XSizeYSizeMatchRowsAndCols) {
    // This is the exact convention the Rao nozzle transposition bug got backwards:
    // the first coordinate binds to matrix rows, the second to columns.
    BicubicInterpolator interp = make_interpolator(ALL_ERROR, make_bilinear_grid());
    EXPECT_EQ(interp.x_size, GRID_ROWS);
    EXPECT_EQ(interp.y_size, GRID_COLS);
}

TEST(BicubicInterpolatorTest, XMaxMatchesTrueLastGridCoordinate) {
    // Regression for the x_max()/y_max() off-by-one: with min=0, step=0.5, and 4
    // grid points (indices 0..3), the true last coordinate is 0+0.5*3=1.5, not
    // 0+0.5*4=2.0.
    BicubicInterpolator interp(0.0, 0.5, 0.0, 1.0, Eigen::MatrixXd::Zero(4, 3), ALL_ERROR);
    EXPECT_NEAR(interp.x_max(), 1.5, 1e-12);
}

TEST(BicubicInterpolatorTest, YMaxMatchesTrueLastGridCoordinate) {
    BicubicInterpolator interp(0.0, 1.0, 0.0, 0.5, Eigen::MatrixXd::Zero(3, 4), ALL_ERROR);
    EXPECT_NEAR(interp.y_max(), 1.5, 1e-12);
}

// ============================================================
// Interpolation accuracy
// ============================================================

TEST(BicubicInterpolatorTest, ExactRecoveryAtGridPoints) {
    Eigen::MatrixXd values = make_bilinear_grid();
    BicubicInterpolator interp = make_interpolator(ALL_ERROR, values);
    for (size_t i = 0; i < GRID_ROWS; i++) {
        for (size_t j = 0; j < GRID_COLS; j++) {
            double expected = GRID_A + GRID_B * static_cast<double>(i) + GRID_C * static_cast<double>(j);
            EXPECT_NEAR(interp.interpolate(static_cast<double>(i), static_cast<double>(j)), expected, 1e-9)
                << "at grid node (" << i << ", " << j << ")";
        }
    }
}

TEST(BicubicInterpolatorTest, RecoversBilinearSurfaceOffGrid) {
    BicubicInterpolator interp = make_interpolator(ALL_ERROR, make_bilinear_grid());
    // Interior points only: at least one full cell of margin from every edge, so
    // the 4-point stencil never clamps into a repeated node. A bilinear surface
    // (degree 1 in both axes) must be reproduced exactly regardless of where the
    // query falls within the interior.
    struct Point { double x, y; };
    for (Point p : {Point{1.3, 2.6}, Point{3.6, 1.4}, Point{2.5, 4.9}, Point{1.9, 5.5}}) {
        double expected = GRID_A + GRID_B * p.x + GRID_C * p.y;
        EXPECT_NEAR(interp.interpolate(p.x, p.y), expected, 1e-9) << "at (" << p.x << ", " << p.y << ")";
    }
}

TEST(BicubicInterpolatorTest, RecoversQuadraticSurfaceInX) {
    BicubicInterpolator interp = make_interpolator(ALL_ERROR, make_quadratic_in_x_grid());
    for (double x : {1.2, 1.5, 1.8, 2.3, 2.6, 3.4}) {
        double expected = x * x;
        EXPECT_NEAR(interp.interpolate(x, 2.0), expected, 1e-9) << "at x=" << x;
    }
}

// ============================================================
// Boundary handling — all 4 (axis, direction) combinations
// ============================================================

TEST(BicubicInterpolatorTest, XLowErrorThrowsBelowMin) {
    BoundsOptions opts{
        .x_low = BoundsHandling::ERROR, .x_high = BoundsHandling::CLAMP,
        .y_low = BoundsHandling::ERROR, .y_high = BoundsHandling::ERROR,
    };
    BicubicInterpolator interp = make_interpolator(opts, make_bilinear_grid());
    EXPECT_THROW(interp.interpolate(-0.1, 2.0), std::invalid_argument);
}

TEST(BicubicInterpolatorTest, XHighClampMatchesValueAtXMax) {
    BoundsOptions opts{
        .x_low = BoundsHandling::ERROR, .x_high = BoundsHandling::CLAMP,
        .y_low = BoundsHandling::ERROR, .y_high = BoundsHandling::ERROR,
    };
    BicubicInterpolator interp = make_interpolator(opts, make_bilinear_grid());
    double at_max = interp.interpolate(interp.x_max(), 2.0);
    double beyond = interp.interpolate(interp.x_max() + 0.5, 2.0);
    EXPECT_NEAR(beyond, at_max, 1e-9);
}

TEST(BicubicInterpolatorTest, XHighErrorThrowsJustBeyondTrueMax) {
    // Regression for the x_max() off-by-one: pre-fix, x_max() was one full grid
    // step too large, so a query just past the TRUE last grid coordinate would
    // silently not throw and would extrapolate instead.
    BicubicInterpolator interp = make_interpolator(ALL_ERROR, make_bilinear_grid());
    EXPECT_THROW(interp.interpolate(interp.x_max() + 0.01, 2.0), std::invalid_argument);
}

TEST(BicubicInterpolatorTest, YLowErrorThrowsBelowMinWithAsymmetricOpts) {
    // Regression for the handle_boundary dispatch bug: the y<y_min branch used
    // to dispatch on y_high instead of y_low. A symmetric bounds config can't
    // distinguish the two, so y_low and y_high are deliberately set differently.
    BoundsOptions opts{
        .x_low = BoundsHandling::ERROR, .x_high = BoundsHandling::ERROR,
        .y_low = BoundsHandling::ERROR, .y_high = BoundsHandling::CLAMP,
    };
    BicubicInterpolator interp = make_interpolator(opts, make_bilinear_grid());
    EXPECT_THROW(interp.interpolate(1.0, -0.1), std::invalid_argument);
}

TEST(BicubicInterpolatorTest, YLowClampClampsBelowMinWithAsymmetricOpts) {
    BoundsOptions opts{
        .x_low = BoundsHandling::ERROR, .x_high = BoundsHandling::ERROR,
        .y_low = BoundsHandling::CLAMP, .y_high = BoundsHandling::ERROR,
    };
    BicubicInterpolator interp = make_interpolator(opts, make_bilinear_grid());
    double at_min = interp.interpolate(1.0, 0.0);
    double below = interp.interpolate(1.0, -0.3);
    EXPECT_NEAR(below, at_min, 1e-9);
}

TEST(BicubicInterpolatorTest, YHighClampMatchesValueAtYMax) {
    BoundsOptions opts{
        .x_low = BoundsHandling::ERROR, .x_high = BoundsHandling::ERROR,
        .y_low = BoundsHandling::ERROR, .y_high = BoundsHandling::CLAMP,
    };
    BicubicInterpolator interp = make_interpolator(opts, make_bilinear_grid());
    double at_max = interp.interpolate(1.0, interp.y_max());
    double beyond = interp.interpolate(1.0, interp.y_max() + 0.7);
    EXPECT_NEAR(beyond, at_max, 1e-9);
}

TEST(BicubicInterpolatorTest, YHighErrorThrowsAboveMax) {
    BicubicInterpolator interp = make_interpolator(ALL_ERROR, make_bilinear_grid());
    EXPECT_THROW(interp.interpolate(1.0, interp.y_max() + 0.01), std::invalid_argument);
}
