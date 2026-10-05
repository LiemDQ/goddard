#include <algorithm>
#include <cmath>
#include <vector>
#include "goddard/newton.hpp"
#include "gtest/gtest.h"

using namespace Goddard;

namespace {

// f(x) = x^2 - 2, root sqrt(2)
NewtonFunction square_minus_two(double x) {
    return NewtonFunction{.value = x*x - 2.0, .derivative = 2.0*x};
}

// f(x) = atan(x), root 0. Full Newton steps diverge from |x0| > 1.3917.
NewtonFunction arctangent(double x) {
    return NewtonFunction{.value = std::atan(x), .derivative = 1.0/(1.0 + x*x)};
}

// Circle x^2 + y^2 = 4 and hyperbola x y = 1. Adding and subtracting 2xy gives (x + y)^2 = 6 and
// (x - y)^2 = 2, so the root with x > y > 0 is ((sqrt(6) + sqrt(2))/2, (sqrt(6) - sqrt(2))/2).
NewtonSystemFunction<2> circle_and_hyperbola(const Eigen::Vector2d& x) {
    NewtonSystemFunction<2> result;
    result.value << x(0)*x(0) + x(1)*x(1) - 4.0, x(0)*x(1) - 1.0;
    result.jacobian << 2.0*x(0), 2.0*x(1),
                       x(1), x(0);
    return result;
}

Eigen::Vector2d circle_and_hyperbola_root() {
    return Eigen::Vector2d((std::sqrt(6.0) + std::sqrt(2.0))/2.0, (std::sqrt(6.0) - std::sqrt(2.0))/2.0);
}

} // anonymous namespace

TEST(NewtonTests, ScalarConvergesToAnalyticRoot) {
    NewtonOptions options;
    options.step_abstol = 1e-10;
    const NewtonResult result = newton_solve(1.0, square_minus_two, options);

    EXPECT_EQ(result.status, NewtonStatus::CONVERGED);
    // The final step is below 1e-10 and is applied, so the error is at round-off.
    EXPECT_NEAR(result.x, std::sqrt(2.0), 1e-14);
}

TEST(NewtonTests, FinalStepIsAppliedWhenTheStepTestPasses) {
    std::vector<double> evaluated_points;
    NewtonOptions options;
    options.step_abstol = 1e-3;
    const NewtonResult result = newton_solve(1.0, [&](double x) {
        evaluated_points.push_back(x);
        return square_minus_two(x);
    }, options);

    ASSERT_EQ(result.status, NewtonStatus::CONVERGED);
    EXPECT_LE(std::abs(result.step), options.step_abstol);
    EXPECT_DOUBLE_EQ(result.x, evaluated_points.back() + result.step);
    // Quadratic convergence: with the last step applied the error is of the order of the square
    // of that step, far below the step tolerance.
    EXPECT_NEAR(result.x, std::sqrt(2.0), 1e-9);
    EXPECT_GT(std::abs(evaluated_points.back() - std::sqrt(2.0)), 1e-9);
}

TEST(NewtonTests, ResidualConvergenceStopsAtTheEvaluatedPoint) {
    std::vector<double> evaluated_points;
    NewtonOptions options;
    options.residual_abstol = 1e-3;
    const NewtonResult result = newton_solve(1.0, [&](double x) {
        evaluated_points.push_back(x);
        return square_minus_two(x);
    }, options);

    ASSERT_EQ(result.status, NewtonStatus::CONVERGED);
    EXPECT_EQ(result.x, evaluated_points.back());
    EXPECT_LE(std::abs(result.residual), options.residual_abstol);
    EXPECT_DOUBLE_EQ(result.residual, result.x*result.x - 2.0);
}

TEST(NewtonTests, RelativeStepToleranceScalesWithTheIterate) {
    // f(x) = x^2 - 1e12, root 1e6. An absolute step of 1e-3 passes a relative tolerance of 1e-8
    // here, where it would fail the same value used as an absolute tolerance.
    NewtonOptions options;
    options.step_reltol = 1e-8;
    const NewtonResult result = newton_solve(3e6, [](double x) {
        return NewtonFunction{.value = x*x - 1e12, .derivative = 2.0*x};
    }, options);

    ASSERT_EQ(result.status, NewtonStatus::CONVERGED);
    EXPECT_LE(std::abs(result.step), options.step_reltol*1e6*1.01);
    EXPECT_NEAR(result.x, 1e6, 1e-8*1e6);
}

TEST(NewtonTests, StepLimiterTurnsADivergingIterationIntoAConvergingOne) {
    NewtonOptions options;
    options.step_abstol = 1e-12;

    const NewtonResult unlimited = newton_solve(2.0, arctangent, options);
    EXPECT_NE(unlimited.status, NewtonStatus::CONVERGED);

    std::vector<double> evaluated_points;
    const double largest_step = 1.0;
    const NewtonResult limited = newton_solve(2.0, [&](double x) {
        evaluated_points.push_back(x);
        return arctangent(x);
    }, [&](int, double, double step) {
        return std::clamp(step, -largest_step, largest_step);
    }, options);

    ASSERT_EQ(limited.status, NewtonStatus::CONVERGED);
    EXPECT_NEAR(limited.x, 0.0, 1e-12);
    for (size_t i = 1; i < evaluated_points.size(); i++) {
        EXPECT_LE(std::abs(evaluated_points[i] - evaluated_points[i - 1]), largest_step*(1.0 + 1e-12));
    }
}

TEST(NewtonTests, StepLimiterReceivesTheIterationIndexAndCurrentPoint) {
    std::vector<double> evaluated_points;
    std::vector<int> limiter_iterations;
    std::vector<double> limiter_points;
    NewtonOptions options;
    options.step_abstol = 1e-10;
    const NewtonResult result = newton_solve(1.0, [&](double x) {
        evaluated_points.push_back(x);
        return square_minus_two(x);
    }, [&](int iteration, double x, double step) {
        limiter_iterations.push_back(iteration);
        limiter_points.push_back(x);
        return step;
    }, options);

    ASSERT_EQ(result.status, NewtonStatus::CONVERGED);
    // The limiter sees every step but the final, converged one.
    ASSERT_EQ(limiter_points.size() + 1, evaluated_points.size());
    for (size_t i = 0; i < limiter_points.size(); i++) {
        EXPECT_EQ(limiter_iterations[i], static_cast<int>(i));
        EXPECT_EQ(limiter_points[i], evaluated_points[i]);
    }
}

TEST(NewtonTests, EvaluationBudgetIsReported) {
    int evaluations = 0;
    NewtonOptions options;
    options.step_abstol = 1e-12;
    options.max_iterations = 3;
    const NewtonResult result = newton_solve(100.0, [&](double x) {
        evaluations++;
        return square_minus_two(x);
    }, options);

    EXPECT_EQ(result.status, NewtonStatus::MAX_ITERATIONS);
    EXPECT_EQ(evaluations, options.max_iterations);
    EXPECT_EQ(result.iterations, options.max_iterations);
}

TEST(NewtonTests, ZeroDerivativeIsReportedAsNonFiniteStep) {
    NewtonOptions options;
    options.step_abstol = 1e-12;
    const NewtonResult result = newton_solve(0.0, square_minus_two, options);

    EXPECT_EQ(result.status, NewtonStatus::NON_FINITE_STEP);
    EXPECT_EQ(result.x, 0.0);
}

TEST(NewtonTests, SystemConvergesToAnalyticRoot) {
    NewtonOptions options;
    options.step_abstol = 1e-10;
    const NewtonSystemResult<2> result =
        newton_solve(Eigen::Vector2d(2.0, 0.3), circle_and_hyperbola, options);

    ASSERT_EQ(result.status, NewtonStatus::CONVERGED);
    EXPECT_NEAR(result.x(0), circle_and_hyperbola_root()(0), 1e-13);
    EXPECT_NEAR(result.x(1), circle_and_hyperbola_root()(1), 1e-13);
}

TEST(NewtonTests, SystemStepLimiterBoundsEveryAppliedStep) {
    std::vector<Eigen::Vector2d> evaluated_points;
    const double largest_step = 0.05;
    NewtonOptions options;
    options.step_abstol = 1e-10;
    options.max_iterations = 200;
    const NewtonSystemResult<2> result = newton_solve(Eigen::Vector2d(2.0, 0.3),
        [&](const Eigen::Vector2d& x) {
            evaluated_points.push_back(x);
            return circle_and_hyperbola(x);
        },
        [&](int, const Eigen::Vector2d&, const Eigen::Vector2d& step) -> Eigen::Vector2d {
            // one factor for both components, so the direction of the step is kept
            return step*std::min(1.0, largest_step/step.cwiseAbs().maxCoeff());
        }, options);

    ASSERT_EQ(result.status, NewtonStatus::CONVERGED);
    EXPECT_NEAR(result.x(0), circle_and_hyperbola_root()(0), 1e-13);
    EXPECT_NEAR(result.x(1), circle_and_hyperbola_root()(1), 1e-13);
    for (size_t i = 1; i < evaluated_points.size(); i++) {
        const double applied = (evaluated_points[i] - evaluated_points[i - 1]).cwiseAbs().maxCoeff();
        EXPECT_LE(applied, largest_step*(1.0 + 1e-12));
    }
}

TEST(NewtonTests, SingularJacobianIsReportedAsNonFiniteStep) {
    // At the origin the Jacobian of the circle and hyperbola system is the zero matrix.
    NewtonOptions options;
    options.step_abstol = 1e-10;
    const NewtonSystemResult<2> result =
        newton_solve(Eigen::Vector2d(0.0, 0.0), circle_and_hyperbola, options);

    EXPECT_EQ(result.status, NewtonStatus::NON_FINITE_STEP);
}
