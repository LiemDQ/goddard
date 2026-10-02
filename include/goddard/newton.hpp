#pragma once
#include <cmath>
#include <type_traits>
#include "eigen3/Eigen/Dense"

/**
 * Newton's method for a scalar equation and for small systems of fixed size, shared by the
 * library's solvers. Internal: not part of the Python API.
 *
 * One iteration, from the current point x:
 *   1. Evaluate the function: f and its derivative. The Newton step is dx = -f/f'.
 *   2. If |f| <= residual_abstol, stop at x.
 *   3. If |dx| <= step_abstol + step_reltol*|x|, or the step no longer changes x, apply the
 *      step and stop.
 *   4. Otherwise pass the step through the limiter and apply the result.
 *
 * The function may set state as a side effect (e.g. a `Gas`). When the iteration stops on the
 * step test, that state belongs to the last evaluated point, one converged step before the
 * returned x. A caller that needs the state at x sets it once more.
 *
 * `newton_solve` does not throw: the caller turns a status other than CONVERGED into its own
 * error.
 */

namespace Goddard {

/** f(x) and f'(x) at one point. */
struct NewtonFunction {
    double value;
    double derivative;
};

/** f(x) and the Jacobian df_i/dx_j at one point. */
template <int N>
struct NewtonSystemFunction {
    Eigen::Matrix<double, N, 1> value;
    Eigen::Matrix<double, N, N> jacobian;
};

enum class NewtonStatus {
    CONVERGED,
    /** `NewtonOptions::max_iterations` evaluations were made without meeting a tolerance. */
    MAX_ITERATIONS,
    /** The step to apply was not finite, e.g. because of a zero derivative. */
    NON_FINITE_STEP
};

/**
 * Tolerances are in the units of x (steps) and of f (residual). A tolerance of zero disables
 * its test, except that an exactly zero step or residual always counts as converged. For a
 * system, |.| is the largest absolute component.
 */
struct NewtonOptions {
    double step_abstol = 0.0;
    double step_reltol = 0.0;
    double residual_abstol = 0.0;
    /** Maximum number of function evaluations. */
    int max_iterations = 50;
};

struct NewtonResult {
    NewtonStatus status;
    double x;
    /** Newton step at the last evaluated point. Already included in x if the step test passed. */
    double step;
    /** f at the last evaluated point. */
    double residual;
    /** Number of function evaluations. */
    int iterations;
};

template <int N>
struct NewtonSystemResult {
    NewtonStatus status;
    Eigen::Matrix<double, N, 1> x;
    /** Newton step at the last evaluated point. Already included in x if the step test passed. */
    Eigen::Matrix<double, N, 1> step;
    /** f at the last evaluated point. */
    Eigen::Matrix<double, N, 1> residual;
    /** Number of function evaluations. */
    int iterations;
};

/**
 * Solve f(x) = 0 from the initial guess `x0`.
 *
 * @param function    Callable as `NewtonFunction(double x)`.
 * @param limit_step  Callable as `double(int iteration, double x, double step)`, with the
 *                    zero-based iteration, the current point and its Newton step. Returns the
 *                    step to apply. Not called for the final, converged step.
 */
template <typename Function, typename LimitStep>
NewtonResult newton_solve(double x0, Function&& function, LimitStep&& limit_step,
                          const NewtonOptions& options)
{
    static_assert(std::is_invocable_r_v<NewtonFunction, Function&, double>,
                  "newton_solve: function must be callable as NewtonFunction(double x)");
    static_assert(std::is_invocable_r_v<double, LimitStep&, int, double, double>,
                  "newton_solve: limit_step must be callable as double(int iteration, double x, double step)");

    NewtonResult result{NewtonStatus::MAX_ITERATIONS, x0, 0.0, 0.0, 0};
    for (int k = 0; k < options.max_iterations; k++) {
        const NewtonFunction evaluated = function(result.x);
        result.iterations = k + 1;
        result.residual = evaluated.value;
        result.step = -evaluated.value/evaluated.derivative;

        if (std::abs(result.residual) <= options.residual_abstol) {
            result.status = NewtonStatus::CONVERGED;
            return result;
        }
        const double step_tolerance = options.step_abstol + options.step_reltol*std::abs(result.x);
        if (std::abs(result.step) <= step_tolerance || result.x + result.step == result.x) {
            result.x += result.step;
            result.status = NewtonStatus::CONVERGED;
            return result;
        }

        const double limited_step = limit_step(k, result.x, result.step);
        if (!std::isfinite(limited_step)) {
            result.status = NewtonStatus::NON_FINITE_STEP;
            return result;
        }
        result.x += limited_step;
    }
    return result;
}

/** Solve f(x) = 0 from the initial guess `x0`, taking full Newton steps. */
template <typename Function>
NewtonResult newton_solve(double x0, Function&& function, const NewtonOptions& options)
{
    return newton_solve(x0, function, [](int, double, double step) { return step; }, options);
}

/**
 * Solve the system f(x) = 0 of N equations from the initial guess `x0`.
 *
 * @param function    Callable as `NewtonSystemFunction<N>(const Eigen::Matrix<double, N, 1>& x)`.
 * @param limit_step  Callable as `Eigen::Matrix<double, N, 1>(int iteration, const
 *                    Eigen::Matrix<double, N, 1>& x, const Eigen::Matrix<double, N, 1>& step)`.
 *                    Returns the step to apply. Not called for the final, converged step.
 */
template <int N, typename Function, typename LimitStep>
NewtonSystemResult<N> newton_solve(const Eigen::Matrix<double, N, 1>& x0, Function&& function,
                                   LimitStep&& limit_step, const NewtonOptions& options)
{
    using Vector = Eigen::Matrix<double, N, 1>;
    static_assert(std::is_invocable_r_v<NewtonSystemFunction<N>, Function&, const Vector&>,
                  "newton_solve: function must be callable as "
                  "NewtonSystemFunction<N>(const Eigen::Matrix<double, N, 1>& x)");
    static_assert(std::is_invocable_r_v<Vector, LimitStep&, int, const Vector&, const Vector&>,
                  "newton_solve: limit_step must be callable as Eigen::Matrix<double, N, 1>("
                  "int iteration, const Eigen::Matrix<double, N, 1>& x, "
                  "const Eigen::Matrix<double, N, 1>& step)");

    NewtonSystemResult<N> result{NewtonStatus::MAX_ITERATIONS, x0, Vector::Zero(), Vector::Zero(), 0};
    for (int k = 0; k < options.max_iterations; k++) {
        const NewtonSystemFunction<N> evaluated = function(result.x);
        result.iterations = k + 1;
        result.residual = evaluated.value;
        result.step = evaluated.jacobian.partialPivLu().solve(-evaluated.value);

        if (result.residual.cwiseAbs().maxCoeff() <= options.residual_abstol) {
            result.status = NewtonStatus::CONVERGED;
            return result;
        }
        const double step_tolerance =
            options.step_abstol + options.step_reltol*result.x.cwiseAbs().maxCoeff();
        if (result.step.cwiseAbs().maxCoeff() <= step_tolerance || result.x + result.step == result.x) {
            result.x += result.step;
            result.status = NewtonStatus::CONVERGED;
            return result;
        }

        const Vector limited_step = limit_step(k, result.x, result.step);
        if (!limited_step.allFinite()) {
            result.status = NewtonStatus::NON_FINITE_STEP;
            return result;
        }
        result.x += limited_step;
    }
    return result;
}

/** Solve the system f(x) = 0 of N equations from the initial guess `x0`, taking full Newton steps. */
template <int N, typename Function>
NewtonSystemResult<N> newton_solve(const Eigen::Matrix<double, N, 1>& x0, Function&& function,
                                   const NewtonOptions& options)
{
    using Vector = Eigen::Matrix<double, N, 1>;
    return newton_solve(
        x0, function, [](int, const Vector&, const Vector& step) -> Vector { return step; }, options);
}

} // namespace Goddard
