#pragma once
#include <optional>
#include <vector>
#include "goddard/moc.hpp"

namespace Goddard {

/**
 * Main class for performing 2D nozzle supersonic flow simulations using
 * the method of characteristics (MoC). Simulation specification is done with `MocOptions`.
 *
 * There are two modes of usage:
 * - Design mode, where a specific mach number or other end goal is specified
 * and the nozzle profile is determined by the solver.
 * - Analysis mode, where an arbitrary nozzle profile is provided and the solver determines
 * the flow field and performance metrics.
 *
 * solve() is the orchestrator: it validates `options`, resolves the throat and thermo model
 * (MocThermo), builds and measures the start line (moc_initialization.hpp), then dispatches
 * to one of the two kernels by MocMode -- DirectMarch for DESIGN_MIN_LENGTH, InverseMarch for
 * ANALYSIS and DESIGN_RAO (moc_direct_march.hpp, moc_inverse_march.hpp) -- and assembles the
 * result.
*/
class MocNozzle {
public:

    MocNozzle(MocOptions options): options(options) { validate_moc_options(options); }

    MocNozzle(Gas gas, MocOptions options): options(options), m_gas(gas) { validate_moc_options(options); }

    /**
     * March the characteristic net and return the solved flow field.
     *
     * Never throws for numerical failures: a solve that breaks down returns a MocResult with
     * `converged == false` and `failure` identifying what went wrong and where. It does throw
     * std::invalid_argument if `options` is self-inconsistent (checked here as well as in
     * the constructor, since `options` is public and callers do adjust it in between), and
     * NotImplementedError for MocMode::DESIGN_CENTERLINE.
     */
    MocResult solve();
    /** True once solve() has run to completion on this instance. */
    bool is_solved() const;

    /// Solver configuration. Public and re-read by solve(), so it can be adjusted between solves.
    MocOptions options;

protected:

    /**
     * Test-only hook (kept protected; set only by a test subclass): when present, solve()
     * seeds the inverse kernel's first front directly from these points instead of running
     * build_start_line, so a hand-built front with a known exact solution (e.g. uniform flow,
     * or a manufactured source-flow field) can be marched without requiring a throat/KL/fan
     * construction consistent with it. See InverseMarchUniformFlow.StaysUniform and
     * InverseMarchSourceFlow.SecondOrderConvergence (test/test_moc_inverse_march.cpp).
     */
    std::optional<std::vector<CharacteristicPoint>> m_inverse_front_override;

private:
    std::optional<Gas> m_gas;
    bool m_is_solved = false;
};

} // namespace Goddard
