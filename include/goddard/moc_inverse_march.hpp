#pragma once
#include <optional>
#include <vector>
#include "goddard/characteristic_net.hpp"
#include "goddard/moc.hpp"
#include "goddard/moc_context.hpp"
#include "goddard/moc_initialization.hpp"
#include "goddard/moc_unit_processes.hpp"

namespace Goddard {

/**
 * Reference-plane kernel for analysis and Rao design: prescribes each marching front and
 * traces its points' characteristics back to the previous front.
 *
 * Contrast with DirectMarch: that kernel discovers where characteristics next intersect by
 * pairing chain leading edges, so its two families keep whatever density they were seeded
 * with. This kernel instead fixes the front's shape and point distribution every pass, so
 * both families stay resolved at the same density everywhere and the wall is sampled at
 * every step.
 */
class InverseMarch {
public:
    InverseMarch(const MocSolveContext& ctx, CharacteristicNet& net) : m_ctx(ctx), m_net(net) {}

    /**
     * Build the inverse kernel's first marching front F_0 from `start`.
     *
     * A Kliegel-Levine transonic line (start.family empty) already spans axis to wall and is
     * returned unchanged. A centered-fan line (start.family == PLUS) is a C+ characteristic
     * and cannot be marched from directly; the front is built instead on the plane through
     * the fan's first axis point, each point carrying the state of the simple-wave ray
     * through it, with the uniform state beyond the last ray above and the contour's angle
     * at the wall point.
     *
     * Extension points are given their full thermodynamic state via
     * MocThermo::set_state_from_nu, so FROZEN/EQUILIBRIUM chemistry stay consistent here
     * exactly as everywhere else in the kernel.
     *
     * @throws ConvergenceError if a thermodynamic update for an extension point fails;
     *         caught by solve()'s existing initialization exception boundary.
     */
    std::vector<CharacteristicPoint> initial_front(const StartLine& start);

    /**
     * Seed the net with the first marching front.
     * @throws ConvergenceError if `front0` has fewer than two points.
     */
    void seed(const std::vector<CharacteristicPoint>& front0);

    /**
     * Advance the march, repeatedly building the next front by prescribing its geometry and
     * tracing each new point's two characteristics back to the previous front, until a step
     * lands exactly on the exit plane (MocStepLimiter::EXIT).
     *
     * @return the failure that aborted the march (a unit-process error, a degenerate step
     *         length, or the pass safety cap being reached), or std::nullopt once a pass
     *         lands on the exit plane.
     */
    std::optional<MocFailure> run();

    /// Per-pass front geometry, copied into MocResult::pass_diagnostics.
    std::vector<MocPassDiagnostics> pass_diagnostics;

private:
    /** The step chosen for one pass: its axial length and which bound set it. */
    struct Step {
        double dx;
        MocStepLimiter limiter;
    };

    /** Pre-loop check that the seeded front is a valid Cauchy surface for the march: every
     *  segment increasing in y and steeper than the characteristics through its endpoints. */
    std::optional<MocFailure> check_front_spacelike(const std::vector<size_t>& front) const;

    /**
     * "Step 1" of the march: the axial step length for this pass and what bounds it --
     * the domain-of-dependence CFL condition, the wall-foot condition (the top interior
     * point's C- must reach the previous front below its wall point), the wall-turn
     * condition (MocOptions::max_wall_turn_per_step on a curving contour), or the distance
     * remaining to the exit plane.
     */
    Step choose_step(const std::vector<size_t>& front, const std::vector<double>& offset,
                     double max_abs_offset, int pass) const;

    /**
     * "Step 2": the new front's geometry -- the previous front's shape translated
     * downstream by `dx`, stretched in y with the wall radius, and relaxed toward a
     * vertical plane by removing a fraction of every point's axial offset from the wall
     * point (capped so the relaxation moves no point by more than half a step).
     */
    void prescribe_front(const std::vector<size_t>& front, const std::vector<double>& s,
                         const std::vector<double>& offset, double max_abs_offset, double dx,
                         std::vector<double>& x_new, std::vector<double>& y_new) const;

    /** The state at a marching front, plus whether it was found by mirroring the trace
     *  ray across the axis (see trace_foot). */
    struct Foot {
        CharacteristicPoint state;
        bool mirrored;
    };

    /**
     * ONE copy of "trace a ray back to `front` -> fit its Riemann invariants there with the
     * ENO-clamped quadratic -> set the resulting thermodynamic state" (previously
     * triplicated across the three unit processes below).
     *
     * @param allow_axis_mirror When the direct trace misses `front`, retry by mirroring the
     *        ray (and negating its angle) across the axis -- used for the C+ family (and
     *        the wall's own C+ trace); the C- family never needs it (see trace_plain vs.
     *        trace_with_axis_mirror in the .cpp).
     * @param err Set to the thermodynamic chokepoint's result when a foot is found
     *        (MocErrorCode::NONE on success); untouched when no foot is found at all --
     *        callers distinguish that case by the returned std::nullopt.
     */
    std::optional<Foot> trace_foot(double x0, double y0, double angle,
                                   const std::vector<size_t>& front, bool allow_axis_mirror,
                                   MocErrorCode& err) const;

    /**
     * Interior unit process: solve for the flow state at the prescribed point (x, y), whose
     * two characteristics are traced back to `front` and interpolated there, then
     * transported forward with the same axisymmetric source terms as the axisymmetric
     * interior process (moc_unit_processes.hpp).
     *
     * Handles the near-axis case where the C+ foot's trace would cross the axis before
     * meeting `front`: the foot is found by mirroring the ray (and negating its
     * interpolated theta) about the axis, and the C+ source term's sin(theta)/y factor is
     * evaluated at the new point itself rather than averaged with the (now negative-y)
     * mirrored foot, whose average with the new point's small positive y would otherwise
     * pass near zero.
     */
    PointResult solve_interior_point(double x, double y, const std::vector<size_t>& front);

    /**
     * Axis unit process: solve for the flow state at the prescribed axis point (x, 0). theta
     * is pinned to 0; the single C- foot is traced back to `front` and the axis-limit source
     * term is applied via axis_source_correction (moc_unit_processes.hpp), the same algebra
     * solve_axis_point's corrector uses.
     */
    PointResult solve_axis_point(double x, const std::vector<size_t>& front);

    /**
     * Wall unit process: solve for the flow state at the prescribed wall point (x, y). theta
     * is fixed by the contour; the single C+ foot is traced back to `front` and K+ is
     * transported with cplus_source_term(foot, wall_point), iterated (retracing the foot
     * each pass) until nu stops moving.
     */
    PointResult solve_wall_point(double x, double y, const std::vector<size_t>& front);

    /** "Step 7" bookkeeping: reduce the just-built front to one MocPassDiagnostics record. */
    static MocPassDiagnostics measure_front(
        const std::vector<CharacteristicPoint>& front, int pass, Step step);

    const MocSolveContext& m_ctx;
    CharacteristicNet& m_net;
};

} // namespace Goddard
