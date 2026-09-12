#pragma once
#include "goddard/moc.hpp"
#include "goddard/moc_context.hpp"
#include "goddard/moc_unit_processes.hpp"
#include "goddard/moc_thermo.hpp"

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
*/
class MocNozzle {
public:

    MocNozzle(MocOptions options): m_options(options) { validate_moc_options(m_options); }

    MocNozzle(Gas gas, MocOptions options): m_options(options), m_gas(gas) { validate_moc_options(m_options); }

    /**
     * March the characteristic net and return the solved flow field.
     *
     * Never throws for numerical failures: a solve that breaks down returns a MocResult with
     * `converged == false` and `failure` identifying what went wrong and where. It does throw
     * std::invalid_argument if `m_options` is self-inconsistent (checked here as well as in
     * the constructor, since `m_options` is public and callers do adjust it in between), and
     * NotImplementedError for MocMode::DESIGN_CENTERLINE.
     */
    MocResult solve();
    /** True once solve() has run to completion on this instance. */
    bool is_solved() const;

    MocOptions m_options;


protected:

    auto setup_nozzle_profile(const NozzleGeometry& geometry) -> NozzleProfile;

    // generate initial data line (both the perfect-gas and Cantera-backed paths;
    // for perfect gas the throat is populated with dummy values)
    std::vector<CharacteristicPoint> generate_initial_data_line(
        const ThroatCondition& throat,
        const NozzleGeometry& geometry,
        size_t num_points);

    // propagate kernel region (C+/C- intersections). Returns the failure that
    // aborted the march (a unit-process error, or the iteration safety cap being
    // reached), or std::nullopt if every chain terminated cleanly.
    std::optional<MocFailure> solve_characteristic_kernel(CharacteristicNet& net);

    // Minimum-length design only. Check PointResult::error for a numerical failure.
    // Resolves the wall angle from the theta schedule, then delegates to
    // solve_wall_point_design (moc_unit_processes.hpp).
    PointResult solve_wall_point(
        const CharacteristicPoint& interior_parent,
        const CharacteristicPoint& previous_wall_point,
        int wall_point_index);

    struct LeadingEdgeView {
        CharacteristicFamily family;
        std::vector<double> y_values;
        std::vector<size_t> chain_indices;
        std::vector<size_t> leading_pt_indices;
    };
    // Generate struct-of-array of leading edge values for a specified family.
    LeadingEdgeView leading_edges(const CharacteristicNet& net, CharacteristicFamily family) const;

    void update_leading_edges(LeadingEdgeView& view, const CharacteristicNet& net, CharacteristicFamily family) const;

    /**
     * Reorder a C+ leading-edge view by descending y (closest to the wall first), tie-broken
     * by ascending x.
     *
     * When several individual C+ chains compete for the same partner in a single kernel
     * pass (e.g. many chains seeded from a Kliegel-Levine transonic line, all wanting the
     * one C- freshly born from a wall reflection), the per-pass search must resolve the
     * geometrically closest competitor first. Iterating in chain-creation order instead
     * lets a far-away C+ (e.g. the transonic line's axis point) claim a partner meant for
     * a much closer chain, producing a physically invalid, oversized jump.
     *
     * Multiple C+ chains can also sit at exactly y=0 simultaneously (e.g. several axis
     * reflections in flight at once, which dual-family KL-init seeding makes common): the
     * y-only comparator leaves their relative order unspecified under `std::sort`
     * (unstable). Ties are broken by ascending x -- the chain whose leading point is
     * further upstream reaches its next partner first, so it must claim before a chain
     * that is already further downstream.
     */
    void sort_plus_edges_by_proximity(LeadingEdgeView& view, const CharacteristicNet& net) const;

    /** Outcome of searching for a C+ leading edge's C- pairing partner. */
    struct PairSearch {
        bool any_cminus_above = false; ///< A C- exists above this C+ at all. False means the C+ is wall-bound.
        std::optional<size_t> chain_idx = std::nullopt; ///< Partner C- chain; empty means SKIP (or wall-bound, when any_cminus_above is false).
        size_t pt_idx = 0;      ///< Partner's leading point index.
        size_t edgevec_idx = 0; ///< Partner's index within the C- leading-edge view (the claim key).
        double dy = 0.0;        ///< Height from the C+ to the partner.
    };

    /**
     * Find the C- leading edge that a C+ at `plus_y` pairs with this pass: the nearest one
     * above it that no closer competitor has already claimed.
     *
     * When the truly-nearest C- above is already claimed, the result is empty rather than
     * the next-nearest unclaimed one -- settling for a farther partner would violate lattice
     * adjacency. A claim implies someone else progressed this pass, so waiting cannot
     * deadlock.
     */
    PairSearch find_pair_partner(
        double plus_y,
        const LeadingEdgeView& minus_edges,
        const std::vector<bool>& claimed,
        const CharacteristicNet& net) const;

    /**
     * Measure the initial data line: its fit to the prescribed wall, its Mach and Mach-angle
     * spread, the grading of its characteristics' axis arrivals, its own spacelike margin,
     * and the mass flow it carries against the 1-D critical value.
     *
     * Called once, after the line is built and the mesh-control anchors are set, for every
     * mode and both initializers. See MocInitDiagnostics for why each quantity is there.
     *
     * @param data_line The initial data line, ordered axis to wall.
     * @return The populated diagnostics record.
     */
    MocInitDiagnostics record_init_diagnostics(
        const std::vector<CharacteristicPoint>& data_line) const;

    /**
     * Relative error of the mass flow carried across the initial data line, against the
     * 1-D critical mass flow through the throat.
     *
     * This is the only check on the start line that appeals to physics rather than to the
     * line's own construction: a start line whose series has been evaluated outside its
     * range, or whose state variables have been misinterpreted, will not carry the right
     * mass however self-consistent it looks.
     *
     * For a centered fan the line is a C+ characteristic that stops one point short of the
     * throat lip (the lip is seeded separately), so it under-counts by O(1/N) by
     * construction; judge that path by whether the deficit shrinks with N.
     *
     * @param data_line The initial data line, ordered axis to wall.
     * @return Signed relative error, or NaN for frozen/equilibrium chemistry, where the
     *         density is not recoverable from the stored point state.
     */
    double start_line_mass_flow_error(
        const std::vector<CharacteristicPoint>& data_line) const;

    // -- Inverse (reference-plane) marching kernel --
    // Used by MocMode::ANALYSIS and MocMode::DESIGN_RAO. The unit processes below are
    // siblings of solve_interior_point/solve_axis_point (moc_unit_processes.hpp) for this
    // kernel, not replacements -- the minimum-length kernel is unchanged.

    /**
     * Build the inverse kernel's first marching front F_0 from the initial data line
     * solve() already constructed (generate_initial_data_line).
     *
     * A Kliegel-Levine transonic line (m_initial_line_family empty) already spans axis to
     * wall and is returned unchanged. A centered-fan line (m_initial_line_family == PLUS) is
     * a C+ characteristic and cannot be marched from directly; the front is built instead on
     * the plane through the fan's first axis point, each point carrying the state of the
     * simple-wave ray through it, with the uniform state beyond the last ray above and the
     * contour's angle at the wall point.
     *
     * Extension points are given their full thermodynamic state via
     * MocThermo::set_state_from_nu, so FROZEN/EQUILIBRIUM chemistry stay consistent here
     * exactly as everywhere else in the kernel.
     *
     * @throws ConvergenceError if a thermodynamic update for an extension point fails;
     *         caught by solve()'s existing initialization exception boundary.
     */
    std::vector<CharacteristicPoint> build_inverse_initial_front(
        const std::vector<CharacteristicPoint>& data_line);

    /**
     * Inverse reference-plane marching kernel: starting from the front already
     * seeded as `net.fronts.back()`, repeatedly builds the next front by prescribing its
     * geometry and tracing each new point's two characteristics back to the previous front,
     * until a step lands exactly on the exit plane (MocStepLimiter::EXIT).
     *
     * Contrast with solve_characteristic_kernel: that kernel discovers where characteristics
     * next intersect by pairing chain leading edges, so its two families keep whatever
     * density they were seeded with. This kernel instead fixes the front's shape and point
     * distribution every pass, so both families stay resolved at the same density
     * everywhere and the wall is sampled at every step.
     *
     * @return the failure that aborted the march (a unit-process error, a degenerate step
     *         length, or the pass safety cap being reached), or std::nullopt once a pass
     *         lands on the exit plane.
     */
    std::optional<MocFailure> solve_inverse_characteristic_kernel(CharacteristicNet& net);

    /**
     * Inverse-march interior unit process: solve for the flow state at the prescribed
     * point (x_new, y_new), whose two characteristics are traced back to `front` (the
     * previous marching front) and interpolated there (Sec. 5 of B.md), then transported
     * forward with the same axisymmetric source terms as the axisymmetric interior process
     * (moc_unit_processes.hpp).
     *
     * Handles the near-axis case where the C+ foot's trace would cross the axis before
     * meeting `front`: the foot is found by mirroring the ray (and negating its
     * interpolated theta) about the axis, and the C+ source term's sin(theta)/y factor is
     * evaluated at the new point itself rather than averaged with the (now negative-y)
     * mirrored foot, whose average with the new point's small positive y would otherwise
     * pass near zero.
     */
    PointResult solve_inverse_march_interior_point(
        double x_new, double y_new,
        const CharacteristicNet& net,
        const std::vector<size_t>& front);

    /**
     * Inverse-march axis unit process: solve for the flow state at the prescribed axis
     * point (x_new, 0). theta is pinned to 0; the single C- foot is traced back to `front`
     * and the axis-limit source term is applied via axis_source_correction (moc_unit_processes.hpp),
     * the same algebra solve_axis_point's corrector uses (dy/y_avg = -2 there is an algebraic
     * identity whenever the far point sits at y=0, so it generalizes unchanged to a
     * traced-back foot that is not an actual net parent).
     */
    PointResult solve_inverse_march_axis_point(
        double x_new,
        const CharacteristicNet& net,
        const std::vector<size_t>& front);

    /**
     * Inverse-march wall unit process: solve for the flow state at the prescribed wall
     * point (x_new, y_new). theta is fixed by the contour (NozzleProfile::theta_at); the
     * single C+ foot is traced back to `front` and K+ is transported with
     * cplus_source_term(foot, wall_point), iterated (retracing the foot each pass) until nu
     * stops moving -- the same fixed-point structure the axisymmetric wall solve uses, but
     * converged on nu directly rather than on the source-term residual, since here the
     * position is prescribed and only the foot (and hence the source term) is unknown.
     */
    PointResult solve_inverse_march_wall_point(
        double x_new, double y_new,
        const CharacteristicNet& net,
        const std::vector<size_t>& front);

    /**
     * Test-only hook (kept protected; set only by a test subclass): when present, solve()
     * seeds the inverse kernel's first front directly from these points instead of running
     * generate_initial_data_line, so a hand-built front with a known exact solution (e.g.
     * uniform flow, or a manufactured source-flow field) can be marched without requiring a
     * throat/KL/fan construction consistent with it. See InverseMarch.UniformFlowStaysUniform
     * and InverseMarch.SourceFlowSecondOrder (test/test_moc_inverse_march.cpp).
     */
    std::optional<std::vector<CharacteristicPoint>> m_inverse_front_override;

    bool m_is_solved = false;
    // Per-pass front geometry, copied into MocResult::pass_diagnostics.
    std::vector<MocPassDiagnostics> m_pass_diagnostics;
    // Target front spacing at the throat, and the throat radius in net units. Both are set
    // once from the seeded initial data line, so mesh control is anchored to a physical
    // length rather than to the mesh's own statistics.
    double m_reference_spacing = 0.0;
    double m_throat_radius = 1.0;
    // Raw wall-boundary-condition residual of the Kliegel-Levine start line, copied from
    // the initializer for MocInitDiagnostics::wall_bc_residual. Reset by solve().
    double m_init_wall_bc_residual = 0.0;
    std::optional<Gas> m_gas;
    std::optional<MocThermo> m_thermo;
    std::vector<double> m_theta_schedule;

    // Characteristic family the initial data line lies along, if any. A centered
    // expansion fan is collinear along one characteristic (the canonical C+ case); a
    // transonic start line crosses many and leaves this empty. Set by the initial-data-line
    // generators (which know the strategy they used) and consumed by solve() when seeding
    // the net. This cannot be recovered by inspecting the points: a Riemann invariant is
    // only constant along a characteristic for planar flow, not axisymmetric.
    std::optional<CharacteristicFamily> m_initial_line_family;

    // Messages collected during solve() (see MocLog); copied into MocResult::messages at
    // the end of solve(). Reset at the top of every solve() call.
    MocLog m_log;

    // Everything a kernel or unit process needs for one solve, besides the points it works
    // on: emplaced in solve() once m_thermo and the resolved wall contour exist. Unit
    // processes are called as free functions (moc_unit_processes.hpp) with *m_context.
    std::optional<MocSolveContext> m_context;
};

} // namespace Goddard
