#pragma once
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

    // unit processes

    PointResult solve_interior_point(
        const CharacteristicPoint& c_minus_parent,
        const CharacteristicPoint& c_plus_parent);

    // Planar algebraic special case
    PointResult solve_interior_point_planar(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2);

    // Axisymmetric flow
    PointResult solve_interior_point_axisymmetric(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2);

    // iterative path (generalized compatibility equation with arbitrary source term)
    PointResult solve_interior_point_iterative(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2);

    // nullopt: no wall hit within the profile bounds (legitimate outflow).
    // Populated with a nonzero PointResult::error: a numerical failure occurred.
    std::optional<PointResult> solve_wall_point(
        const CharacteristicPoint& interior_parent,
        const CharacteristicPoint& previous_wall_point,
        int wall_point_index);

    PointResult solve_wall_flow(
        const CharacteristicPoint& interior_parent,
        double theta_wall);

    // compute flow + position from previous wall point
    PointResult solve_wall_point_design(
        const CharacteristicPoint& interior_parent,
        const CharacteristicPoint& previous_wall_point,
        double theta_wall);

    // Compute flow at known wall position. nullopt: no wall hit within the
    // profile bounds (legitimate outflow), as opposed to a populated PointResult
    // with a nonzero error, which is a numerical failure.
    std::optional<PointResult> solve_wall_point_analysis(
        const CharacteristicPoint& interior_parent);

    /** The first point is a special case, as it lies on the axis but is assigned a
     * nonzero theta. This is because the calculations are started on the characteristic line
     * along which theta is known.
     *
     * This leads to a small physical inconsistency, but it is necessary to bootstrap the downstream marching.
     */
    PointResult solve_initial_axis_point_centered_exp(
        const CharacteristicPoint& expansion_point);

    /**
     * Compute flow properties at centerline for axisymmetric flow.
     *
     * A special method is needed because the axisymmetric compatibility
     * equations have a singularity on the axis of rotation.
     * */
    PointResult solve_axis_point(
        const CharacteristicPoint& off_axis_parent);

    // Determine where parent characteristic intersects with arbitrary wall profile.
    std::pair<double, double> intersect_characteristic_with_wall(
        const CharacteristicPoint& interior_parent,
        const NozzleProfile& wall);
    
    /**
     * Returns the geometric length of the nozzle measured from the throat.
     */
    double maximum_nozzle_length() const;
    
    double gamma_s_from_mach(double mach) const;
    double gamma_s_from_nu(double nu) const;
    /**
     * Get mach number from characteristic Prandtl-Meyer expansion angle.
     */
    double mach_from_nu(const CharacteristicPoint& point, double mach_guess = 0.0) const;

    /**
     * Get Prandtl-Meyer expansion angle from mach number.
     */
    double nu_from_mach(double mach) const;
    
    /**
     * Iteratively find the mach number of the intersecting node.
     */
    double find_node_mach(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2,
        double source_delta,
        double mach_guess = 0.0
    );

    /* Find where wall intersects with line extending outwards from a characteristic point, for a specified angle */
    std::pair<double,double> find_wall_hit(const CharacteristicPoint& p, const NozzleProfile& wall, double char_angle) const;

    ThermodynamicContext build_thermo_context();

    void update_thermodynamic_state(CharacteristicPoint& point);

    // These three are the critical chokepoint for the perfect-gas Prandtl-Meyer
    // inversion and the Cantera/table lookups: a failure (PM inversion
    // non-convergence, or an out-of-range table query) is reported via the
    // returned MocErrorCode instead of being laundered into point.mach as a
    // sentinel value. Callers MUST check the return value before using point's
    // newly-set fields.
    MocErrorCode update_thermodynamic_state_from_nu(CharacteristicPoint& point, double nu, double mach_guess = 0.0);
    MocErrorCode update_thermodynamic_state_from_mach(CharacteristicPoint& point, double mach);
    MocErrorCode update_thermodynamic_state_from_V(CharacteristicPoint& point, double V);

    /**
     * Source term for axisymmetric flow along C+ characteristic.
     */
    double cplus_source_term(
        const CharacteristicPoint& p1, 
        double new_y) const;
    
    /**
     * Source term for axisymmetric flow along C+ characteristic.
     */
    double cplus_source_term(
        const CharacteristicPoint& p1, 
        const CharacteristicPoint& p3) const;
    /**
     * Source term for axisymmetric flow along C- characteristic.
     */
    double cminus_source_term(
        const CharacteristicPoint& p1, 
        double new_y) const;


    struct LeadingEdgeView {
        ChainMetadata::Family family;
        std::vector<double> y_values;
        std::vector<size_t> chain_indices;
        std::vector<size_t> leading_pt_indices;
    };
    // Generate struct-of-array of leading edge values for a specified family.
    LeadingEdgeView leading_edges(const CharacteristicNet& net, ChainMetadata::Family family) const;

    void update_leading_edges(LeadingEdgeView& view, const CharacteristicNet& net, ChainMetadata::Family family) const;

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
     * deadlock. Shared by the marching kernel and build_front so mesh control measures
     * exactly the steps the kernel will take; they diverged when the criterion had its own
     * copy of this search.
     */
    PairSearch find_pair_partner(
        double plus_y,
        const LeadingEdgeView& minus_edges,
        const std::vector<bool>& claimed,
        const CharacteristicNet& net) const;

    /** One segment of the marching front: the C+/C- parent pair that will form a mesh cell. */
    struct FrontSegment {
        size_t lower_pt_idx;      ///< C+ leading point (below).
        size_t upper_pt_idx;      ///< C- leading point (above).
        size_t lower_plus_chain;  ///< C+ chain led by lower_pt_idx.
        size_t upper_minus_chain; ///< C- chain led by upper_pt_idx.
        double arc;               ///< Straight-line length of the segment.
    };

    /**
     * The marching front, as the ordered list of C+/C- parent pairs the next pass will march.
     *
     * Replays find_pair_partner over the same sorted C+ edges with the same claim
     * bookkeeping the kernel uses, so there is exactly one segment per pairing that will
     * actually happen. Wall-bound and skipped C+ edges contribute no segment.
     */
    std::vector<FrontSegment> build_front(
        const CharacteristicNet& net,
        const LeadingEdgeView& plus_edges,
        const LeadingEdgeView& minus_edges) const;

    /**
     * Hold the marching front's point spacing between bounds, then rebuild both leading-edge
     * views. Runs once at the top of each kernel pass, before pairing.
     *
     * Segments longer than max_front_spacing_factor x the local target are subdivided down
     * to the target; rungs whose neighbours have crowded closer than
     * min_front_spacing_factor x the target are retired. This is the classical marched-line
     * mesh control that keeps a characteristic mesh from being merely whatever the
     * characteristics happen to sample -- see MocOptions::max_front_spacing_factor for why
     * the axisymmetric net cannot manage its own density.
     *
     * Segment length is measured as arc length, not height. The failure this exists to
     * prevent is a front segment stretching downstream while *shrinking* in height as it
     * turns tangent to the C- family, which a height-only measure reads as finer than
     * typical exactly when it is worst.
     *
     * State for an inserted rung is interpolated across the front between the bracketing
     * points: the region being refined contains no characteristic of either family to
     * interpolate *along* (that absence is exactly the defect), so this reconstructs missing
     * data rather than resampling existing data. The error is O(spacing) and refinement
     * holds the spacing near the target, keeping it the same order as the scheme's own
     * truncation error -- an argument only as good as the order-of-accuracy test that
     * checks it.
     *
     * Excluded for DESIGN_MIN_LENGTH, where the C- count *defines* the contour: inserting
     * would change the nozzle being designed rather than the mesh resolving it.
     *
     * @param diag  Filled with this pass's front geometry and insert/retire counts.
     * @return number of rungs inserted this pass
     */
    size_t control_front_spacing(
        CharacteristicNet& net,
        LeadingEdgeView& plus_edges,
        LeadingEdgeView& minus_edges,
        MocPassDiagnostics& diag);

    /**
     * The C+ and C- side lengths (s, t) of the mesh cell a front segment will form.
     *
     * The segment is the cell's diagonal: (above - below) = s * u_plus - t * u_minus. Empty
     * when that system is singular or either side is non-positive, which means the segment
     * is no longer spacelike and no cell exists.
     */
    std::optional<std::pair<double, double>> cell_sides(
        const CharacteristicNet& net, const FrontSegment& segment) const;

/**
     * The points of the previous marching front local to `segment`, ordered by ascending y.
     *
     * Each front point was produced by pairing two points one step upstream, so the
     * predecessors of the bracketing points' chains span the region a new rung's
     * characteristics reach back into. Returns fewer than two points when a bracketing chain
     * has no history -- freshly born at an axis or wall reflection, or itself just inserted --
     * in which case no inverse solve is possible.
     */
    std::vector<size_t> previous_front_points(
        const CharacteristicNet& net, const FrontSegment& segment) const;

    /**
     * Solve for the flow state at a prescribed location from the previous marching front.
     *
     * The inverse of the ordinary interior unit process: instead of two parents fixing where
     * the new point lands, the location is prescribed and the two characteristics through it
     * are traced *back* to the previous front, where their feet are interpolated. The
     * compatibility relations then transport the invariants forward to the point.
     *
     * This is what makes an inserted rung a genuine solution rather than a fabricated one.
     * Interpolating (theta, nu) directly across the current front -- the obvious thing to do --
     * produces a point satisfying *neither* compatibility relation: its K+ and K- are blends
     * taken from characteristics that do not pass through it. Here both invariants are
     * evaluated at the actual feet of the actual characteristics through the point. Some
     * cross-family reconstruction remains unavoidable, because the characteristic being
     * created does not yet exist -- that absence is the defect being repaired -- but it is
     * confined to the feet, and the resulting state is on the solution manifold.
     *
     * @param seed  Initial guess (a cross-front interpolation is adequate); only seeds the iteration.
     * @return the solved point, or an error code when a foot falls outside the previous front
     *         (extrapolation is refused rather than clamped).
     */
    PointResult solve_inverse_interior_point(
        double x_new, double y_new,
        const CharacteristicNet& net,
        const std::vector<size_t>& previous_front,
        const CharacteristicPoint& seed);

        /** Record front spacing/aspect/margin statistics for `front` into `diag`. */
    void record_front_diagnostics(
        const CharacteristicNet& net,
        const std::vector<FrontSegment>& front,
        double target_spacing,
        MocPassDiagnostics& diag) const;

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

    // -- Inverse (reference-plane) marching kernel (Package B) --
    // See instructions/moc_fix/B.md for the algorithm this group of methods implements;
    // the unit processes below are new siblings of solve_interior_point_axisymmetric /
    // solve_axis_point / solve_wall_point_analysis for the INVERSE kernel, not
    // replacements -- the DIRECT kernel and its unit processes above are unchanged.

    /**
     * Resolve MocOptions::march_scheme's AUTO value into a concrete scheme, once per
     * solve(). AUTO selects INVERSE for axisymmetric ANALYSIS and DESIGN_RAO, DIRECT
     * otherwise; MocMode::DESIGN_MIN_LENGTH always resolves to DIRECT (validate_moc_options
     * rejects INVERSE explicitly requested together with DESIGN_MIN_LENGTH before this is
     * reached).
     */
    MocMarchScheme resolve_march_scheme() const;

    /**
     * Build the inverse kernel's first marching front F_0 from the initial data line
     * solve() already constructed (generate_initial_data_line).
     *
     * A Kliegel-Levine transonic line (m_initial_line_family empty) already spans axis to
     * wall and is returned unchanged. A centered-fan line (m_initial_line_family == PLUS)
     * stops short of the wall at its last ray's endpoint P; this extends the front from P
     * straight up to the wall with the fan's uniform post-last-ray state, then corrects the
     * new wall point's nu from the planar C+ compatibility relation with P (see B.md,
     * "Initial front").
     *
     * Not const: extension points are given their full thermodynamic state via
     * update_thermodynamic_state_from_nu (a non-const chokepoint), so FROZEN/EQUILIBRIUM
     * chemistry stay consistent here exactly as everywhere else in the kernel.
     *
     * @throws ConvergenceError if a thermodynamic update for an extension point fails;
     *         caught by solve()'s existing initialization exception boundary.
     */
    std::vector<CharacteristicPoint> build_inverse_initial_front(
        const std::vector<CharacteristicPoint>& data_line);

    /**
     * Append a front's points to `net` (no chain bookkeeping -- the inverse kernel does not
     * use CharacteristicNet::c_chains/membership beyond an empty placeholder) and register
     * CharacteristicNet::fronts, wall_point_indices/wall_x/wall_y and axis_point_indices for
     * it. Used both to seed F_0 and, every pass, to record the front the pass just built.
     *
     * @return the point indices of the seeded front, axis to wall.
     */
    std::vector<size_t> seed_inverse_front(
        CharacteristicNet& net, const std::vector<CharacteristicPoint>& front_points) const;

    /**
     * Inverse reference-plane marching kernel (Package B): starting from the front already
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
     * forward with the same axisymmetric source terms as solve_inverse_interior_point.
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
     * and the axis-limit source term is applied with the same algebra as solve_axis_point's
     * corrector (dy/y_avg = -2 there is an algebraic identity whenever the far point sits at
     * y=0, so it generalizes unchanged to a traced-back foot that is not an actual net
     * parent).
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
     * stops moving -- the same fixed-point structure solve_wall_point_analysis uses, but
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

    // Logging helpers
    template <typename... Args>
    void log_warning(std::string_view fmt, Args&&... args) {
        m_messages.push_back("Warning: " + std::vformat(fmt, std::make_format_args(args...)));
    }
    void log_warning(const std::string& msg);
    
    template <typename... Args>
    void log_info(std::string_view fmt, Args&&... args) {
        m_messages.push_back("Info: " + std::vformat(fmt, std::make_format_args(args...)));
    }
    void log_info(const std::string& msg);

    // Verbose kernel/initialization trace, only active when m_options.log_level ==
    // MocLogLevel::DEBUG. Collected in m_messages like log_warning/log_info (so it
    // surfaces via MocResult::messages and the Python binding with no extra
    // plumbing), and additionally echoed live to stderr immediately as each call
    // happens -- useful for a hang or a solve that never returns (maxiter reached),
    // where messages collected only in the returned MocResult would never be seen.
    template <typename... Args>
    void log_debug(std::string_view fmt, Args&&... args) {
        if (m_options.log_level == MocLogLevel::DEBUG) {
            log_debug(std::vformat(fmt, std::make_format_args(args...)));
        }
    }
    void log_debug(const std::string& msg);

    bool m_is_solved = false;
    // Running totals for the whole march, so the front-size cap applies across passes
    // rather than per pass. Reset by solve().
    size_t m_inserted_characteristics = 0;
    size_t m_retired_characteristics = 0;
    // Per-pass front geometry, copied into MocResult::pass_diagnostics.
    std::vector<MocPassDiagnostics> m_pass_diagnostics;
    // Target front spacing at the throat, and the throat radius in net units. Both are set
    // once from the seeded initial data line, so mesh control is anchored to a physical
    // length rather than to the mesh's own statistics.
    double m_reference_spacing = 0.0;
    double m_throat_radius = 1.0;
    // Set once per march when a front segment is found already non-spacelike, so the
    // warning that mesh control engaged too late is reported but not repeated per pass.
    bool m_warned_non_spacelike = false;
    // Raw wall-boundary-condition residual of the Kliegel-Levine start line, copied from
    // the initializer for MocInitDiagnostics::wall_bc_residual. Reset by solve().
    double m_init_wall_bc_residual = 0.0;
    std::optional<Gas> m_gas;
    std::vector<double> m_theta_schedule;

    // Characteristic family the initial data line lies along, if any. A centered
    // expansion fan is collinear along one characteristic (the canonical C+ case); a
    // transonic start line crosses many and leaves this empty. Set by the initial-data-line
    // generators (which know the strategy they used) and consumed by solve() when seeding
    // the net. This cannot be recovered by inspecting the points: a Riemann invariant is
    // only constant along a characteristic for planar flow, not axisymmetric.
    std::optional<ChainMetadata::Family> m_initial_line_family;

    std::vector<std::string> m_messages;

    PrandtlMeyerTable pm_table;

    double m_L_ref; //reference length for dimensionalization
    double m_P_ref; //reference pressure for dimensionalization
    double m_T_ref; //reference temperature for dimensionalization
    double m_S_ref; //reference entropy
};

} // namespace Goddard