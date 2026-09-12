#pragma once
#include <vector>
#include <string>
#include <string_view>
#include <format>
#include <memory>
#include <utility>
#include <optional>

#include "cantera/core.h" 
#include "goddard/characteristics.hpp"
#include "goddard/characteristic_net.hpp"
#include "goddard/chemistry.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/profile.hpp"
#include "goddard/gas.hpp"

namespace Goddard {

// Forward declaration (defined in nozzle.hpp)
struct ThroatCondition;

/** Dimensionality of the flow field the characteristics are marched through. */
enum class MocFlowKind {
    PLANAR,         ///< Two-dimensional flow; the compatibility relations have no source term.
    AXISYMMETRIC    ///< Flow about a centerline; the source term is singular on the axis.
};

/**
 * What the solver is being asked to do: derive a wall contour, or resolve the flow through
 * a contour that is already known.
 */
enum class MocMode {
    DESIGN_MIN_LENGTH,  ///< Minimum length nozzle with uniform exit flow.
    DESIGN_RAO,         ///< Rao-type length-optimized thrust nozzle.
    DESIGN_CENTERLINE,  ///< Optimal nozzle from prescribed centerline values. @warning Currently unimplemented; throws NotImplementedError.
    ANALYSIS            ///< Wall contour is input; the flow field is solved through it.
};

/** How much diagnostic output a solve collects into MocResult::messages. */
enum class MocLogLevel {
    NORMAL, ///< Default: only log_warning()/log_info() messages collected.
    DEBUG   ///< also collect (and echo live to stderr) a verbose kernel trace --
            // initial-data-line construction, every interior-point pairing, wall
            // hits/outflow terminations, and axis reflections. Verbose; meant for
            // diagnosing a non-converging or misbehaving solve, not routine use.
};

/** What limited the length of the last inverse-march step (MocPassDiagnostics::step_limiter). */
enum class MocStepLimiter {
    NONE,      ///< No step has been taken yet.
    CFL,       ///< Bounded by the domain-of-dependence CFL condition (MocOptions::inverse_cfl).
    WALL_FOOT, ///< Bounded so the top interior point's C- foot stays below the previous wall point.
    WALL_TURN, ///< Bounded by MocOptions::max_wall_turn_per_step on a curving contour.
    EXIT       ///< Bounded by the distance remaining to the exit plane; the front this step builds is the last one.
};

/** Human-readable name for a MocStepLimiter, for log/diagnostic messages. */
std::string_view to_string(MocStepLimiter limiter);

/**
 * Throat and contour geometry the solver is anchored to.
 *
 * Lengths are expressed in whatever unit `throat_radius` is given in; the two curvature radii
 * are dimensionless multiples of it. `length_fraction` and `expansion_ratio` are consumed only
 * by MocMode::DESIGN_RAO, which uses them to generate the Rao contour it then marches.
 */
struct NozzleGeometry {
    /// Throat radius, in length units. All other lengths in the solve are scaled by it.
    double throat_radius;
    /// Wall radius of curvature upstream of the throat, as a multiple of `throat_radius`.
    double upstream_wall_curvature_radius = 1.5;
    /**
     * Wall radius of curvature downstream of the throat, as a multiple of `throat_radius`.
     * Set to a negative value to request centered-fan initialization instead of a
     * Kliegel-Levine transonic start line.
     */
    double downstream_wall_curvature_radius = 0.382;
    /// Rao design only: length as a fraction of a comparable 15-degree conical nozzle.
    double length_fraction = 0.8;
    /// Rao design only: exit-to-throat area ratio of the contour to generate.
    double expansion_ratio = 5.0;
};

/** Which initial data line to build. AUTO keeps today's rule (negative
 *  NozzleGeometry::downstream_wall_curvature_radius or planar flow selects the fan) and
 *  additionally falls back to the fan when the Kliegel-Levine series misses the wall
 *  angle by more than kl_max_wall_angle_error. */
enum class MocStartLine { AUTO, KLIEGEL_LEVINE, CENTERED_FAN };

/**
 * Options for method of characteristics simulations.
 */
struct MocOptions {
    MocFlowKind flow_type = MocFlowKind::PLANAR;
    GasChemistry chemistry = GasChemistry::PERFECT_GAS;
    MocMode mode = MocMode::DESIGN_MIN_LENGTH;

    /// Number of C+ lines seeded from the initial expansion fan.
    int num_characteristics;
    /// Ratio of specific heats. Used only when `chemistry` is GasChemistry::PERFECT_GAS.
    double gamma;
    /// Tolerances for the iterative unit processes.
    SolverOptions solver_options{.abstol = 1e-10, .reltol = 1e-5};
    /// Throat and contour geometry.
    NozzleGeometry geometry;
    /// Set to MocLogLevel::DEBUG for a verbose kernel trace.
    MocLogLevel log_level = MocLogLevel::NORMAL;

    /// Maximum wall angle (radians), for MocMode::DESIGN_MIN_LENGTH.
    double theta_max;
    /// Target exit Mach number, for the design modes. Not read by the solver; the design
    /// target is set through theta_max.
    double exit_mach;

    /**
     * Optional user-supplied theta schedule for the expansion fan.
     *
     * When provided, overrides the auto-generated schedule and num_characteristics
     * is inferred from the schedule size.
     * Each entry is the flow angle (radians) of a C- characteristic from the expansion.
     * Must be monotonically increasing, with the last entry equal to theta_max.
     */
    std::vector<double> theta_schedule;

    /**
     * Downstream shift (dimensionless, in throat radii) applied to every station of the
     * Kliegel-Levine transonic start line, used only for axisymmetric ANALYSIS/DESIGN_RAO
     * (the KL-init path; see NozzleGeometry::downstream_wall_curvature_radius). The raw
     * sonic (zero-radial-velocity) locus that Kliegel-Levine solves for is not usable as a
     * dual-family (C+ and C-) seeding line: near the axis its Mach angle mu approaches
     * 90 deg, and rigidly translating every station downstream by this amount raises the
     * Mach number (lowering mu) everywhere while preserving the locus's own near-axis
     * curvature -- which a per-station constant-Mach lift does not (it was tried and
     * empirically produces invalid, behind-parent seeding; see
     * instructions/moc_convergence_roadmap.md Sec 2 Step 0). Default 0.1 throat radii is a
     * moderate lift validated against the default geometry; a larger shift trades
     * numerical margin for accuracy, since it extrapolates the KL series further from the
     * throat plane it is expanded about.
     */
    double initial_line_axial_shift = 0.1;

    /** Which initial data line to build; see MocStartLine. */
    MocStartLine start_line = MocStartLine::AUTO;

    /**
     * Largest |theta_series - theta_wall| (rad) at the Kliegel-Levine line's wall end for
     * which the line is still corrected to the contour (a multiplicative rescaling of the
     * flow angle, see MocInitialization::initialize_kliegel_levine). Above it AUTO falls
     * back to the centered fan and a forced KLIEGEL_LEVINE fails with INITIALIZATION_FAILED.
     *
     * Default 0.25 rad (14 deg). The default throat (r_arc = 0.382) misses by 0.145 rad, and
     * the corrected line is measurably better there than the fan: monotone wall Mach, 1%
     * start-line mass-flow error against the fan's 46%, and higher coverage on the Rao
     * contour. The threshold exists to catch a series that is not describing the throat at
     * all, not to reject the correction where it works.
     */
    double kl_max_wall_angle_error = 0.25;

    NozzleProfile nozzle_profile; // wall geometry -- for analysis mode

    /**
     * Inverse march only: fraction of the domain-of-dependence step taken each pass, in
     * (0, 1]. The full domain-of-dependence step (cfl = 1) is the largest step for which
     * every new front point's characteristics still trace back to a point strictly inside
     * the previous front; a fraction below 1 leaves margin against the linearization error
     * in that estimate.
     */
    double inverse_cfl = 0.8;

    /**
     * Inverse march only: largest change of wall angle tolerated per step, in radians.
     * Caps the step length on a curving contour (e.g. a throat expansion arc) so the wall
     * point sampling stays fine enough to resolve the turn.
     */
    double max_wall_turn_per_step = 0.0175;

    /**
     * Inverse march only: fraction of the marching front's shape retained per pass, in
     * [0, 1]. Each new front is the previous one translated downstream, with every point's
     * axial offset from the wall point multiplied by this factor, so a Kliegel-Levine start
     * line (axis end downstream of its wall end) relaxes toward a vertical plane as the
     * march proceeds. The relaxation is additionally capped so it moves no point by more
     * than half a step, which keeps the step bounds valid; on a front with large offsets
     * that cap, not this factor, is what limits the relaxation.
     */
    double front_tilt_decay = 0.9;
};

/**
 * Check a MocOptions for self-consistency, throwing std::invalid_argument on any value
 * outside its usable range.
 *
 * Called from the MocNozzle constructors and again at the top of solve(), since
 * MocNozzle::m_options is public and can be changed after construction. This is programmer
 * error, not a numerical failure, so it throws rather than returning a MocFailure --
 * solve()'s "never throws for numerical failures" contract is unaffected.
 *
 * Mode-dependent fields (theta_max, exit_mach) are deliberately not checked: they are
 * legitimately left unset in ANALYSIS mode.
 */
void validate_moc_options(const MocOptions& options);

/**
 * Result of scanning a characteristic net for intersections between characteristics of the
 * same family.
 *
 * Two characteristics of one family meeting is the discrete signature of an oblique shock:
 * the wave family is coalescing rather than diverging. An isentropic MoC solution is only
 * valid where this does not happen, so a converged solve reporting a nonzero count has
 * produced a field the method cannot represent -- regardless of whether every unit process
 * succeeded.
 */
struct MocCrossings {
    size_t count = 0;          ///< Number of same-family characteristic crossings found.
    double first_x = 0.0;      ///< Position of the most upstream crossing.
    double first_y = 0.0;
    CharacteristicFamily first_family = CharacteristicFamily::UNSPECIFIED; ///< Family that crossed there.
};

/**
 * Scan a net for same-family characteristic crossings (see MocCrossings).
 *
 * Compares every pair of same-family characteristic segments whose x-ranges overlap, rather
 * than only chains that are neighbours on the front: the net carries no chain-adjacency
 * relation, and coalescence is a statement about any two characteristics of a family, not
 * only about neighbours. Segments sharing an endpoint are excluded -- chains legitimately
 * meet at reflection points.
 */
MocCrossings find_like_characteristic_crossings(const CharacteristicNet& net);

/**
 * Flow state sampled across the nozzle exit plane, as parallel arrays ordered from the axis
 * outward. This is what compute_thrust_coefficient() integrates over.
 */
struct ExitPlane {
    std::vector<double> y;           ///< Radial (or transverse) station, in length units.
    std::vector<double> mach;        ///< Local Mach number.
    std::vector<double> theta;       ///< Local flow angle, in radians.
    std::vector<double> pressure;    ///< Static pressure; Pa for frozen/equilibrium, normalized by the stagnation pressure for perfect gas.
    std::vector<double> temperature; ///< Static temperature; K for frozen/equilibrium, normalized by the stagnation temperature for perfect gas.
    std::vector<double> gamma_s;     ///< Local isentropic exponent.
    std::vector<double> velocity;    ///< Dimensional V (m/s) for frozen/equilibrium; equal to the Mach number for perfect gas.
};

/**
 * Per-pass record of the marching front's geometry, recorded every kernel pass.
 *
 * The front is where every known axisymmetric failure mode shows up first, and none of it
 * is visible in the finished net: a sampling void contains no cells, so per-cell statistics
 * look healthy while it grows. These are the quantities that expose it.
 */
struct MocPassDiagnostics {
    int pass = 0;                       ///< Kernel pass this record describes.
    size_t front_points = 0;            ///< Number of rungs on the front at the start of the pass.
    double min_spacing = 0.0;           ///< Shortest front-segment arc length.
    double max_spacing = 0.0;           ///< Longest front-segment arc length. The void shows up here.
    double mean_spacing = 0.0;          ///< Mean front-segment arc length.

    // The front's two ends, recorded separately. The axisymmetric failure is that they
    // advance at very different axial rates -- the front shears until a near-axis point
    // sits downstream of every available C- partner -- and that is invisible to every
    // aggregate above, all of which average over the front. See summarize_front_shear.
    double front_axis_x = 0.0;          ///< x of the front's lowest (nearest-axis) point.
    double front_wall_x = 0.0;          ///< x of the front's highest (nearest-wall) point.
    double front_axis_spacing = 0.0;    ///< Arc length of the bottom-most front segment.
    double front_wall_spacing = 0.0;    ///< Arc length of the top-most front segment.

    /// Inverse march only: axial step length taken this pass. 0 for minimum-length design.
    double step_dx = 0.0;
    /// Inverse march only: which limiter bound step_dx this pass; NONE for minimum-length design.
    MocStepLimiter step_limiter = MocStepLimiter::NONE;
};

/**
 * How differently the two ends of the marching front advanced over a solve.
 *
 * Both ends of an axisymmetric front accelerate; the failure is that they do so at very
 * different rates, so the front shears. Comparing the mean axial step over the first few
 * passes with the mean over the last few reduces that history to one number per end.
 *
 * Measured on conical AR=8, N=15: the axis end's step grew x20.6 while the wall end's
 * *shrank* to x0.84 -- a shear ratio near 24. A converging planar solve on the same
 * geometry gives x7.9 and x44.6, a ratio near 5.6. The sign differs because planar's wall
 * end outruns its axis end; what matters is the magnitude of the disparity.
 */
struct MocFrontShear {
    double axis_growth = 0.0;   ///< Late mean axial step at the axis end / early mean.
    double wall_growth = 0.0;   ///< Same at the wall end.
    double shear_ratio = 1.0;   ///< max(axis_growth, wall_growth) / min(...); 1 means the ends kept pace.
    bool valid = false;         ///< False when there were too few passes to form both windows.
};

/**
 * Reduce a pass-diagnostics history to the front-shear summary above.
 *
 * @param pass_diagnostics Per-pass records, in order, from MocResult.
 * @param window Number of passes averaged at each end of the march.
 * @return The summary; `valid` is false when fewer than 2*window passes were recorded.
 */
MocFrontShear summarize_front_shear(
    const std::vector<MocPassDiagnostics>& pass_diagnostics, size_t window = 4);

/**
 * Properties of the initial data line, measured once before the march begins.
 *
 * The initialization decides whether an axisymmetric analysis converges at all: with a
 * centered-fan start line the solver grid-converges, and with the Kliegel-Levine start
 * line it gets *worse* under refinement. Anti-convergence means an error that is fixed in
 * absolute terms while the cell size shrinks, so the quantities here are reported both
 * raw and normalized by the characteristic spacing -- a raw value that stays put while its
 * normalized twin grows with N is the signature.
 *
 * All of these are cheap and are filled in for every mode and both initializers, so a
 * healthy solve's numbers are available as the reference for a sick one's.
 */
struct MocInitDiagnostics {
    size_t points = 0;                  ///< Points on the initial data line.

    /**
     * Distance from the data line's wall end to the prescribed wall contour, in length
     * units. The line's top point is registered as *the* wall point, so a nonzero value
     * means the wall march begins from a point that is not on the wall. Zero in design
     * modes, which have no prescribed contour.
     */
    double wall_gap = 0.0;
    /// wall_gap divided by the characteristic spacing. Growth with N is the anti-convergence signature.
    double wall_gap_over_spacing = 0.0;
    /// |theta at the line's wall end - the contour's own angle there|, in radians.
    double wall_theta_mismatch = 0.0;

    double mach_axis = 0.0;             ///< Mach at the data line's axis end.
    double mach_wall = 0.0;             ///< Mach at its wall end.
    double mach_ratio = 1.0;            ///< mach_wall / mach_axis; 1 for a constant-Mach line.

    double mu_axis = 0.0;               ///< Mach angle at the axis end, radians.
    double mu_wall = 0.0;               ///< Mach angle at the wall end, radians.
    /**
     * cot(mu_axis) / cot(mu_wall). Marching step length scales with cot(mu), so this is
     * how much faster one end of the front is licensed to advance than the other at pass 0.
     */
    double cot_mu_ratio = 1.0;

    /**
     * Largest / smallest gap between the points where the data line's C- characteristics
     * reach the axis, on a straight-ray estimate. The start line is uniform in y but its
     * characteristics need not be uniform in where they land; a large value means the axis
     * is fed by a burst of arrivals and then starved.
     */
    double axis_arrival_grading = 1.0;

    /// Smallest normalized spacelike margin over the data line's own segments, defined as
    /// in MocPassDiagnostics. A value <= 0 means the line is crossed by its own
    /// characteristics, i.e. it is not a valid Cauchy surface for the march.
    double min_spacelike_margin = 0.0;

    /**
     * Relative error of the mass flow integrated across the data line against the 1-D
     * critical mass flow through the throat. This is the one measurement that checks the
     * start line against physics rather than against its own construction.
     *
     * NaN for frozen and equilibrium chemistry, where the density needed for the integrand
     * is not recoverable from the stored point state.
     */
    double mass_flow_error = 0.0;

    /**
     * MocOptions::initial_line_axial_shift divided by sqrt(R), the transonic length scale.
     * The shift is specified as an absolute offset in throat radii, but the axial extent of
     * the transonic region scales as sqrt(r_throat * R_curvature), so a fixed shift means
     * different things at different throat curvatures. Zero when the shift is unused.
     */
    double shift_over_transonic_length = 0.0;

    /**
     * Signed distance, in characteristic spacings, from the data line's wall end to the
     * nearest curvature discontinuity on the contour (negative = upstream of it). Near zero
     * means the wall march is seeded exactly on a corner, where the contour's angle query is
     * a one-sided difference across a slope jump. Large sentinel value for a smooth contour.
     */
    double wall_station_to_tangency = 0.0;

    /**
     * How far the raw (uncorrected) Kliegel-Levine series missed its own wall boundary
     * condition at the data line's wall end: 1 - theta_series/theta_wall, measured before
     * initialize_kliegel_levine's wall-consistency correction is applied. Zero for the
     * centered fan (no series to miss).
     */
    double wall_bc_residual = 0.0;

    /**
     * K+ = theta - nu of the topmost interior point (the data line point just below the
     * wall end), after the wall-consistency correction. This is the quantity diagnosis.md
     * Sec A1 tracks: a large negative value here over-expands the first wall solve.
     */
    double kplus_wall_end = 0.0;

    /** Which start line the solve actually used -- meaningful when MocOptions::start_line
     *  is AUTO and the Kliegel-Levine series was rejected in favor of the fan. */
    MocStartLine start_line_used = MocStartLine::AUTO;
};

/**
 * Classification of the numerical/physical failure modes that can occur while
 * resolving a single unit process (interior/wall/axis point solve) or marching
 * the method-of-characteristics net as a whole.
 */
enum class MocErrorCode {
    NONE,                      ///< No error; the point/solve is valid.
    NEGATIVE_NU,               ///< Prandtl-Meyer angle (or generalized PM function) nu is negative beyond tolerance.
    NEGATIVE_THETA,            ///< Flow angle theta is negative beyond tolerance.
    SUBSONIC_MACH,             ///< Mach number is below 1.0; the supersonic compatibility relations no longer apply.
    NONFINITE_VALUE,           ///< One of the point's numeric fields (x, y, theta, nu, mach, mu) is NaN or infinite.
    PM_INVERSION_FAILED,       ///< The Prandtl-Meyer inversion (nu -> Mach), or an equivalent Mach rootfind, failed to converge.
    TABLE_RANGE_EXCEEDED,      ///< A PrandtlMeyerTable lookup (by nu, Mach, or velocity) fell outside the tabulated range.
    NON_DOWNSTREAM_POINT,      ///< The computed intersection lies at or behind (upstream of) one of its parent points.
    WALL_QUERY_OUT_OF_BOUNDS,  ///< A wall-profile query (e.g. theta_at) fell outside the profile's domain.
    INITIALIZATION_FAILED,     ///< Construction of the initial data line (transonic start line) failed to converge.
    MAX_ITERATIONS_REACHED     ///< The characteristic kernel reached its iteration safety cap before all chains terminated.
};

/**
 * Human-readable name for a MocErrorCode, for log/diagnostic messages.
 */
std::string_view to_string(MocErrorCode code);

/**
 * Describes a numerical or physical failure encountered while solving a MocNozzle
 * problem. When MocResult::converged is false, this identifies what went wrong and
 * where, instead of the failure being silently swallowed or thrown as an exception.
 */
struct MocFailure {
    MocErrorCode code = MocErrorCode::NONE; ///< Failure classification; NONE if no failure occurred.
    std::string message;                    ///< Human-readable description of the failure.
    double x = 0.0;                         ///< x-coordinate of the failing point (or a parent's, if the point itself could not be computed).
    double y = 0.0;                         ///< y-coordinate of the failing point (or a parent's, if the point itself could not be computed).
    int kernel_pass = -1;                   ///< Kernel marching pass on which the failure occurred; -1 for pre-kernel (initialization) failures.
};

/** Everything a MocNozzle::solve() produces: the flow field, the contour, and the diagnostics. */
struct MocResult {
    bool converged; ///< True iff failure.code == MocErrorCode::NONE and the kernel completed without hitting the iteration cap.
    CharacteristicNet net; ///< The solved characteristic mesh. Populated even when the march fails partway.
    NozzleProfile profile; ///< Wall contour: computed in design mode, echoed back in analysis mode.

    std::vector<std::string> messages; ///< Warnings and error information collected during the solve.
    MocFailure failure; ///< Populated when converged is false; MocErrorCode::NONE otherwise.

    double exit_mach;       ///< Mach number at the exit plane.
    double nozzle_length;   ///< Distance from throat to exit plane, in length units.
    double area_ratio;      ///< Exit area divided by throat area.

    /// Per-pass marching-front geometry; one entry per kernel pass, in order. Recorded by the
    /// inverse march only; empty for minimum-length design.
    std::vector<MocPassDiagnostics> pass_diagnostics;

    /// Properties of the initial data line, measured before the march begins.
    MocInitDiagnostics init_diagnostics;

    /**
     * Fraction of the target exit radius that the net's outflow boundary actually reached.
     *
     * The outflow boundary is a ragged staircase of independently terminated chains, so
     * even a healthy march ends up to about one characteristic spacing short of the exit
     * lip; this is O(1/num_characteristics) and shrinks under refinement. Judge a solve by
     * whether the shortfall *shrinks with N*, not by its value on a single grid -- a coarse
     * healthy planar solve and a genuinely truncated axisymmetric one are indistinguishable
     * here (both ~0.92 at N=8). Reported rather than folded into `converged` for exactly
     * that reason.
     */
    double exit_coverage = 0.0;
    /// Same-family characteristic crossings in the finished net; nonzero means coalescence.
    /// Meaningful for chain nets, i.e. minimum-length design; always zero for the
    /// front-based nets of analysis/Rao, which prescribe every front directly and build no
    /// characteristic chains (CharacteristicNet::c_chains) for this scan to see; see
    /// CharacteristicNet::fronts for their own mesh record.
    MocCrossings crossings;

    /// True when exit_coverage is within the (deliberately loose) staircase allowance.
    bool reached_exit_plane = false;

    /**
     * Smallest flow angle anywhere in the net (radians) and where it occurs. In a diverging
     * nozzle a markedly negative value marks a compression converging on the axis, i.e. a
     * forming shock, which an isentropic march can only pass through approximately; the
     * solution downstream of it is not to be trusted to better than the size of the dip.
     * Conical nozzles with a circular-arc throat are known to form such a shock (Darwell &
     * Badham 1963; Migdal & Kosson 1965). Reported for every scheme.
     */
    double min_theta = 0.0;
    double min_theta_x = 0.0;
    double min_theta_y = 0.0;

    ExitPlane exit_plane;
};

/** Thrust performance obtained by integrating over a solved exit plane. */
struct ThrustCoefficient {
    double Cf_vacuum;       ///< Vacuum thrust coefficient.
    double Cf;              ///< Thrust coefficient at the specified ambient pressure.
    double momentum_thrust; ///< Momentum component, normalized by p0 * A_throat.
    double pressure_thrust; ///< Pressure component, normalized by p0 * A_throat.
};

/**
 * Compute thrust coefficient by integrating over the MoC exit plane.
 *
 * Uses the relation rho*V^2 = gamma_s * p * M^2 which holds for all chemistry types.
 * The integrand at each exit plane point is: f = p * (gamma_s * M^2 * cos^2(theta) + 1).
 *
 * @param result   MoC solution result (must have populated exit_plane with gamma_s)
 * @param flow_type  PLANAR or AXISYMMETRIC (determines integration measure)
 * @param ambient_pressure_ratio  p_amb / p0 (0 for vacuum)
 */
ThrustCoefficient compute_thrust_coefficient(
    const MocResult& result,
    MocFlowKind flow_type,
    double ambient_pressure_ratio = 0.0);


} // namespace Goddard