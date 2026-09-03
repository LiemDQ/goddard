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
    /// Target exit Mach number, for the design modes.
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

    /**
     * Shifts the Kliegel-Levine start line's points toward the axis (negative) or toward the
     * wall (positive), between uniform spacing in y (0) and uniform spacing in where their
     * C- characteristics reach the axis (+/-1).
     *
     * Those two distributions differ by a factor of ~7, because the near-axis region is a
     * double zero: y -> 0 and cot(mu) -> 0 together, the start line being near-sonic on the
     * axis. Uniform in y therefore lands the C- from the lower third of the line within a few
     * percent of a throat radius of each other. That grading is set by the flow rather than
     * the mesh, so it is independent of num_characteristics (measured 6.2 to 7.2 for N = 8 to
     * 61) -- which is why refining the grid never cured the axisymmetric breakdown.
     *
     * Both directions have been swept, and neither fixes it. Equalizing the arrivals (+1)
     * empties the near-axis mesh and is strictly worse at every N. Clustering toward the axis
     * (negative) helps, but only slightly and only in exit coverage: measured 0.697 -> 0.721
     * at N=61 for AR=4, against a target of 1.0. The apparent cure at N=15 and clustering
     * -0.2 (coverage 0.730 -> 0.981, area ratio 2.13 -> 3.85) is a coarse-grid threshold
     * artifact -- the march happens to clear a barrier that N=31 and N=61 do not -- and must
     * not be read as a fix. Default 0 because the benefit does not survive refinement.
     */
    double initial_line_clustering = 0.0;

    /** Which initial data line to build; see MocStartLine. */
    MocStartLine start_line = MocStartLine::AUTO;

    /** Largest |theta_series - theta_wall| (rad) at the KL line's wall end before the
     *  series is judged untrustworthy there. Default 0.035 (2 deg). */
    double kl_max_wall_angle_error = 0.035;

    /**
     * Upper bound on marching-front point spacing, as a multiple of the local target
     * spacing (see front_spacing_growth). A front segment longer than this is subdivided
     * into as many pieces as it takes to bring it back to the target.
     *
     * The kernel's front is a ladder whose rungs each own one C+ and one C- chain, and a
     * unit process always places its result strictly between its two parents -- so the net
     * cannot add characteristics on its own. In axisymmetric flow the wall-born C- family
     * descends far more slowly than the near-sonic family seeded below it, opening a
     * sampling void at the throat-arc expansion fan that stretches every pass and never
     * heals. It eventually elongates a mesh cell far enough that a front segment becomes
     * tangent to a characteristic and the unit process degenerates. Bounding the spacing
     * is the classical marched-line cure and is applied unconditionally; there is no
     * disable switch, because convergence above area ratio ~4 depends on it.
     *
     * Must be > 1. A factor at or below 1 triggers on essentially every segment and the
     * refinement runs away.
     */
    double max_front_spacing_factor = 1.5;

    /**
     * Lower bound on marching-front point spacing, as a multiple of the local target
     * spacing. A rung whose neighbours have crowded closer than this is retired.
     *
     * Compression regions (notably a Rao contour's turn-back) drive same-family
     * characteristics together; without deletion the front's point count only ever grows
     * and the cells become ill-conditioned. Must be well below max_front_spacing_factor,
     * or refinement and coarsening thrash against each other.
     */
    double min_front_spacing_factor = 0.35;

    /**
     * How much the target spacing is allowed to grow with the local nozzle radius, from
     * 0 (constant spacing everywhere) to 1 (spacing proportional to radius).
     *
     * A nozzle's front legitimately coarsens as it expands: the front's arc length grows
     * roughly as the local radius while the rung count is conserved, and the flow
     * downstream is smoother, so the same truncation error tolerates a longer step. Tying
     * the target to the radius permits exactly that much coarsening and no more. The
     * alternative of measuring against the front's own median spacing was tried and fails:
     * a threshold computed from the spacings it judges is satisfied by any uniformly
     * coarsening mesh, so it never sees a void whose width is set by the flow rather than
     * by the mesh. At the default 1.0 the front holds roughly num_characteristics rungs
     * from throat to exit.
     */
    double front_spacing_growth = 1.0;

    /**
     * Largest tolerated ratio of the C- to the C+ side of a mesh cell, before the rung that
     * forms it is retired.
     *
     * A front segment is a cell diagonal, so bounding its length bounds the long side but
     * leaves the short side free to collapse: a C+ whose partner has marched on while it
     * stalled ends up sitting almost exactly on that partner's C- characteristic, and the
     * cell degenerates while its diagonal still looks perfectly healthy. That is a
     * same-family convergence, and the classical remedy is to delete one of the two
     * redundant points -- they carry nearly the same C- information. Measured healthy fronts
     * wander between 1 and 5; the runaway this bounds climbs geometrically past 200.
     */
    double max_cell_aspect_ratio = 6.0;

    /** Safety cap on marching-front size, so a bad criterion cannot run away. */
    size_t max_front_points = 0; ///< 0 selects the default, 4 * num_characteristics.

    NozzleProfile nozzle_profile; // wall geometry -- for analysis mode
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
    ChainMetadata::Family first_family = ChainMetadata::Family::UNSPECIFIED; ///< Family that crossed there.
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
    double target_spacing = 0.0;        ///< Target spacing this pass (see MocOptions::front_spacing_growth).
    double max_cell_aspect = 0.0;       ///< Largest C-/C+ cell side ratio t/s; diverges as the front turns tangent to the C- family.
    double min_spacelike_margin = 0.0;  ///< Smallest normalized spacelike margin over front segments; reaching 0 is the NON_DOWNSTREAM_POINT failure.
    size_t inserted = 0;                ///< Rungs inserted this pass.
    size_t retired = 0;                 ///< Rungs retired this pass.

    // The front's two ends, recorded separately. The axisymmetric failure is that they
    // advance at very different axial rates -- the front shears until a near-axis point
    // sits downstream of every available C- partner -- and that is invisible to every
    // aggregate above, all of which average over the front. See summarize_front_shear.
    double front_axis_x = 0.0;          ///< x of the front's lowest (nearest-axis) point.
    double front_wall_x = 0.0;          ///< x of the front's highest (nearest-wall) point.
    double front_axis_spacing = 0.0;    ///< Arc length of the bottom-most front segment.
    double front_wall_spacing = 0.0;    ///< Arc length of the top-most front segment.
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
     * centered fan (no series to miss) and when the chemistry is not PERFECT_GAS (the
     * series needs a scalar gamma that only the perfect-gas path carries independently
     * of the throat solve).
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

    /// Rungs added to the marching front by mesh control (MocOptions::max_front_spacing_factor).
    size_t inserted_characteristics = 0;
    /// Rungs retired from the marching front by mesh control (MocOptions::min_front_spacing_factor).
    size_t retired_characteristics = 0;

    /// Per-pass marching-front geometry; one entry per kernel pass, in order.
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
    MocCrossings crossings;

    /// True when exit_coverage is within the (deliberately loose) staircase allowance.
    bool reached_exit_plane = false;

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