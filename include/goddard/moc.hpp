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

enum class MocFlowKind {
    PLANAR,
    AXISYMMETRIC
};

enum class MocMode {
    DESIGN_MIN_LENGTH,  // Minimum length nozzle with uniform exit flow
    DESIGN_RAO, // Rao-type length-optimized thrust nozzle
    DESIGN_CENTERLINE, // Designs an optimal nozzle based on prescribed centerline values. @warning Currently unimplemented!
    ANALYSIS    // wall contour is input
};

enum class MocLogLevel {
    NORMAL, // default: only log_warning()/log_info() messages collected
    DEBUG   // also collect (and echo live to stderr) a verbose kernel trace --
            // initial-data-line construction, every interior-point pairing, wall
            // hits/outflow terminations, and axis reflections. Verbose; meant for
            // diagnosing a non-converging or misbehaving solve, not routine use.
};

struct NozzleGeometry {
    double throat_radius;
    double upstream_wall_curvature_radius = 1.5;
    double downstream_wall_curvature_radius = 0.382; // set to negative value for centered-fan initialization
    double length_fraction = 0.8;
    double expansion_ratio = 5.0;
};


/**
 * Options for method of characteristics simulations.
 */
struct MocOptions {
    MocFlowKind flow_type = MocFlowKind::PLANAR;
    GasChemistry chemistry = GasChemistry::PERFECT_GAS;
    MocMode mode = MocMode::DESIGN_MIN_LENGTH;

    int num_characteristics;    // number of C+ lines from initial expansion fan
    double gamma;               // used only for PERFECT_GAS
    SolverOptions solver_options{.abstol = 1e-10, .reltol = 1e-5};
    NozzleGeometry geometry;     // throat geometry
    MocLogLevel log_level = MocLogLevel::NORMAL; // set to DEBUG for verbose kernel tracing

    double theta_max;           // max wall angle (radians) for minimum length nozzle design mode
    double exit_mach;           // exit mach number -- for design mode

    // Optional user-supplied theta schedule for the expansion fan.
    // When provided, overrides the auto-generated schedule and num_characteristics
    // is inferred from the schedule size.
    // Each entry is the flow angle (radians) of a C- characteristic from the expansion.
    // Must be monotonically increasing, with the last entry equal to theta_max.
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

    NozzleProfile nozzle_profile; // wall geometry -- for analysis mode
};

struct ExitPlane {
    std::vector<double> y;
    std::vector<double> mach;
    std::vector<double> theta;
    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> gamma_s;
    std::vector<double> velocity;  // dimensional V (m/s) for frozen/equil, =Mach for perfect gas
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

struct MocResult {
    bool converged; ///< True iff failure.code == MocErrorCode::NONE and the kernel completed without hitting the iteration cap.
    CharacteristicNet net;
    NozzleProfile profile; //computed in design mode, echoed in analysis

    std::vector<std::string> messages; // warnings and error information
    MocFailure failure; ///< Populated when converged is false; MocErrorCode::NONE otherwise.

    double exit_mach;
    double nozzle_length;   // throat to exit plane
    double area_ratio;      // exit area / throat area

    ExitPlane exit_plane;
};

struct ThrustCoefficient {
    double Cf_vacuum;       // vacuum thrust coefficient
    double Cf;              // thrust coefficient at specified ambient pressure
    double momentum_thrust; // momentum component (normalized by p0 * A_throat)
    double pressure_thrust; // pressure component (normalized by p0 * A_throat)
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