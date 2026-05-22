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

struct ThroatGeometry {
    double throat_radius;
    double upstream_wall_curvature_radius = 1.5;
    double downstream_wall_curvature_radius = 0.382;
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
    ThroatGeometry geometry;     // throat geometry

    double theta_max;           // max wall angle (radians) for minimum length nozzle design mode
    double exit_mach;           // exit mach number -- for design mode

    // Optional user-supplied theta schedule for the expansion fan.
    // When provided, overrides the auto-generated schedule and num_characteristics
    // is inferred from the schedule size.
    // Each entry is the flow angle (radians) of a C- characteristic from the expansion.
    // Must be monotonically increasing, with the last entry equal to theta_max.
    std::vector<double> theta_schedule;

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

struct MocResult {
    bool converged;
    CharacteristicNet net;
    NozzleProfile profile; //computed in design mode, echoed in analysis

    std::vector<std::string> messages; // warnings and error information

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