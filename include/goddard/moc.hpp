#pragma once
#include <vector>
#include <string>
#include <memory>
#include <utility>
#include <optional>

#include "cantera/core.h" 
#include "goddard/characteristics.hpp"
#include "goddard/prandtlmeyer.hpp"
namespace Goddard {

// Forward declaration (defined in nozzle.hpp)
struct ThroatCondition;

enum class MocFlowKind {
    PLANAR,
    AXISYMMETRIC
};

enum class MocChemistry {
    PERFECT_GAS,    // constant gamma, algebraic PRandtl-Meyer
    FROZEN,         // variable gamma, but no composition change
    EQUILIBRIUM     // full chemical equilibrium at each point
};


enum class MocMode {
    DESIGN_MIN_LENGTH,  // Minimum length nozzle with uniform exit flow
    DESIGN_RAO, // Rao-type length-optimized thrust nozzle 
    ANALYSIS    // wall contour is input
};

struct ThroatGeometry {
    double throat_radius;
    double upstream_wall_curvature_radius;
    double downstream_wall_curvature_radius;
};

/** 
 * Geometric representation of a nozzle wall profile. 
 * Used as both input (analysis) or output (design). 
 * 
 * The first coordinate pair (at index 0) should be the coordinates of the throat,
 * which should have an x-coordinate of 0.0 by convention.
 * 
 * */
class NozzleProfile {
public:
    std::vector<double> x;
    std::vector<double> y;
    
    // Interpolate wall slope at a given x-position.
    // Only used in analysis mode.
    double slope_at(double x_query) const;
    
    // Interpolate wall angle at a given x-position.
    double theta_at(double x_query) const;
    
    std::pair<double, double> at(size_t idx) const;

    void push_back(std::pair<double, double>&& coords);

    size_t size() const;

    static NozzleProfile load_profile_csv(const std::string& filename);

    void save_profile_csv(const std::string& filename);

private: 
    size_t find_index(double x_query) const;
};

// For designing nozzle contours where the initial expansion shape is not prescribed
// it is necessary to specific the centerline mach/pressure distribution instead.
// TODO for a future feature.
// enum class ExpansionSource {
//     WALL_CONTOUR,           // prescribe wall shape (circular arc, spline, etc.)
//     CENTERLINE_DISTRIBUTION // prescribe Mach or pressure along axis
// };

// struct ExpansionSpec {
//     ExpansionSource source;
//     // For WALL_CONTOUR
//     NozzleProfile wall_profile; 

//     // for CENTERLINE_DISTRIBUTION
//     std::vector<double> x_stations;
//     std::vector<double> mach_values;
// };

/**
 * Options for method of characteristics simulations.
 */
struct MocOptions {
    MocFlowKind flow_type = MocFlowKind::PLANAR;
    MocChemistry chemistry = MocChemistry::PERFECT_GAS;
    MocMode mode = MocMode::DESIGN_MIN_LENGTH;

    int num_characteristics;    // number of C+ lines from initial expansion fan
    double gamma;               // used only for PERFECT_GAS
    double reltol = 1e-5;
    double abstol = 1e-10;
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

struct MocResult {
    bool converged;
    CharacteristicNet net;
    NozzleProfile profile; //computed in design mode, echoed in analysis

    std::vector<std::string> messages; // warnings and error information 

    double exit_mach;
    double nozzle_length;   // throat to exit plane
    double area_ratio;      // exit area / throat area
};

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

    MocNozzle(MocOptions options): m_options(options) {}

    MocNozzle(std::shared_ptr<Cantera::Solution> gas, MocOptions options): m_options(options), m_gas(gas) {}

    MocResult solve();

    MocOptions m_options;

private:
    // generate initial data line (Cantera-backed path)
    std::vector<CharacteristicPoint> generate_initial_data_line(
        const ThroatCondition& throat,
        const ThroatGeometry& geometry,
        size_t num_points);

    // generate initial data line for perfect gas (no Cantera dependency)
    std::vector<CharacteristicPoint> generate_initial_data_line_perfect_gas(
        size_t num_points);
    // propagate kernel region (C+/C- intersections)
    void solve_kernel_region(CharacteristicNet& net);
    // compute wall points (straightening section)
    void solve_wall_region(CharacteristicNet& net);

    // unit processes

    CharacteristicPoint solve_interior_point(
        const CharacteristicPoint& c_minus_parent,
        const CharacteristicPoint& c_plus_parent);
    
    // algebraic special case (planar + perfect gas)
    CharacteristicPoint solve_interior_point_algebraic(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2);
    
    // iterative path (axisymmetric and/or variable gamma)
    CharacteristicPoint solve_interior_point_iterative(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2);

    CharacteristicPoint solve_wall_point(
        const CharacteristicPoint& interior_parent,
        const CharacteristicPoint& previous_wall_point,
        int wall_point_index);

    CharacteristicPoint solve_wall_flow(
        const CharacteristicPoint& interior_parent,
        double theta_wall);
    
    // compute flow + position from previous wall point
    CharacteristicPoint solve_wall_point_design(
        const CharacteristicPoint& interior_parent,
        const CharacteristicPoint& previous_wall_point,
        double theta_wall);
    
    // Compute flow at known wall position.
    CharacteristicPoint solve_wall_point_analysis(
        const CharacteristicPoint& interior_parent,
        double theta_wall, double x_wall, double y_wall);

    // The first point is a special case, as it lies on the axis 
    // but is assigned a nonzero theta. This is because the calculations 
    // are started on the characteristic line along which theta is known.
    // This leads to a small physical inconsistency, but it is necessary to 
    // bootstrap the downstream marching.
    CharacteristicPoint solve_initial_axis_point(
        const CharacteristicPoint& expansion_point);

    /**
     * Compute flow properties at centerline for axisymmetric flow.
     *  
     * A special method is needed because the axisymmetric compatibility 
     * equations have a singularity on the axis of rotation.
     * */ 
    CharacteristicPoint solve_axis_point(
        const CharacteristicPoint& off_axis_parent);

    // Determine where parent characteristic intersects with arbitrary wall profile.
    std::pair<double, double> intersect_characteristic_with_wall(
        const CharacteristicPoint& interior_parent,
        const NozzleProfile& wall);
    
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
    ) const;

    void update_thermodynamic_state(CharacteristicPoint& point);
    void update_thermodynamic_state_from_nu(CharacteristicPoint& point, double nu, double mach_guess = 0.0);
    void update_thermodynamic_state_from_mach(CharacteristicPoint& point, double mach);

    double cplus_source_term(
        const CharacteristicPoint& p1, 
        double new_y) const;

    double cminus_source_term(
        const CharacteristicPoint& p1, 
        double new_y) const;


    std::shared_ptr<Cantera::Solution> m_gas;
    std::vector<double> m_theta_schedule;
    std::vector<std::string> m_messages;

    PrandtlMeyerTable pm_table;

    double m_L_ref; //reference length for dimensionalization
    double m_P_ref; //reference pressure for dimensionalization
    double m_T_ref; //reference temperature for dimensionalization
    double m_S_ref; //reference entropy
};

} // namespace Goddard