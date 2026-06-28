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

    MocNozzle(MocOptions options): m_options(options) {}

    MocNozzle(Gas gas, MocOptions options): m_options(options), m_gas(gas) {}

    MocResult solve();
    bool is_solved() const;

    MocOptions m_options;
    

protected:
   
    // generate initial data line (Cantera-backed path)
    std::vector<CharacteristicPoint> generate_initial_data_line(
        const ThroatCondition& throat,
        const ThroatGeometry& geometry,
        size_t num_points);

    // generate initial data line for perfect gas (no Cantera dependency)
    std::vector<CharacteristicPoint> generate_initial_data_line_perfect_gas(
        size_t num_points);


    // propagate kernel region (C+/C- intersections)
    void solve_characteristic_kernel(CharacteristicNet& net);

    // unit processes

    CharacteristicPoint solve_interior_point(
        const CharacteristicPoint& c_minus_parent,
        const CharacteristicPoint& c_plus_parent);
    
    // Planar algebraic special case
    CharacteristicPoint solve_interior_point_planar(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2);
    
    // Axisymmetric flow
    CharacteristicPoint solve_interior_point_axisymmetric(
        const CharacteristicPoint& p1,
        const CharacteristicPoint& p2);
    
    // iterative path (generalized compatibility equation with arbitrary source term)
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
        const CharacteristicPoint& interior_parent);

    /** The first point is a special case, as it lies on the axis but is assigned a 
     * nonzero theta. This is because the calculations are started on the characteristic line 
     * along which theta is known.
     * 
     * This leads to a small physical inconsistency, but it is necessary to bootstrap the downstream marching.
     */
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
    void update_thermodynamic_state_from_nu(CharacteristicPoint& point, double nu, double mach_guess = 0.0);
    void update_thermodynamic_state_from_mach(CharacteristicPoint& point, double mach);
    void update_thermodynamic_state_from_V(CharacteristicPoint& point, double V);

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

    bool m_is_solved;
    std::optional<Gas> m_gas;
    std::vector<double> m_theta_schedule;
    std::vector<std::string> m_messages;

    PrandtlMeyerTable pm_table;

    double m_L_ref; //reference length for dimensionalization
    double m_P_ref; //reference pressure for dimensionalization
    double m_T_ref; //reference temperature for dimensionalization
    double m_S_ref; //reference entropy
};

} // namespace Goddard