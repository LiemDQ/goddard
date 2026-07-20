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
     * Reorder a C+ leading-edge view by descending y (closest to the wall first).
     *
     * When several individual C+ chains compete for the same partner in a single kernel
     * pass (e.g. many chains seeded from a Kliegel-Levine transonic line, all wanting the
     * one C- freshly born from a wall reflection), the per-pass search must resolve the
     * geometrically closest competitor first. Iterating in chain-creation order instead
     * lets a far-away C+ (e.g. the transonic line's axis point) claim a partner meant for
     * a much closer chain, producing a physically invalid, oversized jump.
     */
    void sort_plus_edges_by_proximity(LeadingEdgeView& view) const;

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

    bool m_is_solved;
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