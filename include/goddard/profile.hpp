#pragma once
#include <vector>
#include <utility>
#include <string>
#include <array>
#include <Eigen/Dense>

namespace Goddard {
/** 
 * Geometric representation of a nozzle wall profile. 
 * Used as both input (analysis) or output (design). 
 * 
 * 
 * */
class NozzleProfile {
public:
    std::vector<double> x;
    std::vector<double> y;
    size_t throat_index = 0;
    
    // Interpolate wall slope at a given x-position.
    // Only used in analysis mode.
    double slope_at(double x_query) const;

    
    // Interpolate wall angle at a given x-position.
    double theta_at(double x_query) const;
    
    // Interpolate profile height at a given x-position.
    double radius_at(double x_query) const;
    
    /** Calculate the cross-sectional area of nozzle at x-position assuming the nozzle is axisymmetric.
     * */
    double area_at(double x_query) const;
    // Get slope at a given index.
    double slope_at_idx(size_t idx) const;
    // Get wall angle at a given index.
    double theta_at_idx(size_t idx) const;

    double x_max() const;

    double x_min() const;
    
    // Get the length of the nozzle, measured from the throat.
    double length() const;

    std::pair<size_t, double> radius_max() const;

    // Find the maximum wall angle across all segments.
    double max_theta() const;

    std::pair<double, double> at(size_t idx) const;

    void push_back(std::pair<double, double>&& coords);

    size_t size() const;

    static NozzleProfile load_profile_csv(const std::string& filename);
    /**
     * Generate a conical nozzle profile. 
     * 
     * @param area_ratio Ratio of the exit area to the throat area. 
     * @param r_throat Throat radius. 
     * @param angle Conical expansion angle from the centerline, in degrees.
     * @param n_points Number of points in the profile. 
     */
    static NozzleProfile generate_conical_nozzle(
        double area_ratio, double r_throat = 1.0,  
        double angle = 15.0, size_t n_points = 50);

    /**
     * Generate a thrust-optimized parabolic (TOP) nozzle based on the approximations 
     * by Rao.
     * 
     * @note The Rao approximation generates parameters for a Bézier curve construction. For more
     * control over the geometry, use `generate_bezier_nozzle` instead.
     * 
     * ## References
     * 
     * 1. G. V. R. Rao, “Exhaust Nozzle Contour for Optimum Thrust,” Journal of Jet Propulsion, vol. 28, no. 6, pp. 377–382, Jun. 1958, doi: 10.2514/8.7324.
     * 
     * 2. G. V. R. Rao, “Approximation of optimum thrust nozzle contour,” Ars Journal, vol. 30, no. 6, p. 561, 1960.
     * 
     * @param area_ratio Ratio of the exit area to the throat area. 
     * @param r_throat Throat radius. 
     * @param length_frac Fraction of length of comparable 15-degree conical nozzle. 
     * @param n_points Number of points in the profile. 
     */
    static NozzleProfile generate_Rao_TOP_nozzle(
        double area_ratio, double r_throat = 1.0,
        double length_frac = 0.8, size_t n_points = 50);
    /**
     * Generate a parabolic nozzle parametrized by Bezier curves. 
     * 
     * @param area_ratio Ratio of the exit area to the throat area.
     * @param theta_n Maximum expansion angle from centerline, in degrees.
     * @param theta_e Expansion angle from centerline at the nozzle exit, in degrees.
     * @param r_expansion_curve Radius of curvature of the expansion region 
     * as a fraction of the throat radius.
     * @param r_throat Throat radius. 
     * @param length_frac Fraction of length of comparable 15-degree conical nozzle. 
     * @param n_points Number of points in the profile. 
     */
    static NozzleProfile generate_bezier_nozzle(
        double area_ratio,
        double theta_n, double theta_e, 
        double r_expansion_curve = 0.382,
        double r_throat = 1.0, double length_frac = 0.8,
        size_t n_points = 50);

    void save_profile_csv(const std::string& filename);

private: 
    size_t find_index(double x_query) const;
};

} // namespace Goddard
