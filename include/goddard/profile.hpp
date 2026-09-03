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
 * The contour is stored as two parallel coordinate arrays in a meridional plane: `x` runs
 * downstream along the nozzle axis and `y` is the wall radius (axisymmetric) or half-height
 * (planar). Both carry whatever length unit the caller supplied for the throat radius, so a
 * profile built with the default `r_throat = 1.0` is expressed in throat radii.
 * 
 * Angles returned by the query methods are in **radians**; the static generators take their
 * shape angles in **degrees**, matching the convention used to tabulate them.
 * */
class NozzleProfile {
public:
    /** Axial coordinates of the wall points, ascending (length units). */
    std::vector<double> x;
    /** Wall radius (axisymmetric) or half-height (planar) at each point in `x` (length units). */
    std::vector<double> y;
    /** Index into `x`/`y` of the throat, i.e. the origin `length()` is measured from. */
    size_t throat_index = 0;
    
    /**
     * Interpolate wall slope dy/dx at a given axial position. Only used in analysis mode.
     *
     * @param x_query Axial position (length units); must lie within the profile's domain.
     * @return Dimensionless slope dy/dx.
     */
    double slope_at(double x_query) const;

    /**
     * Interpolate the wall angle at a given axial position.
     *
     * @param x_query Axial position (length units); must lie within the profile's domain.
     * @return Wall angle atan(dy/dx), in radians.
     */
    double theta_at(double x_query) const;
    
    /**
     * Interpolate the wall radius at a given axial position.
     *
     * @param x_query Axial position (length units); must lie within the profile's domain.
     * @return Wall radius (length units).
     */
    double radius_at(double x_query) const;
    
    /** Calculate the cross-sectional area of nozzle at x-position assuming the nozzle is axisymmetric.
     *
     * @param x_query Axial position (length units); must lie within the profile's domain.
     * @return Cross-sectional area (length units squared).
     * */
    double area_at(double x_query) const;

    /**
     * Wall slope dy/dx at a given point index, by finite difference.
     *
     * Second-order (weighted two-sided) wherever a downstream neighbour exists, first-order
     * (backward) at the last point.
     *
     * @param idx Index into `x`/`y`; must be at least 1.
     * @return Dimensionless slope dy/dx.
     */
    double slope_at_idx(size_t idx) const;

    /**
     * Wall angle at a given point index.
     *
     * @param idx Index into `x`/`y`; must be at least 1.
     * @return Wall angle atan(dy/dx), in radians.
     */
    double theta_at_idx(size_t idx) const;

    /** Axial coordinate of the last (most downstream) point, in length units. */
    double x_max() const;

    /** Axial coordinate of the first (most upstream) point, in length units. */
    double x_min() const;
    
    /** Get the length of the nozzle, measured from the throat, in length units. */
    double length() const;

    /**
     * Largest wall radius on the profile.
     *
     * @return Index of the widest point and its radius (length units).
     */
    std::pair<size_t, double> radius_max() const;

    /** Find the maximum wall angle across all segments, in radians. */
    double max_theta() const;

    /**
     * Coordinates of one wall point.
     *
     * @param idx Index into `x`/`y`.
     * @return The (x, y) pair at that index, in length units.
     */
    std::pair<double, double> at(size_t idx) const;

    /** Append a wall point, given as an (x, y) pair in length units. */
    void push_back(std::pair<double, double>&& coords);

    /** Number of points on the profile. */
    size_t size() const;

    /**
     * Read a profile from a two-column (x, y) CSV file.
     *
     * @param filename Path to the CSV file.
     */
    static NozzleProfile load_profile_csv(const std::string& filename);
    /**
     * Generate a conical nozzle profile. 
     * 
     * The contour is a circular expansion arc off the throat, followed by a straight cone
     * tangent to it. When `r_expansion_curve` is zero the arc is omitted and the cone starts
     * at a sharp throat corner.
     * 
     * @param area_ratio Ratio of the exit area to the throat area. 
     * @param r_expansion_curve Radius of curvature of the throat expansion arc, in the same
     * length units as `r_throat`. Zero gives a sharp throat corner.
     * @param r_throat Throat radius, in length units. 
     * @param theta_n Conical expansion angle from the centerline, in degrees.
     * @param n_points Number of points in the profile. 
     * @throws std::invalid_argument if any argument is outside its usable range, or if the
     * expansion arc alone already exceeds the exit radius.
     */
    static NozzleProfile generate_conical_nozzle(
        double area_ratio, double r_expansion_curve,
        double r_throat = 1.0, double theta_n = 15.0, 
        size_t n_points = 50);

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
     * @param area_ratio Ratio of the exit area to the throat area. Must be at least 3.55,
     * the lower bound of Rao's tabulated data.
     * @param r_throat Throat radius, in length units. 
     * @param length_frac Fraction of length of comparable 15-degree conical nozzle. Must be
     * between 0.6 and 1.0.
     * @param n_points Number of points in the profile. 
     * @throws std::invalid_argument if any argument is outside its usable range.
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
     * @param r_throat Throat radius, in length units. 
     * @param length_frac Fraction of length of comparable 15-degree conical nozzle. 
     * @param n_points Number of points in the profile. 
     * @throws std::invalid_argument if any argument is outside its usable range.
     */
    static NozzleProfile generate_bezier_nozzle(
        double area_ratio,
        double theta_n, double theta_e, 
        double r_expansion_curve = 0.382,
        double r_throat = 1.0, double length_frac = 0.8,
        size_t n_points = 50);

    /**
     * Generate just the circular expansion arc that turns the flow off the throat.
     *
     * This is the leading segment shared by the conical and Bézier contours, useful on its own
     * as the wall geometry a method-of-characteristics solve is seeded against.
     *
     * @param theta_n Wall angle the arc turns to, measured from the centerline, in degrees.
     * @param r_expansion_curve Radius of curvature of the arc, as a fraction of `r_throat`.
     * @param r_throat Throat radius, in length units.
     * @param n_points Number of points in the arc.
     */
    static NozzleProfile generate_throat_expansion_curve(
        double theta_n, double r_expansion_curve, 
        double r_throat, size_t n_points = 50);
        

    /**
     * Write the profile to a two-column (x, y) CSV file.
     *
     * @param filename Path to write to.
     */
    void save_profile_csv(const std::string& filename);

private: 
    size_t find_index(double x_query) const;
    
    /**
     * Populate profile with expansion curve from the throat, taking into account
     * curvature radius at the throat.
     * 
     * @warning Appends to the existing profile!
     */
    void populate_throat_expansion_curve(
        double theta_n, double r_expansion_curve, 
        double r_throat, size_t n_points = 50);
};

} // namespace Goddard
