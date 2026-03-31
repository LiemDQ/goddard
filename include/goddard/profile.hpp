#pragma once
#include <vector>
#include <utility>
#include <string>

namespace Goddard {
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

    // Find the maximum wall angle across all segments.
    double max_theta() const;

    std::pair<double, double> at(size_t idx) const;

    void push_back(std::pair<double, double>&& coords);

    size_t size() const;

    static NozzleProfile load_profile_csv(const std::string& filename);

    void save_profile_csv(const std::string& filename);

private: 
    size_t find_index(double x_query) const;
};

} // namespace Goddard