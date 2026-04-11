#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "goddard/profile.hpp"
#include "goddard/error.hpp"

namespace Goddard {

double to_radians(double deg) {
    return deg/360.0 * 2 * M_PI;
}
    
// -- NozzleProfile --

NozzleProfile NozzleProfile::generate_conical_nozzle(
    double area_ratio, 
    double r_throat, 
    double angle, 
    size_t n_points)
{
    NozzleProfile profile;
    if (angle <= 0.0 || angle >= 90.0 ) {
        throw std::invalid_argument("Specified angle must be between 0 and 90 degrees.");
    }
    if (area_ratio <= 1.0) {
        throw std::invalid_argument("Area ratio must be greater than 1.");
    }
    if (r_throat <= 0.0) {
        throw std::invalid_argument("Throat radius must be positive.");
    }
    double rad = to_radians(angle);
    double r_exit = std::sqrt(area_ratio*r_throat*r_throat);
    double length = (r_exit-r_throat) / tan(rad);
    for (size_t i = 0; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        double x = length * frac;
        double r = r_throat + (r_exit - r_throat) * frac;
        profile.push_back({x,r});
    }
    return profile;
}

NozzleProfile NozzleProfile::generate_TOP_nozzle(
    double area_ratio, 
    double r_throat, 
    double length_frac, 
    size_t n_points)
{
    if (length_frac <= 0.0 || length_frac >= 1.0 ) {
        throw std::invalid_argument("Specified length fraction must be between 0 and 1.");
    }
    if (area_ratio <= 1.0) {
        throw std::invalid_argument("Area ratio must be greater than 1.");
    }
    if (r_throat <= 0.0) {
        throw std::invalid_argument("Throat radius must be positive.");
    }

    throw NotImplementedError("TOP nozzle profiles are not implemented.");

    double theta_n = to_radians(30.0); //TODO: find equation for theta_n
    double theta_e = to_radians(15.0);
   
    return NozzleProfile::generate_bezier_nozzle(
        area_ratio, theta_n, theta_e, 0.382, r_throat, length_frac, n_points);
}

NozzleProfile NozzleProfile::generate_bezier_nozzle(
    double area_ratio,
    double theta_n, double theta_e, 
    double r_expansion_curve,
    double r_throat, double length_frac,
    size_t n_points)
{
    NozzleProfile profile;
    if (length_frac <= 0.0 || length_frac >= 1.0 ) {
        throw std::invalid_argument("Specified length fraction must be between 0 and 1.");
    }
    if (area_ratio <= 1.0) {
        throw std::invalid_argument("Area ratio must be greater than 1.");
    }
    if (r_throat <= 0.0) {
        throw std::invalid_argument("Throat radius must be positive.");
    }

    double r_exit = std::sqrt(area_ratio*r_throat*r_throat);
    double length = length_frac*(r_exit-r_throat) / tan(15);
    double theta_exp_rad = to_radians(theta_n - 90.0);
    double n_x = r_expansion_curve*r_throat*cos(theta_exp_rad);
    double n_y = r_expansion_curve*r_throat*sin(theta_exp_rad) + r_expansion_curve*r_throat + r_throat;
    
    double theta_n_rad = to_radians(theta_n);
    double theta_e_rad = to_radians(theta_e);
    double slope_n = tan(theta_n_rad);
    double slope_e = tan(theta_e_rad);

    double mid_n_intercept = n_y - slope_n * n_x;
    double mid_e_intercept = r_exit - slope_e * length;
    double mid_x = (mid_e_intercept- mid_n_intercept)/(slope_n-slope_e);
    double mid_y = (slope_n*mid_e_intercept - slope_e*mid_n_intercept)/(slope_n - slope_e);
    size_t i = 0;
    //Generate expansion curve
    for (size_t i = 0; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        double theta = frac*(theta_exp_rad - to_radians(-90.0)) + to_radians(-90.0);
        double x = r_expansion_curve*r_throat*cos(theta);
        double r = r_expansion_curve*r_throat*sin(theta) + r_expansion_curve*r_throat + r_throat;
        profile.push_back({x,r});
    }
    for (size_t i = 0; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        double a = (1.0 - frac)*(1.0 - frac);
        double x = a * n_x + 2 * (1.0 - frac) * frac * mid_x + frac*frac*length;
        double r = a * n_y + 2 * (1.0 - frac) * frac * mid_y + frac*frac*r_exit;
        profile.push_back({x,r});
    }
    return profile;
}


double NozzleProfile::slope_at(double x_query) const {
    return slope_at_idx(find_index(x_query));
}

double NozzleProfile::theta_at(double x_query) const {
    return theta_at_idx(find_index(x_query));
}

double NozzleProfile::radius_at(double x_query) const {
    size_t idx = find_index(x_query); // first index with value larger than x
    double x2 = x[idx];
    double x1 = x[idx-1];
    
    double y2 = y[idx];
    double y1 = y[idx-1];

    double weight = (x_query - x1)/(x2-x1);

    return y1 + weight * (y2-y1);
}

double NozzleProfile::area_at(double x_query) const {
    double radius = radius_at(x_query);
    return M_PI * radius * radius;
}

double NozzleProfile::max_theta() const {
    double max_val = 0.0;
    for (size_t i = 1; i < x.size(); i++) {
        double segment_theta = theta_at_idx(i);
        if (segment_theta > max_val) {
            max_val = segment_theta;
        }
    }
    return max_val;
}

double NozzleProfile::slope_at_idx(size_t idx) const {
    
    double x2 = x[idx];
    double x1 = x[idx-1];
    double y2 = y[idx];
    double y1 = y[idx-1];
    double left_slope = (y2-y1)/(x2-x1);
    // use second order finite difference if possible
    if (idx < x.size()-1) {
        double x3 = x[idx+1];
        double y3 = y[idx+1];
        double right_slope = (y3-y2)/(x3-x2);
        // the dx's can be different, so use a weighted average of slopes
        // to compute final slope. 
        double dx_ratio = (x2-x1)/(x3-x2);
        double total_weight = 1/(1.0 + dx_ratio);
        return (dx_ratio*left_slope+right_slope)*total_weight;
    }
    else {
        return left_slope;
    }
}

double NozzleProfile::theta_at_idx(size_t idx) const {
    return atan(slope_at_idx(idx));
}

double NozzleProfile::x_min() const {
    return x.front();
}

double NozzleProfile::x_max() const {
    return x.back();
}

std::pair<size_t, double> NozzleProfile::radius_max() const {
    double ymax = 0.0;
    size_t idx = 0;
    for (size_t i = 0; i < y.size(); i++) {
        double val = y[i];
        if (val > ymax) {
            ymax = val;
            idx = i;
        }
    }
    return {idx,ymax};
}

std::pair<double, double> NozzleProfile::at(size_t idx) const {
    return {x[idx], y[idx]};
}

void NozzleProfile::push_back(std::pair<double, double>&& coords) {
    x.push_back(coords.first);
    y.push_back(coords.second);
}

size_t NozzleProfile::size() const {
    return x.size();
}

NozzleProfile NozzleProfile::load_profile_csv(const std::string& filename) {
    NozzleProfile profile;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    std::string line;
    // skip header line
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x_val, y_val;
        char comma;
        if (iss >> x_val >> comma >> y_val) {
            profile.x.push_back(x_val);
            profile.y.push_back(y_val);
        }
    }
    return profile;
}

void NozzleProfile::save_profile_csv(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    file << std::setprecision(15);
    file << "x,y\n";
    for (size_t i = 0; i < x.size(); i++) {
        file << x[i] << "," << y[i] << "\n";
    }
}

size_t NozzleProfile::find_index(double x_query) const {
    // linear scan acceptable for small grids
    size_t i = 0;
    while (i < x.size()){
        if (x[i] > x_query) return i;
        i++;
    }
    if (i == x.size()) 
        throw std::runtime_error("Queried x: "+ std::to_string(x_query) + " larger than nozzle profile.");
    // query beyond profile: return last segment
    return x.size() - 1;
}

} // namespace Goddard