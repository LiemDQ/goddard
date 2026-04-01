#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "goddard/profile.hpp"

namespace Goddard {


    
// -- NozzleProfile --

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