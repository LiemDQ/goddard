#include <stdexcept>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "goddard/profile.hpp"

namespace Goddard {


    
// -- NozzleProfile --

double NozzleProfile::slope_at(double x_query) const {
    size_t idx = find_index(x_query); // first index with value larger than x
    double x2 = x[idx];
    double x1 = x[idx-1];
    
    double y2 = y[idx];
    double y1 = y[idx-1];

    return (y2-y1)/(x2-x1);
}

double NozzleProfile::theta_at(double x_query) const {

    return atan(slope_at(x_query));
}

double NozzleProfile::max_theta() const {
    double max_val = 0.0;
    for (size_t i = 1; i < x.size(); i++) {
        double segment_theta = std::atan2(y[i] - y[i-1], x[i] - x[i-1]);
        if (segment_theta > max_val) {
            max_val = segment_theta;
        }
    }
    return max_val;
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
    for (size_t i = 0; i < x.size(); i++){
        if (x[i] > x_query) return i;
    }
    // query beyond profile: return last segment
    return x.size() - 1;
}

} // namespace Goddard