#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <format>
#include <sstream>
#include <Eigen/Dense>
#include "goddard/interpolation.hpp"
#include "goddard/profile.hpp"
#include "goddard/error.hpp"

namespace Goddard {

const Eigen::MatrixXd RAO_PARABOLIC_NOZZLE_THETA_N{
    {25.5617,26.6139,27.5803,28.5369,29.3686,30.1838,30.8414,31.5322,32.1663,32.7779,33.3807,33.9012,34.4370,34.9017,35.3841,35.8830,36.3727,36.8407,37.2330,37.6700,38.0577,38.4905,38.9730,39.4520,39.9037,40.3115}, //60%
    {22.8254,23.8411,24.6748,25.4125,26.1280,26.8243,27.4350,28.0045,28.6146,29.2159,29.7727,30.3021,30.8309,31.2661,31.7118,32.2057,32.7292,33.1384,33.5965,34.1088,34.5292,34.9760,35.4224,35.8232,36.1887,36.5406}, //70%
    {20.9399,21.7800,22.5985,23.2548,24.0176,24.7080,25.3028,25.8744,26.4541,26.9658,27.4755,27.9452,28.4016,28.8714,29.3149,29.7242,30.1721,30.5426,30.9198,31.3664,31.7780,32.2108,32.6261,32.9640,33.3351,33.6040}, //80%
    {19.2380,20.0681,20.6769,21.2634,21.9860,22.6712,23.2586,23.8985,24.4258,24.9870,25.5769,26.0429,26.5634,27.0042,27.5186,28.0027,28.5189,28.9264,29.3999,29.8701,30.3479,30.8301,31.2764,31.7049,32.1411,32.5528}, //90%
    {18.5873,19.1132,19.5830,20.0919,20.6083,21.1377,21.6386,22.2131,22.7335,23.2192,23.7676,24.3059,24.8520,25.3438,25.9327,26.4553,26.9886,27.5296,28.0655,28.6458,29.2084,29.7106,30.2760,30.8161,31.3828,31.9021}  //100%
};

const Eigen::MatrixXd RAO_PARABOLIC_NOZZLE_THETA_E{
    {21.7742,20.6971,19.7317,18.9429,18.2514,17.6146,17.0486,16.5806,16.1623,15.8061,15.4317,15.1560,14.8899,14.6116,14.3607,14.1216,13.8289,13.6752,13.4691,13.2873,13.1424,12.9452,12.7829,12.6543,12.4565,12.2803}, //60%
    {18.0200,17.1219,16.3683,15.7497,15.1154,14.6263,14.1833,13.7300,13.3549,13.0672,12.7323,12.4568,12.2224,12.0127,11.7702,11.5461,11.3283,11.1062,10.8934,10.6548,10.4808,10.2730,10.0371,9.8749,9.7167,9.5273}, //70%
    {14.7164,13.9073,13.2844,12.6088,12.0013,11.5078,10.9991,10.6785,10.3975,10.0820,9.7821,9.5159,9.2774,9.0158,8.7743,8.5555,8.3464,8.1572,7.9658,7.7867,7.6207,7.4152,7.2495,7.1448,7.0238,6.9744}, //80%
    {12.1468,11.3280,10.6192,9.9891,9.5019,9.0344,8.6403,8.3318,8.0808,7.8244,7.6194,7.4299,7.2194,7.0767,6.9314,6.7878,6.6398,6.4987,6.3788,6.3404,6.2423,6.1407,6.1005,6.0785,6.0037,5.9733}, //90%
    {9.9110,9.2236,8.5201,7.9191,7.4036,7.0047,6.6976,6.4119,6.1369,5.8704,5.6635,5.5323,5.4085,5.2457,5.0939,4.9696,4.9031,4.8233,4.7652,4.6989,4.5978,4.5539,4.5371,4.5025,4.4380,4.4216} //100%
};

constexpr double RAO_PARABOLIC_EXP_RATIO_LOG10_STEP = 0.05799;
constexpr double RAO_PARABOLIC_MIN_EXP_RATIO_LOG10 =  0.550258197;
constexpr double RAO_PARABOLIC_MIN_LENGTH_RATIO = 0.60;
constexpr double RAO_PARABOLIC_LENGTH_RATIO_STEP = 0.10;

double to_radians(double deg) {
    return deg/360.0 * 2 * M_PI;
}
    
// -- NozzleProfile --

NozzleProfile NozzleProfile::generate_conical_nozzle(
    double area_ratio, 
    double r_expansion_curve,
    double r_throat, 
    double theta_n, 
    size_t n_points)
{
    NozzleProfile profile;
    if (theta_n <= 0.0 || theta_n >= 90.0 ) {
        throw std::invalid_argument("Specified angle must be between 0 and 90 degrees.");
    }
    if (area_ratio <= 1.0) {
        throw std::invalid_argument("Area ratio must be greater than 1.");
    }
    if (r_throat <= 0.0) {
        throw std::invalid_argument("Throat radius must be positive.");
    }
    double rad = to_radians(theta_n);
    double r_exit = std::sqrt(area_ratio*r_throat*r_throat);
    double length = (r_exit-r_throat) / tan(rad);
    profile.populate_throat_expansion_curve(theta_n, r_expansion_curve, r_throat, n_points);
    for (size_t i = 0; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        double x = length * frac;
        double r = r_throat + (r_exit - r_throat) * frac;
        profile.push_back({x,r});
    }
    return profile;
}

NozzleProfile NozzleProfile::generate_Rao_TOP_nozzle(
    double area_ratio, 
    double r_throat, 
    double length_frac, 
    size_t n_points)
{
    if (length_frac < RAO_PARABOLIC_MIN_LENGTH_RATIO || length_frac > 1.0 ) {
        throw std::invalid_argument("Specified length fraction must be between 0.6 and 1.");
    }
    if (area_ratio < pow(10,RAO_PARABOLIC_MIN_EXP_RATIO_LOG10)) {
        throw std::invalid_argument(std::format("Area ratio must be greater than {}.", pow(10,RAO_PARABOLIC_MIN_EXP_RATIO_LOG10)));
    }
    if (r_throat <= 0.0) {
        throw std::invalid_argument("Throat radius must be positive.");
    }

    // Data is stored in 5x26 matrices: rows are length_frac
    // (5 values, 60%-100%), columns are log10(area_ratio) (26 values). BicubicInterpolator
    // binds its first coordinate to matrix rows and its second to columns.
    BoundsOptions bound_opts{
        .x_low = BoundsHandling::ERROR,
        .x_high = BoundsHandling::ERROR,
        .y_low = BoundsHandling::ERROR,
        .y_high = BoundsHandling::CLAMP,
    };

    BicubicInterpolator theta_n_intp{
        RAO_PARABOLIC_MIN_LENGTH_RATIO,
        RAO_PARABOLIC_LENGTH_RATIO_STEP,
        RAO_PARABOLIC_MIN_EXP_RATIO_LOG10,
        RAO_PARABOLIC_EXP_RATIO_LOG10_STEP,
        RAO_PARABOLIC_NOZZLE_THETA_N,
        bound_opts
    };

    BicubicInterpolator theta_e_intp{
        RAO_PARABOLIC_MIN_LENGTH_RATIO,
        RAO_PARABOLIC_LENGTH_RATIO_STEP,
        RAO_PARABOLIC_MIN_EXP_RATIO_LOG10,
        RAO_PARABOLIC_EXP_RATIO_LOG10_STEP,
        RAO_PARABOLIC_NOZZLE_THETA_E,
        bound_opts
    };

    double theta_n = theta_n_intp.interpolate(length_frac, log10(area_ratio));
    double theta_e = theta_e_intp.interpolate(length_frac, log10(area_ratio));
   
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
    if (length_frac <= 0.0 || length_frac > 1.0 ) {
        throw std::invalid_argument("Specified length fraction must be between 0 and 1.");
    }
    if (area_ratio <= 1.0) {
        throw std::invalid_argument("Area ratio must be greater than 1.");
    }
    if (r_throat <= 0.0) {
        throw std::invalid_argument("Throat radius must be positive.");
    }

    double r_exit = std::sqrt(area_ratio*r_throat*r_throat);
    double length = length_frac*(r_exit-r_throat) / tan(to_radians(15.0));
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
    // size_t i = 0;
    //Generate expansion curve
    profile.populate_throat_expansion_curve(theta_n, r_expansion_curve, r_throat, n_points/2);

    // Start at i = 1: frac = 0 (i = 0) reproduces (n_x, n_y) exactly, which the
    // throat-arc loop above already pushed as its last point. Including it again
    // would create a zero-length segment at the junction (NaN slopes in slope_at_idx).
    for (size_t i = 1; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        double a = (1.0 - frac)*(1.0 - frac);
        double x = a * n_x + 2 * (1.0 - frac) * frac * mid_x + frac*frac*length;
        double r = a * n_y + 2 * (1.0 - frac) * frac * mid_y + frac*frac*r_exit;
        profile.push_back({x,r});
    }
    return profile;
}

NozzleProfile NozzleProfile::generate_throat_expansion_curve(
        double theta_n, double r_expansion_curve, 
        double r_throat, size_t n_points = 50)
{
    NozzleProfile nozzle;
    nozzle.populate_throat_expansion_curve(theta_n, r_expansion_curve, r_throat, n_points);
    return nozzle;
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

double NozzleProfile::length() const {
    return x_max() - x[throat_index];
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
    if (x.size() < 2) {
        throw std::runtime_error("NozzleProfile::find_index requires at least two points.");
    }
    // linear scan acceptable for small grids
    size_t i = 0;
    while (i < x.size()){
        // Return the upper index of the bracketing segment [i-1, i]. Callers index
        // x[i-1], so for a query at or below the first point clamp to the first segment
        // (i == 1) rather than returning 0 and underflowing.
        if (x[i] > x_query) return std::max<size_t>(i, 1);
        i++;
    }
    // No element strictly greater than x_query was found, i.e. x_query >= x.back().
    // A query exactly at the last grid point (e.g. the nozzle exit plane) is valid
    // and should use the last segment; only a query strictly beyond it is an error.
    if (x_query > x.back())
        throw std::runtime_error(
            std::format("Queried x: {} larger than nozzle profile (xmax = {})", x_query,  x_max())
        );
    return x.size() - 1;
}

void NozzleProfile::populate_throat_expansion_curve(
    double theta_n, double r_expansion_curve, 
    double r_throat, size_t n_points) {

    double theta_exp_rad = to_radians(theta_n - 90.0);
    //Generate expansion curve
    for (size_t i = 0; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        
        double theta = frac*(theta_exp_rad - to_radians(-90.0)) + to_radians(-90.0);
        double x = r_expansion_curve*r_throat*cos(theta);
        double r = r_expansion_curve*r_throat*sin(theta) + r_expansion_curve*r_throat + r_throat;
        push_back({x,r});
    }
}

} // namespace Goddard