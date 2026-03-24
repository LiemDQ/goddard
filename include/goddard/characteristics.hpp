#pragma once
#include <vector>
#include <utility>

namespace Goddard {

struct CharacteristicPoint {
    double mach;
    double theta; // flow angle
    double nu; // Prandtl-Meyer angle (or generalized PM function)
    double pressure;
    double temperature;
    double gamma_s; // local isentropic gamma

    // Riemann invariants
    double K_plus; 
    double K_minus;

    // geometry
    double x;
    double y;

    // Characteristic slopes
    double mu; // Mach angle = asin(1/M)

    // for chemistry
    std::vector<double> cantera_state; 
};

/**
 * Get the coordinates of a downstream characteristic, 
 * based on the intersection of the characteristics of two upstream parent points.
 */
std::pair<double, double> characteristic_intersection_coordinates(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2,
    double angle1,
    double angle2);

enum class NetTopology {
    TRIANGULAR,
    WAVEFRONT
};

class CharacteristicNet {
    
    public:
    
    int num_c_plus;
    int num_c_minus;
    
    using Wavefront = std::vector<CharacteristicPoint>;
    // access point at (i_plus, j_minus)
    std::vector<Wavefront> wavefronts;

    // TODO: this implementation is only valid for triangular indexing!
    CharacteristicPoint& at(int i, int j){
        return wavefronts.front()[row_offset(j) + i];
    }

    const CharacteristicPoint& at(int i, int j) const {
        return wavefronts.front()[row_offset(j) + i];
    }

    // Wall contour (profile output)
    std::vector<double> wall_x;
    std::vector<double> wall_y;

    NetTopology topology;

private:
    int row_offset(int j) const; //triangular indexing
};

} // namespace Goddard