#pragma once
#include <vector>
#include <optional>

#include "goddard/prandtlmeyer.hpp"
#include "goddard/gas.hpp"
namespace Goddard {

struct ThermodynamicContext {
    std::optional<Gas> gas;
    const PrandtlMeyerTable& table;
    double T_ref = 273.15;
    double P_ref = 101325.0;
    double gamma_s = 1.4;
};
class CharacteristicPoint {
    public:
    double theta; // flow angle
    double nu; // Prandtl-Meyer angle (or generalized PM function)
    double pressure;
    double temperature;
    double gamma_s; // local isentropic gamma
    double mach;
    double V; //velocity
    // Characteristic slopes
    double mu; // Mach angle = asin(1/M)
    // Riemann invariants
    double K_plus; 
    double K_minus;
    // geometry
    double x;
    double y;
    // for chemistry
    std::vector<double> cantera_state; 

    void update_thermodynamic_state_from_nu(ThermodynamicContext& ctxt, double nu, double mach_guess = 1.0);
    void update_thermodynamic_state_from_mach(ThermodynamicContext& ctxt, double mach);
    void update_thermodynamic_state_from_V(ThermodynamicContext& ctxt, double V);
    void update_Ks();

    private:
    void update_thermodynamic_state(ThermodynamicContext& ctxt);
};

constexpr double average_angle(double angle1, double angle2) {
    return 0.5*(angle1+angle2);
}

constexpr double average_cminus_angle(
    const CharacteristicPoint& p1, const CharacteristicPoint& p2)
{
    return average_angle(p1.theta-p1.mu, p2.theta-p2.mu);
}

constexpr double average_cplus_angle(
    const CharacteristicPoint& p1, const CharacteristicPoint& p2) 
{
    return average_angle(p1.theta + p1.mu, p2.theta + p2.mu);
}

/**
 * Set temperature and pressure based on isentropic relations and thermodynamic
 * state of upstream characteristic node.
 */
void characteristic_isentropic_PT_from_parent(
    CharacteristicPoint& point, const CharacteristicPoint& parent);

/**
 * Get the coordinates of a downstream characteristic,
 * based on the intersection of the characteristics of two upstream parent points
 * and specified characteristic angles.
 * 
 * @param p1 Point along C- characteristic
 * @param p2 Point along C+ characteristic
 * @param angle1 C- characteristic angle
 * @param angle2 C+ characteristic angle
 */
std::pair<double, double> characteristic_intersection_with_angle(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2,
    double angle1,
    double angle2);


} // namespace Goddard