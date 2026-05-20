#pragma once
#include <cmath>
#include <vector>
#include "goddard/error.hpp"
#include "goddard/characteristics.hpp"
#include "goddard/moc.hpp"
#include "goddard/nozzle.hpp"

namespace Goddard {



class MocInitialization {
    public:

    std::vector<CharacteristicPoint> initialize_sauer(const ThroatGeometry& geometry, const ThroatCondition& throat);
    std::vector<CharacteristicPoint> initialize_kliegel_levine();

    protected:
    inline double delta() const {
        if (m_options.flow_type == MocFlowKind::AXISYMMETRIC) return 1.0;
        else return 0.0;
    }

    inline double sauer_alpha(double gamma, const ThroatGeometry& geometry) const {
        return std::sqrt((1 + delta())/((gamma+1)*geometry.downstream_wall_curvature_radius * geometry.throat_radius));
    }

    // Kliegel-Levine utility functions
    
    double KL_axial_coordinate(double x, const ThroatGeometry& geometry, double gamma) {
        double r = geometry.throat_radius;
        double R = geometry.downstream_wall_curvature_radius/r;
        return std::sqrt(2*R/(gamma+1))*x/r;
    }

    double KL_u1(double r, double z) {
        return 0.5*r*r - 0.25 + z;
    }

    double KL_v1(double r, double z) {
        return 0.25*r*r*r - 0.25*r + r*z;
    }

    double KL_u2(double r, double z, double gamma) {
        double r2 = r*r;
        return (2*gamma+9)/24*r2*r2 - (4*gamma+15)/24*r2 + (10*gamma + 57)/288 + z*(r2 - 5.0/8.0)-(2*gamma-3)/6*z*z;
    }

    double KL_v2(double r, double z, double gamma) {
        throw NotImplementedError("Kliegel-Levine v2 not implemented");
        return 0.0;
    }

    double KL_u3(double r, double z, double gamma) {
        double r2 = r*r;
        throw NotImplementedError("Kliegel-Levine u3 not implemented");
        return 0.0;
    }

    double KL_v3(double r, double z, double gamma) {
        throw NotImplementedError("Kliegel-Levine v3 not implemented");
        return 0.0;
    }

    double KL_xMach(double r, double z, double gamma, double R) {
        double denom = 1.0/(R+1);
        double u1 = KL_u1(r, z);
        double u2 = KL_u2(r, z, gamma);
        double u3 = KL_u3(r, z, gamma);
        return 1 + denom*u1 + denom*denom * (u1 + u2) + denom*denom*denom*(u1 + 2*u2 + u3);
    }

    double KL_axisymmetric_coeff(double gamma, double R) {
        return std::sqrt((gamma + 1)/(2*(R+1)));
    }

    double KL_yMach(double r, double z, double gamma, double R) {
        double denom = 1.0/(R+1);
        double v1 = KL_v1(r, z);
        double v2 = KL_v2(r, z, gamma);
        double v3 = KL_v3(r, z, gamma);
        //TODO: only valid for axisymmetric geometries
        
        return 1 + denom*v1 + denom*denom * (1.5*v1 + v2) + denom*denom*denom*(15.0/8*v1 + 5.0/2*v2 + v3);
    }

    MocOptions m_options;
    PrandtlMeyerTable pm_table;

};

} // namespace Goddard