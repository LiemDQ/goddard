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

    std::vector<CharacteristicPoint> initialize_sauer(const ThroatCondition& throat);
    std::vector<CharacteristicPoint> initialize_kliegel_levine(const ThroatCondition& throat);
    std::vector<CharacteristicPoint> initialize_centered_expansion(const ThroatCondition& throat);

    ThroatGeometry geometry;
    protected:
    inline double delta() const {
        if (m_options.flow_type == MocFlowKind::AXISYMMETRIC) return 1.0;
        else return 0.0;
    }

    inline double sauer_alpha(double gamma) const {
        return std::sqrt((1 + delta())/((gamma+1)*geometry.downstream_wall_curvature_radius * geometry.throat_radius));
    }

    // Kliegel-Levine utility functions
    
    double KL_z_coordinate(double x, double gamma) const;
    double KL_dzdx(double x, double gamma) const;
    double KL_u1(double r, double z) const;
    double KL_u2(double r, double z, double gamma) const;
    double KL_u3(double r, double z, double gamma) const;
    double KL_v1(double r, double z) const;
    double KL_v2(double r, double z, double gamma) const;
    double KL_v3(double r, double z, double gamma) const;
    double KL_dv1dz(double r, double z) const;
    double KL_dv2dz(double r, double z, double gamma) const;
    double KL_dv3dz(double r, double z, double gamma) const;
    double KL_xMach(double x, double y, double gamma, double R) const;
    double KL_axisymmetric_coeff(double gamma, double R) const;
    double KL_yMach(double x, double y, double gamma, double R) const;
    double KL_dyMachdx(double x, double y, double gamma, double R) const;
    double KL_solve_transonic_x(double y, double gamma, double R, double x_guess = 0.0) const;

    MocOptions m_options;
    PrandtlMeyerTable table;
    ThermodynamicContext thermo;

};

} // namespace Goddard