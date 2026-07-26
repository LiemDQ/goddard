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
    MocInitialization(NozzleGeometry geom, ThermodynamicContext& thermo, const MocOptions& options);
    
    /**
     * Create an initial dataline using the Sauer method. The data points form a parabola shape
     * near and downstream of the nozzle throat. 
     * 
     * @note This method is not recommended due to poor accuracy and numerical stability when the throat curvature is 
     * large relative to the throat radius. It is mostly used for testing against other methods
     * as the implementation is relatively simple.
     */
    std::vector<CharacteristicPoint> initialize_sauer(const ThroatCondition& throat);
    /**
     * Create an initial dataline using the Kliegel-Levine expansion method.
     * This is the recommended method for most MoC solver modes except in the
     * case of prescribed centerline conditions.
     *
     * The line is not evaluated on the raw sonic (zero-radial-velocity) locus: every
     * station is shifted downstream by `MocOptions::initial_line_axial_shift` (throat
     * radii) first. This keeps the Mach angle away from 90 deg near the axis, which is
     * required for the caller to seed the line with both characteristic families
     * (see CharacteristicNet::add_initial_data_line).
     */
    std::vector<CharacteristicPoint> initialize_kliegel_levine(const ThroatCondition& throat);
    /**
     * Create an initial dataline based on a centered expansion at the nozzle throat.
     * Used in minimum length nozzle design mode.
     * 
     * @warning Unlike most initialization methods, the points provided by this 
     * method all lie on the same C+ characteristic. Special handling is needed 
     * in the characteristic net construction. 
     */
    std::vector<CharacteristicPoint> initialize_centered_expansion(const ThroatCondition& throat);

    NozzleGeometry geometry;

    protected:
    inline double delta() const {
        if (m_options.flow_type == MocFlowKind::AXISYMMETRIC) return 1.0;
        else return 0.0;
    }

    inline double sauer_alpha(double gamma) const {
        // Dimensionless (throat-radius units): uses the curvature ratio R = R_c / r_t.
        return std::sqrt((1 + delta())/((gamma+1)*KL_R()));
    }

    // Kliegel-Levine utility functions
    
    double KL_z_coordinate(double x, double gamma) const;
    double KL_R() const;
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

    ThermodynamicContext m_thermo;
    MocOptions m_options;

};

} // namespace Goddard