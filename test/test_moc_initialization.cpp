#include "goddard/moc_initialization.hpp"
#include <cmath>
#include <vector>
#include "gtest/gtest.h"

using namespace Goddard;

// ============================================================
// Test-only accessor: exposes MocInitialization's protected Kliegel-Levine
// kernel functions so they can be checked directly against the published
// series (Kliegel & Levine, "Transonic Flow in Small Throat Radius of
// Curvature Nozzles", AIAA Journal Vol. 7, No. 7, July 1969, pp. 1375-1378).
// The public API alone can't isolate these: initialize_kliegel_levine()
// only ever evaluates the series on the (solved) transonic line, never at
// the throat plane z=0 where the paper's Eqs. (9)-(14) are stated.
// ============================================================
class MocInitializationTestAccess : public MocInitialization {
public:
    using MocInitialization::MocInitialization;
    using MocInitialization::KL_xMach;
    using MocInitialization::KL_yMach;
    using MocInitialization::KL_solve_transonic_x;
    using MocInitialization::sauer_alpha;
    using MocInitialization::delta;
};

namespace {

ThermodynamicContext make_perfect_gas_context(double gamma) {
    static PrandtlMeyerTable dummy_table; // unused: PERFECT_GAS never touches it
    return ThermodynamicContext{
        .gas = std::nullopt,
        .table = dummy_table,
        .T_ref = 1.0,
        .P_ref = 1.0,
        .gamma_s = gamma,
    };
}

ThroatCondition make_perfect_gas_throat(double gamma) {
    return ThroatCondition{
        .converged = true,
        .speed_of_sound = 1.0,
        .H_stagnation = 1.0,
        .P_inlet = 1.0,
        .S_inlet = 1.0,
        .gamma_s = gamma,
        .dlV_dlP_T = -1.0,
        .dlV_dlT_P = 1.0,
        .state = {},
    };
}

MocOptions make_options(double gamma, int num_characteristics, MocFlowKind flow_type) {
    MocOptions opts;
    opts.flow_type = flow_type;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::ANALYSIS;
    opts.gamma = gamma;
    opts.num_characteristics = num_characteristics;
    opts.geometry.throat_radius = 1.0;
    return opts;
}

// Kliegel & Levine (1969) Eq. (10): third-order throat axis velocity u0,
// evaluated at the throat plane (z=0), expanded in powers of 1/(R+1).
double reference_kl_u0(double gamma, double R) {
    double Rp1 = R + 1.0;
    double term2 = -1.0 / (4.0 * Rp1);
    double term3 = (10.0 * gamma - 15.0) / (288.0 * Rp1 * Rp1);
    double term4 = -(2708.0 * gamma * gamma + 2079.0 * gamma + 2115.0) / (82944.0 * Rp1 * Rp1 * Rp1);
    return 1.0 + term2 + term3 + term4;
}

// Kliegel & Levine (1969) Eq. (12): third-order throat wall velocity uw,
// evaluated at the throat plane (z=0), expanded in powers of 1/(R+1).
double reference_kl_uw(double gamma, double R) {
    double Rp1 = R + 1.0;
    double term2 = 1.0 / (4.0 * Rp1);
    double term3 = -(14.0 * gamma - 57.0) / (288.0 * Rp1 * Rp1);
    double term4 = (2364.0 * gamma * gamma - 3915.0 * gamma + 14337.0) / (82944.0 * Rp1 * Rp1 * Rp1);
    return 1.0 + term2 + term3 + term4;
}

} // namespace

// ============================================================
// Kliegel-Levine kernel vs. published closed-form equations
// ============================================================

TEST(KliegelLevineClosedForm, AxisAndWallVelocityMatchPublishedEquations) {
    // KL_xMach(r, z, gamma, R) implements Eq. (7) with u1/u2/u3 from the
    // paper's Appendix (Eqs. A1, A3, A5). At r=0, z=0 and r=1, z=0 these
    // reduce algebraically (by hand-derivation, not just numerically) to
    // exactly Eq. (10) and Eq. (12). Sweep gamma and R broadly since the
    // match is exact, not an approximation.
    struct Case { double gamma; double R; };
    std::vector<Case> cases = {
        {1.4, 0.625}, {1.4, 0.1}, {1.4, 1.0}, {1.4, 5.0},
        {1.2, 0.625}, {1.667, 2.0}, {1.25, 0.3},
    };

    NozzleGeometry geom;
    geom.throat_radius = 1.0;
    geom.downstream_wall_curvature_radius = 1.0;
    ThermodynamicContext thermo = make_perfect_gas_context(1.4);
    MocOptions opts = make_options(1.4, 5, MocFlowKind::AXISYMMETRIC);
    MocInitializationTestAccess init(geom, thermo, opts);

    for (const auto& c : cases) {
        double u0 = init.KL_xMach(0.0, 0.0, c.gamma, c.R);
        double uw = init.KL_xMach(1.0, 0.0, c.gamma, c.R);

        EXPECT_NEAR(u0, reference_kl_u0(c.gamma, c.R), 1e-10)
            << "axis throat velocity vs Eq.(10), gamma=" << c.gamma << " R=" << c.R;
        EXPECT_NEAR(uw, reference_kl_uw(c.gamma, c.R), 1e-10)
            << "wall throat velocity vs Eq.(12), gamma=" << c.gamma << " R=" << c.R;
    }
}

TEST(KliegelLevineClosedForm, MatchesTable1ExperimentalCase) {
    // Kliegel & Levine (1969) Table 1, "Theory" column: R=0.625 nozzle
    // (Cuffel et al. experiment, air, gamma=1.4). Published to 3 sig figs.
    double gamma = 1.4;
    double R = 0.625;

    NozzleGeometry geom;
    geom.throat_radius = 1.0;
    geom.downstream_wall_curvature_radius = R;
    ThermodynamicContext thermo = make_perfect_gas_context(gamma);
    MocOptions opts = make_options(gamma, 5, MocFlowKind::AXISYMMETRIC);
    MocInitializationTestAccess init(geom, thermo, opts);

    double u0 = init.KL_xMach(0.0, 0.0, gamma, R);
    double uw = init.KL_xMach(1.0, 0.0, gamma, R);

    EXPECT_NEAR(u0, 0.816, 2e-3) << "axis throat velocity u0 vs published Table 1 value";
    EXPECT_NEAR(uw, 1.24, 2e-3) << "wall throat velocity uw vs published Table 1 value";
}

TEST(KliegelLevineClosedForm, WallLeadsAxisThroughTransonicRegion) {
    // Physical content of Table 1 / Figs. 2-3: for a convex downstream throat
    // wall (R>0), the wall reaches sonic velocity upstream of the geometric
    // throat while the axis is still subsonic there -- at the throat plane
    // z=0, u_wall > 1 > u_axis. Holds for both small and large R since the
    // leading-order term (+-1/(4(R+1))) dominates in both limits.
    double gamma = 1.4;
    NozzleGeometry geom;
    geom.throat_radius = 1.0;
    geom.downstream_wall_curvature_radius = 1.0;
    ThermodynamicContext thermo = make_perfect_gas_context(gamma);
    MocOptions opts = make_options(gamma, 5, MocFlowKind::AXISYMMETRIC);
    MocInitializationTestAccess init(geom, thermo, opts);

    for (double R : {0.1, 0.3, 0.625, 1.0, 2.0, 5.0, 20.0}) {
        double u0 = init.KL_xMach(0.0, 0.0, gamma, R);
        double uw = init.KL_xMach(1.0, 0.0, gamma, R);
        EXPECT_LT(u0, 1.0) << "axis should be subsonic at throat plane, R=" << R;
        EXPECT_GT(uw, 1.0) << "wall should be supersonic at throat plane, R=" << R;
    }
}

// ============================================================
// Newton solve for the transonic (zero radial-velocity) line
// ============================================================

TEST(KliegelLevineTransonicLine, SolvedLineHasZeroRadialVelocity) {
    // KL_solve_transonic_x finds x such that the radial velocity component
    // (KL_yMach, the series of Eq. (8)) vanishes. Verify the returned points
    // actually satisfy their own defining equation -- a regression guard on
    // both the Newton iteration and the KL_v1/v2/v3 polynomials.
    double gamma = 1.4;
    double R = 0.8;
    NozzleGeometry geom;
    geom.throat_radius = 1.0;
    geom.downstream_wall_curvature_radius = R;
    ThermodynamicContext thermo = make_perfect_gas_context(gamma);
    MocOptions opts = make_options(gamma, 9, MocFlowKind::AXISYMMETRIC);
    MocInitializationTestAccess init(geom, thermo, opts);

    double alpha = init.sauer_alpha(gamma);
    for (double y : {0.2, 0.4, 0.6, 0.8, 1.0}) {
        double x_guess = (gamma + 1.0) * alpha / (2.0 * (3.0 + init.delta())) * (1.0 - y * y);
        double x = init.KL_solve_transonic_x(y, gamma, R, x_guess);
        // The solver converges the residual normalized by y (see KL_solve_transonic_x)
        // to tol=1e-6, so check the same normalized quantity here.
        double normalized_residual = init.KL_yMach(x, y, gamma, R) / y;
        EXPECT_NEAR(normalized_residual, 0.0, 1e-6) << "radial velocity should vanish on transonic line at y=" << y;
    }
}

// ============================================================
// Kliegel-Levine -> Sauer asymptotic consistency
// ============================================================

TEST(KliegelLevineVsSauer, ConvergesToSauerAsCurvatureRadiusGrows) {
    // Kliegel & Levine (1969), p. 1376: "Since the two solutions must be
    // identical in the limit of large radii of curvature, one can transform
    // Hall's results..." -- the Sauer parabola is exactly the leading
    // (1/(R+1)) term of the KL series (KL_v1(r,z)=0 gives z=0.25*(1-r^2),
    // which is the Sauer x-coordinate after the z<->x change of variables).
    // So initialize_kliegel_levine's transonic line must converge to
    // initialize_sauer's as R grows -- once the constant downstream shift
    // (MocOptions::initial_line_axial_shift, applied so the line is usable for
    // dual-family seeding; see instructions/moc_convergence_roadmap.md Sec 2
    // Step 0) is disabled, since that shift is an intentional, R-independent
    // offset that this asymptotic identity was never about.
    double gamma = 1.4;
    MocOptions opts = make_options(gamma, 9, MocFlowKind::AXISYMMETRIC);
    opts.initial_line_axial_shift = 0.0;
    ThermodynamicContext thermo = make_perfect_gas_context(gamma);
    ThroatCondition throat = make_perfect_gas_throat(gamma);

    auto max_x_diff = [&](double R) {
        NozzleGeometry geom;
        geom.throat_radius = 1.0;
        geom.downstream_wall_curvature_radius = R;
        MocInitialization init(geom, thermo, opts);
        auto sauer_pts = init.initialize_sauer(throat);
        auto kl_pts = init.initialize_kliegel_levine(throat);
        double max_diff = 0.0;
        for (size_t i = 1; i < sauer_pts.size(); i++) { // skip y=0: KL uses an epsilon fallback there
            max_diff = std::max(max_diff, std::abs(sauer_pts[i].x - kl_pts[i].x));
        }
        return max_diff;
    };

    double diff_moderate_R = max_x_diff(2.0);
    double diff_large_R = max_x_diff(1.0e4);

    EXPECT_LT(diff_large_R, diff_moderate_R / 100.0)
        << "KL/Sauer transonic-line discrepancy should shrink well faster than linearly as R grows";
    EXPECT_LT(diff_large_R, 1e-4);
}

// ============================================================
// Sauer transonic line vs. closed form
// ============================================================

TEST(SauerInitialization, TransonicLineMatchesClosedForm) {
    // Sauer's parabola x(y) = (gamma+1)*alpha/(2*(3+delta)) * (1-y^2), with
    // alpha = sqrt((1+delta)/((gamma+1)*R)). Cross-checked independently
    // above (KliegelLevineVsSauer) as the R->infinity limit of the
    // Kliegel-Levine series, not merely copied from the implementation.
    for (double gamma : {1.2, 1.4, 1.667}) {
        for (double R : {0.5, 1.0, 3.0}) {
            for (MocFlowKind flow : {MocFlowKind::PLANAR, MocFlowKind::AXISYMMETRIC}) {
                double delta = (flow == MocFlowKind::AXISYMMETRIC) ? 1.0 : 0.0;

                NozzleGeometry geom;
                geom.throat_radius = 1.0;
                geom.downstream_wall_curvature_radius = R;
                ThermodynamicContext thermo = make_perfect_gas_context(gamma);
                MocOptions opts = make_options(gamma, 5, flow);
                MocInitialization init(geom, thermo, opts);
                ThroatCondition throat = make_perfect_gas_throat(gamma);

                auto points = init.initialize_sauer(throat);
                double alpha = std::sqrt((1.0 + delta) / ((gamma + 1.0) * R));

                for (const auto& pt : points) {
                    double x_ref = (gamma + 1.0) * alpha / (2.0 * (3.0 + delta)) * (1.0 - pt.y * pt.y);
                    EXPECT_NEAR(pt.x, x_ref, 1e-12)
                        << "Sauer x(y) vs closed form, gamma=" << gamma << " R=" << R << " y=" << pt.y;
                    EXPECT_NEAR(pt.theta, 0.0, 1e-15);
                }
            }
        }
    }
}

TEST(SauerInitialization, MachAtAxisAndWallMatchesClosedForm) {
    // At y=1 (wall), x=0 by construction, so the Mach relation reduces to
    // M = 1 + (gamma+1)*alpha^2/(2*(1+delta)). At y=0 (axis), the y^2 term
    // vanishes, so M = 1 + alpha*x_axis.
    double gamma = 1.4;
    double R = 1.5;
    double delta = 1.0; // axisymmetric

    NozzleGeometry geom;
    geom.throat_radius = 1.0;
    geom.downstream_wall_curvature_radius = R;
    ThermodynamicContext thermo = make_perfect_gas_context(gamma);
    MocOptions opts = make_options(gamma, 5, MocFlowKind::AXISYMMETRIC);
    MocInitialization init(geom, thermo, opts);
    ThroatCondition throat = make_perfect_gas_throat(gamma);

    auto points = init.initialize_sauer(throat);
    ASSERT_GE(points.size(), 2u);

    double alpha = std::sqrt((1.0 + delta) / ((gamma + 1.0) * R));
    double x_axis = (gamma + 1.0) * alpha / (2.0 * (3.0 + delta));

    const auto& axis_pt = points.front();
    const auto& wall_pt = points.back();
    ASSERT_NEAR(axis_pt.y, 0.0, 1e-15);
    ASSERT_NEAR(wall_pt.y, 1.0, 1e-15);

    double expected_mach_axis = 1.0 + alpha * x_axis;
    double expected_mach_wall = 1.0 + (gamma + 1.0) * alpha * alpha / (2.0 * (1.0 + delta));

    EXPECT_NEAR(axis_pt.mach, expected_mach_axis, 1e-9);
    EXPECT_NEAR(wall_pt.mach, expected_mach_wall, 1e-9);
}

// ============================================================
// Public-API structural contract
// ============================================================

TEST(KliegelLevineInitialization, ProducesRequestedNumberOfPointsWithMonotonicY) {
    // The line is rigidly shifted downstream of the raw sonic (v=0) locus (default
    // MocOptions::initial_line_axial_shift), so theta is pinned to exactly 0 only at
    // the axis (r=0, where the radial-velocity series vanishes identically) -- it is
    // strictly positive everywhere else, and increases monotonically off-axis for a
    // downstream-shifted diverging locus.
    double gamma = 1.4;
    double R = 1.0;
    int n = 11;

    NozzleGeometry geom;
    geom.throat_radius = 1.0;
    geom.downstream_wall_curvature_radius = R;
    ThermodynamicContext thermo = make_perfect_gas_context(gamma);
    MocOptions opts = make_options(gamma, n, MocFlowKind::AXISYMMETRIC);
    MocInitialization init(geom, thermo, opts);
    ThroatCondition throat = make_perfect_gas_throat(gamma);

    auto points = init.initialize_kliegel_levine(throat);
    ASSERT_EQ(points.size(), static_cast<size_t>(n));
    EXPECT_NEAR(points[0].theta, 0.0, 1e-15) << "axis theta should be pinned to 0";
    for (size_t i = 0; i < points.size(); i++) {
        EXPECT_NEAR(points[i].y, static_cast<double>(i) / (n - 1), 1e-12) << "index " << i;
        EXPECT_GT(points[i].mach, 0.0) << "index " << i;
        if (i > 0) {
            EXPECT_GT(points[i].theta, 0.0) << "index " << i;
            EXPECT_GE(points[i].theta, points[i - 1].theta) << "theta should increase monotonically off-axis, index " << i;
        }
    }
}
