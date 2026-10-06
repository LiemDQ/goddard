#include "goddard/shocks.hpp"
#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/error.hpp"
#include "goddard/global.hpp"
#include "goddard/numerics.hpp"
#include <cmath>
#include <string>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

// ============================================================
//  Perfect-gas normal shock
// ============================================================

TEST(NormalShock, SonicLimit) {
    // M1 = 1 is a vanishingly weak shock (sound wave): no change
    ShockResult r = normal_shock(1.0, 1.4);
    ASSERT_TRUE(r.valid);
    EXPECT_NEAR(r.mach_out, 1.0, 1e-12);
    EXPECT_NEAR(r.static_pressure_ratio, 1.0, 1e-12);
    EXPECT_NEAR(r.static_temperature_ratio, 1.0, 1e-12);
    EXPECT_NEAR(r.total_pressure_ratio, 1.0, 1e-12);
}

TEST(NormalShock, DownstreamSubsonic) {
    // Post-shock Mach must be subsonic for any supersonic inflow
    for (double M1 : {1.01, 1.5, 2.0, 3.0, 5.0, 10.0}) {
        ShockResult r = normal_shock(M1, 1.4);
        ASSERT_TRUE(r.valid) << "M1=" << M1;
        EXPECT_LT(r.mach_out, 1.0) << "M2 must be < 1 for M1=" << M1;
        EXPECT_GT(r.mach_out, 0.0) << "M2 must be > 0 for M1=" << M1;
    }
}

TEST(NormalShock, StaticPressureIncreases) {
    for (double M1 : {1.5, 2.0, 5.0}) {
        ShockResult r = normal_shock(M1, 1.4);
        ASSERT_TRUE(r.valid);
        EXPECT_GT(r.static_pressure_ratio, 1.0)
            << "Static pressure must increase across shock at M1=" << M1;
    }
}

TEST(NormalShock, StaticTemperatureIncreases) {
    for (double M1 : {1.5, 2.0, 5.0}) {
        ShockResult r = normal_shock(M1, 1.4);
        ASSERT_TRUE(r.valid);
        EXPECT_GT(r.static_temperature_ratio, 1.0)
            << "Static temperature must increase across shock at M1=" << M1;
    }
}

TEST(NormalShock, TotalPressureDecreases) {
    // Entropy increases across a shock, so total pressure must decrease
    for (double M1 : {1.01, 1.5, 2.0, 5.0}) {
        ShockResult r = normal_shock(M1, 1.4);
        ASSERT_TRUE(r.valid);
        EXPECT_LE(r.total_pressure_ratio, 1.0)
            << "Total pressure must not increase across shock at M1=" << M1;
    }
}

TEST(NormalShock, DensityRatioConsistency) {
    // The density ratio can be computed two independent ways:
    //   rho2/rho1 = (P2/P1) / (T2/T1)          (ideal gas law)
    //   rho2/rho1 = (gamma+1)*M1^2 / ((gamma-1)*M1^2 + 2)   (Rankine-Hugoniot)
    double gamma = 1.4;
    for (double M1 : {1.5, 2.0, 3.0, 5.0, 10.0}) {
        ShockResult r = normal_shock(M1, gamma);
        ASSERT_TRUE(r.valid);

        double rho_ratio_from_state = r.static_pressure_ratio / r.static_temperature_ratio;
        double rho_ratio_RH = (gamma + 1) * M1 * M1
            / ((gamma - 1) * M1 * M1 + 2);

        EXPECT_NEAR(rho_ratio_from_state, rho_ratio_RH,
                    max_fp_error(rho_ratio_RH, 1e-10, 1e-12))
            << "Density ratio inconsistency at M1=" << M1;
        EXPECT_NEAR(r.density_ratio, rho_ratio_RH, max_fp_error(rho_ratio_RH, 1e-12, 1e-12));
    }
}

TEST(NormalShock, InvalidSubsonicMach) {
    ShockResult r = normal_shock(0.5, 1.4);
    EXPECT_FALSE(r.valid);
}

TEST(NormalShock, InvalidGamma) {
    ShockResult r = normal_shock(2.0, 0.5);
    EXPECT_FALSE(r.valid);
}

TEST(NormalShock, DifferentGamma) {
    // Monatomic gas (gamma = 5/3): same invariants must hold
    double gamma = 5.0 / 3.0;
    for (double M1 : {1.5, 3.0, 5.0}) {
        ShockResult r = normal_shock(M1, gamma);
        ASSERT_TRUE(r.valid);
        EXPECT_LT(r.mach_out, 1.0);
        EXPECT_GT(r.static_pressure_ratio, 1.0);
        EXPECT_GT(r.static_temperature_ratio, 1.0);
        EXPECT_LE(r.total_pressure_ratio, 1.0);
    }
}

struct NormalShockData {
    double gamma;
    double mach_in;
    double mach_out;
    double static_pressure_ratio;
    double static_temperature_ratio;
    double total_pressure_ratio;
};

TEST(NormalShock, AndersonExercise3p5) {
    // Values from Anderson example 3.5
    NormalShockData data{1.4, 3.0, 0.4752, 10.333, 2.679, -1};
    ShockResult r = normal_shock(data.mach_in, data.gamma);
    ASSERT_TRUE(r.valid);
    EXPECT_NEAR(r.mach_out, data.mach_out,
                max_fp_error(data.mach_out, 1e-4, 1e-6));
    EXPECT_NEAR(r.static_pressure_ratio, data.static_pressure_ratio,
                max_fp_error(data.static_pressure_ratio, 1e-4, 1e-6));
    EXPECT_NEAR(r.static_temperature_ratio, data.static_temperature_ratio,
                max_fp_error(data.static_temperature_ratio, 1e-4, 1e-6));
}

TEST(NormalShock, AndersonExercise3p6) {
    // Values from Anderson example 3.6
    double mach1 = 2.0;
    double gamma = 1.4;

    ShockResult r = normal_shock(mach1, gamma);
    double T1 = 519.0; // rankine
    double P1 = 2116.0; // lb/ft^2
    double stag_factor = stagnation_factor(r.mach_out, gamma);
    double P1_stag = perfect_gas_stagnation_pressure(P1, mach1, gamma);
    double T2_stag = stag_factor*r.static_temperature_ratio*T1;
    double P2_stag = r.total_pressure_ratio*P1_stag;
    ASSERT_TRUE(r.valid);

    EXPECT_NEAR(T2_stag, 934.2,
                max_fp_error(934.2, 1e-4, 1e-6));
    EXPECT_NEAR(P2_stag, 11935.0,
                max_fp_error(11935.0, 1e-4, 1e-6));
}

class NormalShockAnderson : public ::testing::TestWithParam<NormalShockData> {};

TEST_P(NormalShockAnderson, TableA2) {
    // normal shock tables (Appendix B or Table A.2). Each entry should list
    // M1, gamma, and the expected M2, P2/P1, T2/T1, P02/P01.
    auto data = GetParam();
    ShockResult r = normal_shock(data.mach_in, data.gamma);
    ASSERT_TRUE(r.valid);
    EXPECT_NEAR(r.mach_out, data.mach_out,
                max_fp_error(data.mach_out, 1e-3, 1e-4));
    EXPECT_NEAR(r.static_pressure_ratio, data.static_pressure_ratio,
                max_fp_error(data.static_pressure_ratio, 1e-3, 1e-4));
    EXPECT_NEAR(r.static_temperature_ratio, data.static_temperature_ratio,
                max_fp_error(data.static_temperature_ratio, 1e-3, 1e-4));
    EXPECT_NEAR(r.total_pressure_ratio, data.total_pressure_ratio,
                max_fp_error(data.total_pressure_ratio, 1e-3, 1e-4));
}


// Format: {M1, gamma, M2, P2/P1, T2/T1, P02/P01}
INSTANTIATE_TEST_SUITE_P(
    Anderson,
    NormalShockAnderson,
    ::testing::Values(
        NormalShockData{1.4, 1.060, 0.9444, 1.144, 1.039, 0.9998},
        NormalShockData{1.4, 1.180, 0.8549, 1.458, 1.115, 0.9946},
        NormalShockData{1.4, 1.380, 0.7483, 2.055, 1.242, 0.9630},
        NormalShockData{1.4, 1.580, 0.6746, 2.746, 1.374, 0.9026},
        NormalShockData{1.4, 1.780, 0.6210, 3.530, 1.517, 0.8215},
        NormalShockData{1.4, 1.980, 0.5808, 4.407, 1.671, 0.7302},
        NormalShockData{1.4, 2.500, 0.5130, 7.125, 2.137, 0.4990},
        NormalShockData{1.4, 3.500, 0.4512, 14.12, 3.315, 0.2129},
        NormalShockData{1.4, 4.500, 0.4236, 23.46, 4.875, 0.0917},
        NormalShockData{1.4, 6.000, 0.4042, 41.83, 7.941, 0.02965},
        NormalShockData{1.4, 10.00, 0.3876, 116.5, 20.39, 0.003045},
        NormalShockData{1.4, 18.00, 0.3810, 377.8, 63.94, 0.1807e-3},
        NormalShockData{1.4, 36.00, 0.3787, 1512., 252.9, 0.5874e-5},
        NormalShockData{1.4, 50.00, 0.3784, 2916., 487.1, 0.1144e-5}
    )
);

// ============================================================
//  Perfect-gas oblique shock
// ============================================================

TEST(ObliqueShock, DeflectionAngleRoundTrip) {
    // wave_angle → deflection_angle should recover the original deflection
    double gamma = 1.4;
    for (double M : {2.0, 3.0, 5.0}) {
        for (double theta_deg : {5.0, 10.0, 15.0}) {
            double theta = theta_deg * DEG;
            auto [weak_beta, strong_beta] = oblique_shock_wave_angle(M, theta, gamma);

            double theta_recovered_weak = oblique_shock_deflection_angle(M, weak_beta, gamma);
            double theta_recovered_strong = oblique_shock_deflection_angle(M, strong_beta, gamma);

            EXPECT_NEAR(theta_recovered_weak, theta, 1e-8)
                << "Weak round-trip failed at M=" << M << " theta=" << theta_deg << " deg";
            EXPECT_NEAR(theta_recovered_strong, theta, 1e-8)
                << "Strong round-trip failed at M=" << M << " theta=" << theta_deg << " deg";
        }
    }
}

TEST(ObliqueShock, WeakAngleBounds) {
    // Weak shock wave angle must be > Mach angle and < 90 degrees
    double gamma = 1.4;
    for (double M : {1.5, 2.0, 3.0, 5.0}) {
        double mu = asin(1.0 / M);  // Mach angle
        for (double theta_deg : {5.0, 10.0}) {
            double theta = theta_deg * DEG;
            auto [weak_beta, strong_beta] = oblique_shock_wave_angle(M, theta, gamma);

            EXPECT_GT(weak_beta, mu)
                << "Weak beta must exceed Mach angle at M=" << M;
            EXPECT_LT(weak_beta, M_PI / 2.0)
                << "Weak beta must be < 90 deg at M=" << M;
        }
    }
}

TEST(ObliqueShock, StrongAngleBounds) {
    // Strong shock wave angle must be > weak and <= 90 degrees
    double gamma = 1.4;
    for (double M : {2.0, 3.0, 5.0}) {
        double theta = 10.0 * DEG;
        auto [weak_beta, strong_beta] = oblique_shock_wave_angle(M, theta, gamma);

        EXPECT_GT(strong_beta, weak_beta)
            << "Strong beta must exceed weak beta at M=" << M;
        EXPECT_LE(strong_beta, M_PI / 2.0)
            << "Strong beta must be <= 90 deg at M=" << M;
    }
}

TEST(ObliqueShock, NormalComponentConsistency) {
    // The normal shock embedded in an oblique shock should match
    // a standalone normal shock at M_n1 = M * sin(beta)
    double gamma = 1.4;
    double M = 3.0;
    double beta = 40.0 * DEG;

    ObliqueShockResult obl = oblique_shock_from_wave_angle(M, beta, gamma);
    double M_n1 = M * sin(beta);
    ShockResult ns = normal_shock(M_n1, gamma);

    ASSERT_TRUE(obl.valid);
    ASSERT_TRUE(ns.valid);
    EXPECT_NEAR(obl.shock.mach_out, ns.mach_out, 1e-12)
        << "Normal component M2 should match standalone normal shock";
    EXPECT_NEAR(obl.shock.static_pressure_ratio, ns.static_pressure_ratio, 1e-12);
    EXPECT_NEAR(obl.shock.static_temperature_ratio, ns.static_temperature_ratio, 1e-12);
    EXPECT_NEAR(obl.shock.total_pressure_ratio, ns.total_pressure_ratio, 1e-12);
}

TEST(ObliqueShock, FromDeflectionConsistency) {
    // oblique_shock_from_deflection should produce a result whose theta
    // matches the input, and whose beta matches oblique_shock_wave_angle
    double gamma = 1.4;
    double M = 3.0;
    double theta = 15.0 * DEG;

    ObliqueShockResult weak = oblique_shock_from_deflection(M, theta, gamma, true);
    ObliqueShockResult strong = oblique_shock_from_deflection(M, theta, gamma, false);
    auto [weak_beta, strong_beta] = oblique_shock_wave_angle(M, theta, gamma);

    ASSERT_TRUE(weak.valid);
    ASSERT_TRUE(strong.valid);

    EXPECT_NEAR(weak.theta, theta, 1e-12);
    EXPECT_NEAR(strong.theta, theta, 1e-12);
    EXPECT_NEAR(weak.beta, weak_beta, 1e-10)
        << "Weak beta should match oblique_shock_wave_angle";
    EXPECT_NEAR(strong.beta, strong_beta, 1e-10)
        << "Strong beta should match oblique_shock_wave_angle";
}

TEST(ObliqueShock, WeakDownstreamSupersonicStrongSubsonic) {
    // For moderate deflection angles, the weak solution is supersonic
    // downstream and the strong solution is subsonic downstream
    double gamma = 1.4;
    double M = 3.0;
    double theta = 10.0 * DEG;

    ObliqueShockResult weak = oblique_shock_from_deflection(M, theta, gamma, true);
    ObliqueShockResult strong = oblique_shock_from_deflection(M, theta, gamma, false);

    ASSERT_TRUE(weak.valid);
    ASSERT_TRUE(strong.valid);

    EXPECT_GT(weak.mach_out, 1.0)
        << "Weak oblique shock should be supersonic downstream";
    EXPECT_LT(strong.mach_out, 1.0)
        << "Strong oblique shock should be subsonic downstream";
}

// --- Anderson oblique shock textbook data stubs ---

struct ObliqueShockData {
    double mach_in;
    double gamma;
    double deflection_angle_deg;
    double weak_wave_angle_deg;
    double strong_wave_angle_deg;
};

class ObliqueShockAnderson : public ::testing::TestWithParam<ObliqueShockData> {};

TEST_P(ObliqueShockAnderson, TableValues) {
    // TODO: Fill in expected values from Anderson's Modern Compressible Flow
    // oblique shock charts/tables. Each entry should list M1, gamma,
    // deflection angle (deg), and the expected weak/strong wave angles (deg).
    auto data = GetParam();
    auto [weak_beta, strong_beta] = oblique_shock_wave_angle(
        data.mach_in, data.deflection_angle_deg * DEG, data.gamma);

    EXPECT_NEAR(weak_beta / DEG, data.weak_wave_angle_deg, 0.05)
        << "Weak wave angle mismatch at M=" << data.mach_in
        << " theta=" << data.deflection_angle_deg << " deg";
    EXPECT_NEAR(strong_beta / DEG, data.strong_wave_angle_deg, 0.05)
        << "Strong wave angle mismatch at M=" << data.mach_in
        << " theta=" << data.deflection_angle_deg << " deg";
}

// TODO: Replace placeholder data below with values from Anderson tables.
// Format: {M1, gamma, theta_deg, weak_beta_deg, strong_beta_deg}
INSTANTIATE_TEST_SUITE_P(
    Anderson,
    ObliqueShockAnderson,
    ::testing::Values(
        ObliqueShockData{2.0, 1.4, 10.0, 39.3139, 83.7124}
    )
);


// ============================================================
//  Perfect-gas reflected shock
// ============================================================

TEST(ReflectedShock, EndWallGasIsAtRest) {
    // Gas 2 moves at u_p = u1 (1 - rho1/rho2) toward the wall. The reflected shock must bring it
    // to rest: in the reflected-shock frame gas 2 enters at M_R a2 = u_p + W_R and gas 5 leaves at
    // W_R = u_p/(rho5/rho2 - 1). Check that M_R from the perfect-gas relation satisfies this.
    for (double gamma : {1.4, 5.0/3.0}) {
        for (double M1 : {1.5, 2.0, 4.0, 8.0}) {
            ReflectedShockResult r = reflected_shock(M1, gamma);
            ASSERT_TRUE(r.valid);
            const double a1 = 1.0;
            const double a2 = a1*std::sqrt(r.incident.static_temperature_ratio);
            const double u1 = M1*a1;
            const double particle_velocity = u1*(1.0 - 1.0/r.incident.density_ratio);
            const double wave_speed = particle_velocity/(r.reflected.density_ratio - 1.0);
            EXPECT_NEAR(r.reflected.mach_in*a2, particle_velocity + wave_speed,
                        max_fp_error(particle_velocity + wave_speed, 1e-10, 1e-12))
                << "gamma=" << gamma << " M1=" << M1;
            EXPECT_NEAR(r.incident.mach_in, M1, 1e-14);
        }
    }
}

TEST(ReflectedShock, SonicLimit) {
    // A sound wave reflects as a sound wave
    ReflectedShockResult r = reflected_shock(1.0, 1.4);
    ASSERT_TRUE(r.valid);
    EXPECT_NEAR(r.reflected.mach_in, 1.0, 1e-12);
    EXPECT_NEAR(r.reflected.static_pressure_ratio, 1.0, 1e-12);
}

TEST(ReflectedShock, Invalid) {
    EXPECT_FALSE(reflected_shock(0.5, 1.4).valid);
    EXPECT_FALSE(reflected_shock(2.0, 1.0).valid);
}

// ============================================================
//  Perfect-gas maximum deflection
// ============================================================

TEST(ObliqueShockMaxDeflection, ClosedFormValues) {
    // NACA 1135 eq. 168. At M = 2, gamma = 1.4: beta_max = 64.67 deg and theta_max = 22.97 deg
    // (NACA 1135 chart 2). As M -> infinity, sin^2(beta_max) -> (gamma+1)/(2 gamma), so
    // beta_max = 67.79 deg and theta_max = 45.58 deg.
    EXPECT_NEAR(oblique_shock_max_deflection_wave_angle(2.0, 1.4)/DEG, 64.67, 0.01);
    EXPECT_NEAR(oblique_shock_max_deflection(2.0, 1.4)/DEG, 22.97, 0.01);

    const double beta_limit = std::asin(std::sqrt(2.4/2.8));
    EXPECT_NEAR(oblique_shock_max_deflection_wave_angle(1e4, 1.4), beta_limit, 1e-6);
    EXPECT_NEAR(oblique_shock_max_deflection(1e4, 1.4)/DEG, 45.58, 0.01);

    // A sonic flow cannot be deflected by a shock
    EXPECT_NEAR(oblique_shock_max_deflection(1.0, 1.4), 0.0, 1e-12);
}

TEST(ObliqueShockMaxDeflection, IsMaximumOfThetaBetaM) {
    for (double M : {1.5, 3.0, 10.0}) {
        const double beta_max = oblique_shock_max_deflection_wave_angle(M, 1.4);
        const double theta_max = oblique_shock_max_deflection(M, 1.4);
        for (double offset : {-1e-3, 1e-3}) {
            EXPECT_LT(oblique_shock_deflection_angle(M, beta_max + offset, 1.4), theta_max) << "M=" << M;
        }
        // Weak and strong branches merge at the maximum deflection
        auto [weak_beta, strong_beta] = oblique_shock_wave_angle(M, theta_max, 1.4);
        EXPECT_NEAR(weak_beta, beta_max, 1e-4) << "M=" << M;
        EXPECT_NEAR(strong_beta, beta_max, 1e-4) << "M=" << M;
    }
}

TEST(ObliqueShockMaxDeflection, DetachedShockIsInvalid) {
    const double theta = oblique_shock_max_deflection(3.0, 1.4) + 1e-3;
    auto [weak_beta, strong_beta] = oblique_shock_wave_angle(3.0, theta, 1.4);
    EXPECT_TRUE(std::isnan(weak_beta));
    EXPECT_TRUE(std::isnan(strong_beta));
    EXPECT_FALSE(oblique_shock_from_deflection(3.0, theta, 1.4, true).valid);
    EXPECT_FALSE(oblique_shock_from_deflection(3.0, theta, 1.4, false).valid);
}

TEST(ObliqueShock, ZeroDeflectionIsMachWaveOrNormalShock) {
    const double M = 3.0;
    auto [weak_beta, strong_beta] = oblique_shock_wave_angle(M, 0.0, 1.4);
    EXPECT_NEAR(weak_beta, std::asin(1.0/M), 1e-14);
    EXPECT_NEAR(strong_beta, M_PI/2, 1e-14);

    ObliqueShockResult weak = oblique_shock_from_deflection(M, 0.0, 1.4, true);
    ASSERT_TRUE(weak.valid);
    EXPECT_NEAR(weak.shock.static_pressure_ratio, 1.0, 1e-12);
    EXPECT_NEAR(weak.mach_out, M, 1e-10);
}

TEST(ObliqueShock, BelowMachAngleIsInvalid) {
    EXPECT_FALSE(oblique_shock_from_wave_angle(2.0, 0.9*std::asin(0.5), 1.4).valid);
}

// ============================================================
//  ShockSolver: real-gas fixtures
// ============================================================

namespace {

/** Gas from h2o2.yaml at T [K], P [Pa] and the given mole fractions. */
Gas make_gas(const std::string& composition, GasChemistry chemistry,
             double T = 300.0, double P = Cantera::OneAtm) {
    Goddard::setup_defaults();
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    sol->thermo()->setState_TPX(T, P, composition);
    return Gas(sol, chemistry);
}

// The shock Newton iteration stops once the step in ln(P2/P1) and ln(T2/T1) is below
// SolverOptions::abstol = 1e-6, so ratios carry relative errors of ~1e-6. Comparisons use 1e-5.
constexpr double RATIO_TOL = 1e-5;

void expect_shock_near(const ShockResult& actual, const ShockResult& expected, double tol,
                       const std::string& label) {
    ASSERT_TRUE(actual.valid) << label;
    ASSERT_TRUE(expected.valid) << label;
    EXPECT_NEAR(actual.mach_in, expected.mach_in, max_fp_error(expected.mach_in, tol, 1e-12)) << label;
    EXPECT_NEAR(actual.mach_out, expected.mach_out, max_fp_error(expected.mach_out, tol, 1e-12)) << label;
    EXPECT_NEAR(actual.static_pressure_ratio, expected.static_pressure_ratio,
                max_fp_error(expected.static_pressure_ratio, tol, 1e-12)) << label;
    EXPECT_NEAR(actual.static_temperature_ratio, expected.static_temperature_ratio,
                max_fp_error(expected.static_temperature_ratio, tol, 1e-12)) << label;
    EXPECT_NEAR(actual.density_ratio, expected.density_ratio,
                max_fp_error(expected.density_ratio, tol, 1e-12)) << label;
    EXPECT_NEAR(actual.total_pressure_ratio, expected.total_pressure_ratio,
                max_fp_error(expected.total_pressure_ratio, tol, 1e-12)) << label;
}

const char* chemistry_name(GasChemistry chemistry) {
    return chemistry == GasChemistry::EQUILIBRIUM ? "EQUILIBRIUM" : "FROZEN";
}

} // namespace

// ============================================================
//  ShockSolver on argon: calorically perfect, so it must reproduce the gamma = 5/3 relations
// ============================================================

class ShockSolverArgonTests : public ::testing::TestWithParam<GasChemistry> {};

TEST_P(ShockSolverArgonTests, NormalShockMatchesPerfectGas) {
    // Argon has cp = 5R/2 exactly in h2o2.yaml and no reactions, so frozen and equilibrium
    // shocks both reduce to the perfect-gas relations with gamma = 5/3.
    ShockSolver solver(make_gas("AR:1", GetParam()));
    for (double M : {1.5, 3.0, 6.0}) {
        expect_shock_near(solver.normal_shock(M), normal_shock(M, 5.0/3.0), RATIO_TOL,
                          std::string(chemistry_name(GetParam())) + " M=" + std::to_string(M));
    }
}

TEST_P(ShockSolverArgonTests, ReflectedShockMatchesPerfectGas) {
    // Regression test: the reflected shock must start from state 2 and the particle velocity
    // u_p, not from state 1 and u1.
    ShockSolver solver(make_gas("AR:1", GetParam()));
    for (double M : {1.5, 3.0, 6.0}) {
        const std::string label = std::string(chemistry_name(GetParam())) + " M=" + std::to_string(M);
        ReflectedShockResult actual = solver.reflected_shock(M);
        ReflectedShockResult expected = reflected_shock(M, 5.0/3.0);
        ASSERT_TRUE(actual.valid) << label;
        expect_shock_near(actual.incident, expected.incident, RATIO_TOL, label + " incident");
        expect_shock_near(actual.reflected, expected.reflected, RATIO_TOL, label + " reflected");
    }
}

TEST_P(ShockSolverArgonTests, ObliqueShockMatchesPerfectGas) {
    ShockSolver solver(make_gas("AR:1", GetParam()));
    const double gamma = 5.0/3.0;
    const double M = 3.0;

    ObliqueShockResult from_beta = solver.oblique_shock_from_wave_angle(M, 35.0*DEG);
    ObliqueShockResult from_beta_pg = oblique_shock_from_wave_angle(M, 35.0*DEG, gamma);
    ASSERT_TRUE(from_beta.valid);
    EXPECT_NEAR(from_beta.theta, from_beta_pg.theta, 1e-6);
    EXPECT_NEAR(from_beta.mach_out, from_beta_pg.mach_out, max_fp_error(from_beta_pg.mach_out, RATIO_TOL, 0.0));
    expect_shock_near(from_beta.shock, from_beta_pg.shock, RATIO_TOL, "wave angle");

    // The maximum deflection is a flat maximum, so its value is far more accurate than the
    // golden-section tolerance on the wave angle (1e-6 rad).
    EXPECT_NEAR(solver.max_deflection(M), oblique_shock_max_deflection(M, gamma), 1e-8);

    // Bisection stops at a deflection residual of 1e-6 rad; d(theta)/d(beta) is O(1) away from
    // the maximum, so the wave angle agrees to ~1e-6 rad.
    for (bool weak : {true, false}) {
        ObliqueShockResult actual = solver.oblique_shock_from_deflection(M, 15.0*DEG, weak);
        ObliqueShockResult expected = oblique_shock_from_deflection(M, 15.0*DEG, gamma, weak);
        ASSERT_TRUE(actual.valid) << "weak=" << weak;
        EXPECT_NEAR(actual.beta, expected.beta, 1e-5) << "weak=" << weak;
        EXPECT_NEAR(actual.theta, 15.0*DEG, 1e-6) << "weak=" << weak;
    }
}

INSTANTIATE_TEST_SUITE_P(Chemistry, ShockSolverArgonTests,
    ::testing::Values(GasChemistry::FROZEN, GasChemistry::EQUILIBRIUM),
    [](const ::testing::TestParamInfo<GasChemistry>& param_info) {
        return std::string(chemistry_name(param_info.param));
    });

// ============================================================
//  ShockSolver with frozen chemistry on air
// ============================================================

class ShockSolverFrozenTests : public ::testing::Test {
protected:
    ShockSolverFrozenTests() {
        Goddard::setup_defaults();
        gas = Cantera::newSolution("h2o2.yaml", "ohmech");
        // Set to a simple diatomic-like state at moderate temperature
        // so gamma is approximately constant
        gas->thermo()->setState_TPX(300.0, Cantera::OneAtm, "N2:0.79, O2:0.21");
        gas->thermo()->saveState(initial_state);
    }

    std::shared_ptr<Cantera::Solution> gas;
    std::vector<double> initial_state;
};

TEST_F(ShockSolverFrozenTests, ApproachesPerfectGas) {
    // At 300 K where air behaves as a perfect gas with gamma ~ 1.4,
    // the Cantera result should be close to the perfect-gas solution
    double gamma = gas->thermo()->cp_mass() / gas->thermo()->cv_mass();
    double mach = 2.0;

    ShockResult r_perf = normal_shock(mach, gamma);
    ASSERT_TRUE(r_perf.valid);

    gas->thermo()->restoreState(initial_state);
    Gas g(*gas, GasChemistry::FROZEN);
    ShockSolver solver(g, {.abstol = 5e-5});
    ShockResult r_cant = solver.normal_shock(mach);
    ASSERT_TRUE(r_cant.valid);

    EXPECT_NEAR(r_cant.mach_out, r_perf.mach_out,
                max_fp_error(r_perf.mach_out, 1e-2, 1e-4))
        << "Cantera M2 should approximate perfect gas at low T";
    EXPECT_NEAR(r_cant.static_pressure_ratio, r_perf.static_pressure_ratio,
                max_fp_error(r_perf.static_pressure_ratio, 1e-2, 1e-4))
        << "Cantera P2/P1 should approximate perfect gas at low T";
    EXPECT_NEAR(r_cant.static_temperature_ratio, r_perf.static_temperature_ratio,
                max_fp_error(r_perf.static_temperature_ratio, 1e-2, 1e-4))
        << "Cantera T2/T1 should approximate perfect gas at low T";

    // Vibrational excitation of N2 and O2 at M = 2 changes theta_max by well under 1%
    EXPECT_NEAR(solver.max_deflection(mach), oblique_shock_max_deflection(mach, gamma),
                1e-2*oblique_shock_max_deflection(mach, gamma));
}

TEST_F(ShockSolverFrozenTests, RankineHugoniot) {
    // Verify conservation of mass, momentum, and energy across the shock.
    double mach = 2.5;
    double gamma1 = gas->thermo()->cp_mass() / gas->thermo()->cv_mass();
    double a1 = gas_sonic_velocity(*gas->thermo(), gamma1);
    double u1 = mach * a1;
    double rho1 = gas->thermo()->density();
    double P1 = gas->thermo()->pressure();
    double h1 = gas->thermo()->enthalpy_mass();

    Gas g(*gas, GasChemistry::FROZEN);
    ShockSolver solver(g, {.abstol = 5e-5});
    ShockResult r = solver.normal_shock(mach);
    ASSERT_TRUE(r.valid);

    // Access post-shock state via solver
    const Gas& post = solver.post_shock_state();
    double rho2 = post.density();
    double P2 = post.pressure();
    double h2 = post.enthalpy_mass();
    double u2 = rho1 * u1 / rho2;
    EXPECT_NEAR(r.density_ratio, rho2/rho1, max_fp_error(rho2/rho1, 1e-12, 0.0));

    // Momentum: P1 + rho1*u1^2 = P2 + rho2*u2^2
    double momentum1 = P1 + rho1 * u1 * u1;
    double momentum2 = P2 + rho2 * u2 * u2;
    EXPECT_NEAR(momentum2, momentum1,
                max_fp_error(momentum1, 1e-4, 1.0))
        << "Momentum must be conserved across shock";

    // Energy: h1 + u1^2/2 = h2 + u2^2/2
    double energy1 = h1 + 0.5 * u1 * u1;
    double energy2 = h2 + 0.5 * u2 * u2;
    EXPECT_NEAR(energy2, energy1,
                max_fp_error(energy1, 1e-4, 1.0))
        << "Total enthalpy must be conserved across shock";
}

TEST_F(ShockSolverFrozenTests, InvalidSubsonic) {
    Gas g(*gas, GasChemistry::FROZEN);
    ShockSolver solver(g, {.abstol = 5e-5});
    EXPECT_FALSE(solver.normal_shock(0.5).valid);
    EXPECT_FALSE(solver.reflected_shock(0.5).valid);
    EXPECT_FALSE(solver.oblique_shock_from_wave_angle(2.0, 0.9*std::asin(0.5)).valid);
}

TEST_F(ShockSolverFrozenTests, DetachedShockIsInvalid) {
    Gas g(*gas, GasChemistry::FROZEN);
    ShockSolver solver(g);
    const double theta_max = solver.max_deflection(3.0);
    EXPECT_FALSE(solver.oblique_shock_from_deflection(3.0, theta_max + 1e-3, true).valid);
    EXPECT_FALSE(solver.oblique_shock_from_deflection(3.0, theta_max + 1e-3, false).valid);
}

// ============================================================
//  ShockSolver with equilibrium chemistry on air
// ============================================================

class ShockSolverEquilibriumTests : public ::testing::Test {
protected:
    // At M = 8 from 300 K and 1 atm the frozen post-shock state is ~3400 K and ~75 atm, where
    // O2 is partly dissociated (h2o2.yaml carries O but no N, so N2 is inert).
    static constexpr double MACH = 8.0;
};

TEST_F(ShockSolverEquilibriumTests, ConservesMassMomentumEnergy) {
    Gas pre = make_gas("N2:0.79, O2:0.21", GasChemistry::EQUILIBRIUM);
    const double rho1 = pre.density();
    const double P1 = pre.pressure();
    const double h1 = pre.enthalpy_mass();

    ShockSolver solver(pre);
    const double u1 = MACH*get_frozen_properties(solver.pre_shock_state()).speed_of_sound;
    ShockResult r = solver.normal_shock_from_velocity(u1);
    ASSERT_TRUE(r.valid);
    EXPECT_NEAR(r.mach_in, MACH, 1e-12);

    const Gas& post = solver.post_shock_state();
    const double rho2 = post.density();
    const double u2 = u1/r.density_ratio;
    EXPECT_NEAR(r.density_ratio, rho2/rho1, max_fp_error(rho2/rho1, 1e-12, 0.0));
    EXPECT_NEAR(post.pressure() + rho2*u2*u2, P1 + rho1*u1*u1, RATIO_TOL*(rho1*u1*u1));
    EXPECT_NEAR(post.enthalpy_mass() + 0.5*u2*u2, h1 + 0.5*u1*u1, RATIO_TOL*(0.5*u1*u1));

    // The post-shock state is at equilibrium: re-equilibrating at (T2, P2) does not move it.
    Gas check = post.clone();
    check.equilibrate_TP(post.temperature(), post.pressure());
    EXPECT_NEAR(check.enthalpy_mass(), post.enthalpy_mass(), 1e-6*(0.5*u1*u1));
}

TEST_F(ShockSolverEquilibriumTests, DissociationLowersTemperatureAndRaisesDensity) {
    ShockSolver frozen(make_gas("N2:0.79, O2:0.21", GasChemistry::FROZEN));
    ShockSolver equilibrium(make_gas("N2:0.79, O2:0.21", GasChemistry::EQUILIBRIUM));
    ShockResult r_frozen = frozen.normal_shock(MACH);
    ShockResult r_equilibrium = equilibrium.normal_shock(MACH);
    ASSERT_TRUE(r_frozen.valid);
    ASSERT_TRUE(r_equilibrium.valid);
    // Dissociation absorbs energy: lower temperature (by ~3% here), higher compression. The
    // 1e-2 margin is far above the ~1e-6 solver error.
    EXPECT_LT(r_equilibrium.static_temperature_ratio, (1.0 - 1e-2)*r_frozen.static_temperature_ratio);
    EXPECT_GT(r_equilibrium.density_ratio, r_frozen.density_ratio);
    EXPECT_LT(r_equilibrium.mach_out, 1.0);
}

TEST_F(ShockSolverEquilibriumTests, ReflectedShockBringsGasToRest) {
    // In the lab frame gas 2 moves at u_p toward the wall and gas 5 is at rest. With
    // W_R = u_p/(rho5/rho2 - 1), the jump conditions across the reflected shock are
    //   P5 - P2 = rho2 (u_p + W_R) u_p,   h5 - h2 = u_p (u_p + 2 W_R)/2.
    const std::string air = "N2:0.79, O2:0.21";
    ShockSolver incident(make_gas(air, GasChemistry::EQUILIBRIUM));
    ShockSolver reflected(make_gas(air, GasChemistry::EQUILIBRIUM));
    const double mach = 5.0;

    const double u1 = mach*get_frozen_properties(reflected.pre_shock_state()).speed_of_sound;
    ShockResult r2 = incident.normal_shock(mach);
    ReflectedShockResult r = reflected.reflected_shock(mach);
    ASSERT_TRUE(r.valid);
    expect_shock_near(r.incident, r2, 1e-12, "incident");

    // incident and reflected own separate Gas objects, so both states stay valid
    const Gas& state2 = incident.post_shock_state();
    const Gas& state5 = reflected.post_shock_state();
    const double P2 = state2.pressure();
    const double rho2 = state2.density();
    const double h2 = state2.enthalpy_mass();
    EXPECT_NEAR(r.reflected.density_ratio, state5.density()/rho2, max_fp_error(state5.density()/rho2, 1e-10, 0.0));

    const double particle_velocity = u1*(1.0 - 1.0/r.incident.density_ratio);
    const double wave_speed = particle_velocity/(r.reflected.density_ratio - 1.0);
    EXPECT_NEAR(state5.pressure() - P2, rho2*(particle_velocity + wave_speed)*particle_velocity,
                RATIO_TOL*state5.pressure());
    EXPECT_NEAR(state5.enthalpy_mass() - h2, 0.5*particle_velocity*(particle_velocity + 2*wave_speed),
                RATIO_TOL*0.5*particle_velocity*(particle_velocity + 2*wave_speed));
    EXPECT_NEAR(r.reflected.mach_in*state2.speed_of_sound(), particle_velocity + wave_speed,
                RATIO_TOL*(particle_velocity + wave_speed));
}

TEST_F(ShockSolverEquilibriumTests, MixedChemistryReflectedShock) {
    // Frozen incident shock followed by an equilibrium reflected shock (CEA's incident_frozen)
    const std::string air = "N2:0.79, O2:0.21";
    ShockSolver solver(make_gas(air, GasChemistry::EQUILIBRIUM));
    ShockSolver frozen(make_gas(air, GasChemistry::FROZEN));
    const double mach = 5.0;

    ReflectedShockResult mixed = solver.reflected_shock(mach, GasChemistry::FROZEN, GasChemistry::EQUILIBRIUM);
    ReflectedShockResult all_frozen = frozen.reflected_shock(mach);
    ASSERT_TRUE(mixed.valid);
    ASSERT_TRUE(all_frozen.valid);
    expect_shock_near(mixed.incident, all_frozen.incident, 1e-12, "incident");
    // Behind the reflected shock (~4500 K) dissociation lowers the temperature
    EXPECT_LT(mixed.reflected.static_temperature_ratio, all_frozen.reflected.static_temperature_ratio);
}

TEST_F(ShockSolverEquilibriumTests, ObliqueRoundTrip) {
    // theta(beta) followed by beta(theta) recovers beta on both branches. Bisection stops at a
    // deflection residual of 1e-6 rad; d(theta)/d(beta) is O(1) away from the maximum.
    ShockSolver solver(make_gas("N2:0.79, O2:0.21", GasChemistry::EQUILIBRIUM));
    const double mach = 7.0;
    const double theta_max = solver.max_deflection(mach);
    for (double beta_deg : {20.0, 80.0}) {
        ObliqueShockResult forward = solver.oblique_shock_from_wave_angle(mach, beta_deg*DEG);
        ASSERT_TRUE(forward.valid);
        ASSERT_LT(forward.theta, theta_max);
        const bool weak = beta_deg < 50.0;
        ObliqueShockResult inverse = solver.oblique_shock_from_deflection(mach, forward.theta, weak);
        ASSERT_TRUE(inverse.valid) << "beta=" << beta_deg;
        EXPECT_NEAR(inverse.beta, forward.beta, 1e-5) << "beta=" << beta_deg;
        EXPECT_NEAR(inverse.mach_out, forward.mach_out, max_fp_error(forward.mach_out, 1e-4, 0.0));
    }
}

// ============================================================
//  ShockSolver with perfect gas chemistry
// ============================================================

TEST(ShockSolverPerfectGas, MatchesFreeFunction) {
    // ShockSolver with PERFECT_GAS should produce identical results to the free function
    Goddard::setup_defaults();
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    sol->thermo()->setState_TPX(300.0, Cantera::OneAtm, "N2:0.79, O2:0.21");

    Gas g(*sol, GasChemistry::PERFECT_GAS);
    double gamma = g.gamma_s();
    ShockSolver solver(g);

    for (double mach : {1.5, 2.0, 3.0, 5.0}) {
        expect_shock_near(solver.normal_shock(mach), normal_shock(mach, gamma), 1e-12, "normal");
        ReflectedShockResult r_solver = solver.reflected_shock(mach);
        ReflectedShockResult r_free = reflected_shock(mach, gamma);
        expect_shock_near(r_solver.reflected, r_free.reflected, 1e-12, "reflected");
    }
    EXPECT_NEAR(solver.max_deflection(3.0), oblique_shock_max_deflection(3.0, gamma), 1e-14);
}

TEST(ShockSolverPerfectGas, ObliqueMatchesFreeFunction) {
    Goddard::setup_defaults();
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    sol->thermo()->setState_TPX(300.0, Cantera::OneAtm, "N2:0.79, O2:0.21");

    Gas g(*sol, GasChemistry::PERFECT_GAS);
    double gamma = g.gamma_s();
    ShockSolver solver(g);

    double mach = 3.0;
    double theta = 15.0 * DEG;

    ObliqueShockResult r_solver = solver.oblique_shock_from_deflection(mach, theta, true);
    ObliqueShockResult r_free = oblique_shock_from_deflection(mach, theta, gamma, true);

    ASSERT_TRUE(r_solver.valid);
    ASSERT_TRUE(r_free.valid);
    EXPECT_NEAR(r_solver.beta, r_free.beta, 1e-12);
    EXPECT_NEAR(r_solver.mach_out, r_free.mach_out, 1e-12);
}

// ============================================================
//  ShockSolver state management
// ============================================================

TEST(ShockSolverState, PreShockStateRestored) {
    Goddard::setup_defaults();
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    sol->thermo()->setState_TPX(300.0, Cantera::OneAtm, "N2:0.79, O2:0.21");
    double T_orig = sol->thermo()->temperature();
    double P_orig = sol->thermo()->pressure();

    Gas g(*sol, GasChemistry::EQUILIBRIUM);
    ShockSolver solver(g, {.abstol = 5e-5});
    solver.reflected_shock(3.0, GasChemistry::FROZEN, GasChemistry::EQUILIBRIUM);

    // pre_shock_state() should restore original T and P, and the construction chemistry
    const Gas& pre = solver.pre_shock_state();
    EXPECT_NEAR(pre.temperature(), T_orig, 1e-10);
    EXPECT_NEAR(pre.pressure(), P_orig, 1e-6);
    EXPECT_EQ(pre.chemistry, GasChemistry::EQUILIBRIUM);
}

TEST(ShockSolverState, PostShockStateDiffers) {
    Goddard::setup_defaults();
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    sol->thermo()->setState_TPX(300.0, Cantera::OneAtm, "N2:0.79, O2:0.21");
    double T_orig = sol->thermo()->temperature();

    Gas g(*sol, GasChemistry::FROZEN);
    ShockSolver solver(g, {.abstol = 5e-5});
    solver.normal_shock(2.0);

    // post_shock_state() should have a different temperature
    const Gas& post = solver.post_shock_state();
    EXPECT_GT(post.temperature(), T_orig * 1.1)
        << "Post-shock temperature should be significantly higher";
}

TEST(ShockSolverState, PerfectGasPostShockStateFollowsJump) {
    // For PERFECT_GAS, post_shock_state() must be the pre-shock gas with T2 = T1*(T2/T1) and
    // P2 = P1*(P2/P1) from the returned jump, and an unchanged (frozen) composition. A reflected
    // shock stores state 5: T5/T1 = (T2/T1)(T5/T2).
    Goddard::setup_defaults();
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    sol->thermo()->setState_TPX(300.0, Cantera::OneAtm, "N2:0.79, O2:0.21");
    const double T1 = sol->thermo()->temperature();
    const double P1 = sol->thermo()->pressure();

    Gas g(*sol, GasChemistry::PERFECT_GAS);
    const std::vector<double> Y1 = g.mass_fractions();
    const double a1 = g.speed_of_sound();
    ShockSolver solver(g);

    auto expect_post_shock = [&](double temperature_ratio, double pressure_ratio,
                                 const std::string& label) {
        ASSERT_GT(temperature_ratio, 1.0) << label;
        const Gas& post = solver.post_shock_state();
        EXPECT_NEAR(post.temperature(), T1*temperature_ratio,
                    max_fp_error(T1*temperature_ratio, 1e-12, 1e-9)) << label;
        EXPECT_NEAR(post.pressure(), P1*pressure_ratio,
                    max_fp_error(P1*pressure_ratio, 1e-12, 1e-6)) << label;
        const std::vector<double> Y2 = post.mass_fractions();
        ASSERT_EQ(Y2.size(), Y1.size()) << label;
        for (size_t k = 0; k < Y1.size(); k++) {
            EXPECT_NEAR(Y2[k], Y1[k], 1e-14) << label << ": composition must stay frozen, species " << k;
        }
    };

    ShockResult normal = solver.normal_shock(3.0);
    ASSERT_TRUE(normal.valid);
    expect_post_shock(normal.static_temperature_ratio, normal.static_pressure_ratio, "normal");

    ShockResult from_velocity = solver.normal_shock_from_velocity(2.0*a1);
    ASSERT_TRUE(from_velocity.valid);
    expect_post_shock(from_velocity.static_temperature_ratio, from_velocity.static_pressure_ratio,
                      "normal from velocity");

    ReflectedShockResult reflected = solver.reflected_shock(3.0);
    ASSERT_TRUE(reflected.valid);
    expect_post_shock(
        reflected.incident.static_temperature_ratio*reflected.reflected.static_temperature_ratio,
        reflected.incident.static_pressure_ratio*reflected.reflected.static_pressure_ratio,
        "reflected (state 5)");

    ObliqueShockResult from_wave_angle = solver.oblique_shock_from_wave_angle(3.0, 40.0*DEG);
    ASSERT_TRUE(from_wave_angle.valid);
    expect_post_shock(from_wave_angle.shock.static_temperature_ratio,
                      from_wave_angle.shock.static_pressure_ratio, "oblique from wave angle");

    ObliqueShockResult from_deflection = solver.oblique_shock_from_deflection(3.0, 15.0*DEG, true);
    ASSERT_TRUE(from_deflection.valid);
    expect_post_shock(from_deflection.shock.static_temperature_ratio,
                      from_deflection.shock.static_pressure_ratio, "oblique from deflection");
}

// ============================================================
//  ShockSolver unsupported chemistry
// ============================================================

TEST(ShockSolverUnsupported, KineticRejected) {
    EXPECT_THROW(ShockSolver(make_gas("N2:0.79, O2:0.21", GasChemistry::KINETIC)), std::invalid_argument);
}

TEST(ShockSolverUnsupported, PerShockChemistryMustBeFrozenOrEquilibrium) {
    ShockSolver solver(make_gas("N2:0.79, O2:0.21", GasChemistry::FROZEN));
    EXPECT_THROW(solver.reflected_shock(3.0, GasChemistry::KINETIC, GasChemistry::FROZEN), std::invalid_argument);
    EXPECT_THROW(solver.reflected_shock(3.0, GasChemistry::FROZEN, GasChemistry::PERFECT_GAS), std::invalid_argument);

    ShockSolver perfect(make_gas("N2:0.79, O2:0.21", GasChemistry::PERFECT_GAS));
    EXPECT_THROW(perfect.reflected_shock(3.0, GasChemistry::FROZEN, GasChemistry::FROZEN), std::invalid_argument);
}
