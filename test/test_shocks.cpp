#include "goddard/shocks.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/error.hpp"
#include "goddard/global.hpp"
#include "goddard/numerics.hpp"
#include <cmath>
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
    // TODO: Fill in expected values from Anderson's Modern Compressible Flow
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
//  Cantera-based normal shock
// ============================================================

class CantNormalShockTests : public ::testing::Test {
protected:
    CantNormalShockTests() {
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

TEST_F(CantNormalShockTests, Converges) {
    double mach = 2.0;
    ShockResult r = normal_shock(*gas->thermo(), mach);
    ASSERT_TRUE(r.valid) << "Cantera normal shock should converge at M=2";
    EXPECT_GT(r.mach_out, 0.0);
    EXPECT_LT(r.mach_out, 1.0);
}

TEST_F(CantNormalShockTests, ApproachesPerfectGas) {
    // At 300 K where air behaves as a perfect gas with gamma ~ 1.4,
    // the Cantera result should be close to the perfect-gas solution
    double gamma = gas->thermo()->cp_mass() / gas->thermo()->cv_mass();
    double mach = 2.0;

    ShockResult r_perf = normal_shock(mach, gamma);
    ASSERT_TRUE(r_perf.valid);

    gas->thermo()->restoreState(initial_state);
    ShockResult r_cant = normal_shock(*gas->thermo(), mach);
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
}

TEST_F(CantNormalShockTests, RankineHugoniot) {
    // TODO: Depends on gas_stagnation_pressure convergence fix.
    // Verify conservation of mass, momentum, and energy across the shock.
    // Pre-shock state
    double mach = 2.5;
    double gamma1 = gas->thermo()->cp_mass() / gas->thermo()->cv_mass();
    double a1 = gas_sonic_velocity(*gas->thermo(), gamma1);
    double u1 = mach * a1;
    double rho1 = gas->thermo()->density();
    double P1 = gas->thermo()->pressure();
    double h1 = gas->thermo()->enthalpy_mass();

    ShockResult r = normal_shock(*gas->thermo(), mach);
    ASSERT_TRUE(r.valid);

    // Post-shock state is left in gas->thermo() by the solver
    double rho2 = gas->thermo()->density();
    double P2 = gas->thermo()->pressure();
    double h2 = gas->thermo()->enthalpy_mass();
    // Back out u2 from mass conservation: rho1*u1 = rho2*u2
    double u2 = rho1 * u1 / rho2;

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

TEST_F(CantNormalShockTests, InvalidSubsonic) {
    ShockResult r = normal_shock(*gas->thermo(), 0.5);
    EXPECT_FALSE(r.valid);
}
