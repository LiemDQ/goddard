#include "goddard/detonations.hpp"
#include "goddard/shocks.hpp"
#include "goddard/gas.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/global.hpp"
#include "goddard/numerics.hpp"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include "gtest/gtest.h"

using namespace Goddard;

namespace {

// The Newton iterations stop once the step in ln(P2/P1) and ln(T2/T1) is below
// SolverOptions::abstol = 1e-6, so ratios carry relative errors of ~1e-6. Comparisons use 1e-5.
constexpr double RATIO_TOL = 1e-5;

const DetonationBranch BRANCHES[] = {DetonationBranch::OVERDRIVEN, DetonationBranch::UNDERDRIVEN};

const char* branch_name(DetonationBranch branch) {
    return branch == DetonationBranch::OVERDRIVEN ? "OVERDRIVEN" : "UNDERDRIVEN";
}

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

/** Compare the jump across a real-gas detonation with a perfect-gas one; `velocity` is not compared. */
void expect_detonation_near(const DetonationResult& actual, const DetonationResult& expected, double tol,
                            const std::string& label) {
    ASSERT_TRUE(actual.valid) << label;
    ASSERT_TRUE(expected.valid) << label;
    EXPECT_NEAR(actual.drive_factor, expected.drive_factor, max_fp_error(expected.drive_factor, tol, 0.0)) << label;
    EXPECT_NEAR(actual.mach_in, expected.mach_in, max_fp_error(expected.mach_in, tol, 0.0)) << label;
    EXPECT_NEAR(actual.mach_out, expected.mach_out, max_fp_error(expected.mach_out, tol, 0.0)) << label;
    EXPECT_NEAR(actual.static_pressure_ratio, expected.static_pressure_ratio,
                max_fp_error(expected.static_pressure_ratio, tol, 0.0)) << label;
    EXPECT_NEAR(actual.static_temperature_ratio, expected.static_temperature_ratio,
                max_fp_error(expected.static_temperature_ratio, tol, 0.0)) << label;
    EXPECT_NEAR(actual.density_ratio, expected.density_ratio,
                max_fp_error(expected.density_ratio, tol, 0.0)) << label;
    EXPECT_NEAR(actual.molecular_weight_ratio, expected.molecular_weight_ratio,
                max_fp_error(expected.molecular_weight_ratio, tol, 0.0)) << label;
    EXPECT_NEAR(actual.total_pressure_ratio, expected.total_pressure_ratio,
                max_fp_error(expected.total_pressure_ratio, tol, 0.0)) << label;
    expect_shock_near(actual.von_neumann, expected.von_neumann, tol, label + " von Neumann");
}

/**
 * Residual of the one-gamma Hugoniot with heat release, divided by P1 v1:
 *   gamma/(gamma - 1) (P2 v2/(P1 v1) - 1) - (P2/P1 - 1)(1 + v2/v1)/2 - q/(P1 v1).
 */
double hugoniot_residual(const DetonationResult& r, double gamma, double heat_release) {
    const double volume_ratio = 1.0/r.density_ratio;
    return gamma/(gamma - 1.0)*(r.static_pressure_ratio*volume_ratio - 1.0)
        - 0.5*(r.static_pressure_ratio - 1.0)*(1.0 + volume_ratio) - heat_release;
}

/** Gas from h2o2.yaml at T [K], P [Pa] and the given mole fractions. */
Gas make_gas(const std::string& composition, GasChemistry chemistry,
             double T = 298.15, double P = Cantera::OneAtm) {
    Goddard::setup_defaults();
    auto sol = Cantera::newSolution("h2o2.yaml", "ohmech");
    sol->thermo()->setState_TPX(T, P, composition);
    return Gas(sol, chemistry);
}

// One-gamma reacting mixture: reactant R and product P are isomers (one argon atom each, so the
// same molar mass) with the same constant cp. R carries the heat of reaction and P a large
// entropy, so equilibrium is all P. It is exactly the calorically perfect gas of the free
// functions. A made-up element would be simpler, but Cantera's Solution::clone(), which
// Gas::stagnation_pressure() uses, loses custom element definitions.
constexpr double GAMMA = 1.4;
constexpr double HEAT_RELEASE = 40.0;   // q/(R T1)
constexpr double T1 = 300.0;            // K
constexpr double P1 = 1e5;              // Pa

Gas make_one_gamma_gas(GasChemistry chemistry) {
    Goddard::setup_defaults();
    const double R = Cantera::GasConstant;
    const double cp = GAMMA/(GAMMA - 1.0)*R;                 // J/(kmol K)
    const double heat_of_reaction = HEAT_RELEASE*R*T1;        // q M, J/kmol
    std::ostringstream yaml;
    yaml << std::setprecision(17)
         << "units: {length: m, quantity: kmol, energy: J}\n"
         << "phases:\n"
         << "- name: one-gamma\n"
         << "  thermo: ideal-gas\n"
         << "  elements: [Ar]\n"
         << "  species: [R, P]\n"
         << "  state: {T: " << T1 << ", P: " << P1 << ", X: {R: 1.0}}\n"
         << "species:\n"
         << "- name: R\n"
         << "  composition: {Ar: 1}\n"
         << "  thermo: {model: constant-cp, T0: 298.15, h0: " << heat_of_reaction
         << ", s0: 0.0, cp0: " << cp << ", T-max: 100000.0}\n"
         << "- name: P\n"
         << "  composition: {Ar: 1}\n"
         << "  thermo: {model: constant-cp, T0: 298.15, h0: 0.0, s0: " << 30.0*R
         << ", cp0: " << cp << ", T-max: 100000.0}\n";
    Cantera::AnyMap root = Cantera::AnyMap::fromYamlString(yaml.str());
    Cantera::AnyMap phase = root["phases"].asVector<Cantera::AnyMap>()[0];
    return Gas(Cantera::newSolution(phase, root), chemistry);
}

} // namespace

// ============================================================
//  Perfect-gas free functions
// ============================================================

TEST(PerfectGasDetonation, ChapmanJouguetIsSonicAndOnTheHugoniot) {
    for (double gamma : {1.2, 1.4}) {
        for (double Q : {5.0, 40.0}) {
            const std::string label = "gamma=" + std::to_string(gamma) + " Q=" + std::to_string(Q);
            DetonationResult cj = chapman_jouguet_detonation(gamma, Q);
            ASSERT_TRUE(cj.valid) << label;
            // M_CJ = sqrt(H) + sqrt(H + 1), H = (gamma^2 - 1) Q/(2 gamma)
            const double H = (gamma*gamma - 1.0)*Q/(2.0*gamma);
            EXPECT_NEAR(cj.mach_in, std::sqrt(H) + std::sqrt(H + 1.0), 1e-12) << label;
            EXPECT_NEAR(cj.mach_out, 1.0, 1e-12) << label;
            EXPECT_NEAR(hugoniot_residual(cj, gamma, Q), 0.0, 1e-10*Q) << label;
            // Rayleigh line: P2/P1 - 1 = gamma M1^2 (1 - rho1/rho2)
            EXPECT_NEAR(cj.static_pressure_ratio - 1.0,
                        gamma*cj.mach_in*cj.mach_in*(1.0 - 1.0/cj.density_ratio), 1e-10*cj.static_pressure_ratio)
                << label;
            expect_shock_near(cj.von_neumann, normal_shock(cj.mach_in, gamma), 1e-14, label);
        }
    }
}

TEST(PerfectGasDetonation, BranchesSatisfyJumpConditions) {
    // Both roots lie on the Hugoniot and on the Rayleigh line of the same wave speed; the
    // overdriven root has subsonic products and the under-driven root supersonic ones.
    const double gamma = 1.3;
    const double Q = 20.0;
    const DetonationResult cj = chapman_jouguet_detonation(gamma, Q);
    for (double f : {1.001, 1.3, 3.0}) {
        const DetonationResult strong = detonation(f, gamma, Q, DetonationBranch::OVERDRIVEN);
        const DetonationResult weak = detonation(f, gamma, Q, DetonationBranch::UNDERDRIVEN);
        const std::string label = "f=" + std::to_string(f);
        ASSERT_TRUE(strong.valid) << label;
        ASSERT_TRUE(weak.valid) << label;
        for (const DetonationResult& r : {strong, weak}) {
            EXPECT_NEAR(r.mach_in, f*cj.mach_in, 1e-12*r.mach_in) << label;
            EXPECT_NEAR(hugoniot_residual(r, gamma, Q), 0.0, 1e-10*Q) << label;
            EXPECT_NEAR(r.static_pressure_ratio - 1.0,
                        gamma*r.mach_in*r.mach_in*(1.0 - 1.0/r.density_ratio), 1e-10*r.static_pressure_ratio)
                << label;
        }
        EXPECT_LT(strong.mach_out, 1.0) << label;
        EXPECT_GT(weak.mach_out, 1.0) << label;
        EXPECT_GT(strong.static_pressure_ratio, cj.static_pressure_ratio) << label;
        EXPECT_LT(weak.static_pressure_ratio, cj.static_pressure_ratio) << label;
    }
}

TEST(PerfectGasDetonation, NoHeatReleaseIsNormalShock) {
    // With Q = 0 the CJ speed is the sound speed: the strong root is a normal shock at M = f
    // and the weak root a wave that changes nothing.
    const double gamma = 1.4;
    for (double f : {1.5, 3.0}) {
        const DetonationResult strong = detonation(f, gamma, 0.0, DetonationBranch::OVERDRIVEN);
        const ShockResult shock = normal_shock(f, gamma);
        ASSERT_TRUE(strong.valid);
        EXPECT_NEAR(strong.static_pressure_ratio, shock.static_pressure_ratio, 1e-12*shock.static_pressure_ratio);
        EXPECT_NEAR(strong.static_temperature_ratio, shock.static_temperature_ratio,
                    1e-12*shock.static_temperature_ratio);
        EXPECT_NEAR(strong.density_ratio, shock.density_ratio, 1e-12*shock.density_ratio);
        EXPECT_NEAR(strong.mach_out, shock.mach_out, 1e-12);
        EXPECT_NEAR(strong.total_pressure_ratio, shock.total_pressure_ratio, 1e-12);

        const DetonationResult weak = detonation(f, gamma, 0.0, DetonationBranch::UNDERDRIVEN);
        ASSERT_TRUE(weak.valid);
        EXPECT_NEAR(weak.static_pressure_ratio, 1.0, 1e-12);
        EXPECT_NEAR(weak.density_ratio, 1.0, 1e-12);
    }
}

TEST(PerfectGasDetonation, ReflectedShockBringsProductsToRest) {
    // In the reflected-shock frame the products enter at M_R a2 and leave at W_R = M_R a2 rho2/rho5,
    // so the lab-frame velocity they lose is u_p = M_R a2 (1 - rho2/rho5).
    const double gamma = 1.4;
    const double Q = 40.0;
    for (double f : {1.0, 1.5}) {
        for (DetonationBranch branch : BRANCHES) {
            const ReflectedDetonationResult r = reflected_detonation(f, gamma, Q, branch);
            ASSERT_TRUE(r.valid);
            const double particle_mach = r.incident.mach_out*(r.incident.density_ratio - 1.0);
            EXPECT_NEAR(particle_mach, r.reflected.mach_in*(1.0 - 1.0/r.reflected.density_ratio), 1e-12)
                << "f=" << f << " " << branch_name(branch);
        }
    }
}

TEST(PerfectGasDetonation, InvalidInputs) {
    EXPECT_FALSE(detonation(0.99, 1.4, 10.0).valid);
    EXPECT_FALSE(chapman_jouguet_detonation(1.0, 10.0).valid);
    EXPECT_FALSE(chapman_jouguet_detonation(1.4, -1.0).valid);
    EXPECT_FALSE(reflected_detonation(0.99, 1.4, 10.0).valid);
    // A drive factor of 1 up to round-off is the CJ detonation
    EXPECT_TRUE(detonation(1.0 - 1e-14, 1.4, 10.0, DetonationBranch::UNDERDRIVEN).valid);
}

// ============================================================
//  DetonationSolver on the one-gamma mixture: must reproduce the free functions
// ============================================================

class OneGammaDetonationTests : public ::testing::TestWithParam<GasChemistry> {};

TEST_P(OneGammaDetonationTests, ChapmanJouguetMatchesPerfectGas) {
    DetonationSolver solver(make_one_gamma_gas(GetParam()));
    const DetonationResult cj = solver.chapman_jouguet();
    expect_detonation_near(cj, chapman_jouguet_detonation(GAMMA, HEAT_RELEASE), RATIO_TOL, "CJ");
    const double a1 = std::sqrt(GAMMA*Cantera::GasConstant*T1/solver.pre_detonation_state().molecular_weight());
    EXPECT_NEAR(cj.velocity, cj.mach_in*a1, max_fp_error(cj.velocity, 1e-12, 0.0));
}

TEST_P(OneGammaDetonationTests, DrivenDetonationsMatchPerfectGas) {
    DetonationSolver solver(make_one_gamma_gas(GetParam()));
    for (double f : {1.01, 1.2, 2.0}) {
        for (DetonationBranch branch : BRANCHES) {
            expect_detonation_near(solver.detonation(f, branch), detonation(f, GAMMA, HEAT_RELEASE, branch),
                                   RATIO_TOL, "f=" + std::to_string(f) + " " + branch_name(branch));
        }
    }
}

TEST_P(OneGammaDetonationTests, ReflectedDetonationMatchesPerfectGas) {
    DetonationSolver solver(make_one_gamma_gas(GetParam()));
    for (double f : {1.0, 1.5}) {
        for (DetonationBranch branch : BRANCHES) {
            const std::string label = "f=" + std::to_string(f) + " " + branch_name(branch);
            const ReflectedDetonationResult actual = solver.reflected_detonation(f, branch);
            const ReflectedDetonationResult expected = reflected_detonation(f, GAMMA, HEAT_RELEASE, branch);
            ASSERT_TRUE(actual.valid) << label;
            expect_detonation_near(actual.incident, expected.incident, RATIO_TOL, label);
            expect_shock_near(actual.reflected, expected.reflected, RATIO_TOL, label + " reflected");
        }
    }
}

INSTANTIATE_TEST_SUITE_P(Chemistry, OneGammaDetonationTests,
    ::testing::Values(GasChemistry::FROZEN, GasChemistry::EQUILIBRIUM),
    [](const ::testing::TestParamInfo<GasChemistry>& param_info) {
        return std::string(param_info.param == GasChemistry::EQUILIBRIUM ? "EQUILIBRIUM" : "FROZEN");
    });

// ============================================================
//  DetonationSolver on stoichiometric H2/O2
// ============================================================

class HydrogenOxygenDetonationTests : public ::testing::Test {
protected:
    HydrogenOxygenDetonationTests() : solver(make_gas(MIXTURE, GasChemistry::FROZEN)) {
        const Gas& pre = solver.pre_detonation_state();
        rho1 = pre.density();
        P1 = pre.pressure();
        h1 = pre.enthalpy_mass();
    }

    static constexpr const char* MIXTURE = "H2:2, O2:1";
    DetonationSolver solver;
    double rho1;
    double P1;
    double h1;
};

TEST_F(HydrogenOxygenDetonationTests, ChapmanJouguetConservesMassMomentumEnergy) {
    const DetonationResult cj = solver.chapman_jouguet();
    ASSERT_TRUE(cj.valid);
    const double u1 = cj.velocity;
    EXPECT_NEAR(cj.mach_in*get_frozen_properties(solver.pre_detonation_state()).speed_of_sound, u1, 1e-9*u1);

    const Gas& post = solver.post_detonation_state();
    const double rho2 = post.density();
    const double u2 = u1/cj.density_ratio;
    EXPECT_NEAR(cj.density_ratio, rho2/rho1, max_fp_error(rho2/rho1, 1e-12, 0.0));
    EXPECT_NEAR(cj.static_pressure_ratio, post.pressure()/P1, max_fp_error(post.pressure()/P1, 1e-12, 0.0));
    EXPECT_NEAR(post.pressure() + rho2*u2*u2, P1 + rho1*u1*u1, RATIO_TOL*(rho1*u1*u1));
    EXPECT_NEAR(post.enthalpy_mass() + 0.5*u2*u2, h1 + 0.5*u1*u1, RATIO_TOL*(0.5*u1*u1));
    // Jouguet condition: the products leave at their equilibrium sound speed
    EXPECT_NEAR(u2, post.speed_of_sound(), RATIO_TOL*u2);

    // The products are at equilibrium: re-equilibrating at (T2, P2) does not move them.
    Gas check = post.clone();
    check.equilibrate_TP(post.temperature(), post.pressure());
    EXPECT_NEAR(check.enthalpy_mass(), post.enthalpy_mass(), 1e-6*(0.5*u1*u1));
}

TEST_F(HydrogenOxygenDetonationTests, ChapmanJouguetIsEntropyMinimumOnHugoniot) {
    // Along the Hugoniot, entropy is stationary where the Rayleigh line is tangent (Jouguet's
    // rule), and the CJ point is its minimum: both roots at a faster wave speed have more entropy.
    solver.chapman_jouguet();
    const double s_cj = solver.post_detonation_state().entropy_mass();
    for (DetonationBranch branch : BRANCHES) {
        solver.detonation(1.1, branch);
        EXPECT_GT(solver.post_detonation_state().entropy_mass(), s_cj) << branch_name(branch);
    }
}

TEST_F(HydrogenOxygenDetonationTests, BranchesBracketChapmanJouguet) {
    const DetonationResult cj = solver.chapman_jouguet();
    for (double f : {1.01, 1.5}) {
        const DetonationResult strong = solver.detonation(f, DetonationBranch::OVERDRIVEN);
        const DetonationResult weak = solver.detonation(f, DetonationBranch::UNDERDRIVEN);
        ASSERT_TRUE(strong.valid);
        ASSERT_TRUE(weak.valid);
        EXPECT_LT(strong.mach_out, 1.0) << "f=" << f;
        EXPECT_GT(weak.mach_out, 1.0) << "f=" << f;
        EXPECT_GT(strong.static_pressure_ratio, cj.static_pressure_ratio) << "f=" << f;
        EXPECT_LT(weak.static_pressure_ratio, cj.static_pressure_ratio) << "f=" << f;
        EXPECT_NEAR(strong.velocity, f*cj.velocity, 1e-12*strong.velocity);
    }
}

TEST_F(HydrogenOxygenDetonationTests, DrivenDetonationConservesMassMomentumEnergy) {
    for (DetonationBranch branch : BRANCHES) {
        const DetonationResult r = solver.detonation(1.3, branch);
        ASSERT_TRUE(r.valid) << branch_name(branch);
        const Gas& post = solver.post_detonation_state();
        const double u1 = r.velocity;
        const double rho2 = post.density();
        const double u2 = u1/r.density_ratio;
        EXPECT_NEAR(post.pressure() + rho2*u2*u2, P1 + rho1*u1*u1, RATIO_TOL*(rho1*u1*u1)) << branch_name(branch);
        EXPECT_NEAR(post.enthalpy_mass() + 0.5*u2*u2, h1 + 0.5*u1*u1, RATIO_TOL*(0.5*u1*u1)) << branch_name(branch);
        EXPECT_NEAR(r.mach_out, u2/post.speed_of_sound(), 1e-12*r.mach_out) << branch_name(branch);
    }
}

TEST_F(HydrogenOxygenDetonationTests, VonNeumannStateIsFrozenShock) {
    const DetonationResult cj = solver.chapman_jouguet();
    ShockSolver frozen(make_gas(MIXTURE, GasChemistry::FROZEN));
    expect_shock_near(cj.von_neumann, frozen.normal_shock_from_velocity(cj.velocity), 1e-12, "von Neumann");
    // The spike: the leading shock compresses more than the whole detonation
    EXPECT_GT(cj.von_neumann.static_pressure_ratio, cj.static_pressure_ratio);
    const Gas& spike = solver.von_neumann_state();
    EXPECT_EQ(spike.chemistry, GasChemistry::FROZEN);
    EXPECT_NEAR(spike.pressure(), cj.von_neumann.static_pressure_ratio*P1, 1e-10*spike.pressure());
}

TEST_F(HydrogenOxygenDetonationTests, ReflectedShockBringsProductsToRest) {
    // As for ShockSolver: with W_R = u_p/(rho5/rho2 - 1),
    //   P5 - P2 = rho2 (u_p + W_R) u_p,   h5 - h2 = u_p (u_p + 2 W_R)/2.
    const DetonationResult cj = solver.chapman_jouguet();
    const Gas& state2 = solver.post_detonation_state();
    const double P2 = state2.pressure();
    const double rho2 = state2.density();
    const double h2 = state2.enthalpy_mass();
    const double a2 = state2.speed_of_sound();

    const ReflectedDetonationResult r = solver.reflected_detonation();
    ASSERT_TRUE(r.valid);
    EXPECT_NEAR(r.incident.velocity, cj.velocity, 1e-12*cj.velocity);
    const Gas& state5 = solver.post_detonation_state();
    EXPECT_NEAR(r.reflected.density_ratio, state5.density()/rho2, max_fp_error(state5.density()/rho2, 1e-10, 0.0));

    const double particle_velocity = cj.velocity*(1.0 - 1.0/cj.density_ratio);
    const double wave_speed = particle_velocity/(r.reflected.density_ratio - 1.0);
    EXPECT_NEAR(state5.pressure() - P2, rho2*(particle_velocity + wave_speed)*particle_velocity,
                RATIO_TOL*state5.pressure());
    EXPECT_NEAR(state5.enthalpy_mass() - h2, 0.5*particle_velocity*(particle_velocity + 2*wave_speed),
                RATIO_TOL*0.5*particle_velocity*(particle_velocity + 2*wave_speed));
    EXPECT_NEAR(r.reflected.mach_in*a2, particle_velocity + wave_speed, RATIO_TOL*(particle_velocity + wave_speed));
}

TEST_F(HydrogenOxygenDetonationTests, UnburnedChemistryDoesNotChangeResult) {
    // The unburned gas is taken as given, and the products are always at equilibrium.
    DetonationSolver equilibrium(make_gas(MIXTURE, GasChemistry::EQUILIBRIUM));
    expect_detonation_near(equilibrium.chapman_jouguet(), solver.chapman_jouguet(), 1e-12, "CJ");
}

TEST_F(HydrogenOxygenDetonationTests, SlowerThanChapmanJouguetIsInvalid) {
    const DetonationResult cj = solver.chapman_jouguet();
    EXPECT_FALSE(solver.detonation(0.99).valid);
    EXPECT_FALSE(solver.detonation_from_velocity(0.99*cj.velocity, DetonationBranch::UNDERDRIVEN).valid);
    EXPECT_FALSE(solver.reflected_detonation(0.99).valid);
    // After an invalid result the post state is the unburned gas
    EXPECT_NEAR(solver.post_detonation_state().pressure(), P1, 1e-10*P1);
}

// ============================================================
//  DetonationSolver on the RP-1311 example 7 mixture (5% H2, 5% O2, 90% Ar, 10 mmHg)
// ============================================================

class DiluteDetonationTests : public ::testing::Test {
protected:
    static Gas make_dilute_gas(GasChemistry chemistry) {
        return make_gas("H2:0.05, O2:0.05, AR:0.9", chemistry, 300.0, 10.0*133.322387415);
    }
};

TEST_F(DiluteDetonationTests, OverdrivenMatchesEquilibriumShock) {
    // Above the CJ speed (~1094 m/s) an equilibrium normal shock is an overdriven detonation;
    // ShockSolver reaches it from the frozen shock, DetonationSolver from the CJ-fitted guess.
    DetonationSolver detonations(make_dilute_gas(GasChemistry::FROZEN));
    ShockSolver shocks(make_dilute_gas(GasChemistry::EQUILIBRIUM));
    ASSERT_LT(detonations.chapman_jouguet().velocity, 1100.0);
    for (double u1 : {1100.0, 1200.0, 1400.0}) {
        const DetonationResult r = detonations.detonation_from_velocity(u1, DetonationBranch::OVERDRIVEN);
        const ShockResult shock = shocks.normal_shock_from_velocity(u1);
        ASSERT_TRUE(r.valid) << "u1=" << u1;
        ShockResult as_shock = r.von_neumann;
        as_shock.mach_in = r.mach_in;
        as_shock.mach_out = r.mach_out;
        as_shock.static_pressure_ratio = r.static_pressure_ratio;
        as_shock.static_temperature_ratio = r.static_temperature_ratio;
        as_shock.density_ratio = r.density_ratio;
        as_shock.total_pressure_ratio = r.total_pressure_ratio;
        expect_shock_near(as_shock, shock, RATIO_TOL, "u1=" + std::to_string(u1));
    }
}

TEST_F(DiluteDetonationTests, NearChapmanJouguetNeverReturnsWrongBranch) {
    // Close to the CJ speed the two roots merge; the solver may fail to converge, but it must
    // not return a result from the other branch.
    DetonationSolver solver(make_dilute_gas(GasChemistry::FROZEN));
    for (double excess : {1e-8, 1e-6, 1e-4}) {
        for (DetonationBranch branch : BRANCHES) {
            try {
                const DetonationResult r = solver.detonation(1.0 + excess, branch);
                ASSERT_TRUE(r.valid);
                if (branch == DetonationBranch::OVERDRIVEN) {
                    EXPECT_LT(r.mach_out, 1.0) << "excess=" << excess;
                }
                else {
                    EXPECT_GT(r.mach_out, 1.0) << "excess=" << excess;
                }
            }
            catch (const ConvergenceError&) {
                SUCCEED();
            }
        }
    }
}

// ============================================================
//  DetonationSolver: unsupported inputs and state management
// ============================================================

TEST(DetonationSolverUnsupported, RejectsPerfectGasAndKinetic) {
    EXPECT_THROW(DetonationSolver(make_gas("H2:2, O2:1", GasChemistry::KINETIC)), std::invalid_argument);
    EXPECT_THROW(DetonationSolver(make_gas("H2:2, O2:1", GasChemistry::PERFECT_GAS)), std::invalid_argument);
}

TEST(DetonationSolverUnsupported, MixtureWithoutHeatReleaseHasNoDetonation) {
    DetonationSolver solver(make_gas("N2:0.79, O2:0.21", GasChemistry::FROZEN));
    EXPECT_FALSE(solver.chapman_jouguet().valid);
    EXPECT_FALSE(solver.detonation(1.5).valid);
    EXPECT_FALSE(solver.reflected_detonation().valid);
}

TEST(DetonationSolverState, AccessorsRestoreStatesAndChemistry) {
    Gas gas = make_gas("H2:2, O2:1", GasChemistry::FROZEN);
    const double T_initial = gas.temperature();
    const double P_initial = gas.pressure();
    DetonationSolver solver(gas);
    const DetonationResult cj = solver.chapman_jouguet();
    const DetonationResult overdriven = solver.detonation(1.5);
    ASSERT_TRUE(overdriven.valid);

    const Gas& pre = solver.pre_detonation_state();
    EXPECT_NEAR(pre.temperature(), T_initial, 1e-10);
    EXPECT_NEAR(pre.pressure(), P_initial, 1e-6);
    EXPECT_EQ(pre.chemistry, GasChemistry::FROZEN);

    const Gas& post = solver.post_detonation_state();
    EXPECT_EQ(post.chemistry, GasChemistry::EQUILIBRIUM);
    EXPECT_NEAR(post.temperature(), overdriven.static_temperature_ratio*T_initial, 1e-8*post.temperature());

    const Gas& spike = solver.von_neumann_state();
    EXPECT_EQ(spike.chemistry, GasChemistry::FROZEN);
    EXPECT_NEAR(spike.temperature(), overdriven.von_neumann.static_temperature_ratio*T_initial,
                1e-8*spike.temperature());

    // The CJ result is solved once and is not changed by later solves
    expect_detonation_near(solver.chapman_jouguet(), cj, 0.0, "CJ");
    EXPECT_NEAR(solver.post_detonation_state().temperature(), cj.static_temperature_ratio*T_initial,
                1e-8*cj.static_temperature_ratio*T_initial);
}
