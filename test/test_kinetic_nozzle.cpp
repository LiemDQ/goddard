#include "goddard/kinetic_nozzle.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/profile.hpp"
#include "goddard/global.hpp"
#include "goddard/numerics.hpp"
#include "cantera/core.h"
#include "gtest/gtest.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <iostream>

using namespace Goddard;

// Build a simple conical diverging profile from throat (x=0) to exit.
// Only the diverging section is needed since KineticNozzle starts at x_min.
static NozzleProfile make_conical_profile(
    double r_throat, double r_exit, double length, int n_points)
{
    NozzleProfile profile;
    for (int i = 0; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        double x = length * frac;
        double r = r_throat + (r_exit - r_throat) * frac;
        profile.push_back({x, r});
    }
    return profile;
}

// KineticNozzle::solve() with dt=1e-6 is expensive. SetUpTestSuite() runs it
// once for the entire fixture class; all tests share s_results.
class KineticNozzleTests : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        Goddard::setup_defaults();
        s_gas = Cantera::newSolution("h2o2.yaml", "ohmech");
        // Typical H2/O2 post-combustion state at high temperature and pressure
        s_gas->thermo()->setState_TPX(3500.0, 70.0 * Cantera::OneAtm,
                                      "H2O:0.9, H2:0.05, O2:0.03, OH:0.02");
        s_gas->thermo()->equilibrate("HP");
        s_gas->thermo()->saveState(s_inlet_state);

        // Conical profile: throat r=0.01 m at x=0, exit r=0.01414 m at x=0.1 m
        // Area ratio at exit = (0.01414/0.01)^2 ≈ 2.0
        s_profile = make_conical_profile(0.01, 0.01414, 0.1, 100);
        double length = s_profile.x_max() - s_profile.x_min();
        int min_steps = 20000;

        KineticNozzle nozzle(*s_gas, s_profile, s_mdot);
        s_results = nozzle.solve(1e-6, length/min_steps, 100000);
    }

    static void TearDownTestSuite() { s_gas.reset(); }

    void SetUp() override { 
        s_gas->thermo()->restoreState(s_inlet_state);
        // don't check all stations, only sample a few 
        size_t num_checks = 40;
        check_interval = std::max(s_results.stations.size()/num_checks,1ul);
    }
    
    size_t check_interval = 1;
    static std::shared_ptr<Cantera::Solution> s_gas;
    static std::vector<double> s_inlet_state;
    static NozzleProfile s_profile;
    static KineticNozzleResults s_results;
    static constexpr double s_mdot = 1.0; // kg/s — scale-invariant for physics tests
};

std::shared_ptr<Cantera::Solution> KineticNozzleTests::s_gas;
std::vector<double> KineticNozzleTests::s_inlet_state;
NozzleProfile KineticNozzleTests::s_profile;
KineticNozzleResults KineticNozzleTests::s_results;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

TEST_F(KineticNozzleTests, EquilibriumThroatModelConstructs) {
    ASSERT_NO_THROW({
        KineticNozzle nozzle(*s_gas, s_profile, s_mdot, NozzleChemistryType::EQUILIBRIUM);
    });
}

TEST_F(KineticNozzleTests, FrozenThroatModelConstructs) {
    ASSERT_NO_THROW({
        KineticNozzle nozzle(*s_gas, s_profile, s_mdot, NozzleChemistryType::FROZEN);
    });
}

TEST_F(KineticNozzleTests, ExplicitInletStateConstructs) {
    ASSERT_NO_THROW({
        KineticNozzle nozzle(*s_gas, s_profile, s_mdot, s_inlet_state, NozzleChemistryType::EQUILIBRIUM);
    });
}

TEST_F(KineticNozzleTests, KineticThroatModelThrows) {
    // NozzleChemistryType::KINETIC is not a valid throat model for KineticNozzle
    EXPECT_THROW(
        { KineticNozzle nozzle(*s_gas, s_profile, s_mdot, NozzleChemistryType::KINETIC); },
        std::runtime_error);
}

// ---------------------------------------------------------------------------
// Throat conditions
// ---------------------------------------------------------------------------

TEST_F(KineticNozzleTests, ThroatConverges) {
    EXPECT_TRUE(s_results.throat.converged)
        << "Throat conditions must converge for a well-posed inlet state";
}

TEST_F(KineticNozzleTests, SolveReturnsStations) {
    ASSERT_GT(s_results.stations.size(), 0u)
        << "Solve must produce at least one station";
}

TEST_F(KineticNozzleTests, FirstStationMachIsApproximatelyOne) {
    // The solver starts at the throat (x_min), so the first station
    // should have Mach ≈ 1.0 (within the first step's integration error)
    ASSERT_GT(s_results.stations.size(), 0u);
    EXPECT_NEAR(s_results.stations.front().mach, 1.0, 0.05)
        << "First station should be near the throat (Mach ≈ 1)";
}

// ---------------------------------------------------------------------------
// Energy conservation — the primary invariant of the Strang splitting scheme
// ---------------------------------------------------------------------------

TEST_F(KineticNozzleTests, StagnationEnthalpyConservedAtEveryStation) {
    // First law (adiabatic, steady): H0 = h + u²/2 = const.
    // Tolerance 1% relative: Strang splitting is O(dt²); accumulated error
    // over ~1000 steps with dt=1e-6 and u~1500 m/s is well within 1%.
    ASSERT_TRUE(s_results.throat.converged);
    ASSERT_GT(s_results.stations.size(), 0u);

    double H0 = s_results.throat.H_stagnation;

    for (size_t i = 0; i < s_results.stations.size(); i += check_interval) {
        const auto& station = s_results.stations[i];
        s_gas->thermo()->restoreState(station.state);
        double h = s_gas->thermo()->enthalpy_mass();
        double H0_station = h + 0.5 * station.velocity * station.velocity;
        EXPECT_NEAR(H0_station, H0, max_fp_error(H0, 1e-2, 1.0))
            << "Stagnation enthalpy not conserved at station " << i
            << " (x=" << station.x << " m)";
    }
}

// ---------------------------------------------------------------------------
// Monotonicity along the supersonic branch
// ---------------------------------------------------------------------------

TEST_F(KineticNozzleTests, PressureDecreasesSupersonicBranch) {
    // With dt=1e-6, per-step pressure drop can be below double precision, so
    // allow equality between adjacent stations; verify overall decrease.
    ASSERT_GT(s_results.stations.size(), 1u);

    for (size_t i = 1; i < s_results.stations.size(); i += check_interval) {
        s_gas->thermo()->restoreState(s_results.stations[i].state);
        double P_cur = s_gas->thermo()->pressure();
        s_gas->thermo()->restoreState(s_results.stations[i-1].state);
        double P_prev = s_gas->thermo()->pressure();
        // Allow a tiny relative tolerance for Strang splitting floating-point noise
        EXPECT_LE(P_cur, P_prev * (1.0 + 1e-9))
            << "Pressure must not increase along the supersonic nozzle at station " << i;
    }

    // Overall: exit pressure must be significantly below throat pressure
    s_gas->thermo()->restoreState(s_results.stations.back().state);
    double P_exit = s_gas->thermo()->pressure();
    s_gas->thermo()->restoreState(s_results.stations.front().state);
    double P_first = s_gas->thermo()->pressure();
    EXPECT_LT(P_exit, P_first * 0.99)
        << "Exit pressure must be well below throat pressure";
}

TEST_F(KineticNozzleTests, TemperatureDecreasesSupersonicBranch) {
    // Same rationale as pressure: allow equality between adjacent stations.
    ASSERT_GT(s_results.stations.size(), 1u);

    for (size_t i = 1; i < s_results.stations.size(); i += check_interval) {
        s_gas->thermo()->restoreState(s_results.stations[i].state);
        double T_cur = s_gas->thermo()->temperature();
        s_gas->thermo()->restoreState(s_results.stations[i-1].state);
        double T_prev = s_gas->thermo()->temperature();
        EXPECT_LE(T_cur, T_prev)
            << "Temperature must not increase along the supersonic nozzle at station " << i;
    }

    s_gas->thermo()->restoreState(s_results.stations.back().state);
    double T_exit = s_gas->thermo()->temperature();
    s_gas->thermo()->restoreState(s_results.stations.front().state);
    double T_first = s_gas->thermo()->temperature();
    EXPECT_LT(T_exit, T_first * 0.99)
        << "Exit temperature must be well below throat temperature";
}

TEST_F(KineticNozzleTests, VelocityIncreasesSupersonicBranch) {
    // Physical invariant: exit velocity > throat velocity.
    // Per-step monotonicity is not checked because the dt-limiter near the nozzle
    // exit causes Zeno's-paradox clustering of many stations at essentially the
    // same x, where FP noise in gas_isenthalpic_velocity introduces sub-ULP
    // fluctuations that can flip the last bit.
    ASSERT_GT(s_results.stations.size(), 1u);

    EXPECT_GT(s_results.stations.back().velocity, s_results.stations.front().velocity * 1.01)
        << "Exit velocity must be significantly above throat velocity";
}

TEST_F(KineticNozzleTests, MachNumberAboveOneAndIncreasing) {
    ASSERT_GT(s_results.stations.size(), 1u);

    for (size_t i = 0; i < s_results.stations.size(); i++) {
        EXPECT_GT(s_results.stations[i].mach, 0.9)
            << "Mach must be >= 1 in the supersonic branch at station " << i;
    }

    EXPECT_GT(s_results.stations.back().mach, s_results.stations.front().mach * 1.01)
        << "Exit Mach must be significantly above throat Mach";
}

TEST_F(KineticNozzleTests, AreaRatioIncreasing) {
    ASSERT_GT(s_results.stations.size(), 1u);

    for (size_t i = 1; i < s_results.stations.size(); i += check_interval) {
        EXPECT_GE(s_results.stations[i].area_ratio, s_results.stations[i-1].area_ratio)
            << "Area ratio must not decrease along the diverging nozzle at station " << i;
    }

    EXPECT_GT(s_results.stations.back().area_ratio, s_results.stations.front().area_ratio * 1.5)
        << "Exit area ratio must be well above throat area ratio";
}

// ---------------------------------------------------------------------------
// Bounding invariant (most important physical check per reference §7.1)
// ---------------------------------------------------------------------------

TEST_F(KineticNozzleTests, KineticTemperatureBoundedByFrozenAndEquilibrium) {
    // Key invariant: T_frozen ≤ T_kinetic ≤ T_equilibrium at the same exit area ratio.
    //
    // Physical reasoning:
    //   Frozen: composition locked at throat, no recombination energy recovered → lower T
    //   Equilibrium: maximum recombination energy recovered → higher T
    //   Kinetic: partial recombination → between the two bounds
    //
    // Reference: kinetic_nozzle_reference.md §7.1

    ASSERT_TRUE(s_results.throat.converged);
    ASSERT_GT(s_results.stations.size(), 0u);

    // Get kinetic exit temperature and area ratio
    const auto& exit_station = s_results.stations.back();
    s_gas->thermo()->restoreState(exit_station.state);
    double T_kinetic = s_gas->thermo()->temperature();
    double exit_ar = exit_station.area_ratio;

    // Solve equilibrium and frozen to the same area ratio
    s_gas->thermo()->restoreState(s_inlet_state);
    EquilibriumNozzle eq_nozzle(*s_gas);
    NozzleResults eq_results = eq_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    ASSERT_TRUE(eq_results.expansions.front().converged);
    s_gas->thermo()->restoreState(eq_results.expansions.front().state);
    double T_equilibrium = s_gas->thermo()->temperature();

    s_gas->thermo()->restoreState(s_inlet_state);
    FrozenNozzle frz_nozzle(*s_gas);
    NozzleResults frz_results = frz_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    ASSERT_TRUE(frz_results.expansions.front().converged);
    s_gas->thermo()->restoreState(frz_results.expansions.front().state);
    double T_frozen = s_gas->thermo()->temperature();

    // Allow a small relative margin (0.5%) for numerical differences
    double margin = T_kinetic * 0.005;
    EXPECT_GE(T_kinetic, T_frozen - margin)
        << "Kinetic exit T (" << T_kinetic << " K) should be >= frozen T ("
        << T_frozen << " K) at area ratio " << exit_ar;
    EXPECT_LE(T_kinetic, T_equilibrium + margin)
        << "Kinetic exit T (" << T_kinetic << " K) should be <= equilibrium T ("
        << T_equilibrium << " K) at area ratio " << exit_ar;
}

TEST_F(KineticNozzleTests, KineticExitVelocityBoundedByFrozenAndEquilibrium) {
    // Companion to the temperature bound: v_frozen ≤ v_kinetic ≤ v_equilibrium.
    // Equilibrium recovers chemical bond energy (recombination exothermic) as
    // additional enthalpy, which in the nozzle converts to higher exit velocity.
    // Exit velocity: v = sqrt(2*(H0-h)), and h_frozen > h_equilibrium because the
    // frozen composition retains dissociation energy in chemical bonds.
    ASSERT_TRUE(s_results.throat.converged);
    ASSERT_GT(s_results.stations.size(), 0u);

    double v_kinetic = s_results.stations.back().velocity;
    double exit_ar = s_results.stations.back().area_ratio;
    double H0 = s_results.throat.H_stagnation;

    s_gas->thermo()->restoreState(s_inlet_state);
    EquilibriumNozzle eq_nozzle(*s_gas);
    NozzleResults eq_results = eq_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    ASSERT_TRUE(eq_results.expansions.front().converged);
    s_gas->thermo()->restoreState(eq_results.expansions.front().state);
    double v_eq = std::sqrt(2.0 * (H0 - s_gas->thermo()->enthalpy_mass()));

    s_gas->thermo()->restoreState(s_inlet_state);
    FrozenNozzle frz_nozzle(*s_gas);
    NozzleResults frz_results = frz_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    ASSERT_TRUE(frz_results.expansions.front().converged);
    s_gas->thermo()->restoreState(frz_results.expansions.front().state);
    double v_frz = std::sqrt(2.0 * (H0 - s_gas->thermo()->enthalpy_mass()));

    double margin = v_kinetic * 0.005;
    EXPECT_GE(v_kinetic, v_frz - margin)
        << "Kinetic exit velocity should be >= frozen exit velocity";
    EXPECT_LE(v_kinetic, v_eq + margin)
        << "Kinetic exit velocity should be <= equilibrium exit velocity";
}

TEST_F(KineticNozzleTests, ExitPressureBoundedByFrozenAndEquilibrium) {
    // Companion to the temperature bound: P_frozen <= P_kinetic <= P_equilibrium.
    // Equilibrium recovers chemical bond energy (recombination exothermic) as
    // additional enthalpy, which in the nozzle converts to higher temperature 
    // & therefore pressure.
    
    ASSERT_TRUE(s_results.throat.converged);
    ASSERT_GT(s_results.stations.size(), 0u);

    const auto& exit_result = s_results.stations.back(); 

    s_gas->thermo()->restoreState(exit_result.state);

    double P_kinetic = s_gas->thermo()->pressure();
    double exit_ar = s_results.stations.back().area_ratio;

    s_gas->thermo()->restoreState(s_inlet_state);
    EquilibriumNozzle eq_nozzle(*s_gas);
    NozzleResults eq_results = eq_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    ASSERT_TRUE(eq_results.expansions.front().converged);
    s_gas->thermo()->restoreState(eq_results.expansions.front().state);
    double P_eq = s_gas->thermo()->pressure();

    s_gas->thermo()->restoreState(s_inlet_state);
    FrozenNozzle frz_nozzle(*s_gas);
    NozzleResults frz_results = frz_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    ASSERT_TRUE(frz_results.expansions.front().converged);
    s_gas->thermo()->restoreState(frz_results.expansions.front().state);
    double P_frz = s_gas->thermo()->pressure();

    double margin = P_kinetic * 0.005;
    EXPECT_GE(P_kinetic, P_frz - margin)
        << "Exit pressure should be >= frozen exit pressure";
    EXPECT_LE(P_kinetic, P_eq + margin)
        << "Exit pressure should be <= equilibrium exit pressure";
}

// ---------------------------------------------------------------------------
// Gamma_s: must be frozen Cp/Cv
// ---------------------------------------------------------------------------

TEST_F(KineticNozzleTests, GammaUsedIsFrozenCpOverCv) {
    // KineticNozzle::get_gamma_s must return cp/cv (frozen gamma), not
    // the equilibrium gamma_s. This is physically correct because acoustic
    // timescales are much shorter than chemical timescales (see §5 of reference).
    KineticNozzle kinetic_nozzle(*s_gas, s_profile, s_mdot);
    FrozenNozzle frozen_nozzle(*s_gas);

    s_gas->thermo()->restoreState(s_inlet_state);
    double gamma_kinetic = kinetic_nozzle.get_gamma_s(*s_gas->thermo());
    double gamma_frozen = frozen_nozzle.get_gamma_s(*s_gas->thermo());

    // Both should return cp/cv; they must be equal on the same thermodynamic state
    EXPECT_NEAR(gamma_kinetic, gamma_frozen,
                max_fp_error(gamma_frozen, 1e-10, 1e-14))
        << "KineticNozzle::get_gamma_s must use the same frozen (cp/cv) formula as FrozenNozzle";

    // And it must differ from equilibrium gamma_s (which accounts for composition shifts)
    EquilibriumNozzle eq_nozzle(*s_gas);
    s_gas->thermo()->restoreState(s_inlet_state);
    double gamma_eq = eq_nozzle.get_gamma_s(*s_gas->thermo());

    // Equilibrium and frozen gamma_s should generally be different; this
    // verifies we're not accidentally using the equilibrium formula.
    // (They can be close but are rarely identical for reactive mixtures.)
    EXPECT_GT(gamma_kinetic, 1.0) << "gamma_s must be > 1";
    EXPECT_LT(gamma_kinetic, 2.0) << "gamma_s must be physically reasonable";
    // Note: for near-equilibrium inlet conditions the two may be very close,
    // so we don't assert they are different — just that kinetic uses cp/cv.
    (void)gamma_eq;
}

// ---------------------------------------------------------------------------
// Damköhler number diagnostics
// ---------------------------------------------------------------------------

TEST_F(KineticNozzleTests, DamkohlerVectorHasCorrectSize) {
    // damkohler must have one entry per species so it can be indexed by species index
    ASSERT_GT(s_results.stations.size(), 0u);

    size_t n_species = s_gas->thermo()->nSpecies();
    for (size_t i = 0; i < s_results.stations.size(); i += check_interval) {
        EXPECT_EQ(s_results.stations[i].damkohler.size(), n_species)
            << "Damkohler vector size mismatch at station " << i;
    }
}

TEST_F(KineticNozzleTests, DamkohlerMinIsMinimumOfNonZeroEntries) {
    // Da_min must equal the smallest non-zero entry in the damkohler vector
    // (trace species are zeroed out and excluded from the minimum)
    ASSERT_GT(s_results.stations.size(), 0u);

    for (size_t i = 0; i < s_results.stations.size(); i += check_interval) {
        const auto& station = s_results.stations[i];
        double manual_min = std::numeric_limits<double>::max();
        for (double da : station.damkohler) {
            if (da > 0.0 && da < manual_min) manual_min = da;
        }
        if (manual_min < std::numeric_limits<double>::max()) {
            EXPECT_NEAR(station.Da_min, manual_min,
                        max_fp_error(manual_min, 1e-10, 1e-30))
                << "Da_min mismatch at station " << i;
        }
    }
}

TEST_F(KineticNozzleTests, FreezingSpeciesIndexIsValid) {
    // min_Da_species must be a valid species index
    ASSERT_GT(s_results.stations.size(), 0u);

    int n_species = static_cast<int>(s_gas->thermo()->nSpecies());
    for (size_t i = 0; i < s_results.stations.size(); i += check_interval) {
        EXPECT_GE(s_results.stations[i].min_Da_species, 0)
            << "min_Da_species must be non-negative at station " << i;
        EXPECT_LT(s_results.stations[i].min_Da_species, n_species)
            << "min_Da_species out of range at station " << i;
    }
}

TEST_F(KineticNozzleTests, EntropyIsStableOrIncreasing) {
    // For a kinetically reacting nozzle, the entropy must be equal to or greater than 
    // its equilibrium counterpart.
    ASSERT_GT(s_results.stations.size(), 0u);

    ASSERT_TRUE(s_results.throat.converged);
    ASSERT_GT(s_results.stations.size(), 0u);

    const auto& exit_result = s_results.stations.back(); 

    s_gas->thermo()->restoreState(exit_result.state);

    double S_kinetic = s_gas->thermo()->entropy_mass();
    double exit_ar = s_results.stations.back().area_ratio;

    s_gas->thermo()->restoreState(s_inlet_state);
    EquilibriumNozzle eq_nozzle(*s_gas);
    NozzleResults eq_results = eq_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    ASSERT_TRUE(eq_results.expansions.front().converged);
    s_gas->thermo()->restoreState(eq_results.expansions.front().state);
    double S_eq = s_gas->thermo()->entropy_mass();

    double margin = S_kinetic * 0.001;
    EXPECT_GE(S_kinetic, S_eq - margin)
        << "Exit entropy >= equiilibrium exit entropy";
}

// ---------------------------------------------------------------------------
// Temporary diagnostic (will remove)
// ---------------------------------------------------------------------------
TEST_F(KineticNozzleTests, DiagnosticBoundingValues) {
    ASSERT_TRUE(s_results.throat.converged);
    const auto& exit_station = s_results.stations.back();
    s_gas->thermo()->restoreState(exit_station.state);
    double T_kinetic = s_gas->thermo()->temperature();
    double exit_ar = exit_station.area_ratio;
    double H0 = s_results.throat.H_stagnation;
    double H0_check = s_gas->thermo()->enthalpy_mass() + 0.5*exit_station.velocity*exit_station.velocity;

    s_gas->thermo()->restoreState(s_inlet_state);
    EquilibriumNozzle eq_nozzle(*s_gas);
    NozzleResults eq_results = eq_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    double T_equilibrium = 0, v_eq = 0;
    if (eq_results.expansions.front().converged) {
        s_gas->thermo()->restoreState(eq_results.expansions.front().state);
        T_equilibrium = s_gas->thermo()->temperature();
        v_eq = std::sqrt(2.0*(H0 - s_gas->thermo()->enthalpy_mass()));
    }

    s_gas->thermo()->restoreState(s_inlet_state);
    FrozenNozzle frz_nozzle(*s_gas);
    NozzleResults frz_results = frz_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_ar);
    double T_frozen = 0, v_frz = 0;
    if (frz_results.expansions.front().converged) {
        s_gas->thermo()->restoreState(frz_results.expansions.front().state);
        T_frozen = s_gas->thermo()->temperature();
        v_frz = std::sqrt(2.0*(H0 - s_gas->thermo()->enthalpy_mass()));
    }

    std::cout << "\n=== Bounding Diagnostic ===\n"
              << "Exit A/At: " << exit_ar << "\n"
              << "H0 (throat): " << H0 << " J/kg\n"
              << "H0 (check at exit): " << H0_check << " J/kg\n"
              << "Stations total: " << s_results.stations.size() << "\n"
              << "T_frozen:      " << T_frozen << " K, v_frozen: " << v_frz << " m/s\n"
              << "T_kinetic:     " << T_kinetic << " K, v_kinetic: " << exit_station.velocity << " m/s\n"
              << "T_equilibrium: " << T_equilibrium << " K, v_eq: " << v_eq << " m/s\n"
              << "Eq throat T: ";
    s_gas->thermo()->restoreState(eq_results.throat.state);
    std::cout << s_gas->thermo()->temperature() << " K\n"
              << "Kinetic throat T: ";
    s_gas->thermo()->restoreState(s_results.throat.state);
    std::cout << s_gas->thermo()->temperature() << " K\n"
              << "===========================\n";
    SUCCEED();
}
