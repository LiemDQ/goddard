#include "goddard/prandtlmeyer.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/equilibrium.hpp"
#include "cantera/core.h"
#include <cmath>
#include "gtest/gtest.h"

using namespace Goddard;

// Test fixture: builds a PrandtlMeyerTable with a monatomic ideal gas (Ar)
// where gamma = 5/3 is constant. This lets us compare table values against
// the closed-form Prandtl-Meyer function.
class PrandtlMeyerTableTest : public ::testing::Test {
protected:
    void SetUp() override {
        gas = Cantera::newSolution("h2o2.yaml", "ohmech");
        auto thermo = gas->thermo();

        // Pure argon: gamma = 5/3, calorically perfect
        // Set stagnation conditions first
        thermo->setState_TPX(1000.0, 1e6, "AR:1.0");
        gamma = thermo->cp_mass() / thermo->cv_mass();

        double s0 = thermo->entropy_mass();
        double h0 = thermo->enthalpy_mass();

        // Set to throat conditions (isentropic expansion to M=1)
        // For ideal gas: T* = T0 * 2/(gamma+1), P* = P0 * (2/(gamma+1))^(gamma/(gamma-1))
        double T_throat = 1000.0 * 2.0 / (gamma + 1.0);
        double P_throat = 1e6 * pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0));
        thermo->setState_TP(T_throat, P_throat);

        double a_throat = gas_sonic_velocity(*thermo, gamma);

        table.build_table(*thermo, false, s0, h0, a_throat);
    }

    std::shared_ptr<Cantera::Solution> gas;
    PrandtlMeyerTable table;
    double gamma;
};

// ============================================================
// Construction and invariants
// ============================================================

TEST_F(PrandtlMeyerTableTest, TableBuildsSuccessfully) {
    EXPECT_TRUE(table.is_built());
    EXPECT_GT(table.machs.size(), 0u);
}

TEST_F(PrandtlMeyerTableTest, MonotonicVelocity) {
    for (size_t i = 1; i < table.velocities.size(); i++) {
        EXPECT_GT(table.velocities[i], table.velocities[i - 1])
            << "Velocity not monotonically increasing at index " << i;
    }
}

TEST_F(PrandtlMeyerTableTest, MonotonicMach) {
    for (size_t i = 1; i < table.machs.size(); i++) {
        EXPECT_GT(table.machs[i], table.machs[i - 1])
            << "Mach not monotonically increasing at index " << i;
    }
}

TEST_F(PrandtlMeyerTableTest, MonotonicNu) {
    for (size_t i = 1; i < table.nus.size(); i++) {
        EXPECT_GE(table.nus[i], table.nus[i - 1])
            << "Nu not monotonically increasing at index " << i;
    }
}

TEST_F(PrandtlMeyerTableTest, DecreasingEnthalpy) {
    for (size_t i = 1; i < table.enthalpies.size(); i++) {
        EXPECT_LT(table.enthalpies[i], table.enthalpies[i - 1])
            << "Enthalpy not monotonically decreasing at index " << i;
    }
}

TEST_F(PrandtlMeyerTableTest, ConstantGammaForIdealGas) {
    // For a monatomic ideal gas, gamma_s should be constant = 5/3
    for (size_t i = 0; i < table.gamma_s.size(); i++) {
        EXPECT_NEAR(table.gamma_s[i], gamma, 0.01)
            << "gamma_s deviates from expected value at index " << i;
    }
}

TEST_F(PrandtlMeyerTableTest, FirstPointIsSonic) {
    EXPECT_NEAR(table.machs[0], 1.0, 0.02);
    EXPECT_NEAR(table.nus[0], 0.0, 0.01);
}

TEST_F(PrandtlMeyerTableTest, VelocityEnergyConservation) {
    // h0 = V^2/2 + h at every point (energy conservation)
    // h0 = V_throat^2/2 + h_throat, and V_throat = velocities[0], h_throat = enthalpies[0]
    double h0 = table.velocities[0] * table.velocities[0] / 2.0 + table.enthalpies[0];
    for (size_t i = 0; i < table.velocities.size(); i++) {
        double h_total = table.velocities[i] * table.velocities[i] / 2.0 + table.enthalpies[i];
        EXPECT_NEAR(h_total, h0, h0 * 1e-6)
            << "Energy conservation violated at index " << i;
    }
}

// ============================================================
// Interpolation accuracy (compare to closed-form PM function)
// ============================================================

TEST_F(PrandtlMeyerTableTest, InterpolateMachFromNu) {
    // Test at several Mach numbers: compute exact nu, interpolate M, compare
    for (double mach : {1.2, 1.5, 2.0, 2.5, 3.0}) {
        double nu_exact = prandtl_meyer(mach, gamma);
        // Only test if nu is within table range
        if (nu_exact < table.nus.front() || nu_exact > table.nus.back()) continue;
        double mach_interp = table.interpolate_mach(nu_exact);
        EXPECT_NEAR(mach_interp, mach, 0.01)
            << "Mach interpolation error at M=" << mach;
    }
}

TEST_F(PrandtlMeyerTableTest, InterpolateNuFromMach) {
    for (double mach : {1.2, 1.5, 2.0, 2.5, 3.0}) {
        if (mach < table.machs.front() || mach > table.machs.back()) continue;
        double nu_exact = prandtl_meyer(mach, gamma);
        double nu_interp = table.interpolate_nu_from_mach(mach);
        EXPECT_NEAR(nu_interp, nu_exact, 0.005)
            << "Nu interpolation error at M=" << mach;
    }
}

TEST_F(PrandtlMeyerTableTest, InterpolateGammaSFromMach) {
    // For a constant-gamma gas, interpolated gamma should equal gamma everywhere
    for (double mach : {1.1, 1.5, 2.0, 3.0}) {
        if (mach < table.machs.front() || mach > table.machs.back()) continue;
        double gamma_interp = table.interpolate_gamma_s_from_mach(mach);
        EXPECT_NEAR(gamma_interp, gamma, 0.01)
            << "Gamma interpolation error at M=" << mach;
    }
}

TEST_F(PrandtlMeyerTableTest, RoundTripMachNuMach) {
    for (double mach : {1.3, 1.8, 2.5}) {
        if (mach < table.machs.front() || mach > table.machs.back()) continue;
        double nu = table.interpolate_nu_from_mach(mach);
        double mach_back = table.interpolate_mach(nu);
        EXPECT_NEAR(mach_back, mach, 0.02)
            << "Round trip M->nu->M failed at M=" << mach;
    }
}

TEST_F(PrandtlMeyerTableTest, InterpolationAccuracyImproves) {
    // A finer table should have smaller interpolation error
    auto thermo = gas->thermo();
    thermo->setState_TPX(1000.0, 1e6, "AR:1.0");
    double s0 = thermo->entropy_mass();
    double h0 = thermo->enthalpy_mass();
    double T_throat = 1000.0 * 2.0 / (gamma + 1.0);
    double P_throat = 1e6 * pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0));
    thermo->setState_TP(T_throat, P_throat);
    double a_throat = gas_sonic_velocity(*thermo, gamma);

    PrandtlMeyerTable coarse_table;
    coarse_table.build_table(*thermo, false, s0, h0, a_throat, 1e-4, 50);

    // Reset to throat conditions for fine table
    thermo->setState_TP(T_throat, P_throat);
    PrandtlMeyerTable fine_table;
    fine_table.build_table(*thermo, false, s0, h0, a_throat, 1e-4, 500);

    // Compare at M=2.0
    double mach_test = 2.0;
    double nu_exact = prandtl_meyer(mach_test, gamma);

    double error_coarse = std::abs(coarse_table.interpolate_mach(nu_exact) - mach_test);
    double error_fine = std::abs(fine_table.interpolate_mach(nu_exact) - mach_test);

    EXPECT_LT(error_fine, error_coarse)
        << "Finer table should have smaller error";
}

// ============================================================
// Edge cases and error handling
// ============================================================

TEST_F(PrandtlMeyerTableTest, QueryBelowRange) {
    // Querying below the sonic point should throw
    EXPECT_THROW(table.interpolate_mach(-0.1), std::out_of_range);
}

TEST_F(PrandtlMeyerTableTest, QueryAboveRange) {
    // Querying beyond the table's maximum nu should throw
    double nu_max = table.nus.back() + 1.0;
    EXPECT_THROW(table.interpolate_V(nu_max), std::out_of_range);
}

TEST_F(PrandtlMeyerTableTest, QueryOnUnbuiltTable) {
    PrandtlMeyerTable empty_table;
    EXPECT_THROW(empty_table.interpolate_mach(0.5), std::runtime_error);
}

// ============================================================
// Integration quality
// ============================================================

TEST_F(PrandtlMeyerTableTest, NuMatchesClosedForm) {
    // Compare accumulated nu at the last table point against closed-form
    double mach_final = table.machs.back();
    double nu_table = table.nus.back();
    double nu_exact = prandtl_meyer(mach_final, gamma);

    // With trapezoidal rule and 500 points, expect < 0.5% error
    double rel_error = std::abs(nu_table - nu_exact) / nu_exact;
    EXPECT_LT(rel_error, 0.005)
        << "Table nu=" << nu_table << " vs exact nu=" << nu_exact
        << " (rel error=" << rel_error << ")";
}
