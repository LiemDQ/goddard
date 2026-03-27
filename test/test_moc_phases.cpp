#include "goddard/moc.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/nozzle.hpp"
#include "cantera/core.h"
#include <cmath>
#include <fstream>
#include <filesystem>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

// Helper: create a MocNozzle configured for planar perfect gas
static MocNozzle make_perfect_gas_solver(
    double gamma, double theta_max, int num_chars,
    MocFlowKind flow_type = MocFlowKind::PLANAR)
{
    MocOptions opts;
    opts.flow_type = flow_type;
    opts.chemistry = MocChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = num_chars;
    opts.geometry.throat_radius = 1.0;
    return MocNozzle(opts);
}

// ============================================================
// Phase 2: Wall point storage
// ============================================================

TEST(MocPhase2, WallPointsStored) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7);
    auto result = solver.solve();

    // Wall points should be stored with full flow properties
    EXPECT_GT(result.net.wall_points.size(), 0u);
    EXPECT_EQ(result.net.wall_points.size(), result.net.wall_x.size() - 1);
    // wall_x/wall_y has an extra entry for the throat point at (0,1)

    for (const auto& wp : result.net.wall_points) {
        EXPECT_GT(wp.mach, 1.0) << "Wall point should be supersonic";
        EXPECT_GT(wp.x, 0.0) << "Wall point should be downstream of throat";
        EXPECT_GT(wp.y, 0.0) << "Wall point should be above axis";
        EXPECT_GT(wp.pressure, 0.0) << "Wall pressure should be positive";
        EXPECT_GT(wp.temperature, 0.0) << "Wall temperature should be positive";
    }
}

TEST(MocPhase2, WallPointThetaDecreasing) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7);
    auto result = solver.solve();

    // For a min-length nozzle, wall theta should decrease from theta_max to ~0
    for (size_t i = 1; i < result.net.wall_points.size(); i++) {
        EXPECT_LE(result.net.wall_points[i].theta, result.net.wall_points[i-1].theta)
            << "Wall theta should decrease monotonically";
    }
    // Last wall point should have theta near 0
    EXPECT_NEAR(result.net.wall_points.back().theta, 0.0, 1.0 * DEG);
}

// ============================================================
// Phase 2: Exit plane extraction
// ============================================================

TEST(MocPhase2, ExitPlanePopulated) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7);
    auto result = solver.solve();

    EXPECT_GT(result.exit_plane.mach.size(), 0u);
    EXPECT_EQ(result.exit_plane.mach.size(), result.exit_plane.y.size());
    EXPECT_EQ(result.exit_plane.mach.size(), result.exit_plane.theta.size());
    EXPECT_EQ(result.exit_plane.mach.size(), result.exit_plane.pressure.size());
    EXPECT_EQ(result.exit_plane.mach.size(), result.exit_plane.temperature.size());
}

TEST(MocPhase2, ExitPlaneUniformMach) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 10);
    auto result = solver.solve();

    // For a min-length nozzle, exit plane Mach should be uniform
    double expected_exit_mach = result.exit_mach;
    for (size_t i = 0; i < result.exit_plane.mach.size(); i++) {
        EXPECT_NEAR(result.exit_plane.mach[i], expected_exit_mach, 0.05)
            << "Exit plane Mach not uniform at index " << i;
    }
}

TEST(MocPhase2, ExitPlaneZeroTheta) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 10);
    auto result = solver.solve();

    // For a min-length nozzle, exit plane theta should be ~0 (uniform axial flow)
    for (size_t i = 0; i < result.exit_plane.theta.size(); i++) {
        EXPECT_NEAR(result.exit_plane.theta[i], 0.0, 1.0 * DEG)
            << "Exit plane theta not zero at index " << i;
    }
}

// ============================================================
// Phase 2: CSV export/import
// ============================================================

TEST(MocPhase2, CsvRoundTrip) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7);
    auto result = solver.solve();

    std::string tmp_file = "/tmp/goddard_test_profile.csv";

    // Save
    result.profile.save_profile_csv(tmp_file);

    // Load
    auto loaded = NozzleProfile::load_profile_csv(tmp_file);

    ASSERT_EQ(loaded.x.size(), result.profile.x.size());
    for (size_t i = 0; i < loaded.x.size(); i++) {
        EXPECT_NEAR(loaded.x[i], result.profile.x[i], 1e-10);
        EXPECT_NEAR(loaded.y[i], result.profile.y[i], 1e-10);
    }

    // Cleanup
    std::filesystem::remove(tmp_file);
}

// ============================================================
// Phase 2: Convergence checking
// ============================================================

TEST(MocPhase2, ConvergedFlagTrue) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7);
    auto result = solver.solve();
    EXPECT_TRUE(result.converged);
    EXPECT_TRUE(result.messages.empty());
}

// ============================================================
// Phase 3: Axisymmetric basic solve
// ============================================================

TEST(MocPhase3, AxiPerfectGasSolves) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7, MocFlowKind::AXISYMMETRIC);
    auto result = solver.solve();

    EXPECT_TRUE(result.converged);
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.nozzle_length, 0.0);
    EXPECT_GT(result.area_ratio, 1.0);
}

TEST(MocPhase3, AxiSymmetryOnAxis) {
    auto solver = make_perfect_gas_solver(1.4, 12.0 * DEG, 5, MocFlowKind::AXISYMMETRIC);
    auto result = solver.solve();

    // All axis points (first point in each wavefront after the initial data line)
    // should have theta=0 and y=0
    for (size_t i = 1; i < result.net.wavefronts.size(); i++) {
        const auto& wf = result.net.wavefronts[i];
        ASSERT_FALSE(wf.empty());
        EXPECT_NEAR(wf[0].theta, 0.0, 1e-8)
            << "Axis theta should be 0 in wavefront " << i;
        EXPECT_NEAR(wf[0].y, 0.0, 1e-8)
            << "Axis y should be 0 in wavefront " << i;
    }
}

TEST(MocPhase3, AxiExitMachReasonable) {
    // Axisymmetric nozzle should produce a different exit Mach than planar
    // for the same theta_max, due to the source terms
    double gamma = 1.4;
    double theta_max = 15.0 * DEG;

    auto planar = make_perfect_gas_solver(gamma, theta_max, 10, MocFlowKind::PLANAR);
    auto axi = make_perfect_gas_solver(gamma, theta_max, 10, MocFlowKind::AXISYMMETRIC);

    auto result_planar = planar.solve();
    auto result_axi = axi.solve();

    // Both should produce valid results
    EXPECT_TRUE(result_planar.converged);
    EXPECT_TRUE(result_axi.converged);

    // Axisymmetric exit Mach should differ from planar
    EXPECT_NE(result_planar.exit_mach, result_axi.exit_mach);

    // Both should be supersonic
    EXPECT_GT(result_axi.exit_mach, 1.0);
}

TEST(MocPhase3, AxiMonotonicWall) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7, MocFlowKind::AXISYMMETRIC);
    auto result = solver.solve();

    // Wall x should be monotonically increasing
    for (size_t i = 1; i < result.net.wall_x.size(); i++) {
        EXPECT_GT(result.net.wall_x[i], result.net.wall_x[i-1])
            << "Wall x should increase monotonically";
    }
    // Wall y should be monotonically increasing (nozzle expands)
    for (size_t i = 1; i < result.net.wall_y.size(); i++) {
        EXPECT_GE(result.net.wall_y[i], result.net.wall_y[i-1])
            << "Wall y should increase monotonically";
    }
}

TEST(MocPhase3, AxiAreaRatioSquared) {
    auto solver = make_perfect_gas_solver(1.4, 15.0 * DEG, 7, MocFlowKind::AXISYMMETRIC);
    auto result = solver.solve();

    // For axisymmetric, area_ratio = (y_exit/y_throat)^2
    double y_ratio = result.net.wall_y.back() / result.net.wall_y.front();
    EXPECT_NEAR(result.area_ratio, y_ratio * y_ratio, 0.01);
}

// ============================================================
// Phase 5: Frozen chemistry integration test
// ============================================================

class MocFrozenTest : public ::testing::Test {
protected:
    void SetUp() override {
        gas = Cantera::newSolution("h2o2.yaml", "ohmech");
        auto thermo = gas->thermo();
        // H2/O2 combustion products
        thermo->setState_TPX(3000.0, 3e6, "H2O:0.8, OH:0.1, H2:0.05, O2:0.05");
    }
    std::shared_ptr<Cantera::Solution> gas;
};

TEST_F(MocFrozenTest, PlanarFrozenSolves) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = MocChemistry::FROZEN;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.theta_max = 15.0 * DEG;
    opts.num_characteristics = 7;
    opts.geometry.throat_radius = 0.05;

    MocNozzle nozzle(gas, opts);
    auto result = nozzle.solve();

    EXPECT_TRUE(result.converged);
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.nozzle_length, 0.0);
    EXPECT_GT(result.area_ratio, 1.0);

    // Wall should be monotonically increasing in x
    for (size_t i = 1; i < result.net.wall_x.size(); i++) {
        EXPECT_GT(result.net.wall_x[i], result.net.wall_x[i-1]);
    }
}

TEST_F(MocFrozenTest, FrozenExitMachConsistent) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = MocChemistry::FROZEN;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.theta_max = 12.0 * DEG;
    opts.num_characteristics = 5;
    opts.geometry.throat_radius = 0.05;

    MocNozzle nozzle(gas, opts);
    auto result = nozzle.solve();

    // Exit Mach should be positive and reasonable
    EXPECT_GT(result.exit_mach, 1.5);
    EXPECT_LT(result.exit_mach, 10.0);
}

// ============================================================
// Phase 5: Equilibrium chemistry integration test
// ============================================================

TEST_F(MocFrozenTest, PlanarEquilibriumSolves) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = MocChemistry::EQUILIBRIUM;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.theta_max = 15.0 * DEG;
    opts.num_characteristics = 7;
    opts.geometry.throat_radius = 0.05;

    MocNozzle nozzle(gas, opts);
    auto result = nozzle.solve();

    EXPECT_TRUE(result.converged);
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.nozzle_length, 0.0);
    EXPECT_GT(result.area_ratio, 1.0);
}

TEST_F(MocFrozenTest, EquilibriumVsFrozenDiffers) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.theta_max = 15.0 * DEG;
    opts.num_characteristics = 7;
    opts.geometry.throat_radius = 0.05;

    // Frozen
    opts.chemistry = MocChemistry::FROZEN;
    MocNozzle frozen_nozzle(gas, opts);
    auto frozen_result = frozen_nozzle.solve();

    // Reset gas state for equilibrium run
    gas->thermo()->setState_TPX(3000.0, 3e6, "H2O:0.8, OH:0.1, H2:0.05, O2:0.05");

    // Equilibrium
    opts.chemistry = MocChemistry::EQUILIBRIUM;
    MocNozzle equil_nozzle(gas, opts);
    auto equil_result = equil_nozzle.solve();

    // Both should converge
    EXPECT_TRUE(frozen_result.converged);
    EXPECT_TRUE(equil_result.converged);

    // Exit Mach should differ (equilibrium shifts composition)
    // The difference may be small for this composition, but should not be identical
    EXPECT_GT(std::abs(frozen_result.exit_mach - equil_result.exit_mach), 1e-6);
}
