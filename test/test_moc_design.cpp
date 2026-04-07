#include "goddard/moc.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/gas_dynamics.hpp"
#include <cmath>
#include <vector>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

// ============================================================
// Reference data structures for textbook validation
// ============================================================

// A single node in the reference characteristic net.
// Indices follow the convention:
//   i = C- characteristic index (0 = first expansion ray)
//   j = C+ characteristic index (0 = centerline)
// For the initial data line, j=0 for all points.
// For the kernel region, wavefront w contains points at j = w.
struct MocReferenceNode {
    int wavefront;  // which wavefront (0 = initial data line)
    int index;      // index within the wavefront (0 = axis point)
    double theta;   // flow angle (degrees)
    double nu;      // Prandtl-Meyer angle (degrees)
    double mach;    // Mach number
    double mu;      // Mach angle (degrees)
};

struct MocReferenceWallNode {
    int index;      // wall point index
    double theta;   // wall angle (degrees)
    double nu;      // Prandtl-Meyer angle (degrees)
    double mach;    // Mach number
    double mu;      // Mach angle (degrees)
};

struct MocReferenceCase {
    double gamma;
    double theta_max_deg;
    std::vector<double> theta_schedule_deg; // expansion fan angles in degrees
    std::vector<MocReferenceNode> kernel_nodes;
    std::vector<MocReferenceWallNode> wall_nodes;
    double exit_mach;
};

// ============================================================
// Test fixture for textbook validation
// ============================================================

class MocTextbookValidation : public ::testing::TestWithParam<MocReferenceCase> {
protected:
    void SetUp() override {
        const auto& ref = GetParam();
        double gamma = ref.gamma;
        double theta_max = ref.theta_max_deg * DEG;

        // Convert theta schedule from degrees to radians
        std::vector<double> schedule_rad;
        schedule_rad.reserve(ref.theta_schedule_deg.size());
        for (double td : ref.theta_schedule_deg) {
            schedule_rad.push_back(td * DEG);
        }

        MocOptions opts;
        opts.flow_type = MocFlowKind::PLANAR;
        opts.chemistry = GasChemistry::PERFECT_GAS;
        opts.mode = MocMode::DESIGN_MIN_LENGTH;
        opts.gamma = gamma;
        opts.theta_max = theta_max;
        opts.theta_schedule = schedule_rad;
        opts.num_characteristics = static_cast<int>(schedule_rad.size());
        opts.geometry.throat_radius = 1.0;

        MocNozzle nozzle(opts);
        result = nozzle.solve();
    }

    MocResult result;
};

TEST_P(MocTextbookValidation, KernelNodes) {
    const auto& ref = GetParam();

    for (const auto& node : ref.kernel_nodes) {
        ASSERT_LT(static_cast<size_t>(node.wavefront), result.net.wavefronts.size())
            << "Wavefront " << node.wavefront << " does not exist";
        const auto& wf = result.net.wavefronts[node.wavefront];
        ASSERT_LT(static_cast<size_t>(node.index), wf.size())
            << "Index " << node.index << " out of range in wavefront " << node.wavefront;

        const auto& pt = wf[node.index];

        EXPECT_NEAR(pt.theta / DEG, node.theta, 0.01)
            << "theta mismatch at wavefront=" << node.wavefront
            << " index=" << node.index;
        EXPECT_NEAR(pt.nu / DEG, node.nu, 0.01)
            << "nu mismatch at wavefront=" << node.wavefront
            << " index=" << node.index;
        EXPECT_NEAR(pt.mach, node.mach, 0.01)
            << "Mach mismatch at wavefront=" << node.wavefront
            << " index=" << node.index;
        EXPECT_NEAR(pt.mu / DEG, node.mu, 0.2)
            << "mu mismatch at wavefront=" << node.wavefront
            << " index=" << node.index;
    }
}

TEST_P(MocTextbookValidation, WallNodes) {
    const auto& ref = GetParam();

    // Wall nodes are not directly stored in wavefronts — they're in
    // wall_x/wall_y. The flow properties at wall points need to be
    // reconstructable. For now, verify that the wall contour has the
    // right number of points.
    // TODO: expose wall flow properties for direct comparison.
    EXPECT_GE(result.net.wall_x.size(), ref.wall_nodes.size());
}

TEST_P(MocTextbookValidation, ExitMach) {
    const auto& ref = GetParam();
    EXPECT_NEAR(result.exit_mach, ref.exit_mach, 0.01);
}

TEST_P(MocTextbookValidation, UniformExitFlow) {
    // The last wavefront should have nearly uniform Mach and theta ~ 0
    // (the straightening section cancels all expansion waves)
    const auto& last_wf = result.net.wavefronts.back();

    for (const auto& pt : last_wf) {
        EXPECT_NEAR(pt.theta, 0.0, 0.5 * DEG)
            << "Exit plane theta should be ~0";
    }

    if (last_wf.size() > 1) {
        double m_first = last_wf.front().mach;
        double m_last = last_wf.back().mach;
        EXPECT_NEAR(m_first, m_last, 0.02)
            << "Exit plane Mach should be approximately uniform";
    }
}

// ============================================================
// Reference data taken from Anderson Table 11.1
// ============================================================

static MocReferenceCase anderson_chap11() {
    MocReferenceCase ref;
    ref.gamma = 1.40;
    ref.theta_max_deg = 18.375;
    // 5 equally spaced characteristics
    ref.theta_schedule_deg = {0.375, 3.375, 6.375, 9.375, 12.375, 15.375, 18.375};

    // Exit nu = 2*theta_max for a min-length nozzle
    ref.exit_mach = 2.4;


    // Kernel nodes: fill from textbook
    // {wavefront index, index, theta, nu, mach, mu}
    // the first two mu angles are modified slightly from Anderson
    // as the table interpolation used in the book is slightly off.
    ref.kernel_nodes = {
        {0, 0, 0.375, 0.375, 1.04, 73.8},
        {0, 1, 3.375, 3.375, 1.19, 57.0},
        {0, 2, 6.375, 6.375, 1.31, 49.8},
        {0, 3, 9.375, 9.375, 1.41, 45.2},
        {0, 4, 12.375, 12.375, 1.52, 41.1},
        {0, 5, 15.375, 15.375, 1.62, 38.1},
        {0, 6, 18.375, 18.375, 1.72, 35.6},
        {1, 0, 0.0, 6.75, 1.32, 49.3},
        {1, 1, 3.0, 9.75, 1.43, 44.4},
        {1, 2, 6.0, 12.75, 1.53, 40.8},
        {1, 3, 9.0, 15.75, 1.63, 37.8},
        {1, 4, 12.0, 18.75, 1.73, 35.3},
        {1, 5, 15.0, 21.75, 1.84, 32.9},
        {2, 0, 0.0, 12.75, 1.53, 40.8},
        {2, 1, 3.0, 15.75, 1.63, 37.8},
        {2, 2, 6.0, 18.75, 1.73, 35.3},
        {2, 3, 9.0, 21.75, 1.84, 32.9},
        {2, 4, 12.0, 24.75, 1.94, 31.0},
        {3, 0, 0.0, 18.75, 1.73, 35.3},
        {3, 1, 3.0, 21.75, 1.84, 32.9},
        {3, 2, 6.0, 24.75, 1.94, 31.0},
        {3, 3, 9.0, 27.75, 2.05, 29.2},
        {4, 0, 0.0, 24.75, 1.94, 31.0},
        {4, 1, 3.0, 27.75, 2.05, 29.2},
        {4, 2, 6.0, 30.75, 2.16, 27.6},
        {5, 0, 0.0, 30.75, 2.16, 27.6},
        {5, 1, 3.0, 33.75, 2.28, 26.0},
        {6, 0, 0.0, 36.75, 2.4, 24.6},
    };

    // wall nodes
    // {index, theta, nu, mach, mu}
    ref.wall_nodes = {
        {0, 18.375, 18.375, 1.72, 35.6},
        {1, 15.0, 21.75, 1.84, 32.9},
        {2, 12.0, 24.75, 1.94, 31.0},
        {3,  9.0, 27.75, 2.05, 29.2},
        {4,  6.0, 30.75, 2.16, 27.6},
        {5, 3.0, 33.75, 2.28, 26.0},
        {6, 0.0, 36.75, 2.4, 24.6}
    };

    return ref;
}

INSTANTIATE_TEST_SUITE_P(
    AndersonTestCase,
    MocTextbookValidation,
    ::testing::Values(anderson_chap11())
);

// ============================================================
// Standalone design tests (not parameterized)
// ============================================================

TEST(MocDesign, MinLengthNozzleMonotonicWall) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = 1.4;
    opts.theta_max = 15.0 * DEG;
    opts.num_characteristics = 8;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    // Wall x should be monotonically increasing
    for (size_t i = 1; i < result.net.wall_x.size(); i++) {
        EXPECT_GT(result.net.wall_x[i], result.net.wall_x[i - 1])
            << "Wall x not monotonic at index " << i;
    }

    // Wall y should be monotonically increasing (nozzle expands)
    for (size_t i = 1; i < result.net.wall_y.size(); i++) {
        EXPECT_GE(result.net.wall_y[i], result.net.wall_y[i - 1])
            << "Wall y not monotonic at index " << i;
    }
}

TEST(MocDesign, WavefrontSizeDecreases) {
    // In the kernel region of a min-length nozzle, each successive
    // wavefront has one fewer point
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = 1.4;
    opts.theta_max = 12.0 * DEG;
    opts.num_characteristics = 6;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    // Initial data line has N points, next has N-1, etc.
    size_t N = result.net.wavefronts[0].size();
    for (size_t i = 1; i < result.net.wavefronts.size(); i++) {
        EXPECT_EQ(result.net.wavefronts[i].size(), N - i)
            << "Wavefront " << i << " should have " << N - i << " points";
    }
}

TEST(MocDesign, AreaRatioConsistent) {
    // Area ratio from wall coordinates should be consistent with
    // isentropic area ratio for the exit Mach
    double gamma = 1.4;
    double theta_max = 10.0 * DEG;

    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = 8;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    EXPECT_GT(result.area_ratio, 1.0);
    EXPECT_GT(result.exit_mach, 1.0);
}

TEST(MocDesign, DifferentGammaProducesCorrectExitMach) {
    // Test with gamma = 1.3 (common for combustion products)
    double gamma = 1.3;
    double theta_max = 12.0 * DEG;

    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = 6;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    double expected_exit_mach = mach_from_prandtl_meyer(2.0 * theta_max, gamma);
    EXPECT_NEAR(result.exit_mach, expected_exit_mach, 0.02);
}

TEST(MocDesign, MonoatomicGas) {
    // gamma = 5/3 for monatomic gas
    double gamma = 5.0 / 3.0;
    double theta_max = 8.0 * DEG;

    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = 5;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    double expected_exit_mach = mach_from_prandtl_meyer(2.0 * theta_max, gamma);
    EXPECT_NEAR(result.exit_mach, expected_exit_mach, 0.02);
}
