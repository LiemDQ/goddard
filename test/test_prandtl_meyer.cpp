#include "goddard/prandtlmeyer.hpp"
#include <cmath>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

// --- Known values from gas dynamics tables ---

TEST(PrandtlMeyer, AtMach1IsZero) {
    EXPECT_NEAR(prandtl_meyer(1.0, 1.4), 0.0, 1e-12);
}

TEST(PrandtlMeyer, KnownValueGamma14) {
    // M=2, gamma=1.4: nu = 26.3798 degrees
    double nu = prandtl_meyer(2.0, 1.4);
    EXPECT_NEAR(nu / DEG, 26.3798, 0.001);
}

TEST(PrandtlMeyer, KnownValueMach3) {
    // M=3, gamma=1.4: nu = 49.7573 degrees
    double nu = prandtl_meyer(3.0, 1.4);
    EXPECT_NEAR(nu / DEG, 49.7573, 0.001);
}

TEST(PrandtlMeyer, MaxAngleGamma14) {
    // nu_max = (sqrt((g+1)/(g-1)) - 1) * 90 degrees
    // For gamma=1.4: nu_max = (sqrt(6) - 1) * 90 = 130.45 degrees
    double nu_max = (std::sqrt(6.0) - 1.0) * 90.0;
    // At very high Mach, nu should approach nu_max
    double nu = prandtl_meyer(1000.0, 1.4);
    EXPECT_NEAR(nu / DEG, nu_max, 0.5);
}

// --- Round trip: nu -> M -> nu ---

class PrandtlMeyerRoundTrip : public ::testing::TestWithParam<double> {};

TEST_P(PrandtlMeyerRoundTrip, ForwardInverseConsistency) {
    double mach = GetParam();
    double gamma = 1.4;
    double nu = prandtl_meyer(mach, gamma);
    double mach_recovered = mach_from_prandtl_meyer(nu, gamma, mach);
    EXPECT_NEAR(mach_recovered, mach, 1e-8);
}

INSTANTIATE_TEST_SUITE_P(
    MachValues,
    PrandtlMeyerRoundTrip,
    ::testing::Values(1.01, 1.1, 1.5, 2.0, 3.0, 5.0, 10.0)
);

// --- Round trip with different gamma values ---

class PrandtlMeyerGammaRoundTrip
    : public ::testing::TestWithParam<std::pair<double, double>> {};

TEST_P(PrandtlMeyerGammaRoundTrip, ForwardInverseConsistency) {
    auto [mach, gamma] = GetParam();
    double nu = prandtl_meyer(mach, gamma);
    double mach_recovered = mach_from_prandtl_meyer(nu, gamma, mach);
    EXPECT_NEAR(mach_recovered, mach, 1e-8);
}

INSTANTIATE_TEST_SUITE_P(
    MachGammaValues,
    PrandtlMeyerGammaRoundTrip,
    ::testing::Values(
        std::make_pair(2.0, 1.2),
        std::make_pair(2.0, 1.3),
        std::make_pair(2.0, 1.4),
        std::make_pair(2.0, 1.667),
        std::make_pair(3.0, 1.2),
        std::make_pair(3.0, 1.667)
    )
);

// --- Derivative check against finite differences ---

TEST(PrandtlMeyer, DerivativeFiniteDifference) {
    double gamma = 1.4;
    double eps = 1e-6;

    for (double mach : {1.5, 2.0, 3.0, 5.0}) {
        double analytical = prandtl_meyer_derivative(mach, gamma);
        double numerical = (prandtl_meyer(mach + eps, gamma)
                            - prandtl_meyer(mach - eps, gamma)) / (2.0 * eps);
        EXPECT_NEAR(analytical, numerical, 1e-6)
            << "Derivative mismatch at M=" << mach;
    }
}

// --- Monotonicity with gamma ---

TEST(PrandtlMeyer, MonotonicInGamma) {
    double mach = 2.5;
    double nu_12 = prandtl_meyer(mach, 1.2);
    double nu_14 = prandtl_meyer(mach, 1.4);
    double nu_16 = prandtl_meyer(mach, 1.667);
    // Higher gamma -> lower nu at same Mach
    EXPECT_GT(nu_12, nu_14);
    EXPECT_GT(nu_14, nu_16);
}

// --- Cold start (no guess) ---

TEST(PrandtlMeyer, InverseWithoutGuess) {
    double gamma = 1.4;
    double nu = prandtl_meyer(2.5, gamma);
    // mach_guess = 0.0 triggers the default heuristic
    double mach_recovered = mach_from_prandtl_meyer(nu, gamma);
    EXPECT_NEAR(mach_recovered, 2.5, 1e-8);
}
