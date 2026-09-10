#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/gas_dynamics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

// ============================================================
// Grid-convergence tests: the MoC solution must converge as the
// number of characteristics N increases. These tests protect the
// axisymmetric source terms and the analysis marching, whose errors
// are invisible at a single fixed N.
// ============================================================

static MocOptions make_options(MocFlowKind kind, double gamma, double theta_max, int n) {
    MocOptions opts;
    opts.flow_type = kind;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = n;
    opts.geometry.throat_radius = 1.0;
    return opts;
}

static MocResult solve_design(MocFlowKind kind, double gamma, double theta_max, int n) {
    MocNozzle solver(make_options(kind, gamma, theta_max, n));
    return solver.solve();
}

// Analyze a designed contour with centered-fan initialization (the design
// contour has a sharp throat corner, so the fan is the consistent start line).
// march_scheme defaults to AUTO (INVERSE for analysis, planar and axisymmetric); test_moc_inverse_march.cpp's InverseMarch.PlanarDesignRoundTrip and
// diagnosis.md A8/B.md's "Corrections after implementation" motivate forcing it
// explicitly where the two schemes are being compared.
static MocResult solve_analysis_of(const MocResult& design_result,
                                   MocFlowKind kind, double gamma, int n,
                                   MocMarchScheme scheme = MocMarchScheme::AUTO) {
    MocOptions opts = make_options(kind, gamma, 0.0, n);
    opts.mode = MocMode::ANALYSIS;
    opts.geometry.downstream_wall_curvature_radius = -1.0; // centered-fan init
    opts.nozzle_profile = design_result.profile;
    opts.march_scheme = scheme;
    MocNozzle solver(opts);
    return solver.solve();
}

// Isentropic 1D area-Mach relation A/A* (quasi-1D mass conservation).
static double area_ratio_1d(double mach, double gamma) {
    double t = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
    return std::pow(t, (gamma + 1.0) / (2.0 * (gamma - 1.0))) / mach;
}

// Supersonic Mach from A/A* (bisection; relation is monotone for M > 1).
static double mach_from_area_ratio_1d(double area_ratio, double gamma) {
    double lo = 1.0 + 1e-9, hi = 50.0;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (area_ratio_1d(mid, gamma) > area_ratio) hi = mid;
        else lo = mid;
    }
    return 0.5 * (lo + hi);
}

// Mirrors test_moc_inverse_march.cpp's helper of the same name (anonymous-namespace
// helpers are not shared across translation units): non-dimensional perfect-gas mass
// flux rho*V/(rho*_crit*a*_crit) as a function of Mach.
static double perfect_gas_mass_flux(double mach, double gamma) {
    const double exponent = -(gamma + 1.0) / (2.0 * (gamma - 1.0));
    return mach * std::pow(1.0 + 0.5 * (gamma - 1.0) * mach * mach, exponent);
}

// Relative error of the mass flow integrated across the exit plane against the 1-D
// critical mass flow through the throat, mirroring
// InverseMarch::exit_plane_mass_flow_error (test_moc_inverse_march.cpp) and
// MocNozzle::start_line_mass_flow_error's construction (moc_nozzle.cpp) applied to the
// exit plane instead of the start line.
static double exit_plane_mass_flow_error(const MocResult& result, double gamma, double r_throat) {
    const ExitPlane& ep = result.exit_plane;
    if (ep.y.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    double mdot = 0.0;
    for (size_t i = 1; i < ep.y.size(); i++) {
        const double dy = ep.y[i] - ep.y[i - 1];
        const double theta_avg = 0.5 * (ep.theta[i] + ep.theta[i - 1]);
        const double y_avg = 0.5 * (ep.y[i] + ep.y[i - 1]);
        const double flux = 0.5 * (perfect_gas_mass_flux(ep.mach[i], gamma)
                                   + perfect_gas_mass_flux(ep.mach[i - 1], gamma));
        // Exit plane is a vertical cut (constant x): dx=0, so the downstream normal flux
        // reduces to cos(theta_avg)*dy.
        const double normal_flux = std::cos(theta_avg) * dy;
        mdot += flux * normal_flux * (2.0 * M_PI * y_avg);
    }
    const double mdot_reference = perfect_gas_mass_flux(1.0, gamma) * (M_PI * r_throat * r_throat);
    return (mdot - mdot_reference) / mdot_reference;
}

// Area-averaged exit Mach against the 1-D value for the achieved area ratio, mirroring
// the same construction used throughout this file and test_moc_inverse_march.cpp.
static double area_mean_exit_mach(const MocResult& result) {
    const ExitPlane& ep = result.exit_plane;
    double mach_area = 0.0, area = 0.0;
    for (size_t k = 1; k < ep.y.size(); k++) {
        const double dA = M_PI * (ep.y[k] * ep.y[k] - ep.y[k - 1] * ep.y[k - 1]);
        mach_area += 0.5 * (ep.mach[k] + ep.mach[k - 1]) * dA;
        area += dA;
    }
    return mach_area / area;
}

// Wall Mach non-decreasing along the contour (monotone up to `tolerance`), and the
// smallest flow angle in the net bounded below by min_theta_lo (guards against a folded,
// runaway solution) -- the physical checks D.md item 1 / Addendum 2026-09-09 require of
// every default-throat conical solve. Shared by the item-1 default-options test and the
// flipped MocKlInitConvergence tripwires (item 2), which exercise the same throat
// geometry through a different options path (forced KLIEGEL_LEVINE start line).
//
// expect_dip additionally asserts min_theta < min_theta_hi_deg, i.e. that the axis
// compression (diagnosis.md A8) has actually formed by min_theta_hi_deg; pass false at
// grids coarse enough that it need not have formed yet (measured: exactly 0 deg -- the
// axis points' own pinned theta, not a real dip -- at N=8 for this throat).
static void check_conical_default_throat_physics(
    const MocResult& result, const std::string& tag,
    double wall_mach_tol, double min_theta_lo_deg, double min_theta_hi_deg,
    bool expect_dip = true)
{
    std::vector<double> wall_mach;
    for (size_t idx : result.net.wall_point_indices) {
        wall_mach.push_back(result.net.points[idx].mach);
    }
    for (size_t i = 1; i < wall_mach.size(); i++) {
        EXPECT_GE(wall_mach[i], wall_mach[i - 1] - wall_mach_tol)
            << tag << ": wall Mach dropped at wall point " << i;
    }

    ::testing::Test::RecordProperty("min_theta_deg" + tag, std::to_string(result.min_theta / DEG));
    ::testing::Test::RecordProperty("min_theta_x" + tag, std::to_string(result.min_theta_x));
    EXPECT_GT(result.min_theta, min_theta_lo_deg * DEG)
        << tag << ": flow angle dipped to " << result.min_theta / DEG
        << " deg at x=" << result.min_theta_x << " (folded solution?)";
    if (expect_dip) {
        EXPECT_LT(result.min_theta, min_theta_hi_deg * DEG)
            << tag << ": the r_arc=0.382 axis compression (diagnosis.md A8) is expected here "
            << "(got " << result.min_theta / DEG << " deg); if it has vanished, the geometry "
            << "or the physics changed -- update this test";
    }
}

// ------------------------------------------------------------
// Planar design: the Riemann invariants are exact, so the exit Mach
// follows algebraically from nu_exit = 2*theta_max at ANY resolution.
// Tolerance: the kernel's root solves (nu <-> M) converge to abstol=1e-10;
// 1e-6 leaves margin for accumulation over the march.
// ------------------------------------------------------------
TEST(MocConvergence, PlanarDesignExitMachExactAtAllN) {
    double gamma = 1.4;
    double theta_max = 15.0 * DEG;
    double expected = mach_from_prandtl_meyer(2.0 * theta_max, gamma);

    for (int n : {8, 16, 32, 64}) {
        auto result = solve_design(MocFlowKind::PLANAR, gamma, theta_max, n);
        ASSERT_TRUE(result.converged) << "N=" << n;
        EXPECT_NEAR(result.exit_mach, expected, 1e-6)
            << "Planar exit Mach must equal PM^-1(2*theta_max) at N=" << n;
    }
}

// ------------------------------------------------------------
// Planar design: uniform sonic throat + uniform exit flow means mass
// conservation forces the computed area ratio toward the 1D isentropic
// area-Mach relation as N grows. Measured baselines (this codebase):
// |AR - AR_1D| = 5.3e-3 at N=8, ~6e-4 at N=64 (first-order decay).
// ------------------------------------------------------------
TEST(MocConvergence, PlanarDesignAreaRatioConvergesTo1D) {
    double gamma = 1.4;
    double theta_max = 15.0 * DEG;

    auto coarse = solve_design(MocFlowKind::PLANAR, gamma, theta_max, 8);
    auto fine = solve_design(MocFlowKind::PLANAR, gamma, theta_max, 64);
    ASSERT_TRUE(coarse.converged);
    ASSERT_TRUE(fine.converged);

    double err_coarse = std::abs(coarse.area_ratio - area_ratio_1d(coarse.exit_mach, gamma));
    double err_fine = std::abs(fine.area_ratio - area_ratio_1d(fine.exit_mach, gamma));

    EXPECT_LT(err_fine, err_coarse / 3.0)
        << "Area-ratio error vs the 1D relation must shrink with N";
    EXPECT_LT(err_fine, 2e-3)
        << "At N=64 the area ratio should satisfy 1D mass conservation closely";
}

// ------------------------------------------------------------
// Planar design->analysis round trip, DIRECT vs INVERSE (D.md item 4 /
// Addendum 2026-09-09). Both schemes reanalyze the same designed contour with
// a centered-fan start line; DIRECT's wall solve (solve_wall_point_analysis,
// src/moc_nozzle.cpp) queries NozzleProfile::theta_at, which is piecewise
// constant per facet, while the inverse kernel's own wall solve interpolates
// the vertex angles linearly (wall_angle_at, src/moc_inverse_march.cpp -- see
// B.md "Corrections after implementation" item 4 and diagnosis.md A8). That
// difference is reported here, not worked around: fixing NozzleProfile::theta_at
// itself is a library change and out of this package's scope (a Package E
// candidate per the addendum).
//
// Measured on this tree: DIRECT fails to converge at both N=8 and N=32 (NEGATIVE_THETA,
// a per-facet kink accumulating into a small negative theta late in the march -- see
// solve_wall_point_analysis/NozzleProfile::theta_at above); the facet-quantization
// artifact is set by the wall polyline's own resolution (60-ish facets from the
// design), not by the characteristic spacing, so refining N does not cure it. INVERSE
// converges at both levels with the error shrinking (1.4e-2 -> 1.4e-3), consistent with
// InverseMarch.PlanarDesignRoundTrip's 4.4e-3 / 1.5e-3 at N=16/32
// (test_moc_inverse_march.cpp, same contour/round-trip construction).
// ------------------------------------------------------------
TEST(MocConvergence, PlanarRoundTripErrorShrinksWithN) {
    double gamma = 1.4;
    double theta_max = 15.0 * DEG;
    const int levels[2] = {8, 32};

    for (MocMarchScheme scheme : {MocMarchScheme::DIRECT, MocMarchScheme::INVERSE}) {
        const std::string scheme_name = (scheme == MocMarchScheme::DIRECT) ? "Direct" : "Inverse";
        double err[2] = {0.0, 0.0};
        bool ok[2] = {false, false};
        for (int i = 0; i < 2; i++) {
            auto design = solve_design(MocFlowKind::PLANAR, gamma, theta_max, levels[i]);
            ASSERT_TRUE(design.converged) << scheme_name << " design N=" << levels[i];
            auto analysis = solve_analysis_of(design, MocFlowKind::PLANAR, gamma, levels[i], scheme);
            ok[i] = analysis.converged;
            RecordProperty(scheme_name + "_converged_N" + std::to_string(levels[i]),
                           ok[i] ? "true" : "false");
            if (ok[i]) {
                err[i] = std::abs(analysis.exit_mach - design.exit_mach);
                RecordProperty(scheme_name + "_exit_mach_error_N" + std::to_string(levels[i]),
                               std::to_string(err[i]));
            } else {
                // Failure honesty holds for either scheme: a non-converging solve must
                // never silently succeed.
                EXPECT_NE(analysis.failure.code, MocErrorCode::NONE)
                    << scheme_name << " N=" << levels[i]
                    << ": a non-converging solve must carry a specific failure code";
                RecordProperty(scheme_name + "_failure_code_N" + std::to_string(levels[i]),
                               std::string(to_string(analysis.failure.code)));
            }
        }

        if (scheme == MocMarchScheme::INVERSE) {
            // The fix under test: the inverse kernel's wall solve does not carry the
            // facet-quantization artifact, so the round trip must actually converge and
            // improve with N -- there is no known limitation to carve out here.
            ASSERT_TRUE(ok[0]) << "INVERSE N=8 must converge";
            ASSERT_TRUE(ok[1]) << "INVERSE N=32 must converge";
            EXPECT_LT(err[1], err[0])
                << "INVERSE round-trip error must not grow with N (N=8 " << err[0]
                << ", N=32 " << err[1] << ")";
            EXPECT_LT(err[1], 5e-3) << "INVERSE N=32 round-trip error " << err[1]
                << " (see InverseMarch.PlanarDesignRoundTrip for the same check at N=16/32)";
        }
        // DIRECT is documented above, not gated: the facet-quantized wall angle is a
        // pre-existing library defect (NozzleProfile::theta_at) this package does not
        // fix. The RecordProperty entries above carry the measured numbers for the report.
    }
}

// ------------------------------------------------------------
// Axisymmetric design: no closed-form exit Mach exists, but the scheme must
// be self-convergent: successive grid refinements must produce shrinking
// increments (first-order marching gives ~halving per doubling of N).
// Measured baselines: |M16-M8| = 2.4e-2, |M32-M16| = 1.2e-2, |M64-M32| = 5.7e-3.
// A broken source term shows up here as stalled or erratic increments.
// ------------------------------------------------------------
TEST(MocConvergence, AxiDesignExitMachCauchyConvergence) {
    double gamma = 1.4;
    double theta_max = 12.0 * DEG;

    double mach[4];
    int levels[4] = {8, 16, 32, 64};
    for (int i = 0; i < 4; i++) {
        auto result = solve_design(MocFlowKind::AXISYMMETRIC, gamma, theta_max, levels[i]);
        ASSERT_TRUE(result.converged) << "N=" << levels[i];
        mach[i] = result.exit_mach;
    }

    double d1 = std::abs(mach[1] - mach[0]);
    double d2 = std::abs(mach[2] - mach[1]);
    double d3 = std::abs(mach[3] - mach[2]);

    EXPECT_LT(d2, 0.75 * d1) << "Refinement increments must shrink (got "
                             << d1 << " -> " << d2 << ")";
    EXPECT_LT(d3, 0.75 * d2) << "Refinement increments must shrink (got "
                             << d2 << " -> " << d3 << ")";
    EXPECT_LT(d3, 1e-2) << "Exit Mach should be nearly grid-independent by N=64";
}

// ------------------------------------------------------------
// Axisymmetric design vs 1D mass conservation: as for the planar case, the
// uniform exit flow must satisfy the 1D area-Mach relation in the converged
// limit.
// TODO: the solution currently carries a known ~+0.08 Mach bias vs the 1D
// relation that saturates with N (suspected N-independent error in the
// centered-fan initial-line construction: O(1) marching steps along the fan
// rays and the dropped source term on the first ray). Tighten this bound to
// ~1e-2 once that is fixed. The bound below still catches gross breakage
// (sign/scaling errors in the source terms shift it by >0.05).
// ------------------------------------------------------------
TEST(MocConvergence, AxiDesign1DConsistencyBounded) {
    double gamma = 1.4;
    double theta_max = 12.0 * DEG;

    auto result = solve_design(MocFlowKind::AXISYMMETRIC, gamma, theta_max, 64);
    ASSERT_TRUE(result.converged);

    double mach_1d = mach_from_area_ratio_1d(result.area_ratio, gamma);
    EXPECT_NEAR(result.exit_mach, mach_1d, 0.12)
        << "Axi design exit Mach vs 1D area-Mach relation (known bias ~+0.08)";
}

// ------------------------------------------------------------
// Axisymmetric design->analysis round trip (D.md item 3 / Addendum 2026-09-09): design
// with DESIGN_MIN_LENGTH (always DIRECT; unchanged) at N = 8, 16, 32, then analyze the
// design's own contour with the default (AUTO -> INVERSE for axisymmetric ANALYSIS)
// scheme. The DIRECT design's exit Mach carries a known +0.08 bias against the 1-D
// area-Mach relation (AxiDesign1DConsistencyBounded) and the inverse-march analysis of
// its contour gives a non-uniform exit plane, so the two exit Machs are not directly
// comparable (the addendum: "do not compare against the DIRECT design's exit Mach").
// What must hold for a correct isentropic analysis instead: the exit-plane mass flow
// matches the throat's, and the area-averaged exit Mach matches the 1-D value for the
// contour's area ratio, both increasingly well with N -- the same construction
// InverseMarch.AxiDesignRoundTripConservesMass (test_moc_inverse_march.cpp) uses for a
// Rao-schedule design; this test exercises the DESIGN_MIN_LENGTH schedule instead.
// ------------------------------------------------------------
TEST(MocConvergence, AxiRoundTripErrorShrinksWithN) {
    double gamma = 1.4;
    double theta_max = 12.0 * DEG;

    double mach_err[3] = {0, 0, 0};
    double mdot_err[3] = {0, 0, 0};
    const int levels[3] = {8, 16, 32};
    for (int i = 0; i < 3; i++) {
        auto design = solve_design(MocFlowKind::AXISYMMETRIC, gamma, theta_max, levels[i]);
        ASSERT_TRUE(design.converged) << "design N=" << levels[i];
        auto analysis = solve_analysis_of(design, MocFlowKind::AXISYMMETRIC, gamma, levels[i]);
        ASSERT_TRUE(analysis.converged) << "N=" << levels[i] << ": "
            << to_string(analysis.failure.code) << " -- " << analysis.failure.message;

        const ExitPlane& ep = analysis.exit_plane;
        ASSERT_GE(ep.y.size(), 2u) << "N=" << levels[i];
        const double mach_mean = area_mean_exit_mach(analysis);
        const double mach_1d = mach_from_area_ratio_1d(analysis.area_ratio, gamma);
        mach_err[i] = std::abs(mach_mean - mach_1d);
        mdot_err[i] = exit_plane_mass_flow_error(analysis, gamma, 1.0);
        RecordProperty("area_mean_exit_mach_N" + std::to_string(levels[i]), std::to_string(mach_mean));
        RecordProperty("mach_1d_N" + std::to_string(levels[i]), std::to_string(mach_1d));
        RecordProperty("mdot_err_N" + std::to_string(levels[i]), std::to_string(mdot_err[i]));
        EXPECT_LT(std::abs(mdot_err[i]), 0.03)
            << "N=" << levels[i] << ": exit-plane mass-flow error " << mdot_err[i];
    }

    EXPECT_LT(mach_err[1], mach_err[0])
        << "Axi round-trip error (area-mean exit Mach vs 1-D) must not grow N=8->16 (coarse "
        << mach_err[0] << ", mid " << mach_err[1] << ")";
    EXPECT_LT(mach_err[2], mach_err[1])
        << "Axi round-trip error (area-mean exit Mach vs 1-D) must not grow N=16->32 (mid "
        << mach_err[1] << ", fine " << mach_err[2] << ")";
    EXPECT_LT(mach_err[2], 5e-2);
    EXPECT_LT(std::abs(mdot_err[2]), std::abs(mdot_err[1]))
        << "Mass-flow error must shrink N=16->32 (mid " << mdot_err[1] << ", fine " << mdot_err[2] << ")";
}

// ------------------------------------------------------------
// Default-options conical convergence (D.md item 1 / Addendum 2026-09-09). Every
// option left at its default: MocMarchScheme::AUTO (resolves to INVERSE for
// axisymmetric ANALYSIS), MocStartLine::AUTO (resolves to KLIEGEL_LEVINE at this
// throat -- the wall-angle mismatch of 0.145 rad is within kl_max_wall_angle_error's
// default 0.25 rad), and NozzleGeometry::downstream_wall_curvature_radius's default
// (0.382), matching the arc radius passed to generate_conical_nozzle. This is what a
// caller gets by only setting num_characteristics/gamma/geometry and a contour --
// distinct from InverseMarch.ConicalConvergesAtCleanThroat and
// .ConicalDefaultThroatReachesExit (test_moc_inverse_march.cpp), which force the
// scheme and start line explicitly at r_arc = 2.0 and 0.382 respectively.
//
// diagnosis.md A8 documents a compression converging on the axis at this throat
// (r_arc = 0.382) near x ~ 3.4, which steepens under refinement and both AR = 4 and
// AR = 8 pass through (the arc/throat geometry is identical upstream of the exit cut
// for both; only where the exit is cut differs). That rules out a pointwise exit-Mach
// Cauchy check as the accuracy metric at AR = 4, whose exit sits just past the
// compression -- the addendum's replacement is used instead: exit-plane mass flow
// within 2% of the throat value at every N, closer at N=61 than N=31, exit-Mach
// Cauchy convergence only for AR = 8 (exit well beyond the compression), and
// min_theta recorded and bounded to (-8, -1) deg -- the compression is expected, and
// its disappearance or an unbounded runaway both indicate a regression.
// ------------------------------------------------------------
TEST(MocDefaultOptionsConvergence, ConicalDefaultThroatConvergesAcrossN) {
    const double gamma = 1.4;
    const int levels[3] = {15, 31, 61};

    for (double ar : {4.0, 8.0}) {
        double exit_mach[3] = {0, 0, 0};
        double mdot_err[3] = {0, 0, 0};
        bool ok[3] = {false, false, false};

        for (int li = 0; li < 3; li++) {
            const int n = levels[li];
            MocOptions opts;
            opts.flow_type = MocFlowKind::AXISYMMETRIC;
            opts.chemistry = GasChemistry::PERFECT_GAS;
            opts.mode = MocMode::ANALYSIS;
            opts.gamma = gamma;
            opts.num_characteristics = n;
            opts.geometry.throat_radius = 1.0;
            opts.geometry.downstream_wall_curvature_radius = 0.382; // == generate_conical_nozzle's arc below
            opts.nozzle_profile = NozzleProfile::generate_conical_nozzle(ar, 0.382, 1.0, 15.0, 60);
            // march_scheme, start_line: left at MocOptions defaults (AUTO, AUTO).

            MocNozzle solver(opts);
            MocResult result = solver.solve();
            ok[li] = result.converged;

            const std::string tag = "_AR" + std::to_string(static_cast<int>(ar)) + "_N" + std::to_string(n);
            RecordProperty("start_line_used" + tag,
                           result.init_diagnostics.start_line_used == MocStartLine::KLIEGEL_LEVINE
                               ? "KLIEGEL_LEVINE" : "CENTERED_FAN");

            EXPECT_TRUE(result.converged) << "AR=" << ar << " N=" << n << ": "
                << to_string(result.failure.code) << " -- " << result.failure.message;
            if (!result.converged) continue;
            EXPECT_TRUE(result.reached_exit_plane) << "AR=" << ar << " N=" << n;

            // 1e-2 matches MocKlInitConvergence.FirstWallHitsHaveMonotoneMach's tolerance
            // at this same throat (r_arc=0.382): the wall-consistency correction (Package A)
            // shrinks the pre-fix dip to a peak-to-trough of 0/0.007/0.014 at N=15/31/61 with
            // the largest single step -0.0062, comfortably under this bound. r_arc=2.0's
            // InverseMarch.ConicalConvergesAtCleanThroat is a clean throat with no such
            // wall-Mach dip and uses a tighter 1e-3, which does not apply here.
            check_conical_default_throat_physics(result, tag, /*wall_mach_tol=*/1e-2,
                                                  /*min_theta_lo_deg=*/-8.0, /*min_theta_hi_deg=*/-1.0);

            mdot_err[li] = exit_plane_mass_flow_error(result, gamma, 1.0);
            EXPECT_LT(std::abs(mdot_err[li]), 0.02)
                << "AR=" << ar << " N=" << n << ": exit-plane mass-flow error " << mdot_err[li];

            const double mach_1d = mach_from_area_ratio_1d(result.area_ratio, gamma);
            const double mach_mean = area_mean_exit_mach(result);
            EXPECT_NEAR(mach_mean, mach_1d, 0.05 * mach_1d) << "AR=" << ar << " N=" << n
                << ": area-mean exit Mach " << mach_mean << " vs 1-D " << mach_1d;

            exit_mach[li] = result.exit_mach;
            RecordProperty("exit_mach" + tag, std::to_string(exit_mach[li]));
            RecordProperty("mdot_err" + tag, std::to_string(mdot_err[li]));
        }

        if (ok[0] && ok[1] && ok[2]) {
            EXPECT_LT(std::abs(mdot_err[2]), std::abs(mdot_err[1]))
                << "AR=" << ar << ": mass-flow error must shrink N=31->61 (mid " << mdot_err[1]
                << ", fine " << mdot_err[2] << ")";
        }
        if (ar > 4.0 && ok[0] && ok[1] && ok[2]) {
            // Exit Mach Cauchy convergence only for AR=8, whose exit lies well beyond the
            // axis compression (diagnosis.md A8); AR=4's exit sits just past it and is not
            // a convergence metric (mass flow is, above).
            EXPECT_LT(std::abs(exit_mach[2] - exit_mach[1]), 0.75 * std::abs(exit_mach[1] - exit_mach[0]))
                << "AR=" << ar << ": exit Mach increments must shrink (M15=" << exit_mach[0]
                << ", M31=" << exit_mach[1] << ", M61=" << exit_mach[2] << ")";
        }
    }
}

// ============================================================
// Kliegel-Levine dual-family seeding grid-convergence tests (Phase 2).
//
// These exercise the KL-init path (positive downstream_wall_curvature_radius,
// the default, selected for axisymmetric ANALYSIS/DESIGN_RAO), distinct from the
// fan-init path (downstream_wall_curvature_radius <= 0) the tests above use.
//
// Before this phase, the KL start line seeded only a C+ per interior point (the
// wall point's C- was the only characteristic reaching the near-axis region),
// leaving that region of the initial mesh empty: the first C- to reach the axis
// descended the entire radius in one step, which blew up the axisymmetric source
// term and produced a Prandtl-Meyer-inversion or non-downstream-intersection
// failure for most (area_ratio, N) combinations (see
// instructions/moc_convergence_roadmap.md Sec 1 for the pre-fix baseline and
// trace evidence). The fix: MocInitialization::initialize_kliegel_levine now
// rigidly shifts the start line downstream of the raw sonic locus
// (MocOptions::initial_line_axial_shift) so it can be seeded with both
// characteristic families (CharacteristicNet::add_initial_data_line), plus two
// pairing-hardening fixes in MocNozzle::solve_characteristic_kernel and
// sort_plus_edges_by_proximity.
// ============================================================

// disable_mesh_control: set the three mesh-control factors so they never trigger (as
// tools/moc_sweep.cpp does for its "mesh_control off" rows), to isolate the initial-data-line
// behavior from the marching front's own refinement/coarsening. Defaults to false so
// existing callers keep the solver's default mesh control unchanged.
static MocResult solve_conical_kl_analysis(
    double area_ratio, int n, double gamma = 1.4, bool disable_mesh_control = false)
{
    MocOptions opts;
    opts.flow_type = MocFlowKind::AXISYMMETRIC;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::ANALYSIS;
    opts.gamma = gamma;
    opts.num_characteristics = n;
    opts.geometry.throat_radius = 1.0;
    // downstream_wall_curvature_radius left at its positive default: KL-init path.
    //
    // start_line is forced rather than left at AUTO: r_arc = 0.382 is the exact throat this
    // helper's own name calls out (the flagship case diagnosis.md Sec A2 documents the KL
    // series recovering under half the wall angle for), so under AUTO's new wall-angle-miss
    // check (MocOptions::kl_max_wall_angle_error, package A) every call here would silently
    // fall back to the centered fan -- which converges worse for this contour at low N than
    // the wall-consistency-corrected KL line does (measured: AR=2 N=15 fan hits
    // NEGATIVE_THETA where corrected KL reaches the exit plane). Forcing KLIEGEL_LEVINE
    // keeps this helper doing what its name and every caller's comments say it does, and the
    // generous threshold below is wide enough for r_arc = 0.382's ~0.145 rad raw miss so the
    // wall-consistency correction is applied instead of the forced path failing honestly
    // (see KliegelLevineInitialization.ForcedSeriesFailsHonestlyBelowThreshold for that
    // failure mode with a threshold tight enough to trigger it).
    opts.start_line = MocStartLine::KLIEGEL_LEVINE;
    opts.kl_max_wall_angle_error = 0.2; // below the 0.25 default; kept explicit so the helper is independent of it
    if (disable_mesh_control) {
        opts.max_front_spacing_factor = 1e9;
        opts.min_front_spacing_factor = 1e-9;
        opts.max_cell_aspect_ratio = 1e9;
    }
    opts.nozzle_profile = NozzleProfile::generate_conical_nozzle(area_ratio, 0.382, 1.0, 15.0, 60);
    MocNozzle solver(opts);
    return solver.solve();
}

// AR=2 converged at every grid level even before this phase; it is a regression
// guard that dual-family seeding and the pairing-hardening fixes did not break
// the already-working case, and the computed area ratio should track the
// requested contour AR increasingly closely (mass conservation through a
// correctly-marched supersonic flow field) as N grows.
// Judged on exit_coverage, not area_ratio.
//
// area_ratio is read off the last point of the outflow staircase, a ragged boundary of
// independently terminated chains, so it lands short of the contour even for a march that
// reached the exit plane. Mesh control moves area_ratio by up to 4% while exit_mach is
// bit-identical (instructions/moc_convergence_roadmap.md Sec 4), which is what disqualifies
// it as the accuracy metric. It is kept below as a loose secondary check with a tolerance
// that reflects what the staircase can actually deliver; exit_coverage carries the
// assertion that matters.
//
// Coverage floors and area_ratio tolerances updated for Package A's wall-consistency
// correction (initialize_kliegel_levine, multiplicative variant): forcing every KL-seeded
// point's flow angle to be consistent with the contour's own wall boundary condition
// changes the near-wall theta distribution enough to move this coarse-grid AR=2 case's
// coverage down a little, even though the correction's target is a different regime (the
// N=15/31/61 wall-Mach-monotonicity and mass-flow defects the correction fixes -- see
// FirstWallHitsHaveMonotoneMach and StartLineMassFlowWithinTwoPercent below). Measured
// 2026-09-02, after the correction: coverage 0.9486/0.9756/0.9942 and area_ratio
// 1.761/1.863/1.934 at N=8/15/31 (was 0.99+/1.93-1.94 uniformly before). The trend that
// actually matters -- coverage improving monotonically with N, toward the AR=2 target --
// still holds and is checked explicitly below.
TEST(MocKlInitConvergence, ConicalAR2ReachesExitPlaneAcrossN) {
    const int levels[3] = {8, 15, 31};
    const double coverage_floor[3] = {0.94, 0.97, 0.99};
    const double area_ratio_tol[3] = {0.25, 0.15, 0.10};
    double previous_coverage = 0.0;
    for (int i = 0; i < 3; i++) {
        const int n = levels[i];
        auto result = solve_conical_kl_analysis(2.0, n);
        ASSERT_TRUE(result.converged) << "N=" << n << ": " << result.failure.message;
        EXPECT_GT(result.exit_mach, 1.0) << "N=" << n;
        EXPECT_TRUE(result.reached_exit_plane) << "N=" << n;
        EXPECT_GE(result.exit_coverage, coverage_floor[i])
            << "N=" << n << ": the march must reach the requested AR=2 exit radius";
        EXPECT_NEAR(result.area_ratio, 2.0, area_ratio_tol[i])
            << "N=" << n << ": staircase readout of the achieved area ratio";
        EXPECT_GE(result.exit_coverage, previous_coverage)
            << "N=" << n << ": coverage must not get worse under refinement";
        previous_coverage = result.exit_coverage;
    }
}

// AR=4 N=8, forced KLIEGEL_LEVINE (see solve_conical_kl_analysis): before Package A's
// wall-consistency correction this stopped at exit_coverage ~0.81 with NEGATIVE_THETA; a
// coverage floor was tracked instead of `converged` because the case sat on the
// convergence boundary of the near-axis void (roadmap Sec 4). Package A (wall-consistent
// KL line) and Package B (inverse march) together move this configuration off that
// boundary. Flipped 2026-09-09 (D.md item 2 / Addendum): asserted on convergence and the
// physical checks of ConicalDefaultThroatConvergesAcrossN (mass conservation, monotone
// wall Mach, the bounded axis compression) rather than a coverage floor, which "exists
// only because nothing converged" (D.md item 2).
TEST(MocKlInitConvergence, ConicalAR4AtCoarseNReachesRecordedCoverage) {
    auto result = solve_conical_kl_analysis(4.0, 8);

    RecordProperty("exit_coverage", std::to_string(result.exit_coverage));
    RecordProperty("exit_mach", std::to_string(result.exit_mach));
    RecordProperty("failure_code", std::string(to_string(result.failure.code)));

    ASSERT_TRUE(result.converged) << to_string(result.failure.code)
        << " -- " << result.failure.message;
    EXPECT_TRUE(result.reached_exit_plane);
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.area_ratio, 2.0);
    // N=8 is coarser than the N=15/31/61 grids diagnosis.md A8 measured the compression
    // on (-3.1/-5.0 deg at N=31/61); measured here: min_theta is exactly 0 deg (the axis
    // points' own pinned theta -- no dip has formed yet), so expect_dip=false skips that
    // assertion while the lower bound still guards against a runaway fold.
    check_conical_default_throat_physics(result, "_AR4_N8", /*wall_mach_tol=*/1e-2,
                                          /*min_theta_lo_deg=*/-90.0, /*min_theta_hi_deg=*/-1.0,
                                          /*expect_dip=*/false);
    EXPECT_LT(std::abs(exit_plane_mass_flow_error(result, 1.4, 1.0)), 0.03)
        << "AR=4 N=8 exit-plane mass-flow error";
}

// Flipped 2026-09-09 (D.md item 2 / Addendum): with the wall-consistent KL line
// (Package A) and the inverse march (Package B), both AR=4 and AR=8 reach the exit plane
// at every N and satisfy the same physical checks as ConicalDefaultThroatConvergesAcrossN
// -- mass conservation across the exit plane and a bounded axis compression (not a
// coverage floor, which "exists only because nothing converged", D.md item 2). This
// mirrors that test but through solve_conical_kl_analysis's forced-KLIEGEL_LEVINE,
// default-mesh-control options path rather than fully-default options.
TEST(MocKlInitConvergence, ConicalAR4FinerNAndAR8Converge) {
    for (int n : {15, 31}) {
        auto result = solve_conical_kl_analysis(4.0, n);
        ASSERT_TRUE(result.converged) << "N=" << n << " AR=4: " << result.failure.message;
        EXPECT_TRUE(result.reached_exit_plane) << "N=" << n;
        check_conical_default_throat_physics(result, "_AR4_N" + std::to_string(n),
                                              /*wall_mach_tol=*/1e-2,
                                              /*min_theta_lo_deg=*/-8.0, /*min_theta_hi_deg=*/-1.0);
        EXPECT_LT(std::abs(exit_plane_mass_flow_error(result, 1.4, 1.0)), 0.02)
            << "AR=4 N=" << n << " exit-plane mass-flow error";
    }
    for (int n : {8, 15, 31}) {
        auto result = solve_conical_kl_analysis(8.0, n);
        ASSERT_TRUE(result.converged) << "N=" << n << " AR=8: " << result.failure.message;
        EXPECT_TRUE(result.reached_exit_plane) << "N=" << n;
        // N=8 is coarser than diagnosis.md A8's measured grid (N=31/61); the compression
        // dip need not have fully formed there (see the AR=4 N=8 comment above), so its
        // lower bound is left open at N=8 and tightened at N=15/31.
        const bool coarse = (n == 8);
        check_conical_default_throat_physics(result, "_AR8_N" + std::to_string(n),
                                              /*wall_mach_tol=*/1e-2,
                                              /*min_theta_lo_deg=*/coarse ? -90.0 : -8.0,
                                              /*min_theta_hi_deg=*/-1.0,
                                              /*expect_dip=*/!coarse);
        EXPECT_LT(std::abs(exit_plane_mass_flow_error(result, 1.4, 1.0)), 0.02)
            << "AR=8 N=" << n << " exit-plane mass-flow error";
    }
}

// ------------------------------------------------------------
// Wall-consistency correction (Package A): the KL start line's wall end is now forced to
// match the contour's own wall angle (initialize_kliegel_levine), instead of the series'
// raw (and, at this throat, badly wrong -- see KliegelLevineClosedForm.SeriesMissesWallAngleAtSmallR)
// value. diagnosis.md Sec A1 traces the pre-fix defect: K+ on the near-wall interior points
// is 8-10 deg too negative, the first wall solve over-expands, and the wall Mach *decreases*
// over the next several wall points before recovering (dips of 3/7/12 wall points at
// N=15/31/61, growing with N -- refining the grid does not cure it, since the defect is a
// property of the series, not the mesh).
// ------------------------------------------------------------

// Sorted by x: net.wall_points() is built by iterating chain_metadata in chain-creation
// order, not by x, and the KL line's dual-family seeding creates one C+ chain per interior
// data-line point in axis-to-wall order -- the axis-seeded chain has to travel the farthest
// (highest mu, near-sonic) and so is the *last* to reach the wall, while the near-wall chain
// arrives almost immediately. So the raw vector is closer to reverse-x order than to
// march order; every comparison here re-sorts by x first, which is what "wall Mach along
// the nozzle" means physically and what diagnosis.md's own point-by-point table reports.
static std::vector<CharacteristicPoint> wall_points_by_x(const MocResult& result) {
    // wall_point_indices is filled by both kernels (the INVERSE march builds no chains, so
    // net.wall_points(), which walks chain terminations, is empty there).
    std::vector<CharacteristicPoint> pts;
    for (size_t idx : result.net.wall_point_indices) pts.push_back(result.net.points[idx]);
    std::sort(pts.begin(), pts.end(),
        [](const CharacteristicPoint& a, const CharacteristicPoint& b) { return a.x < b.x; });
    return pts;
}

// This is the test that matters (A.md): before this package, at r_arc = 0.382, AR = 4, mesh
// control non-binding (isolating the initial-data-line defect from the marching front's own
// refinement), the first wall-adjacent points' Mach *decreases* for several points before
// turning around. Measured on this tree pre-fix: dips 3/7/13 wall points deep at N=15/31/61,
// bottoming 0.045/0.064/0.093 in Mach below the first point. The multiplicative
// wall-consistency correction (see initialize_kliegel_levine) shrinks that to 0/2/4 points
// (peak-to-trough 0/0.007/0.014), and the largest single-step regression across all three N
// is -0.0062 -- two orders of magnitude below the pre-fix largest step of roughly -0.03.
// tolerance is set to 0.01, comfortably above that residual (attributable to the line still
// being a finite-order series, not to the K+ deficit this package targets) and comfortably
// below both the pre-fix defect and the p=1 additive variant's largest single step (-0.0146,
// see the package report), so a regression back toward either would fail this test.
TEST(MocKlInitConvergence, FirstWallHitsHaveMonotoneMach) {
    constexpr double tolerance = 0.01;
    for (int n : {15, 31, 61}) {
        auto result = solve_conical_kl_analysis(4.0, n, 1.4, /*disable_mesh_control=*/true);
        std::vector<CharacteristicPoint> wall_pts = wall_points_by_x(result);
        const size_t window = std::min<size_t>(12, wall_pts.size());
        ASSERT_GE(wall_pts.size(), 8u) << "N=" << n << ": too few wall points to judge monotonicity";
        for (size_t i = 1; i < window; i++) {
            EXPECT_GE(wall_pts[i].mach, wall_pts[i - 1].mach - tolerance)
                << "N=" << n << ": wall Mach dropped from " << wall_pts[i - 1].mach
                << " to " << wall_pts[i].mach << " between wall points " << (i - 1)
                << " and " << i << " (x=" << wall_pts[i - 1].x << " -> " << wall_pts[i].x << ")";
        }
    }
}

// The start line's own mass flow should match the 1-D critical mass flow through the throat
// to within a couple of percent; before this package the raw (uncorrected) series carried
// -4.3% to -4.6% at this throat (the wall-angle deficit biases the near-wall velocity
// components), comfortably outside any reasonable tolerance.
TEST(MocKlInitConvergence, StartLineMassFlowWithinTwoPercent) {
    auto result = solve_conical_kl_analysis(4.0, 31, 1.4, /*disable_mesh_control=*/true);
    EXPECT_LT(std::abs(result.init_diagnostics.mass_flow_error), 0.02)
        << "mass_flow_error = " << result.init_diagnostics.mass_flow_error;
}
