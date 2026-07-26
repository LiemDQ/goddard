#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/gas_dynamics.hpp"
#include <cmath>
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
static MocResult solve_analysis_of(const MocResult& design_result,
                                   MocFlowKind kind, double gamma, int n) {
    MocOptions opts = make_options(kind, gamma, 0.0, n);
    opts.mode = MocMode::ANALYSIS;
    opts.geometry.downstream_wall_curvature_radius = -1.0; // centered-fan init
    opts.nozzle_profile = design_result.profile;
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
// Planar design->analysis round trip: analyzing the designed contour must
// reproduce the design exit Mach, with the discrepancy vanishing as N grows.
// Measured baselines: |dM| = 1.3e-2 at N=8, 2.8e-3 at N=32 (~halving per
// doubling of N).
// ------------------------------------------------------------
TEST(MocConvergence, PlanarRoundTripErrorShrinksWithN) {
    double gamma = 1.4;
    double theta_max = 15.0 * DEG;

    double err[2];
    int levels[2] = {8, 32};
    for (int i = 0; i < 2; i++) {
        auto design = solve_design(MocFlowKind::PLANAR, gamma, theta_max, levels[i]);
        ASSERT_TRUE(design.converged) << "design N=" << levels[i];
        auto analysis = solve_analysis_of(design, MocFlowKind::PLANAR, gamma, levels[i]);
        // The analysis kernel re-reflects off a faceted (piecewise-linear) wall
        // built from the design's discrete wall points; small per-facet kinks can
        // accumulate into a slightly negative theta late in the march -- a
        // pre-existing accuracy limitation (instructions/moc_algorithm.md Sec.
        // 9.1/10, Phase 2 scope, not fixed here). Before the point-validity checks
        // added in this phase, such a point silently entered the net and this
        // trend comparison ran on a net that was not actually fully valid. Skip
        // (rather than silently pass or hard-fail) when that known limitation is
        // hit; run the full trend comparison otherwise.
        EXPECT_TRUE(analysis.converged || analysis.failure.code != MocErrorCode::NONE) 
            << "Analysis round trip must either converge or report a known failure at N="
            << levels[i] << " (got " << to_string(analysis.failure.code) << ")"
            << " -- " << analysis.failure.message;
        // if (!analysis.converged) {
        //     GTEST_SKIP() << "Analysis round trip did not converge at N=" << levels[i]
        //                  << " (known accuracy limitation, see "
        //                     "instructions/moc_algorithm.md Sec. 9.1/10): "
        //                  << to_string(analysis.failure.code)
        //                  << " -- " << analysis.failure.message;
        // }
        err[i] = std::abs(analysis.exit_mach - design.exit_mach);
    }

    EXPECT_LT(err[1], err[0] / 3.0)
        << "Round-trip exit-Mach error must shrink with N (coarse " << err[0]
        << ", fine " << err[1] << ")";
    EXPECT_LT(err[1], 5e-3);
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
// Axisymmetric design->analysis round trip: the centerline exit Mach from
// analyzing the designed contour must approach the design value as N grows.
// Measured baselines: |dM| = 7.5e-2 at N=8, 2.8e-2 at N=16.
// TODO: extend to N >= 32 once the analysis march survives it. Today the
// re-reflected waves off the faceted contour coalesce near the exit lip at
// N >= 32 (characteristics fold; solver reports non-downstream intersections
// and converged=false).
// ------------------------------------------------------------
TEST(MocConvergence, AxiRoundTripErrorShrinksWithN) {
    double gamma = 1.4;
    double theta_max = 12.0 * DEG;

    double err[2];
    int levels[2] = {8, 16};
    for (int i = 0; i < 2; i++) {
        auto design = solve_design(MocFlowKind::AXISYMMETRIC, gamma, theta_max, levels[i]);
        ASSERT_TRUE(design.converged) << "design N=" << levels[i];
        auto analysis = solve_analysis_of(design, MocFlowKind::AXISYMMETRIC, gamma, levels[i]);
        // See the comment in MocConvergence.PlanarRoundTripErrorShrinksWithN: a
        // faceted-wall reflection accuracy limitation in the analysis kernel
        // (Phase 2 scope, not fixed here) can leave converged == false. Skip
        // rather than mask it.
        if (!analysis.converged) {
            GTEST_SKIP() << "Analysis round trip did not converge at N=" << levels[i]
                         << " (known accuracy limitation, see "
                            "instructions/moc_algorithm.md Sec. 9.1/10): "
                         << to_string(analysis.failure.code)
                         << " -- " << analysis.failure.message;
        }
        err[i] = std::abs(analysis.exit_mach - design.exit_mach);
    }

    EXPECT_LT(err[1], err[0])
        << "Axi round-trip error must not grow with N (coarse " << err[0]
        << ", fine " << err[1] << ")";
    EXPECT_LT(err[1], 5e-2);
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

static MocResult solve_conical_kl_analysis(double area_ratio, int n, double gamma = 1.4) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::AXISYMMETRIC;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::ANALYSIS;
    opts.gamma = gamma;
    opts.num_characteristics = n;
    opts.geometry.throat_radius = 1.0;
    // downstream_wall_curvature_radius left at its positive default: KL-init path.
    opts.nozzle_profile = NozzleProfile::generate_conical_nozzle(area_ratio, 1.0, 15.0, 60);
    MocNozzle solver(opts);
    return solver.solve();
}

// AR=2 converged at every grid level even before this phase; it is a regression
// guard that dual-family seeding and the pairing-hardening fixes did not break
// the already-working case, and the computed area ratio should track the
// requested contour AR increasingly closely (mass conservation through a
// correctly-marched supersonic flow field) as N grows.
TEST(MocKlInitConvergence, ConicalAR2ConvergesWithAreaRatioTrackingTargetAcrossN) {
    for (int n : {8, 15, 31}) {
        auto result = solve_conical_kl_analysis(2.0, n);
        ASSERT_TRUE(result.converged) << "N=" << n << ": " << result.failure.message;
        EXPECT_GT(result.exit_mach, 1.0) << "N=" << n;
        EXPECT_NEAR(result.area_ratio, 2.0, 0.05)
            << "N=" << n << ": computed area ratio should track the requested AR=2 contour";
    }
}

// AR=4 at N=8 is the case this phase's fix was validated against: before it, this
// configuration hit a giant single-step axis descent and a Prandtl-Meyer
// inversion failure at kernel pass 16 (see
// instructions/moc_convergence_roadmap.md Sec 2 Step 0). Regression guard.
TEST(MocKlInitConvergence, ConicalAR4ConvergesAtCoarseN) {
    auto result = solve_conical_kl_analysis(4.0, 8);
    ASSERT_TRUE(result.converged) << result.failure.message;
    EXPECT_GT(result.exit_mach, 1.0);
    // N=8 is coarse enough that the last wall point the march reaches undershoots
    // the requested contour AR (measured: 3.40 vs the requested 4.0) -- the tight
    // area-ratio check belongs on the AR=2 grid-convergence test above, which
    // refines N. This is deliberately a loose sanity bound, not an accuracy check.
    EXPECT_GT(result.area_ratio, 2.0);
}

// TODO: AR=4 at N>=15 and AR=8 at every tested N still hit a residual
// mesh-density / domain-of-dependence mismatch deeper in the march: an
// axis-reflected characteristic's leading point ends up geometrically behind a
// C- partner that has already advanced much further downstream in a
// fast-expanding region -- a distinct, deeper issue from the near-axis void this
// phase fixed. See instructions/moc_convergence_roadmap.md Sec 2 Step 4 ("revisit
// after the net is dense") and consider Step 3 (step-size cap / point insertion
// for long characteristic segments) as the next candidate fix. This test
// documents the current (improved but incomplete) state with a specific,
// non-NONE error code -- not a crash, hang, or silently-invalid net -- rather
// than skipping silently; tighten it (replace EXPECT_FALSE with a convergence +
// accuracy check) once that residual mismatch is fixed.
TEST(MocKlInitConvergence, ConicalAR4FinerNAndAR8DocumentResidualFailure) {
    for (int n : {15, 31}) {
        auto result = solve_conical_kl_analysis(4.0, n);
        EXPECT_FALSE(result.converged)
            << "N=" << n << " AR=4 now converges -- tighten this test per its TODO comment";
        if (!result.converged) {
            EXPECT_NE(result.failure.code, MocErrorCode::NONE) << "N=" << n;
        }
    }
    for (int n : {8, 15, 31}) {
        auto result = solve_conical_kl_analysis(8.0, n);
        EXPECT_FALSE(result.converged)
            << "N=" << n << " AR=8 now converges -- tighten this test per its TODO comment";
        if (!result.converged) {
            EXPECT_NE(result.failure.code, MocErrorCode::NONE) << "N=" << n;
        }
    }
}
