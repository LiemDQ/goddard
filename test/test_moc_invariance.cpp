#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include <cmath>
#include "gtest/gtest.h"

using namespace Goddard;

// ============================================================
// Bit-level invariance of the minimum-length design path.
//
// instructions/moc_convergence_roadmap.md Sec 7 requires DESIGN_MIN_LENGTH
// planar and axisymmetric output to be unchanged by any work on the
// Kliegel-Levine start line: min-length always initializes with a centered
// expansion fan (moc_nozzle.cpp, generate_initial_data_line) and is excluded
// by mode from mesh control (control_front_spacing). That invariance was not
// enforced by any test, so a KL-path change could silently move the one
// configuration validated against Anderson Table 11.1.
//
// The references were captured from the tree at the baseline commit with
// %.17g and are exact to the last bit there. The comparison tolerance is
// 1e-12 relative rather than exact equality so the test survives a compiler
// or optimization-flag change without going quiet about a real regression;
// any actual change to the march moves these numbers far more than that.
// ============================================================

namespace {

struct WallPoint { double x, y; };

struct GoldenCase {
    const char* name;
    MocFlowKind flow;
    double gamma;
    double theta_max_deg;
    int num_characteristics;
    double exit_mach;
    double area_ratio;
    double nozzle_length;
    std::vector<WallPoint> wall;
};

MocResult solve_golden(const GoldenCase& c) {
    MocOptions o;
    o.flow_type = c.flow;
    o.chemistry = GasChemistry::PERFECT_GAS;
    o.mode = MocMode::DESIGN_MIN_LENGTH;
    o.gamma = c.gamma;
    o.theta_max = c.theta_max_deg * M_PI / 180.0;
    o.num_characteristics = c.num_characteristics;
    o.geometry.throat_radius = 1.0;
    MocNozzle solver(o);
    return solver.solve();
}

// Relative comparison, falling back to absolute near zero (wall[0] is exactly
// the throat lip at x = 0).
void expect_matches(double actual, double expected, const char* what, size_t idx) {
    const double scale = std::max(1.0, std::abs(expected));
    EXPECT_NEAR(actual, expected, 1e-12 * scale)
        << what << " at index " << idx << " moved: got " << actual
        << ", golden " << expected;
}

std::vector<GoldenCase> golden_cases() {
    return {
        GoldenCase{
            "planar", MocFlowKind::PLANAR, 1.4, 15.0, 20,
            2.1339050332318474, 1.8894505839966758, 5.7525398692310965,
            {
                {0.0,                  1.0},
                {0.84848684094848181,  1.2267562636210831},
                {1.4393456899580401,   1.3799174783727617},
                {1.7157072809023792,   1.4475259131248843},
                {1.9457014930234722,   1.5004601690089807},
                {2.1581412210345166,   1.5462967986351306},
                {2.3633269123960465,   1.5876327152960039},
                {2.5663458100988668,   1.62564355910336},
                {2.7702127376501844,   1.6609274072012852},
                {2.97695386701223,     1.6937961022812793},
                {3.1880680299060695,   1.7243987739468074},
                {3.4047521753479462,   1.7527815177819939},
                {3.6280231724077194,   1.7789188686949713},
                {3.8587889513427034,   1.8027313780510814},
                {4.097892821716087,    1.8240957364637855},
                {4.3461426830770007,   1.8428506051261013},
                {4.6043313015788492,   1.8587998167277315},
                {4.873251102039065,    1.8717138660588886},
                {5.1537054986910533,   1.8813302209203042},
                {5.4465179990578525,   1.8873527680692519},
                {5.7525398692310965,   1.8894505839966758},
            }},
        GoldenCase{
            "axisymmetric", MocFlowKind::AXISYMMETRIC, 1.4, 12.0, 20,
            3.0322104501599916, 4.0581414538115386, 8.3750502491154535,
            {
                {0.0,                  1.0},
                {0.91502911251641772,  1.1939947432958415},
                {1.8244493139662288,   1.3810985288530089},
                {2.2817368434035665,   1.4699644237231251},
                {2.642309415529887,    1.5359398920539205},
                {2.9762877751102863,   1.5932710189157973},
                {3.2955438466260936,   1.6444767407150247},
                {3.611497178272292,    1.6916043046251925},
                {3.92810533556885,     1.7352853749458101},
                {4.2488171773160017,   1.7759533575284077},
                {4.5758241760473801,   1.8137802259827738},
                {4.9108522294067667,   1.8488157449385167},
                {5.2552234457721951,   1.8810141710246369},
                {5.6099884914708635,   1.9102631906388341},
                {5.9759709108738948,   1.9363991426644018},
                {6.3538013410468865,   1.9592191070366725},
                {6.7439452489488216,   1.9784907298895873},
                {7.1467530654196851,   1.9939613396819718},
                {7.5625585713003733,   2.0053656020351727},
                {7.9875177651361913,   2.0123576720398484},
                {8.3750502491154535,   2.0144829246760914},
            }},
    };
}

} // namespace

TEST(MocInvariance, DesignMinLengthUnchanged) {
    for (const GoldenCase& c : golden_cases()) {
        SCOPED_TRACE(c.name);
        MocResult r = solve_golden(c);

        ASSERT_TRUE(r.converged) << c.name << ": " << r.failure.message;
        expect_matches(r.exit_mach, c.exit_mach, "exit_mach", 0);
        expect_matches(r.area_ratio, c.area_ratio, "area_ratio", 0);
        expect_matches(r.nozzle_length, c.nozzle_length, "nozzle_length", 0);

        ASSERT_EQ(r.net.wall_x.size(), c.wall.size()) << c.name << ": wall point count changed";
        ASSERT_EQ(r.net.wall_y.size(), c.wall.size()) << c.name;
        for (size_t i = 0; i < c.wall.size(); i++) {
            expect_matches(r.net.wall_x[i], c.wall[i].x, "wall_x", i);
            expect_matches(r.net.wall_y[i], c.wall[i].y, "wall_y", i);
        }
    }
}

// Mesh control is excluded by mode for DESIGN_MIN_LENGTH (control_front_spacing
// returns early), so changing its knobs must not move the design at all. This is
// the cheap direct check of that exclusion, which the golden values above would
// only catch indirectly.
TEST(MocInvariance, DesignMinLengthIgnoresMeshControlOptions) {
    for (const GoldenCase& c : golden_cases()) {
        SCOPED_TRACE(c.name);
        MocResult baseline = solve_golden(c);

        MocOptions o;
        o.flow_type = c.flow;
        o.chemistry = GasChemistry::PERFECT_GAS;
        o.mode = MocMode::DESIGN_MIN_LENGTH;
        o.gamma = c.gamma;
        o.theta_max = c.theta_max_deg * M_PI / 180.0;
        o.num_characteristics = c.num_characteristics;
        o.geometry.throat_radius = 1.0;
        // Deliberately extreme, but still inside validate_moc_options' bounds.
        o.max_front_spacing_factor = 1.01;
        o.min_front_spacing_factor = 0.49;
        o.max_cell_aspect_ratio = 1.01;
        MocNozzle solver(o);
        MocResult stressed = solver.solve();

        ASSERT_TRUE(stressed.converged) << c.name << ": " << stressed.failure.message;
        EXPECT_EQ(stressed.inserted_characteristics, 0u) << c.name;
        EXPECT_EQ(stressed.retired_characteristics, 0u) << c.name;
        EXPECT_DOUBLE_EQ(stressed.exit_mach, baseline.exit_mach) << c.name;
        EXPECT_DOUBLE_EQ(stressed.area_ratio, baseline.area_ratio) << c.name;
        ASSERT_EQ(stressed.net.wall_x.size(), baseline.net.wall_x.size()) << c.name;
        for (size_t i = 0; i < baseline.net.wall_x.size(); i++) {
            EXPECT_DOUBLE_EQ(stressed.net.wall_x[i], baseline.net.wall_x[i]) << c.name << " i=" << i;
            EXPECT_DOUBLE_EQ(stressed.net.wall_y[i], baseline.net.wall_y[i]) << c.name << " i=" << i;
        }
    }
}
