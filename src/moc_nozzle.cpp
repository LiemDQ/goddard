#include <cmath>
#include <optional>
#include <format>
#include <stdexcept>
#include "goddard/equilibrium.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/moc_context.hpp"
#include "goddard/moc_initialization.hpp"
#include "goddard/moc_direct_march.hpp"
#include "goddard/moc_inverse_march.hpp"
#include "goddard/moc_thermo.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/error.hpp"

namespace Goddard {

namespace {

/**
 * The options a solve actually runs on: a copy of the caller's, in throat radii.
 *
 * The start lines are built at a throat of unit radius (the centered fan's lip at (0, 1), the
 * Kliegel-Levine series in r / r_t) and the Rao contour is generated at r_throat = 1, so the
 * analysis contour is divided by the caller's throat radius to match, and
 * geometry.throat_radius becomes 1. Everything downstream of solve() then sees one unit, and
 * solve() converts the result's lengths back (scale_result_lengths).
 */
MocOptions options_in_throat_radii(const MocOptions& options) {
    MocOptions normalized = options;
    normalized.geometry.throat_radius = 1.0;
    normalized.nozzle_profile = options.nozzle_profile.scaled(1.0 / options.geometry.throat_radius);
    return normalized;
}

/**
 * Convert every length in a result assembled in throat radii to the caller's unit. The
 * dimensionless fields (Mach numbers, angles, area ratio, coverage, the *_over_spacing and
 * grading/margin/error diagnostics) are left alone.
 */
void scale_result_lengths(MocResult& result, double factor) {
    for (CharacteristicPoint& pt : result.net.points) {
        pt.x *= factor;
        pt.y *= factor;
    }
    for (double& value : result.net.wall_x) value *= factor;
    for (double& value : result.net.wall_y) value *= factor;
    result.profile = result.profile.scaled(factor);
    for (double& value : result.exit_plane.y) value *= factor;

    result.nozzle_length *= factor;
    result.min_theta_x *= factor;
    result.min_theta_y *= factor;
    result.crossings.first_x *= factor;
    result.crossings.first_y *= factor;
    result.failure.x *= factor;
    result.failure.y *= factor;

    for (MocPassDiagnostics& diag : result.pass_diagnostics) {
        diag.min_spacing *= factor;
        diag.max_spacing *= factor;
        diag.mean_spacing *= factor;
        diag.front_axis_x *= factor;
        diag.front_wall_x *= factor;
        diag.front_axis_spacing *= factor;
        diag.front_wall_spacing *= factor;
        diag.step_dx *= factor;
    }
    result.init_diagnostics.wall_gap *= factor;
}

// Resolved wall contour for a solve, in throat radii: options.nozzle_profile as given
// (ANALYSIS, DESIGN_MIN_LENGTH; already normalized by options_in_throat_radii), or the Rao
// contour generated for a unit throat (DESIGN_RAO). Returned by value so the caller holds it
// as a solve() local -- never written back into options.nozzle_profile: doing so would leave
// a later solve on the same MocNozzle instance, with mode changed away from DESIGN_RAO,
// reading a stale Rao contour from the prior solve.
NozzleProfile resolve_wall_profile(const MocOptions& options) {
    switch (options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            return options.nozzle_profile;
        }
        case MocMode::DESIGN_RAO: {
            return NozzleProfile::generate_Rao_TOP_nozzle(
                options.geometry.expansion_ratio, 1.0, options.geometry.length_fraction);
        }
        case MocMode::ANALYSIS: {
            return options.nozzle_profile;
        }
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("DESIGN_CENTERLINE is not implemented.");
        }
        default: {
            //unreachable
            throw std::runtime_error("Invalid value for MocMode specified.");
        }
    }
}

// ONE copy of the seven push_backs that appear three times in result assembly below.
void append_exit_plane_point(ExitPlane& exit_plane, const CharacteristicPoint& pt) {
    exit_plane.y.push_back(pt.y);
    exit_plane.mach.push_back(pt.mach);
    exit_plane.theta.push_back(pt.theta);
    exit_plane.pressure.push_back(pt.pressure);
    exit_plane.temperature.push_back(pt.temperature);
    exit_plane.gamma_s.push_back(pt.gamma_s);
    exit_plane.velocity.push_back(pt.V);
}

/**
 * The exit plane of a solved net: for a min-length nozzle, the last characteristic (leading
 * axis point + leading wall point); for the inverse march, the last front it built (axis to
 * wall, as stored), whole even on a partial failure.
 */
ExitPlane exit_plane_of(const CharacteristicNet& net, MocMode mode) {
    ExitPlane exit_plane;
    if (net.points.empty()) return exit_plane;

    if (mode != MocMode::DESIGN_MIN_LENGTH) {
        // The exit plane is exactly the last front the kernel built (or, on a partial
        // failure, whatever front it last completed) -- axis to wall, as stored.
        if (!net.fronts.empty()) {
            for (size_t idx : net.fronts.back()) {
                append_exit_plane_point(exit_plane, net.points[idx]);
            }
        }
        return exit_plane;
    }

    // Add the last characteristic for min length nozzle. Guarded by emptiness: the kernel
    // fails fast, so it can abort before any axis/wall point has been recorded, where
    // leading_axis_point()/leading_wall_point() would otherwise index an empty vector.
    if (!net.axis_point_indices.empty()) {
        append_exit_plane_point(exit_plane, net.leading_axis_point());
    }
    if (!net.wall_point_indices.empty()) {
        append_exit_plane_point(exit_plane, net.leading_wall_point());
    }
    return exit_plane;
}

// Exit-to-throat area ratio, given the throat reference radius r_throat in net units
// (net.wall_y.front() for the ladder, the normalized geometry's 1 for the inverse march --
// see solve()). Callers guard net.wall_y/r_throat emptiness themselves before calling this.
double area_ratio_of(const CharacteristicNet& net, MocFlowKind flow, double r_throat) {
    const double y_ratio = net.wall_y.back() / r_throat;
    return (flow == MocFlowKind::PLANAR) ? y_ratio : y_ratio * y_ratio;
}

} // namespace

bool MocNozzle::is_solved() const {
    return m_is_solved;
}

// -- MocNozzle --
MocResult MocNozzle::solve() {
    // Re-checked here, not just in the constructor: options is public and callers do
    // adjust it between construction and solve().
    validate_moc_options(options);

    MocLog log(options.log_level);
    m_is_solved = false;

    // The solve runs in throat radii (see options_in_throat_radii); `throat_radius` is the
    // caller's unit, used only to convert the result back at the end.
    const double throat_radius = options.geometry.throat_radius;
    const MocOptions normalized = options_in_throat_radii(options);
    NozzleProfile wall = resolve_wall_profile(normalized);

    ThroatCondition throat{};
    std::optional<MocThermo> thermo;

    if (options.chemistry == GasChemistry::PERFECT_GAS) {
        // Pure algebraic path: no Cantera dependency
        // throat populated with dummy values; they are not used in the perfect gas case.
        throat = ThroatCondition{
            .speed_of_sound = 1.0,
            .H_stagnation = 1.0,
            .P_inlet = 1.0, // dimensionless stagnation pressure
            .S_inlet = 1.0,
            .gamma_s = options.gamma,
            .dlV_dlP_T = -1.0,
            .dlV_dlT_P = 1.0,
            .state = {}
        };
        thermo = MocThermo::perfect_gas(options.gamma);
    }
    else {
        // Cantera-backed path for frozen/equilibrium chemistry
        if (!m_gas) {
            throw std::invalid_argument(
                "MocNozzle: FROZEN and EQUILIBRIUM chemistry need a Gas; construct with MocNozzle(gas, options).");
        }
        NozzleOptions nozzle_opts;
        nozzle_opts.chemistry = options.chemistry;
        Nozzle nozzle(*m_gas, nozzle_opts);
        throat = nozzle.solve_throat_conditions();
        thermo = MocThermo::tabulated(*m_gas, options.chemistry, throat);
    }

    // Everything a kernel or unit process needs for this solve, besides the points it works
    // on.
    MocSolveContext ctx{normalized, wall, *thermo, log};

    MocResult result;
    result.throat_radius = throat_radius;

    StartLine line;
    double reference_spacing = 0.0;
    // Throat radius in the net's own unit (throat radii, so 1), for area_ratio and the
    // start-line diagnostics.
    double net_throat_radius = 1.0;
    CharacteristicNet net;

    // Exactly one of these is populated inside the try block below, depending on
    // options.mode; both outlive it so the kernel can be run (and, for the inverse march,
    // its pass diagnostics read) after init_diagnostics is measured.
    std::optional<DirectMarch> direct_march;
    std::optional<InverseMarch> inverse_march;

    // The initial data line (a Newton solve for the Kliegel-Levine transonic line, a
    // centered-fan Prandtl-Meyer inversion, or a PrandtlMeyerTable lookup) can fail
    // numerically before the kernel ever starts marching. Converting those failures
    // here, rather than letting them propagate as exceptions, is what makes solve()
    // never throw for a numerical/convergence failure.
    try {
        if (m_inverse_front_override.has_value()) {
            // Test-only path (see m_inverse_front_override): skip the normal throat/KL/fan
            // construction entirely -- it need not be consistent with a hand-built front --
            // and seed the inverse kernel directly from the override. Treated as a generic
            // (Kliegel-Levine-like) line for InverseMarch::initial_front and
            // measure_start_line: it is not collinear along a single characteristic.
            // Given in the caller's unit, like nozzle_profile.
            line.points = *m_inverse_front_override;
            for (CharacteristicPoint& pt : line.points) {
                pt.x *= 1.0 / throat_radius;
                pt.y *= 1.0 / throat_radius;
            }
            line.family.reset();
            line.used = MocStartLine::KLIEGEL_LEVINE;
        } else {
            line = build_start_line(ctx, throat);
        }

        if (options.mode == MocMode::DESIGN_MIN_LENGTH) {
            direct_march.emplace(ctx, line, net);
            direct_march->seed();

            // Establish net_throat_radius and reference_spacing, once, from the line actually
            // seeded. Both init paths put a wall point at the throat lip first, so
            // wall_y.front() is the throat radius in the units the net is carrying (1).
            if (!net.wall_y.empty() && net.wall_y.front() > 0.0 &&
                options.num_characteristics > 1) {
                net_throat_radius = net.wall_y.front();
                reference_spacing =
                    net_throat_radius / static_cast<double>(options.num_characteristics - 1);
            }
        }
        else {
            // The front-based kernel (InverseMarch): no throat-lip anchor is seeded (see
            // InverseMarch::initial_front / CharacteristicNet::add_front); wall_x.front() is
            // F_0's own wall point, not (0, 1), so the throat radius for area_ratio and
            // diagnostic purposes is the normalized geometry's, i.e. 1.
            inverse_march.emplace(ctx, net);
            inverse_march->seed(m_inverse_front_override.has_value()
                ? line.points : inverse_march->initial_front(line));
            net_throat_radius = normalized.geometry.throat_radius;
            // Without this, reference_spacing would stay 0 for every analysis/Rao solve,
            // which would make MocInitDiagnostics::wall_gap_over_spacing and
            // wall_station_to_tangency silently report spacing = 1; anchoring it the same way
            // the chain ladder does keeps both kernels' diagnostics comparable.
            if (options.num_characteristics > 1) {
                reference_spacing = net_throat_radius / static_cast<double>(options.num_characteristics - 1);
            }
        }
    }
    catch (const NotImplementedError&) {
        // Programmer error (unimplemented mode): not a numerical failure, keep throwing.
        throw;
    }
    catch (const std::out_of_range& e) {
        result.failure = MocFailure{
            MocErrorCode::TABLE_RANGE_EXCEEDED,
            std::string("Initial data line construction failed: ") + e.what(),
            0.0, 0.0, -1
        };
    }
    catch (const std::exception& e) {
        result.failure = MocFailure{
            MocErrorCode::INITIALIZATION_FAILED,
            std::string("Initial data line construction failed: ") + e.what(),
            0.0, 0.0, -1
        };
    }

    if (result.failure.code != MocErrorCode::NONE) {
        result.converged = false;
        result.messages = log.messages;
        m_is_solved = true;
        return result;
    }

    // Measured before the march so it is available even when the kernel fails partway --
    // an initialization defect is exactly the case where the march does not finish.
    result.init_diagnostics = measure_start_line(line, ctx, reference_spacing, net_throat_radius);

    std::optional<MocFailure> kernel_failure = direct_march.has_value()
        ? direct_march->run()
        : inverse_march->run();

    // Minimum flow angle and its location, for every scheme: see MocResult::min_theta.
    if (!net.points.empty()) {
        const CharacteristicPoint* lowest = &net.points.front();
        for (const CharacteristicPoint& pt : net.points) {
            if (pt.theta < lowest->theta) lowest = &pt;
        }
        result.min_theta = lowest->theta;
        result.min_theta_x = lowest->x;
        result.min_theta_y = lowest->y;
        if (result.min_theta < -0.5 * M_PI / 180.0) {
            log.warning("Flow angle reaches {:.2f} deg at (x={:.3f}, y={:.3f} throat radii): a "
                "compression is converging on the axis (a forming shock); the isentropic "
                "solution downstream of it is approximate.", result.min_theta * 180.0 / M_PI,
                result.min_theta_x, result.min_theta_y);
        }
    }

    result.net = net;
    result.messages = log.messages;
    if (kernel_failure.has_value()) {
        result.failure = *kernel_failure;
        result.converged = false;
    } else {
        result.converged = true;
    }

    // Populate performance fields
    result.crossings = find_like_characteristic_crossings(net);
    // Per-pass front geometry: recorded by the inverse march only.
    result.pass_diagnostics = inverse_march.has_value() ? inverse_march->pass_diagnostics
                                                         : std::vector<MocPassDiagnostics>{};

    // Fraction of the target exit radius the net reached (MocResult::exit_coverage). Full for
    // the minimum-length ladder when it completes, since that contour is defined to stop
    // exactly at theta_max; a ladder that fails partway leaves the defaults (0, false). For
    // the front-based kernel the last pass always lands exactly on the exit plane when the
    // march converges; on a solve that fails partway, this instead reports how far the net's
    // last wall point got, which is the diagnostic exit_coverage exists for.
    if (options.mode == MocMode::DESIGN_MIN_LENGTH) {
        if (!kernel_failure.has_value()) {
            result.exit_coverage = 1.0;
            result.reached_exit_plane = true;
        }
    }
    else {
        // The inverse march always lands its last front exactly on the exit plane when it
        // converges (see InverseMarch::run()'s MocStepLimiter::EXIT step). A march that
        // failed partway did not reach it; report how far it got instead.
        if (!kernel_failure.has_value()) {
            result.exit_coverage = 1.0;
            result.reached_exit_plane = true;
        }
        else if (!wall.y.empty() && wall.y.back() > 0.0 && !net.wall_y.empty()) {
            result.exit_coverage = net.wall_y.back() / wall.y.back();
            result.reached_exit_plane = false;
        }
    }
    if (!net.wall_x.empty()) {
        result.nozzle_length = net.wall_x.back();
    }
    if (options.mode != MocMode::DESIGN_MIN_LENGTH) {
        // net.wall_x.front()/wall_y.front() are F_0's own wall point, not the throat lip
        // (0, 1) DirectMarch seeds -- see the front-based seeding branch above -- so the
        // throat radius reference is the normalized geometry's (1) instead.
        if (!net.wall_y.empty()) {
            result.area_ratio = area_ratio_of(net, options.flow_type, net_throat_radius);
        }
    }
    else if (!net.wall_y.empty() && net.wall_y.front() > 0.0) {
        result.area_ratio = area_ratio_of(net, options.flow_type, net.wall_y.front());
    }

    // Exit Mach:
    // For a min-length nozzle the exit flow is uniform by construction, so the last computed
    // point (on the final characteristic) is representative.
    // For analysis/Rao the exit flow is generally not uniform, so no single computed point
    // represents "the" exit Mach; report the centerline value (most downstream axis point)
    // instead, for a number that is consistently defined across modes.
    if (!net.points.empty()) {
        if (options.mode != MocMode::DESIGN_MIN_LENGTH && !net.axis_point_indices.empty()) {
            result.exit_mach = net.leading_axis_point().mach;
        } else {
            result.exit_mach = net.points.back().mach;
        }
    }

    // Build wall profile from computed wall coordinates
    result.profile.x = net.wall_x;
    result.profile.y = net.wall_y;

    // Exit plane extraction: see exit_plane_of() above.
    result.exit_plane = exit_plane_of(net, options.mode);

    // Everything above is in throat radii; report lengths in the caller's unit.
    scale_result_lengths(result, throat_radius);

    m_is_solved = true;
    return result;
}

} // namespace Goddard
