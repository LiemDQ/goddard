#include <cmath>
#include <algorithm>
#include <numeric>
#include <optional>
#include <limits>
#include <format>
#include <iostream>
#include "goddard/equilibrium.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/moc_initialization.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/error.hpp"

namespace Goddard {

namespace {
// Unwraps a PointResult produced while constructing the initial data line. Unlike
// the kernel (which communicates a unit-process failure back to solve() via a
// returned MocFailure), the initial data line runs before the kernel and
// communicates failure through solve()'s exception boundary: throwing here is
// caught by MocNozzle::solve() and converted to MocFailure{INITIALIZATION_FAILED}.
CharacteristicPoint unwrap_or_throw(const PointResult& result, std::string_view context) {
    if (result.error != MocErrorCode::NONE) {
        throw ConvergenceError(std::format(
            "{}: unit process failed at ({}, {}) with error {}.",
            context, result.point.x, result.point.y, to_string(result.error)));
    }
    return result.point;
}
} // namespace

// -- Logging helpers --
void MocNozzle::log_warning(const std::string& msg) {
    m_messages.push_back("Warning: " + msg);
}

void MocNozzle::log_info(const std::string& msg) {
    m_messages.push_back("Info: " + msg);
}

void MocNozzle::log_debug(const std::string& msg) {
    if (m_options.log_level != MocLogLevel::DEBUG) return;
    std::string full = "Debug: " + msg;
    m_messages.push_back(full);
    std::cerr << full << std::endl;
}

bool MocNozzle::is_solved() const {
    return m_is_solved;
}

// -- MocNozzle --
MocResult MocNozzle::solve() {
    // Re-checked here, not just in the constructor: m_options is public and callers do
    // adjust it between construction and solve().
    validate_moc_options(m_options);

    m_messages.clear();
    m_theta_schedule.clear();
    m_initial_line_family.reset();
    m_is_solved = false;
    m_inserted_characteristics = 0;
    m_retired_characteristics = 0;
    m_pass_diagnostics.clear();
    m_reference_spacing = 0.0;
    m_throat_radius = 1.0;
    m_warned_non_spacelike = false;
    m_init_wall_bc_residual = 0.0;

    // Resolved once here so the rest of solve() and the kernel dispatch below work with a
    // concrete scheme; see MocMarchScheme::AUTO.
    const MocMarchScheme scheme = resolve_march_scheme();

    CharacteristicNet net;
    std::vector<CharacteristicPoint> data_line;
    m_options.nozzle_profile = setup_nozzle_profile(m_options.geometry);

    ThroatCondition throat{};

    if (m_options.chemistry == GasChemistry::PERFECT_GAS) {
        // Pure algebraic path: no Cantera dependency
        m_L_ref = m_options.geometry.throat_radius;
        m_P_ref = 1.0; // dimensionless stagnation pressure
        m_T_ref = 1.0; // dimensionless stagnation temperature

        // throat populated with dummy values; they are not used in the perfect gas case.
        throat = ThroatCondition{
            .converged = true,
            .speed_of_sound = 1.0,
            .H_stagnation = 1.0,
            .P_inlet = m_P_ref,
            .S_inlet = 1.0,
            .gamma_s = m_options.gamma,
            .dlV_dlP_T = -1.0,
            .dlV_dlT_P = 1.0,
            .state = {}
        };
    }
    else {
        // Cantera-backed path for frozen/equilibrium chemistry
        NozzleOptions nozzle_opts;
        nozzle_opts.chemistry = m_options.chemistry;
        Nozzle nozzle(*m_gas, nozzle_opts);
        throat = nozzle.solve_throat_conditions();
        m_gas->restore_state(throat.state);

        m_L_ref = m_options.geometry.throat_radius;
        m_P_ref = throat.P_inlet;
        m_T_ref = m_gas->thermo()->temperature();
        m_S_ref = throat.S_inlet;

        double a_throat = m_gas->speed_of_sound();

        pm_table.build_table(
            *m_gas->thermo(),
            m_options.chemistry == GasChemistry::EQUILIBRIUM,
            throat.S_inlet,
            throat.H_stagnation,
            a_throat);
    }

    MocResult result;

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
            // (Kliegel-Levine-like) line for build_inverse_initial_front and
            // record_init_diagnostics: it is not collinear along a single characteristic.
            data_line = *m_inverse_front_override;
            m_initial_line_family.reset();
        } else {
            data_line = generate_initial_data_line(throat, m_options.geometry, m_options.num_characteristics);
        }

        if (scheme == MocMarchScheme::DIRECT) {
            // The initial-data-line generators declared whether the line lies along a single
            // characteristic (m_initial_line_family). A centered expansion fan does, and needs a
            // wall anchor at the throat lip (0, 1): it seeds leading_wall_point() for the first wall
            // solve and gives the area ratio its throat reference. A transonic start line spans
            // axis-to-wall and already carries its own wall point, so no separate anchor is seeded.
            if (m_initial_line_family.has_value()) {
                CharacteristicPoint throat_lip{};
                throat_lip.x = 0.0;
                throat_lip.y = 1.0;
                // The lip angle anchors the first wall solve in DESIGN_MIN_LENGTH, where it equals
                // theta_max. In analysis the wall angles come from the profile and the lip is purely
                // an anchor, so the topmost fan angle is a harmless stand-in.
                throat_lip.theta = m_options.mode == MocMode::DESIGN_MIN_LENGTH
                    ? m_options.theta_max
                    : (data_line.empty() ? 0.0 : data_line.back().theta);
                net.seed_wall_point(throat_lip);
            }
            net.add_initial_data_line(data_line, m_initial_line_family);

            // Anchor mesh control to the geometry, once, from the line actually seeded. Both
            // init paths put a wall point at the throat lip first, so wall_y.front() is the
            // throat radius in whatever units the net is carrying; taking it from the net rather
            // than from geometry.throat_radius keeps the two consistent even when the initial
            // line is built in normalized coordinates.
            if (!net.wall_y.empty() && net.wall_y.front() > 0.0 &&
                m_options.num_characteristics > 1) {
                m_throat_radius = net.wall_y.front();
                m_reference_spacing =
                    m_throat_radius / static_cast<double>(m_options.num_characteristics - 1);
            }
        }
        else {
            // INVERSE: no throat-lip anchor is seeded (see build_inverse_initial_front /
            // seed_inverse_front); wall_x.front() is F_0's own wall point, not (0, 1), so the
            // throat radius for area_ratio/mesh purposes comes from the geometry directly.
            std::vector<CharacteristicPoint> front0 = build_inverse_initial_front(data_line);
            if (front0.size() < 2) {
                throw ConvergenceError(
                    "Inverse march: initial front has fewer than two points.");
            }
            seed_inverse_front(net, front0);
            m_throat_radius = (m_options.geometry.throat_radius > 0.0)
                ? m_options.geometry.throat_radius : 1.0;
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
        result.messages = m_messages;
        m_is_solved = true;
        return result;
    }

    // Measured before the march so it is available even when the kernel fails partway --
    // an initialization defect is exactly the case where the march does not finish.
    result.init_diagnostics = record_init_diagnostics(data_line);

    std::optional<MocFailure> kernel_failure = (scheme == MocMarchScheme::INVERSE)
        ? solve_inverse_characteristic_kernel(net)
        : solve_characteristic_kernel(net);

    result.net = net;
    result.messages = m_messages;
    if (kernel_failure.has_value()) {
        result.failure = *kernel_failure;
        result.converged = false;
    } else {
        result.converged = true;
    }

    // Populate performance fields
    result.crossings = find_like_characteristic_crossings(net);
    result.inserted_characteristics = m_inserted_characteristics;
    result.retired_characteristics = m_retired_characteristics;
    result.pass_diagnostics = m_pass_diagnostics;

    // How much of the requested exit radius the outflow boundary actually reached. The
    // boundary is a ragged staircase of independently terminated chains, so even a healthy
    // march ends up to about one characteristic spacing short; the allowance below is two
    // spacings. This is reported, not folded into `converged`: measured coverage does not
    // separate a coarse-but-healthy solve from a truncated one (planar N=8 and the
    // axisymmetric AR=4 N=8 truncation both sit at 0.92), so the meaningful test is whether
    // the shortfall shrinks under refinement, which only a grid sweep can see.
    if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
        result.exit_coverage = 1.0;
        result.reached_exit_plane = true;
    }
    else if (scheme == MocMarchScheme::INVERSE) {
        // The inverse march always lands its last front exactly on the exit plane (see
        // solve_inverse_characteristic_kernel's MocStepLimiter::EXIT step) -- there is no
        // ragged outflow staircase to fall short of, unlike DIRECT. A march that failed
        // partway did not reach it; report how far it got the same way DIRECT does.
        if (!kernel_failure.has_value()) {
            result.exit_coverage = 1.0;
            result.reached_exit_plane = true;
        }
        else if (!m_options.nozzle_profile.y.empty() && m_options.nozzle_profile.y.back() > 0.0
                 && !net.wall_y.empty()) {
            result.exit_coverage = net.wall_y.back() / m_options.nozzle_profile.y.back();
            result.reached_exit_plane = false;
        }
    }
    else if (!m_options.nozzle_profile.y.empty() && m_options.nozzle_profile.y.back() > 0.0
             && !net.wall_y.empty()) {
        result.exit_coverage = net.wall_y.back() / m_options.nozzle_profile.y.back();
        // Allow two characteristic spacings of staircase shortfall. The allowance is
        // O(1/num_characteristics), so it tightens as the grid refines -- which is the
        // discriminating property: a healthy solve's shortfall shrinks with N, a truncated
        // one's does not. A fixed percentage cannot separate the two (a coarse healthy
        // planar solve and a genuinely truncated axisymmetric one both sit near 0.92).
        const double staircase_allowance =
            2.0 * m_reference_spacing / std::max(m_throat_radius, 1e-30);
        result.reached_exit_plane = result.exit_coverage >= 1.0 - staircase_allowance;

        // A march whose chains all terminated without covering the nozzle has not solved the
        // problem that was posed, however clean each individual unit process was. Reporting
        // that as converged is what let a requested area ratio of 4 be silently answered
        // with 2.4.
        if (!result.reached_exit_plane && result.failure.code == MocErrorCode::NONE) {
            result.failure = MocFailure{
                MocErrorCode::INCOMPLETE_MARCH,
                std::format("Every characteristic terminated but the net reached only {:.1f}% "
                    "of the requested exit radius (wall y={:.6f} of {:.6f}).",
                    100.0 * result.exit_coverage, net.wall_y.back(),
                    m_options.nozzle_profile.y.back()),
                net.wall_x.empty() ? 0.0 : net.wall_x.back(), net.wall_y.back(),
                -1
            };
            result.converged = false;
        }
    }
    if (!net.wall_x.empty()) {
        result.nozzle_length = net.wall_x.back();
    }
    if (scheme == MocMarchScheme::INVERSE) {
        // net.wall_x.front()/wall_y.front() are F_0's own wall point, not the throat lip
        // (0, 1) DIRECT seeds -- see the INVERSE seeding branch above -- so the throat
        // radius reference comes from the geometry directly instead.
        const double r_throat = (m_options.geometry.throat_radius > 0.0)
            ? m_options.geometry.throat_radius : 1.0;
        if (!net.wall_y.empty()) {
            if (m_options.flow_type == MocFlowKind::PLANAR) {
                result.area_ratio = net.wall_y.back() / r_throat;
            } else {
                const double y_ratio = net.wall_y.back() / r_throat;
                result.area_ratio = y_ratio * y_ratio;
            }
        }
    }
    else if (!net.wall_y.empty() && net.wall_y.front() > 0.0) {
        if (m_options.flow_type == MocFlowKind::PLANAR) {
            result.area_ratio = net.wall_y.back() / net.wall_y.front();
        } else {
            double y_ratio = net.wall_y.back() / net.wall_y.front();
            result.area_ratio = y_ratio * y_ratio;
        }
    }

    // Exit Mach:
    // For a min-length nozzle the exit flow is uniform, so the last computed point
    // (on the final characteristic) is representative.
    // In analysis mode the outflow boundary is a ragged staircase of terminated
    // chains and the last computed point is arbitrary; report the centerline exit
    // Mach (most downstream axis point) instead.
    if (!net.points.empty()) {
        if (m_options.mode != MocMode::DESIGN_MIN_LENGTH && !net.axis_point_indices.empty()) {
            result.exit_mach = net.leading_axis_point().mach;
        } else {
            result.exit_mach = net.points.back().mach;
        }
    }

    // Build wall profile from computed wall coordinates
    result.profile.x = net.wall_x;
    result.profile.y = net.wall_y;



    // Exit plane extraction:
    // For a min-length nozzle, the exit plane is the last wavefront + last wall point.
    // The last wavefront contains the axis point, and each preceding wavefront's
    // last point was absorbed into wall calculations.
    if (!net.points.empty()) {
        if (scheme == MocMarchScheme::INVERSE) {
            // The exit plane is exactly the last front the kernel built (or, on a partial
            // failure, whatever front it last completed) -- axis to wall, as stored.
            if (!net.fronts.empty()) {
                for (size_t idx : net.fronts.back()) {
                    const CharacteristicPoint& pt = net.points[idx];
                    result.exit_plane.y.push_back(pt.y);
                    result.exit_plane.mach.push_back(pt.mach);
                    result.exit_plane.theta.push_back(pt.theta);
                    result.exit_plane.pressure.push_back(pt.pressure);
                    result.exit_plane.temperature.push_back(pt.temperature);
                    result.exit_plane.gamma_s.push_back(pt.gamma_s);
                    result.exit_plane.velocity.push_back(pt.V);
                }
            }
            m_is_solved = true;
            return result;
        }
        const auto& outflows = net.outflow_points();
        // Add the last characteristic for min length nozzle. Guarded by
        // emptiness: a kernel that now fails fast can abort before any axis/wall
        // point has been recorded, where leading_axis_point()/leading_wall_point()
        // would otherwise index an empty vector.
        if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
            if (!net.axis_point_indices.empty()) {
                const auto& last_axis_pt = net.leading_axis_point();
                result.exit_plane.y.push_back(last_axis_pt.y);
                result.exit_plane.mach.push_back(last_axis_pt.mach);
                result.exit_plane.theta.push_back(last_axis_pt.theta);
                result.exit_plane.pressure.push_back(last_axis_pt.pressure);
                result.exit_plane.temperature.push_back(last_axis_pt.temperature);
                result.exit_plane.gamma_s.push_back(last_axis_pt.gamma_s);
                result.exit_plane.velocity.push_back(last_axis_pt.V);
            }

            if (!net.wall_point_indices.empty()) {
                const auto& last_wall_pt = net.leading_wall_point();
                result.exit_plane.y.push_back(last_wall_pt.y);
                result.exit_plane.mach.push_back(last_wall_pt.mach);
                result.exit_plane.theta.push_back(last_wall_pt.theta);
                result.exit_plane.pressure.push_back(last_wall_pt.pressure);
                result.exit_plane.temperature.push_back(last_wall_pt.temperature);
                result.exit_plane.gamma_s.push_back(last_wall_pt.gamma_s);
                result.exit_plane.velocity.push_back(last_wall_pt.V);
            }
        }
        else {
            for (const auto& pt : outflows) {
                result.exit_plane.y.push_back(pt.y);
                result.exit_plane.mach.push_back(pt.mach);
                result.exit_plane.theta.push_back(pt.theta);
                result.exit_plane.pressure.push_back(pt.pressure);
                result.exit_plane.temperature.push_back(pt.temperature);
                result.exit_plane.gamma_s.push_back(pt.gamma_s);
                result.exit_plane.velocity.push_back(pt.V);
            }
        }
    }
    m_is_solved = true;

    return result;
}


auto MocNozzle::setup_nozzle_profile(const NozzleGeometry& geometry) -> NozzleProfile {
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            return m_options.nozzle_profile;
        }
        case MocMode::DESIGN_RAO: {
            return NozzleProfile::generate_Rao_TOP_nozzle(geometry.expansion_ratio, 1.0, geometry.length_fraction);
        }
        case MocMode::ANALYSIS: {
            return m_options.nozzle_profile;
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


std::vector<CharacteristicPoint> MocNozzle::generate_initial_data_line(
    const ThroatCondition& throat,
    const NozzleGeometry& geometry,
    size_t num_points)
{
    std::vector<CharacteristicPoint> data_line;
    data_line.reserve(num_points);
    if (!m_options.theta_schedule.empty()) {
        m_theta_schedule = m_options.theta_schedule;
    }
    
    if (m_gas.has_value()) {
        m_gas->thermo()->restoreState(throat.state);
    }

    ThermodynamicContext context = build_thermo_context();

    // Resolve which start line this solve actually uses, once, as a single enum-valued
    // decision that both the fan/KL choice below and the theta_max derivation key off of
    // (previously a bool computed the choice while a separate ad hoc condition guarded the
    // theta_max derivation, and the two could disagree once a fallback was possible).
    // MocStartLine::AUTO keeps today's rule: a positive downstream curvature radius and
    // axisymmetric flow select the Kliegel-Levine series; anything else selects the fan.
    // The series is also not valid for planar flow when forced (validate_moc_options
    // rejects that combination outright).
    MocStartLine resolved_start_line = m_options.start_line;
    if (resolved_start_line == MocStartLine::AUTO) {
        resolved_start_line =
            (m_options.geometry.downstream_wall_curvature_radius <= 0.0 ||
             m_options.flow_type == MocFlowKind::PLANAR)
                ? MocStartLine::CENTERED_FAN
                : MocStartLine::KLIEGEL_LEVINE;
        if (resolved_start_line == MocStartLine::CENTERED_FAN &&
            m_options.flow_type == MocFlowKind::PLANAR &&
            m_options.geometry.downstream_wall_curvature_radius > 0.0 &&
            m_options.mode != MocMode::DESIGN_MIN_LENGTH)
        {
            log_info("Kliegel-Levine initialization is only valid for axisymmetric flow; "
                     "using centered-fan initialization for planar flow.");
        }
    }

    // In analysis mode the expansion is set by the wall contour, not by a design
    // theta_max (which is meaningless there and typically left unset). When the
    // centered-fan initializer will be used, derive theta_max from the profile. Shared
    // between the eager fan path and the AUTO fallback below (see the KL catch clause) so
    // the fallback gets exactly the same theta_max the eager path would have used.
    //
    // DESIGN_RAO is included alongside ANALYSIS: before this package DESIGN_RAO's default
    // positive curvature radius meant it never took the fan branch at all (use_centered_fan
    // was false whenever curvature was positive and axisymmetric), so this branch's
    // ANALYSIS-only guard was never exercised for it. The AUTO wall-angle-miss fallback
    // changes that -- DESIGN_RAO's default r_arc = 0.382 throat misses the threshold just
    // as an ANALYSIS contour at the same throat does -- and DESIGN_RAO's profile is equally
    // available here: solve() calls setup_nozzle_profile() (which generates the Rao
    // contour) before generate_initial_data_line() runs, for every mode.
    MocOptions init_options = m_options;
    auto derive_fan_theta_max = [&]() {
        if (m_options.mode != MocMode::ANALYSIS && m_options.mode != MocMode::DESIGN_RAO) return;
        const NozzleProfile& wall = m_options.nozzle_profile;
        if (wall.size() < 2) {
            throw std::runtime_error(
                "Centered-fan initialization requires a wall profile with at least 2 points.");
        }
        init_options.theta_max = wall.max_theta();
        if (init_options.theta_max <= 0.0) {
            throw std::runtime_error(
                "Wall profile must have a positive expansion angle for centered-fan initialization.");
        }
    };
    if (resolved_start_line == MocStartLine::CENTERED_FAN) {
        derive_fan_theta_max();
    }
    // MocInitialization holds a ThermodynamicContext by value, which carries a reference
    // member (its PrandtlMeyerTable) -- that makes the class copy-constructible but not
    // copy-assignable, so a rebuild after the fallback below uses emplace() (in-place
    // construction) rather than assignment.
    std::optional<MocInitialization> initializer;
    initializer.emplace(geometry, context, init_options);

    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // The minimum length nozzle is a special case as the sonic line is straight, and
            // a centered expansion fan is the "exact" solution.
            // Therefore all initialization methods simplify to straight line initialization.

            // A centered expansion fan is collinear along a single C+ characteristic.
            m_initial_line_family = ChainMetadata::Family::PLUS;

            auto expansion_line = initializer->initialize_centered_expansion(throat);
            // Mirror the fan angles (user-supplied or auto-generated inside the
            // initializer) into m_theta_schedule: the min-length wall solve reads
            // theta_wall = theta_max - m_theta_schedule[k] for each wall point.
            m_theta_schedule.resize(expansion_line.size());
            for (size_t i = 0; i < expansion_line.size(); i++) {
                m_theta_schedule[i] = expansion_line[i].theta;
            }
            CharacteristicPoint upstream_point;
            //for a minimum length nozzle, the sonic line is straight.
            //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
            for (size_t i = 0; i < expansion_line.size(); i++){
                const auto& expansion_point = expansion_line[i];
                // special case: the first point is on the centerline
                if (i == 0) {
                    upstream_point = unwrap_or_throw(
                        solve_initial_axis_point_centered_exp(expansion_point),
                        "Initial data line (min-length axis point)");
                }
                else {
                    upstream_point = unwrap_or_throw(
                        solve_interior_point(expansion_point, upstream_point),
                        "Initial data line (min-length interior point)");
                }
                data_line.push_back(upstream_point);
            }
            break;
        }
        case MocMode::DESIGN_RAO:
        case MocMode::ANALYSIS: {
            // The Kliegel-Levine transonic start line crosses many characteristics, so it is
            // seeded as a generic data line (m_initial_line_family stays empty).
             // the first theta should be small to minimize approximation error.
            CharacteristicPoint upstream_point;
            if (resolved_start_line == MocStartLine::KLIEGEL_LEVINE) {
                try {
                    // The Kliegel-Levine transonic line already spans axis-to-wall with full
                    // thermodynamic state and K+/K- set: it IS the initial data line. Unlike
                    // the centered fan (whose rays all emanate from the throat lip and must
                    // be marched into the flow field), re-marching these points pairwise
                    // would intersect characteristics *behind* their parents. See
                    // CharacteristicNet::add_initial_data_line's nullopt branch for how this
                    // non-collinear line is actually seeded into the net without that
                    // re-marching.
                    data_line = initializer->initialize_kliegel_levine(throat);
                    m_init_wall_bc_residual = initializer->last_wall_bc_residual;
                }
                catch (const KlWallAngleFallback& e) {
                    // MocOptions::start_line was AUTO (a forced KLIEGEL_LEVINE throws
                    // ConvergenceError instead, which is not this type and is left to
                    // propagate to solve()'s exception boundary as INITIALIZATION_FAILED).
                    // Rebuild as the centered fan below, with theta_max derived exactly as
                    // the eager fan path derives it.
                    log_info("Kliegel-Levine series misses the wall angle by {:.6f} rad at "
                             "R = {:.6f}; using centered-fan initialization.", e.miss, e.R);
                    resolved_start_line = MocStartLine::CENTERED_FAN;
                    derive_fan_theta_max();
                    initializer.emplace(geometry, context, init_options);
                }
            }
            if (resolved_start_line == MocStartLine::CENTERED_FAN) {
                m_initial_line_family = ChainMetadata::Family::PLUS;
                auto expansion_line = initializer->initialize_centered_expansion(throat);
                m_theta_schedule.resize(expansion_line.size());
                for (size_t i = 0; i < expansion_line.size(); i++) {
                    m_theta_schedule[i] = expansion_line[i].theta;
                }
                //for a minimum length nozzle, the sonic line is straight.
                //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
                for (size_t i = 0; i < expansion_line.size(); i++){
                    const auto& expansion_point = expansion_line[i];
                    // special case: the first point is on the centerline
                    if (i == 0) {
                        upstream_point = unwrap_or_throw(
                            solve_initial_axis_point_centered_exp(expansion_point),
                            "Initial data line (centered-fan axis point)");
                        data_line.push_back(upstream_point);
                    }
                    else {
                        upstream_point = unwrap_or_throw(
                            solve_interior_point(expansion_point, upstream_point),
                            "Initial data line (centered-fan interior point)");
                        data_line.push_back(upstream_point);
                    }
                }
            }
            break;
        }
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("Centerline nozzle initialization not implemented.");
        }
    }

    for (size_t i = 0; i < data_line.size(); i++) {
        const auto& pt = data_line[i];
        log_debug("INIT[{}] x={:.6f} y={:.6f} theta={:.6f} nu={:.6f} mach={:.6f} mu={:.6f}",
            i, pt.x, pt.y, pt.theta, pt.nu, pt.mach, pt.mu);
        // The Kliegel-Levine and Sauer transonic-line initializers compute Mach
        // directly (not via Prandtl-Meyer inversion) and are not funneled through a
        // unit process that would otherwise validate them; check them here so a
        // breakdown of the small-perturbation series (e.g. a subsonic or non-finite
        // result near the wall) is reported as MocFailure::INITIALIZATION_FAILED
        // instead of an invalid point silently entering the net.
        MocErrorCode validity = check_point_validity(pt, m_options.solver_options.abstol);
        if (validity != MocErrorCode::NONE) {
            throw ConvergenceError(std::format(
                "Initial data line point {} at (x={}, y={}) failed validity check: {}.",
                i, pt.x, pt.y, to_string(validity)));
        }
    }
    return data_line;
}

std::optional<MocFailure> MocNozzle::solve_characteristic_kernel(CharacteristicNet& net) {
    using Family = ChainMetadata::Family;
    LeadingEdgeView plus_edges = leading_edges(net, Family::PLUS);
    sort_plus_edges_by_proximity(plus_edges, net);
    LeadingEdgeView minus_edges = leading_edges(net, Family::MINUS);

    std::vector<std::pair<CharacteristicPoint, PointMembership>> intersections;
    std::vector<size_t> paired_cminus;
    std::vector<bool> cminus_is_intersected(net.c_chains.size(), false);

    intersections.reserve(net.c_chains.size());
    paired_cminus.reserve(net.c_chains.size());


    // A correctly converging 2D MoC kernel needs O(N) marching passes. This cap is a
    // safety bound so a non-converging/runaway net terminates promptly instead of
    // spinning indefinitely; reaching it indicates a kernel that has not converged.
    const int maxiter = 2000;
    int iters = 0;

    // Populated the instant any unit process reports a numerical failure; the march
    // aborts immediately rather than continuing to build on top of an invalid point
    // (which is what let corrupted points silently propagate through the net before).
    std::optional<MocFailure> failure;

    while (net.has_active_chains() && iters < maxiter && !failure.has_value()) {
        // Control spacing before pairing, never during: the pairing loop iterates the
        // leading-edge views below, and inserting or retiring rungs rebuilds them.
        MocPassDiagnostics diag;
        diag.pass = iters;
        control_front_spacing(net, plus_edges, minus_edges, diag);
        m_pass_diagnostics.push_back(diag);
        log_debug("--- kernel pass {}: {} active C+, {} active C- (front {} pts, "
            "spacing {:.6f}/{:.6f}/{:.6f} min/mean/max, target {:.6f}, aspect {:.3f}, "
            "margin {:.3f}, +{} -{}) ---",
            iters, plus_edges.chain_indices.size(), minus_edges.chain_indices.size(),
            diag.front_points, diag.min_spacing, diag.mean_spacing, diag.max_spacing,
            diag.target_spacing, diag.max_cell_aspect, diag.min_spacelike_margin,
            diag.inserted, diag.retired);
        intersections.clear();
        paired_cminus.clear();
        // Sized to the *current* minus-edge view: reflections add chains over time, so a
        // one-shot allocation sized to the initial chain count would be indexed out of range.
        cminus_is_intersected.assign(minus_edges.chain_indices.size(), false);

        /* We use a simple geometric approach to determine which points are intersecting.
        This requires no assumptions about net topology or precomputed traversal maps and is
        fairly efficient except for very large N (100,000+) which are unrealistic for 2D MoC methods.

        1. For each C+ leading point we find the corresponding C- leading point with the
        smallest y (height) that is larger than the C+ y and pair them. These points will intersect
        on the next pass.
        2. After all C+'s are paired, the bottommost unpaired C- pairs with the axis if it exists
        3. Wall interaction the topmost C+: either reflection or absorption depending on circumstances.
        4. Track active-chain count, terminate when it reaches 0 or all leading edges are past the outflow

        Invariant: each C- index should appear at most only once per set of intersections
        */
        for (size_t i = 0; i < plus_edges.chain_indices.size() && !failure.has_value(); i++) {
            PairSearch match = find_pair_partner(
                plus_edges.y_values[i], minus_edges, cminus_is_intersected, net);

            // intersect C+ with closest C- above it.
            if (match.chain_idx.has_value()) [[likely]] {
                const CharacteristicPoint& minus_pt = net.points[match.pt_idx];
                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];

                PointResult result = solve_interior_point(minus_pt, plus_pt);
                if (result.error != MocErrorCode::NONE) {
                    failure = MocFailure{
                        result.error,
                        std::format("Interior point solve failed ({}) pairing "
                            "minus(x={:.6f},y={:.6f}) with plus(x={:.6f},y={:.6f}).",
                            to_string(result.error), minus_pt.x, minus_pt.y, plus_pt.x, plus_pt.y),
                        result.point.x, result.point.y,
                        iters
                    };
                    break;
                }
                log_debug("PAIR minus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) "
                    "plus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) -> (x={:.6f},y={:.6f},mach={:.6f})",
                    minus_pt.x, minus_pt.y, minus_pt.theta, minus_pt.mu,
                    plus_pt.x, plus_pt.y, plus_pt.theta, plus_pt.mu,
                    result.point.x, result.point.y, result.point.mach);
                intersections.push_back({
                    result.point,
                    PointMembership {
                        .c_plus_chain_idx = plus_edges.chain_indices[i],
                        .c_minus_chain_idx = *match.chain_idx
                    }
                });
                paired_cminus.push_back(*match.chain_idx);
                cminus_is_intersected[match.edgevec_idx] = true;
            }
            else if (match.any_cminus_above) {
                // A C- does exist above this C+, but a closer competitor already claimed it
                // this pass (e.g. many individual C+ chains from a Kliegel-Levine transonic
                // line, competing for a single C- freshly born from a wall reflection). Leave
                // this chain active and retry once that C- has advanced on the next pass,
                // rather than wrongly treating it as a wall hit.
                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                log_debug("SKIP plus(x={:.6f},y={:.6f}) waiting for scarce C- partner "
                    "(already claimed this pass)", plus_pt.x, plus_pt.y);
            }
            else [[unlikely]] { // C+ intersects with wall

                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                std::optional<PointResult> maybe_result = solve_wall_point(
                    plus_pt,
                    net.leading_wall_point(),
                    static_cast<int>(net.wall_point_indices.size()) - 1
                );

                if (maybe_result.has_value()) {
                    if (maybe_result->error != MocErrorCode::NONE) {
                        failure = MocFailure{
                            maybe_result->error,
                            std::format("Wall point solve failed ({}) for plus(x={:.6f},y={:.6f}).",
                                to_string(maybe_result->error), plus_pt.x, plus_pt.y),
                            maybe_result->point.x, maybe_result->point.y,
                            iters
                        };
                        break;
                    }
                    log_debug("WALL plus(x={:.6f},y={:.6f}) -> wall hit at (x={:.6f},y={:.6f})",
                        plus_pt.x, plus_pt.y, maybe_result->point.x, maybe_result->point.y);
                    intersections.push_back({
                        maybe_result->point,
                        PointMembership {
                            .c_plus_chain_idx = plus_edges.chain_indices[i],
                            .c_minus_chain_idx = std::nullopt
                        }
                    });
                }
                else {
                    // if no intersection found, this point intersects beyond the profile boundary
                    log_debug("OUTFLOW plus(x={:.6f},y={:.6f}) terminates: "
                        "no wall hit found within profile bounds", plus_pt.x, plus_pt.y);
                    net.terminate_chain(
                        plus_edges.chain_indices[i],
                        ChainMetadata::TerminationType::OUTFLOW
                    );
                }
            }
        }

        if (failure.has_value()) break;

        // find C- characteristic reflecting off axis.
        double min_cminus_y = std::numeric_limits<double>::max();
        std::optional<size_t> min_y_cminus_idx = std::nullopt;
        for (size_t i = 0; i < minus_edges.chain_indices.size(); i++) {
            if (minus_edges.y_values[i] < min_cminus_y) {
                min_cminus_y = minus_edges.y_values[i];
                min_y_cminus_idx = minus_edges.chain_indices[i];
            }
        }
        // if the lowest C- characteristic is unpaired, reflect it off the axis
        if (min_y_cminus_idx.has_value() &&
            std::none_of(
                paired_cminus.cbegin(),
                paired_cminus.cend(),
                [min_y_cminus_idx](size_t x) {return x == *min_y_cminus_idx;})
            )
        {
            const CharacteristicPoint& minus_pt = net.leading_point(*min_y_cminus_idx);
            PointResult axis_result = solve_axis_point(minus_pt);
            if (axis_result.error != MocErrorCode::NONE) {
                failure = MocFailure{
                    axis_result.error,
                    std::format("Axis point solve failed ({}) for minus(x={:.6f},y={:.6f}).",
                        to_string(axis_result.error), minus_pt.x, minus_pt.y),
                    axis_result.point.x, axis_result.point.y,
                    iters
                };
                break;
            }
            log_debug("AXIS minus(x={:.6f},y={:.6f}) reflects off axis", minus_pt.x, minus_pt.y);
            intersections.push_back({
                axis_result.point,
                PointMembership {
                    .c_plus_chain_idx = std::nullopt,
                    .c_minus_chain_idx = min_y_cminus_idx
                }
            });
        }
        // insert all intersections into net
        for (auto&& [pt, mem] : intersections) {
            if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) { //min length nozzle is a special case
                if (!mem.c_minus_chain_idx.has_value()) {
                    // C+ reaches the wall: absorb it (no reflected wave). Cancelling the
                    // reflection here is exactly what produces a minimum-length contour.
                    net.terminate_c_plus_at_wall(*mem.c_plus_chain_idx, pt);
                    net.wall_x.push_back(pt.x);
                    net.wall_y.push_back(pt.y);
                }
                else if (!mem.c_plus_chain_idx.has_value()) {
                    // C- reaches the axis: reflect it into a new upward-marching C+.
                    net.reflect_c_minus_off_axis(*mem.c_minus_chain_idx, pt);
                }
                else {
                    net.add_point(pt, mem);
                }
            }
            else if (maximum_nozzle_length() >= pt.x) { //standard path
                if (!mem.c_minus_chain_idx.has_value()) {
                    // wall intersection
                    net.reflect_c_plus_off_wall(*mem.c_plus_chain_idx, pt);
                    net.wall_x.push_back(pt.x);
                    net.wall_y.push_back(pt.y);
                }
                else if (!mem.c_plus_chain_idx.has_value()) {
                    // axis intersection
                    net.reflect_c_minus_off_axis(*mem.c_minus_chain_idx, pt);
                } else {
                    net.add_point(pt, mem);
                }
            }
            else { // intersection falls outside of nozzle bounds
                log_debug("OUTFLOW pt(x={:.6f},y={:.6f}) terminates: beyond "
                    "maximum_nozzle_length() ({:.6f})", pt.x, pt.y, maximum_nozzle_length());
                if (mem.c_minus_chain_idx.has_value())
                    net.terminate_chain(*mem.c_minus_chain_idx, ChainMetadata::TerminationType::OUTFLOW);
                if (mem.c_plus_chain_idx.has_value())
                    net.terminate_chain(*mem.c_plus_chain_idx, ChainMetadata::TerminationType::OUTFLOW);
            }
        }
        update_leading_edges(plus_edges, net, Family::PLUS);
        sort_plus_edges_by_proximity(plus_edges, net);
        update_leading_edges(minus_edges, net, Family::MINUS);

        iters++;
    }

    if (failure.has_value()) {
        return failure;
    }

    if (iters >= maxiter) {
        return MocFailure{
            MocErrorCode::MAX_ITERATIONS_REACHED,
            std::format("Kernel reached the iteration safety cap ({}) without all "
                "characteristics terminating.", maxiter),
            0.0, 0.0, iters
        };
    }

    return std::nullopt;
}

std::optional<PointResult> MocNozzle::solve_wall_point(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    int wall_point_index)
{
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // for minimum length design: theta is determined by theta schedule.
            // Guard against an empty schedule (the Cantera min-length path does not
            // populate m_theta_schedule) and against indices beyond it (the kernel
            // currently produces more wall hits than scheduled characteristics).
            double theta_wall = m_options.theta_max;
            if (!m_theta_schedule.empty()) {
                size_t k = std::min(static_cast<size_t>(wall_point_index),
                                    m_theta_schedule.size() - 1);
                theta_wall = m_options.theta_max - m_theta_schedule[k];
            }

            return solve_wall_point_design(interior_parent, previous_wall_point, theta_wall);
        }
        case MocMode::DESIGN_RAO:
        case MocMode::ANALYSIS: {

            return solve_wall_point_analysis(interior_parent);
        }
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("DESIGN_CENTERLINE mode is not implemented.");
        }
        default: {
            throw std::runtime_error("Invalid MocMode in solve_wall_point");
        }
    }
    // unreachable
}

PointResult MocNozzle::solve_initial_axis_point_centered_exp(const CharacteristicPoint& expansion_point) {
    CharacteristicPoint point{};
    point.y = 0.0; //point always lies on axis
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {
            point.theta = expansion_point.theta;
            // nu = expansion point nu follows from geometric analysis
            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                point, expansion_point.nu, expansion_point.mach);
            if (thermo_error != MocErrorCode::NONE) {
                point.x = expansion_point.x;
                return {point, thermo_error};
            }
            point.K_plus = point.theta - point.nu;
            point.K_minus = expansion_point.K_minus;

            double c_minus_angle = average_cminus_angle(expansion_point, point);
            point.x = expansion_point.x - expansion_point.y / tan(c_minus_angle);
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {
            // Same approach as planar for the initial axis point:
            // assign theta from the expansion fan to bootstrap marching.
            // The axisymmetric source term is not applied here because the
            // initial expansion fan is modeled as a centered expansion at a point,
            // and the first axis point is a direct consequence of this geometric
            // construction. Source terms enter in subsequent solve_axis_point calls.
            point.theta = expansion_point.theta;
            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                point, expansion_point.nu, expansion_point.mach);
            if (thermo_error != MocErrorCode::NONE) {
                point.x = expansion_point.x;
                return {point, thermo_error};
            }
            point.K_plus = point.theta - point.nu;
            point.K_minus = expansion_point.K_minus;

            double c_minus_angle = average_cminus_angle(expansion_point, point);
            point.x = expansion_point.x - expansion_point.y / tan(c_minus_angle);
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    MocErrorCode validity = check_point_validity(point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {point, validity};
    }
    return {point, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_axis_point(const CharacteristicPoint& off_axis_parent) {
    CharacteristicPoint axis_point{};
    axis_point.y = 0.0;
    axis_point.theta = 0.0; // symmetry condition

    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {

            axis_point.K_minus = off_axis_parent.K_minus;
            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                axis_point,
                axis_point.K_minus - axis_point.theta,
                off_axis_parent.mach);
            if (thermo_error != MocErrorCode::NONE) {
                axis_point.x = off_axis_parent.x;
                return {axis_point, thermo_error};
            }

            axis_point.K_plus = axis_point.theta - axis_point.nu;

            double c_minus_angle = 0.5 * (off_axis_parent.theta + axis_point.theta)
                - 0.5 * (off_axis_parent.mu + axis_point.mu);

            axis_point.x = off_axis_parent.x - off_axis_parent.y/tan(c_minus_angle);
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {

            const double theta = off_axis_parent.theta;
            const double mu = off_axis_parent.mu;
            // Predictor method recommended by Zucrow & Hoffman, ch. 17
            // Use parent point as initial guess for source term.
            // The C- compatibility per dy is d(theta+nu) = [sin(theta)/(M sin(theta-mu))] (dy/y).
            // The descent from the parent to the axis has dy = -y_parent (signed), which makes
            // the source contribution positive (sin(theta-mu) < 0).
            axis_point.K_minus = off_axis_parent.K_minus;

            const double dy = -off_axis_parent.y; // y_axis - y_parent
            double source_pred = sin(theta) / (off_axis_parent.mach * sin(theta - mu))
                * dy / off_axis_parent.y;

            // dtheta + dnu = S => nu_2 - nu_1 = S + theta_1 = S + K_minus_1
            double nu_pred = axis_point.K_minus + source_pred;

            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                axis_point,
                nu_pred,
                off_axis_parent.mach);
            if (thermo_error != MocErrorCode::NONE) {
                axis_point.x = off_axis_parent.x;
                return {axis_point, thermo_error};
            }

            axis_point.K_plus = axis_point.theta - axis_point.nu;

            double c_minus_angle = average_cminus_angle(off_axis_parent, axis_point);
            axis_point.x = off_axis_parent.x - off_axis_parent.y / tan(c_minus_angle);

            // Corrector: apply the averaged source term with the same signed dy:
            // dy/y_avg = (0 - y_parent)/(y_parent/2) = -2, and theta_avg = theta_parent/2,
            // which together reproduce the finite sin(theta)/y limit at the axis.
            double dy_over_y_avg = -2.0;
            double theta_avg = 0.5 * off_axis_parent.theta;
            double mu_avg = 0.5 * (off_axis_parent.mu + axis_point.mu);
            double M_avg = 0.5 * (off_axis_parent.mach + axis_point.mach);
            double source_corrected = sin(theta_avg)/ (M_avg * sin(theta_avg - mu_avg)) * dy_over_y_avg;

            // Corrected nu: K_minus from parent, + source contribution
            double nu_corrected = off_axis_parent.K_minus + source_corrected;
            thermo_error = update_thermodynamic_state_from_nu(axis_point, nu_corrected, axis_point.mach);
            if (thermo_error != MocErrorCode::NONE) {
                return {axis_point, thermo_error};
            }
            axis_point.K_minus = axis_point.nu; // theta=0
            axis_point.K_plus = -axis_point.nu;

            // Recompute position with corrected slope
            c_minus_angle = average_cminus_angle(off_axis_parent, axis_point);
            axis_point.x = off_axis_parent.x - off_axis_parent.y / tan(c_minus_angle);
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    // solve_axis_point has no non-downstream guard today (a documented gap): the
    // resulting axis point's x must not precede its parent's.
    if (axis_point.x < off_axis_parent.x) {
        log_warning("Non-downstream axis intersection at ({}, {}). Parent at x={}.",
            axis_point.x, axis_point.y, off_axis_parent.x);
        return {axis_point, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(axis_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid axis point at ({}, {}): error code {}.",
            axis_point.x, axis_point.y, to_string(validity));
        return {axis_point, validity};
    }

    return {axis_point, MocErrorCode::NONE};
}

std::pair<double, double> MocNozzle::intersect_characteristic_with_wall(
        const CharacteristicPoint& parent,
        const NozzleProfile& wall) {
    
    double char_angle = parent.theta + parent.mu;
    
    // We use a basic predictor-corrector scheme to find the wall intersection
    // Predictor: straight-line intersection
    // Linear scan is appropriate for tens to low hundreds of points. 
    // if npoints grows beyond that, switch to bisection algorithm
    return find_wall_hit(parent, wall, char_angle);
}

PointResult MocNozzle::solve_interior_point(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {
            return solve_interior_point_planar(p1, p2);
        }
        case MocFlowKind::AXISYMMETRIC: {
            return solve_interior_point_axisymmetric(p1, p2);
        }
        default:
            throw std::runtime_error("Invalid MoC flow type specified.");
    }
}

PointResult MocNozzle::solve_interior_point_planar(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3{};
    p3.theta = 0.5*(p1.K_minus + p2.K_plus);
    double mach_guess = 0.5*(p1.mach + p2.mach);

    // Position (x, y) is computed after the thermo update below, so a thermo
    // failure here leaves p3.x/p3.y unset; fall back to a parent's coordinates so
    // the returned point still carries a meaningful location for diagnostics.
    MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
        p3, 0.5*(p1.K_minus - p2.K_plus), mach_guess);
    if (thermo_error != MocErrorCode::NONE) {
        p3.x = p1.x;
        p3.y = p1.y;
        return {p3, thermo_error};
    }
    p3.K_minus = p3.theta + p3.nu;
    p3.K_plus = p3.theta - p3.nu;

    // compute intersection point
    // assumimg characteristics are straight lines.
    double c_minus_angle = average_cminus_angle(p1, p3);
    double c_plus_angle = average_cplus_angle(p2, p3);

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        log_warning("Non-downstream intersection at ({}, {}). Parents at x=({}, {}).",
            p3.x, p3.y, p1.x, p2.x);
        return {p3, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(p3, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid interior point at ({}, {}): error code {}.",
            p3.x, p3.y, to_string(validity));
        return {p3, validity};
    }

    return {p3, MocErrorCode::NONE};
}
PointResult MocNozzle::solve_interior_point_axisymmetric(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3{};

    // Axisymmetric interior point in the (theta, nu) formulation, consistent with the planar,
    // wall, and axis solvers. The Riemann invariants K+ = theta - nu (constant along C+) and
    // K- = theta + nu (constant along C-) pick up source terms in axisymmetric flow:
    //   along C+:  d(theta - nu) = -L dx,   L = sin(mu) sin(theta) / (y cos(theta + mu))
    //   along C-:  d(theta + nu) = +M dx,   M = sin(mu) sin(theta) / (y cos(theta - mu))
    // With L = M = 0 this reduces exactly to the planar solver. Working in nu (rather than the
    // velocity form) keeps it exact for perfect gas, where nu = nu(M); the velocity form would
    // require the *true* velocity for cot(mu) dV/V = dnu to hold, but V is stored as M here.

    // Predictor: straight characteristics, source terms evaluated at the parents.
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // NOTE: watch out for singularities when a parent is on the centerline.
    // The averaging of y should prevent issues.
    double L = sin(p2.mu)*sin(p2.theta)/(0.5*(p2.y+p3.y) * cos(p2.theta + p2.mu));
    double M = sin(p1.mu)*sin(p1.theta)/(0.5*(p1.y+p3.y) * cos(p1.theta - p1.mu));

    double K_plus = (p2.theta - p2.nu) - L*(p3.x - p2.x);
    double K_minus = (p1.theta + p1.nu) + M*(p3.x - p1.x);
    p3.theta = 0.5*(K_minus + K_plus);
    MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
        p3, 0.5*(K_minus - K_plus), 0.5*(p1.mach + p2.mach));
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }

    // Corrector: average the characteristic angles and source terms across the step.
    double c_minus_angle_corrected = average_cminus_angle(p1, p3);
    double c_plus_angle_corrected = average_cplus_angle(p2, p3);
    auto [x_corrected, y_corrected] = characteristic_intersection_with_angle(
        p1, p2,
        c_minus_angle_corrected,
        c_plus_angle_corrected);
    p3.x = x_corrected;
    p3.y = y_corrected;
    double mu1_avg = 0.5*(p1.mu + p3.mu);
    double mu2_avg = 0.5*(p2.mu + p3.mu);
    double theta1_avg = 0.5*(p1.theta + p3.theta);
    double theta2_avg = 0.5*(p2.theta + p3.theta);
    double M_avg = sin(mu1_avg)*sin(theta1_avg)/(0.5*(p1.y+p3.y) * cos(c_minus_angle_corrected));
    double L_avg = sin(mu2_avg)*sin(theta2_avg)/(0.5*(p2.y+p3.y) * cos(c_plus_angle_corrected));

    K_plus = (p2.theta - p2.nu) - L_avg*(p3.x - p2.x);
    K_minus = (p1.theta + p1.nu) + M_avg*(p3.x - p1.x);
    p3.theta = 0.5*(K_minus + K_plus);
    thermo_error = update_thermodynamic_state_from_nu(p3, 0.5*(K_minus - K_plus), p3.mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }
    p3.K_minus = p3.theta + p3.nu;
    p3.K_plus = p3.theta - p3.nu;

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        log_warning("Non-downstream intersection at ({}, {}). Parents at x=({}, {}).",
            p3.x, p3.y, p1.x, p2.x);
        return {p3, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(p3, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid interior point at ({}, {}): error code {}.",
            p3.x, p3.y, to_string(validity));
        return {p3, validity};
    }

    return {p3, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_interior_point_iterative(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    // General procedure:
    // 1. Receive theta and nu from the two incoming characteristics
    // 2. Estimate characteristic slopes at parent nodes
    // 3. Intersect straight line characteristics
    // 4. Compute source terms (if applicable)
    // 5. Solve nonlinear equation in V using upstream states only
    // 6. Get theta from either compatibility equation 
    // 7. Correct characteristic slope by averaging with endpoint 
    // 8. Correct V by averaging with endpoint 
    // 9. Recompute theta and source terms at corrected location
    // 10. [optional] One additional pass if source term changed significantly
    CharacteristicPoint p3{};
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;
    double S1 = cminus_source_term(p1, y);
    double S2 = cplus_source_term(p2, y);

    double mach = find_node_mach(p1, p2, S1 - S2, 0.5 * (p1.mach + p2.mach));
    if (mach < 1.0) {
        return {p3, MocErrorCode::PM_INVERSION_FAILED};
    }
    MocErrorCode thermo_error = update_thermodynamic_state_from_mach(p3, mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }

    // from compatibility equation.
    p3.theta = S1 + p1.theta - (p3.nu - p1.nu);
    p3.K_plus = p3.theta - p3.nu;
    p3.K_minus = p3.theta + p3.nu;

    // corrector step
    const int max_iters = 3;

    for (int i = 0; i < max_iters; i++) {
        c_minus_angle = average_cminus_angle(p1, p3);
        c_plus_angle = average_cplus_angle(p2, p3);
        auto [new_x, new_y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
        p3.x = new_x;
        p3.y = new_y;
        double S1_new = cminus_source_term(p1, new_y);
        double S2_new = cplus_source_term(p2, new_y);
        double mach_new = find_node_mach(p1, p2, S1_new - S2_new, p3.mach);
        if (mach_new < 1.0) {
            return {p3, MocErrorCode::PM_INVERSION_FAILED};
        }
        thermo_error = update_thermodynamic_state_from_mach(p3, mach_new);
        if (thermo_error != MocErrorCode::NONE) {
            return {p3, thermo_error};
        }

        p3.theta = S1_new + p1.theta - (p3.nu - p1.nu);
        p3.K_plus = p3.theta - p3.nu;
        p3.K_minus = p3.theta + p3.nu;
        double residual = std::max(std::abs(S1_new - S1), std::abs(S2_new - S2));
        if (residual < m_options.solver_options.abstol) break;
        S1 = S1_new;
        S2 = S2_new;

        if (i == max_iters - 1 && residual >= m_options.solver_options.abstol) {
            log_warning("Iterative interior solver did not converge."
                "Residual={} after {} iterations at ({}, {}).", residual, max_iters, p3.x, p3.y);
        }
    }

    MocErrorCode validity = check_point_validity(p3, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {p3, validity};
    }
    return {p3, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_wall_flow(
    const CharacteristicPoint& interior_parent,
    double theta_wall)
{
    CharacteristicPoint wall_point{};
    // geometric constraint: the flow at the wall must follow the wall curvature
    wall_point.theta = theta_wall;
    // Predictor: planar K+ preservation
    // Flow solve: compatibility equation along C+ from interior parent
    // K+ is preserved: theta - nu = const along C+
    // this only applies for planar flow.

    // For axisymmetric flow, K+ is not preserved along C+.
    // Use predictor-corrector: start with planar approximation,
    // then correct using the source term.
    // The source term correction requires y, which is computed later
    // in solve_wall_point_design/analysis. We store the predictor result
    // and the corrector is applied in those methods.
    // For now, this gives a first-order approximation.

    wall_point.K_plus = interior_parent.K_plus;
    MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
        wall_point,
        wall_point.theta - wall_point.K_plus,
        interior_parent.mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {wall_point, thermo_error};
    }
    wall_point.K_minus = wall_point.theta + wall_point.nu;

    return {wall_point, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_wall_point_design(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    double theta_wall)
{
    PointResult flow_result = solve_wall_flow(interior_parent, theta_wall);
    if (flow_result.error != MocErrorCode::NONE) {
        flow_result.point.x = interior_parent.x;
        flow_result.point.y = interior_parent.y;
        return flow_result;
    }
    CharacteristicPoint wall_point = flow_result.point;

    // the angle between the interior point and the wall point
    // is given by the C+ characteristic.
    double c_plus_angle = average_cplus_angle(interior_parent, wall_point);

    // the angle between the previous wall point and the current wall point
    // is given by wall_theta.
    double prev_wall_angle = 0.5*(previous_wall_point.theta + wall_point.theta);

    auto [x,y] = characteristic_intersection_with_angle(
        interior_parent, previous_wall_point,
        c_plus_angle, prev_wall_angle);

    wall_point.x = x;
    wall_point.y = y;

    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: { break; }
        case MocFlowKind::AXISYMMETRIC: {
            // Wall corrector step
            // Unlike the case of an interior point, we only have one characteristic line
            // and theta is already known (specified by wall).
            // The general compatibility equation for the C+ characteristic is:
            // $theta_P - theta_B - S_{BP} = nu_P - nu_B$
            // where S is the source term, B is the interior point and P is the wall point.

            // The solution procedure is:
            // 1. Calculate source term from estimated y.
            // 2. Calculate nu from axisymmetric compatibility equation.
            // 3. Obtain thermodynamic values for given value of nu.
            // 4. Recompute new position based on new C+ angle
            // 5. Repeat until source residual is satisfactory.
            // In practice only 1-2 iterations should be needed.
            double old_S = 0.0;
            double residual = 0.0;
            const int max_iter = 4;
            for (int i = 0; i < max_iter; i++){

                double S = cplus_source_term(interior_parent, wall_point);
                residual = std::abs(S - old_S);
                old_S = S;
                if (residual < m_options.solver_options.abstol) break;

                // C+ compatibility nu_P = theta_P - theta_B + nu_B + S, matching the interior
                // solver's K+ = theta - nu decreasing by the source along C+ (d(theta-nu) = -S).
                // The source enters with a +S sign here, NOT -S.
                MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                    wall_point,
                    wall_point.theta-interior_parent.theta + S + interior_parent.nu,
                    interior_parent.mach);
                if (thermo_error != MocErrorCode::NONE) {
                    return {wall_point, thermo_error};
                }

                wall_point.K_plus = wall_point.theta - wall_point.nu;
                wall_point.K_minus = wall_point.theta + wall_point.nu;
                double corrected_cplus_angle = average_cplus_angle(interior_parent, wall_point);
                auto [new_x,new_y] = characteristic_intersection_with_angle(
                    interior_parent,
                    previous_wall_point,
                    corrected_cplus_angle,
                    prev_wall_angle);

                wall_point.x = new_x;
                wall_point.y = new_y;

                if (i == max_iter - 1 && residual >= m_options.solver_options.abstol) {
                    log_warning("Wall design source term iteration did not converge."
                        "Residual={} at ({}, {}).", residual, wall_point.x, wall_point.y);
                }
            }
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    MocErrorCode validity = check_point_validity(wall_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid wall point at ({}, {}): error code {}.",
            wall_point.x, wall_point.y, to_string(validity));
        return {wall_point, validity};
    }

    return {wall_point, MocErrorCode::NONE};
}

std::optional<PointResult> MocNozzle::solve_wall_point_analysis(
        const CharacteristicPoint& interior_parent)
{
    const NozzleProfile& wall = m_options.nozzle_profile;
    auto [x_wall, y_wall] = intersect_characteristic_with_wall(
                interior_parent, wall);

    // no wall hit found
    if (x_wall < 0.0 && y_wall < 0.0) {
        return std::nullopt;
    }

    // theta_at() throws (std::runtime_error, "Queried x ... larger than nozzle
    // profile") when the predicted C+ hits beyond the last profile station -- a
    // known consequence of the analysis-outflow gap (characteristics leaving the
    // contour are not yet terminated as OUTFLOW before the wall solve runs; see
    // instructions/moc_algorithm.md Sec. 9.1, a Phase 2 fix). Converting that here
    // instead of letting it escape is exactly what keeps solve() from throwing.
    double wall_theta;
    try {
        wall_theta = wall.theta_at(x_wall);
    } catch (const std::exception& e) {
        log_warning("Wall profile query out of bounds at x={}: {}", x_wall, e.what());
        CharacteristicPoint failed_point{};
        failed_point.x = interior_parent.x;
        failed_point.y = interior_parent.y;
        return PointResult{failed_point, MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS};
    }

    double char_angle = interior_parent.theta + interior_parent.mu;
    PointResult flow_result = solve_wall_flow(interior_parent, wall_theta);
    if (flow_result.error != MocErrorCode::NONE) {
        flow_result.point.x = interior_parent.x;
        flow_result.point.y = interior_parent.y;
        return flow_result;
    }
    CharacteristicPoint wall_point = flow_result.point;

    wall_point.x = x_wall;
    wall_point.y = y_wall;
    // Corrector (for curved characteristics):
    // Solve the wall point flow using computed theta
    // then recompute characteristic slope as average of parent and wall-point slopes.
    // One correction is usually sufficient.
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {

            //for planar flow, solve_wall_flow equations are exact.
            // Corrector: re-intersect the C+ ray from the interior parent using the average
            // of the parent and predicted wall-point characteristic angles.
            double wall_char_angle = wall_point.theta + wall_point.mu;
            double average_char_angle = 0.5*(wall_char_angle + char_angle);

            auto [x_corrected, y_corrected] = find_wall_hit(interior_parent, wall, average_char_angle);
            try {
                wall_theta = wall.theta_at(x_corrected);
            } catch (const std::exception& e) {
                log_warning("Wall profile query out of bounds at x={}: {}", x_corrected, e.what());
                wall_point.x = interior_parent.x;
                wall_point.y = interior_parent.y;
                return PointResult{wall_point, MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS};
            }
            // Recompute the flow at the corrected wall angle, then restore the position:
            // solve_wall_flow returns flow properties only and leaves x,y at their defaults.
            flow_result = solve_wall_flow(interior_parent, wall_theta);
            if (flow_result.error != MocErrorCode::NONE) {
                flow_result.point.x = interior_parent.x;
                flow_result.point.y = interior_parent.y;
                return flow_result;
            }
            wall_point = flow_result.point;
            wall_point.x = x_corrected;
            wall_point.y = y_corrected;
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {
            // For axisymmetric flow, we must account for the source term.
            double old_S = 0.0;
            const int max_iter = 4;
            for (int i = 0; i < max_iter; i++){

                double c_plus_angle = average_cplus_angle(interior_parent, wall_point);
                auto [x_corrected, y_corrected] = find_wall_hit(interior_parent, wall, c_plus_angle);
                // The corrected C+ angle can miss the contour even when the predictor hit
                // it (e.g. near the exit lip). Terminate the characteristic as outflow rather
                // than marching on with a (-1, -1) miss sentinel.
                if (x_corrected < 0.0 && y_corrected < 0.0) {
                    return std::nullopt;
                }
                wall_point.x = x_corrected;
                wall_point.y = y_corrected;
                // In analysis mode the wall point's flow angle is set by the wall, not
                // prescribed as in design mode. Update it to the wall angle at the corrected
                // position before applying the C+ compatibility relation; otherwise the
                // reflected wave is computed from a stale predictor angle.
                try {
                    wall_theta = wall.theta_at(wall_point.x);
                } catch (const std::exception& e) {
                    log_warning("Wall profile query out of bounds at x={}: {}", wall_point.x, e.what());
                    return PointResult{wall_point, MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS};
                }
                wall_point.theta = wall_theta;

                double S = cplus_source_term(interior_parent, wall_point);
                double residual = std::abs(S - old_S);
                if (residual < m_options.solver_options.abstol) break;
                old_S = S;

                // C+ compatibility nu_P = theta_P - theta_B + nu_B + S (see solve_wall_point_design).
                MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                    wall_point,
                    wall_point.theta-interior_parent.theta + S + interior_parent.nu,
                    interior_parent.mach);
                if (thermo_error != MocErrorCode::NONE) {
                    return PointResult{wall_point, thermo_error};
                }

                if (i == max_iter - 1 && residual >= m_options.solver_options.abstol) {
                    log_warning("Wall analysis source term iteration did not converge."
                        "Residual={} at ({}, {}).", residual, wall_point.x, wall_point.y);
                }
            }
            wall_point.K_plus = wall_point.theta - wall_point.nu;
            wall_point.K_minus = wall_point.theta + wall_point.nu;
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    MocErrorCode validity = check_point_validity(wall_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid analysis wall point at ({}, {}): error code {}.",
            wall_point.x, wall_point.y, to_string(validity));
        return PointResult{wall_point, validity};
    }

    return PointResult{wall_point, MocErrorCode::NONE};
}

double MocNozzle::maximum_nozzle_length() const {
    return m_options.nozzle_profile.length();
}

double MocNozzle::gamma_s_from_mach(double mach) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return m_options.gamma;
        }
        case GasChemistry::FROZEN: 
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_gamma_s_from_mach(mach);
        }
        default:
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    // unreachable
    return m_options.gamma;
}

double MocNozzle::gamma_s_from_nu(double nu) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return m_options.gamma;
        }
        case GasChemistry::FROZEN: 
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_gamma_s_from_nu(nu);
        }
        default:
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    // unreachable
    return m_options.gamma;
}

double MocNozzle::mach_from_nu(const CharacteristicPoint& point, double mach_guess) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return mach_from_prandtl_meyer(point.nu, point.gamma_s, mach_guess);
        }
        case GasChemistry::FROZEN: 
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_mach(point.nu);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

double MocNozzle::nu_from_mach(double mach) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return prandtl_meyer(mach, m_options.gamma);
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_nu_from_mach(mach);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

double MocNozzle::find_node_mach(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2,
    double source_delta,
    double mach_guess)
{

    // initial guess: mach of upstream
    double mach = (mach_guess > 1.0) ? mach_guess : p1.mach;

    double delta_theta = p1.theta - p2.theta;
    const int max_iter = 15;
    // Newton's method to find root
    // Residual function obtained from combining both compatibility equations:
    //
    //$$\epsilon (V_P) = (\theta_{A}- \theta_{B}) + \Delta \nu (V_{A} \rightarrow V_{P})
    // + \Delta \nu (V_{B} \rightarrow V_{P}) - (S_{A}- S_{B})$$
    // where A and B indicate upstream points.
    for (int i = 0; i < max_iter; i++) {
        double nu3, derivative;

        if (m_options.chemistry == GasChemistry::PERFECT_GAS) {
            nu3 = prandtl_meyer(mach, m_options.gamma);
            // dnu/dM for perfect gas
            derivative = 2.0 * prandtl_meyer_derivative(mach, m_options.gamma);
        } else {
            auto [idx, weight] = pm_table.find_mach_index_and_weight(mach);
            nu3 = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
            double V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
            derivative = 2.0 * sqrt(mach * mach - 1.0) / V;
        }

        double residual = delta_theta - source_delta + (nu3 - p1.nu) + (nu3 - p2.nu);
        if (std::abs(residual) < m_options.solver_options.abstol) return mach;

        mach -= residual / derivative;
    }

    // rootfinding has failed
    log_warning("find_node_mach did not converge after {} iterations. "
        "Last Mach={}, delta_theta={}.", max_iter, mach, delta_theta);
    return -1.0;
}

void MocNozzle::update_thermodynamic_state(CharacteristicPoint& point) {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            // the choice of upstream point can be arbitrary due to Crocco's theorem
            // stagnation factor at throat is = 1 by definition
            double stagnation_ratio = stagnation_factor(point.mach, point.gamma_s);
            point.temperature = m_T_ref / stagnation_ratio;
            point.pressure = m_P_ref / pow(stagnation_ratio, point.gamma_s / (point.gamma_s - 1.0));
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            // Guard against points without a stored Cantera state (e.g. the bare
            // throat-lip seed); restoring an empty state vector aborts inside Cantera.
            if (point.cantera_state.empty()) break;
            m_gas->thermo()->restoreState(point.cantera_state);
            point.temperature = m_gas->thermo()->temperature() / m_T_ref;
            point.pressure = m_gas->thermo()->pressure() / m_P_ref;
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

MocErrorCode MocNozzle::update_thermodynamic_state_from_nu(
    CharacteristicPoint& point,
    double nu, double mach_guess)
{
    if (nu < 0.0) {
        log_warning("Negative Prandtl-Meyer angle nu = {}. Subsonic flow or numerical error.", nu);
    }
    point.nu = nu;
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.mach = mach_from_prandtl_meyer(nu, point.gamma_s, mach_guess);
            // mach_from_prandtl_meyer's documented failure sentinel: the Newton
            // solve for the inverse Prandtl-Meyer function did not converge. Report
            // it as data instead of letting -1.0 flow into mach_to_mu() and
            // corrector-step Mach averages downstream.
            if (point.mach < 0.0) {
                log_warning("Prandtl-Meyer inversion failed to converge for nu = {}.", nu);
                return MocErrorCode::PM_INVERSION_FAILED;
            }
            point.V = point.mach;
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            try {
                auto [idx, weight] = pm_table.find_nu_index_and_weight(nu);
                point.V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
                point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
                point.mach = pm_table.interpolate_at_index(idx, weight, pm_table.machs);
                point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            } catch (const std::out_of_range& e) {
                log_warning("Prandtl-Meyer table lookup by nu = {} out of range: {}", nu, e.what());
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

MocErrorCode MocNozzle::update_thermodynamic_state_from_mach(CharacteristicPoint& point, double mach) {
    point.mach = mach;

    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.nu = prandtl_meyer(point.mach, point.gamma_s);
            point.V = mach;
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            try {
                auto [idx, weight] = pm_table.find_mach_index_and_weight(mach);
                point.V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
                point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
                point.nu = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
                point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            } catch (const std::out_of_range& e) {
                log_warning("Prandtl-Meyer table lookup by mach = {} out of range: {}", mach, e.what());
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

MocErrorCode MocNozzle::update_thermodynamic_state_from_V(CharacteristicPoint& point, double V) {
    point.V = V;

    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.mach = V; // for a perfect gas, the velocity is kept dimensionless
            point.nu = prandtl_meyer(point.mach, point.gamma_s);
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            try {
                auto [idx, weight] = pm_table.find_V_index_and_weight(V);
                point.mach = pm_table.interpolate_at_index(idx, weight, pm_table.machs);
                point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
                point.nu = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
                point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            } catch (const std::out_of_range& e) {
                log_warning("Prandtl-Meyer table lookup by V = {} out of range: {}", V, e.what());
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

ThermodynamicContext MocNozzle::build_thermo_context() {
    return ThermodynamicContext{
        .gas = m_gas, 
        .table = pm_table, 
        .T_ref = m_T_ref, 
        .P_ref = m_P_ref, 
        .gamma_s = m_options.gamma};
}

std::pair<double,double> MocNozzle::find_wall_hit(
    const CharacteristicPoint& p, 
    const NozzleProfile& wall, 
    double char_angle) const
{
    double c_plus_slope = tan(char_angle);
    double x,y;
    bool found = false;

    // A wall hit needs at least one segment; a degenerate profile would underflow the
    // loop bound and the extrapolation indices below.
    if (wall.x.size() < 2) {
        throw std::runtime_error("find_wall_hit: nozzle profile has fewer than two points.");
    }

    for (size_t k = 0; k < wall.x.size() - 1; k++) {
        double wall_slope = (wall.y[k+1] - wall.y[k])
            / (wall.x[k+1] - wall.x[k]);
        double wall_intercept = wall.y[k] - wall_slope * wall.x[k];
        double char_intercept = p.y - c_plus_slope * p.x;

        double x_intercept = (wall_intercept - char_intercept) / (c_plus_slope - wall_slope);

        // check if intersection is within the wall segment
        if (x_intercept >= wall.x[k] && x_intercept <= wall.x[k+1]) {
            x = x_intercept;
            y = p.y + c_plus_slope * (x - p.x);
            found = true;
            break;
        }
    }
    if (!found) { 
        // if the intersection cannot be found within the nozzle profile, 
        // trun negative values
        return {-1.0, -1.0};
        // size_t last = wall.x.size() - 1;
        // double wall_slope = (wall.y[last] - wall.y[last-1])/(wall.x[last] - wall.x[last - 1]);
        // double wall_intercept = wall.y[last] - wall_slope * wall.x[last];
        // double char_intercept = p.y - c_plus_slope * p.x;
        // x = (wall_intercept - char_intercept) / (c_plus_slope - wall_slope);
        // y = p.y + c_plus_slope * (x - p.x);
        // log_warning("Extrapolating wall intercept to x = {}, y = {}", x, y);
    }
    return {x,y};
}

double MocNozzle::cplus_source_term(
    const CharacteristicPoint& p, 
    double new_y) const 
{
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return sin(p.theta)/(p.mach * sin(p.theta - p.mu)) * dy/y_avg;
}

double MocNozzle::cplus_source_term(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p3) const 
{
    double y_avg = 0.5 * (p1.y + p3.y);
    double dy = p3.y - p1.y;
    double theta_avg = 0.5 * (p1.theta + p3.theta);
    double theta_plus_mu_avg = 0.5 * ((p1.theta + p1.mu) + (p3.theta + p3.mu));
    double mach_avg = 0.5 * (p1.mach + p3.mach);
    return sin(theta_avg)/(mach_avg * sin(theta_plus_mu_avg)) * dy/y_avg;
}

double MocNozzle::cminus_source_term(
    const CharacteristicPoint& p, 
    double new_y) const 
{
    
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return -sin(p.theta)/(p.mach * sin(p.theta - p.mu)) * dy/y_avg;
}

MocNozzle::LeadingEdgeView MocNozzle::leading_edges(const CharacteristicNet& net, ChainMetadata::Family fam) const {
    LeadingEdgeView view;
    view.family = fam;
    for (size_t i = 0; i < net.chain_metadata.size(); i++) {
        const ChainMetadata& meta = net.chain_metadata[i];
        if (meta.active && meta.family == fam) {
            view.y_values.push_back(net.points[meta.latest_point_idx].y);
            view.chain_indices.push_back(i);
            view.leading_pt_indices.push_back(meta.latest_point_idx);
        }
    }
    return view;
}

void MocNozzle::update_leading_edges(LeadingEdgeView& view, const CharacteristicNet& net, ChainMetadata::Family family) const {
    view.y_values.clear();
    view.chain_indices.clear();
    view.leading_pt_indices.clear();
    for (size_t i = 0; i < net.chain_metadata.size(); i++) {
        const ChainMetadata& meta = net.chain_metadata[i];
        if (meta.active && meta.family == family) {
            view.y_values.push_back(net.points[meta.latest_point_idx].y);
            view.chain_indices.push_back(i);
            view.leading_pt_indices.push_back(meta.latest_point_idx);
        }
    }
}

void MocNozzle::sort_plus_edges_by_proximity(LeadingEdgeView& view, const CharacteristicNet& net) const {
    std::vector<size_t> order(view.chain_indices.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (view.y_values[a] != view.y_values[b]) return view.y_values[a] > view.y_values[b];
        return net.points[view.leading_pt_indices[a]].x < net.points[view.leading_pt_indices[b]].x;
    });

    std::vector<double> y_sorted(order.size());
    std::vector<size_t> chain_sorted(order.size());
    std::vector<size_t> pt_sorted(order.size());
    for (size_t k = 0; k < order.size(); k++) {
        y_sorted[k] = view.y_values[order[k]];
        chain_sorted[k] = view.chain_indices[order[k]];
        pt_sorted[k] = view.leading_pt_indices[order[k]];
    }
    view.y_values = std::move(y_sorted);
    view.chain_indices = std::move(chain_sorted);
    view.leading_pt_indices = std::move(pt_sorted);
}

MocNozzle::PairSearch MocNozzle::find_pair_partner(
    double plus_y,
    const LeadingEdgeView& minus_edges,
    const std::vector<bool>& claimed,
    const CharacteristicNet& net) const
{
    PairSearch search;
    double best_dy = std::numeric_limits<double>::max();

    // Nearest C- above regardless of claim status -- tracked separately from the best
    // *available* partner so a C+ whose true nearest partner was already claimed by a
    // closer competitor is told to wait rather than settling for a farther one.
    double nearest_dy = std::numeric_limits<double>::max();
    bool nearest_is_claimed = false;

    for (size_t j = 0; j < minus_edges.chain_indices.size(); j++) {
        double dy = minus_edges.y_values[j] - plus_y;
        if (dy <= 0) continue;
        search.any_cminus_above = true;

        if (dy < nearest_dy) {
            nearest_dy = dy;
            nearest_is_claimed = claimed[j];
        }
        if (claimed[j]) continue; //skip if already paired

        if (dy < best_dy) {
            best_dy = dy;
            search.chain_idx = minus_edges.chain_indices[j];
            search.pt_idx = minus_edges.leading_pt_indices[j];
            search.edgevec_idx = j;
        }
        // edge case: multiple points at literally the same y due to an expansion fan
        // the tiebreaker is determined by the characteristic angle theta-mu
        if (search.chain_idx.has_value() && best_dy == dy) [[unlikely]] {
            const CharacteristicPoint& old_pt = net.points[search.pt_idx];
            const CharacteristicPoint& new_pt = net.points[minus_edges.leading_pt_indices[j]];
            if (new_pt.theta - new_pt.mu < old_pt.theta - old_pt.mu) {
                search.chain_idx = minus_edges.chain_indices[j];
                search.pt_idx = minus_edges.leading_pt_indices[j];
                search.edgevec_idx = j;
            }
        }
    }

    if (nearest_is_claimed) search.chain_idx = std::nullopt;
    search.dy = best_dy;
    return search;
}

std::vector<MocNozzle::FrontSegment> MocNozzle::build_front(
    const CharacteristicNet& net,
    const LeadingEdgeView& plus_edges,
    const LeadingEdgeView& minus_edges) const
{
    std::vector<FrontSegment> front;
    front.reserve(plus_edges.chain_indices.size());
    std::vector<bool> claimed(minus_edges.chain_indices.size(), false);

    for (size_t i = 0; i < plus_edges.chain_indices.size(); i++) {
        PairSearch match = find_pair_partner(plus_edges.y_values[i], minus_edges, claimed, net);
        // No available partner: this C+ is either wall-bound or waiting. Neither is a
        // front segment, and neither is something mesh control can act on.
        if (!match.chain_idx.has_value()) continue;
        claimed[match.edgevec_idx] = true;

        const CharacteristicPoint& below = net.points[plus_edges.leading_pt_indices[i]];
        const CharacteristicPoint& above = net.points[match.pt_idx];
        front.push_back(FrontSegment{
            plus_edges.leading_pt_indices[i],
            match.pt_idx,
            plus_edges.chain_indices[i],
            *match.chain_idx,
            std::hypot(above.x - below.x, above.y - below.y)
        });
    }
    return front;
}

std::optional<std::pair<double, double>> MocNozzle::cell_sides(
    const CharacteristicNet& net, const FrontSegment& segment) const
{
    const CharacteristicPoint& below = net.points[segment.lower_pt_idx];
    const CharacteristicPoint& above = net.points[segment.upper_pt_idx];
    const double dx = above.x - below.x;
    const double dy = above.y - below.y;

    const double plus_angle = below.theta + below.mu;
    const double minus_angle = above.theta - above.mu;
    const double ux_p = std::cos(plus_angle), uy_p = std::sin(plus_angle);
    const double ux_m = std::cos(minus_angle), uy_m = std::sin(minus_angle);

    const double det = ux_m * uy_p - ux_p * uy_m;
    if (std::abs(det) < 1e-14) return std::nullopt;

    const double side_plus = (ux_m * dy - uy_m * dx) / det;
    const double side_minus = (ux_p * dy - uy_p * dx) / det;
    if (!(side_plus > 0.0) || !(side_minus > 0.0)) return std::nullopt;
    return std::make_pair(side_plus, side_minus);
}

std::vector<size_t> MocNozzle::previous_front_points(
    const CharacteristicNet& net, const FrontSegment& segment) const
{
    std::vector<size_t> found;

    // The predecessor of a point within a chain it currently leads. Anything else (a chain
    // the point does not lead, or one with no history) contributes nothing.
    auto predecessor = [&net, &found](std::optional<size_t> chain_idx, size_t pt_idx) {
        if (!chain_idx.has_value()) return;
        const std::vector<size_t>& chain = net.c_chains[*chain_idx];
        if (chain.size() < 2 || chain.back() != pt_idx) return;
        found.push_back(chain[chain.size() - 2]);
    };

    predecessor(segment.lower_plus_chain, segment.lower_pt_idx);
    predecessor(net.membership[segment.lower_pt_idx].c_minus_chain_idx, segment.lower_pt_idx);
    predecessor(net.membership[segment.upper_pt_idx].c_plus_chain_idx, segment.upper_pt_idx);
    predecessor(segment.upper_minus_chain, segment.upper_pt_idx);

    std::sort(found.begin(), found.end());
    found.erase(std::unique(found.begin(), found.end()), found.end());
    std::sort(found.begin(), found.end(), [&net](size_t a, size_t b) {
        return net.points[a].y < net.points[b].y;
    });
    return found;
}

namespace {

/**
 * Walk upstream from (x, y) along a characteristic of the given angle until the polyline is
 * met. Returns the segment index and the parameter along it, or nothing when the ray misses.
 */
std::optional<std::pair<size_t, double>> trace_back_to_front(
    double x, double y, double angle,
    const CharacteristicNet& net, const std::vector<size_t>& front)
{
    const double dir_x = -std::cos(angle);
    const double dir_y = -std::sin(angle);

    for (size_t i = 0; i + 1 < front.size(); i++) {
        const CharacteristicPoint& a = net.points[front[i]];
        const CharacteristicPoint& b = net.points[front[i + 1]];
        const double seg_x = b.x - a.x;
        const double seg_y = b.y - a.y;

        const double det = dir_x * (-seg_y) - dir_y * (-seg_x);
        if (std::abs(det) < 1e-14) continue; // ray parallel to this segment

        const double rhs_x = a.x - x;
        const double rhs_y = a.y - y;
        const double s = (rhs_x * (-seg_y) - rhs_y * (-seg_x)) / det;
        const double u = (dir_x * rhs_y - dir_y * rhs_x) / det;

        // s > 0 keeps the foot strictly upstream; u within [0,1] keeps it on the segment.
        // Both are refusals rather than clamps: a clamped foot is an extrapolation dressed
        // up as an interpolation.
        if (s > 0.0 && u >= 0.0 && u <= 1.0) return std::make_pair(i, u);
    }
    return std::nullopt;
}

} // namespace

PointResult MocNozzle::solve_inverse_interior_point(
    double x_new, double y_new,
    const CharacteristicNet& net,
    const std::vector<size_t>& previous_front,
    const CharacteristicPoint& seed)
{
    CharacteristicPoint p4{};
    p4.x = x_new;
    p4.y = y_new;
    p4.theta = seed.theta;
    MocErrorCode thermo_error = update_thermodynamic_state_from_nu(p4, seed.nu, seed.mach);
    if (thermo_error != MocErrorCode::NONE) return {p4, thermo_error};

    if (previous_front.size() < 2) return {p4, MocErrorCode::INITIALIZATION_FAILED};

    const bool axisymmetric = (m_options.flow_type == MocFlowKind::AXISYMMETRIC);

    // Feet move as the state at p4 is refined, so iterate. Four passes matches the corrector
    // budget the wall solvers use; in practice this settles in two.
    constexpr int max_iters = 4;
    for (int iter = 0; iter < max_iters; iter++) {
        const double angle_minus = p4.theta - p4.mu;
        const double angle_plus = p4.theta + p4.mu;

        std::optional<std::pair<size_t, double>> hit_minus =
            trace_back_to_front(p4.x, p4.y, angle_minus, net, previous_front);
        std::optional<std::pair<size_t, double>> hit_plus =
            trace_back_to_front(p4.x, p4.y, angle_plus, net, previous_front);
        if (!hit_minus.has_value() || !hit_plus.has_value()) {
            return {p4, MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS};
        }

        // Interpolate each foot along the previous front.
        auto foot_at = [&net, &previous_front](const std::pair<size_t, double>& hit) {
            const CharacteristicPoint& a = net.points[previous_front[hit.first]];
            const CharacteristicPoint& b = net.points[previous_front[hit.first + 1]];
            const double u = hit.second;
            CharacteristicPoint f{};
            f.x = a.x + u * (b.x - a.x);
            f.y = a.y + u * (b.y - a.y);
            f.theta = a.theta + u * (b.theta - a.theta);
            f.nu = a.nu + u * (b.nu - a.nu);
            f.mu = a.mu + u * (b.mu - a.mu);
            f.mach = a.mach + u * (b.mach - a.mach);
            return f;
        };
        const CharacteristicPoint foot_minus = foot_at(*hit_minus);
        const CharacteristicPoint foot_plus = foot_at(*hit_plus);

        // Same source terms as solve_interior_point_axisymmetric, evaluated on the averages
        // between each foot and p4:
        //   along C+:  d(theta - nu) = -L dx,  L = sin(mu) sin(theta) / (y cos(theta + mu))
        //   along C-:  d(theta + nu) = +M dx,  M = sin(mu) sin(theta) / (y cos(theta - mu))
        double source_minus = 0.0;
        double source_plus = 0.0;
        if (axisymmetric) {
            const double mu_m = 0.5 * (foot_minus.mu + p4.mu);
            const double th_m = 0.5 * (foot_minus.theta + p4.theta);
            const double y_m = 0.5 * (foot_minus.y + p4.y);
            const double ang_m = 0.5 * (angle_minus + (foot_minus.theta - foot_minus.mu));
            if (std::abs(y_m) > 1e-12) source_minus = sin(mu_m) * sin(th_m) / (y_m * cos(ang_m));

            const double mu_p = 0.5 * (foot_plus.mu + p4.mu);
            const double th_p = 0.5 * (foot_plus.theta + p4.theta);
            const double y_p = 0.5 * (foot_plus.y + p4.y);
            const double ang_p = 0.5 * (angle_plus + (foot_plus.theta + foot_plus.mu));
            if (std::abs(y_p) > 1e-12) source_plus = sin(mu_p) * sin(th_p) / (y_p * cos(ang_p));
        }

        const double K_minus = (foot_minus.theta + foot_minus.nu)
                             + source_minus * (p4.x - foot_minus.x);
        const double K_plus = (foot_plus.theta - foot_plus.nu)
                            - source_plus * (p4.x - foot_plus.x);

        const double theta_previous = p4.theta;
        p4.theta = 0.5 * (K_minus + K_plus);
        thermo_error = update_thermodynamic_state_from_nu(p4, 0.5 * (K_minus - K_plus), p4.mach);
        if (thermo_error != MocErrorCode::NONE) return {p4, thermo_error};

        if (std::abs(p4.theta - theta_previous) < m_options.solver_options.abstol) break;
    }

    p4.K_minus = p4.theta + p4.nu;
    p4.K_plus = p4.theta - p4.nu;

    MocErrorCode validity = check_point_validity(p4, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) return {p4, validity};
    return {p4, MocErrorCode::NONE};
}

void MocNozzle::record_front_diagnostics(
    const CharacteristicNet& net,
    const std::vector<FrontSegment>& front,
    double target_spacing,
    MocPassDiagnostics& diag) const
{
    diag.target_spacing = target_spacing;
    diag.front_points = front.empty() ? 0 : front.size() + 1;
    if (front.empty()) return;

    double min_spacing = std::numeric_limits<double>::max();
    double max_spacing = 0.0;
    double total = 0.0;
    double max_aspect = 0.0;
    double min_margin = std::numeric_limits<double>::max();

    // The front's two ends, found by scanning rather than by assuming an ordering:
    // build_front inherits the C+ edge sort, which is currently descending in y, and this
    // must not silently invert if that changes.
    double axis_y = std::numeric_limits<double>::max();
    double wall_y = -std::numeric_limits<double>::max();

    for (const FrontSegment& seg : front) {
        const CharacteristicPoint& below = net.points[seg.lower_pt_idx];
        const CharacteristicPoint& above = net.points[seg.upper_pt_idx];
        if (below.y < axis_y) {
            axis_y = below.y;
            diag.front_axis_x = below.x;
            diag.front_axis_spacing = seg.arc;
        }
        if (above.y > wall_y) {
            wall_y = above.y;
            diag.front_wall_x = above.x;
            diag.front_wall_spacing = seg.arc;
        }
        min_spacing = std::min(min_spacing, seg.arc);
        max_spacing = std::max(max_spacing, seg.arc);
        total += seg.arc;

        const double dx = above.x - below.x;
        const double dy = above.y - below.y;
        const double plus_angle = below.theta + below.mu;
        const double minus_angle = above.theta - above.mu;

        // t/s diverges exactly as the front turns tangent to the C- family, which is the
        // observed axisymmetric failure -- the sharpest early warning there is, and one that
        // the segment's own length cannot give: bounding the diagonal bounds t and leaves s
        // free to collapse.
        if (std::optional<std::pair<double, double>> sides = cell_sides(net, seg)) {
            if (sides->first > 1e-14) {
                max_aspect = std::max(max_aspect, sides->second / sides->first);
            }
        }

        // Normalized spacelike margin, matching the definition used by the diagnosis
        // probes: how far the segment sits from being parallel to a characteristic through
        // one of its endpoints, as a fraction of its height. Healthy fronts sit near 1;
        // reaching 0 is the NON_DOWNSTREAM_POINT failure.
        if (dy > 0.0) {
            const double margin_minus = (dy - std::tan(minus_angle) * dx) / dy;
            const double margin_plus = (dy - std::tan(plus_angle) * dx) / dy;
            min_margin = std::min({min_margin, margin_minus, margin_plus});
        }
    }

    diag.min_spacing = min_spacing;
    diag.max_spacing = max_spacing;
    diag.mean_spacing = total / static_cast<double>(front.size());
    diag.max_cell_aspect = max_aspect;
    diag.min_spacelike_margin =
        (min_margin == std::numeric_limits<double>::max()) ? 0.0 : min_margin;

}

double MocNozzle::start_line_mass_flow_error(
    const std::vector<CharacteristicPoint>& data_line) const
{
    // The integrand needs density, which for frozen/equilibrium would take the molar mass
    // that CharacteristicPoint does not carry. Report nothing rather than something wrong.
    if (m_options.chemistry != GasChemistry::PERFECT_GAS) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (data_line.size() < 2) return std::numeric_limits<double>::quiet_NaN();

    const double gamma = m_options.gamma;
    if (!(gamma > 1.0)) return std::numeric_limits<double>::quiet_NaN();

    // rho*V normalized by its stagnation reference rho0*a0, from the isentropic relations:
    //   rho/rho0 = (1 + (g-1)/2 M^2)^(-1/(g-1)),  a/a0 = (1 + (g-1)/2 M^2)^(-1/2)
    // so rho*V/(rho0*a0) = M * (1 + (g-1)/2 M^2)^(-(g+1)/(2(g-1))).
    const double exponent = -(gamma + 1.0) / (2.0 * (gamma - 1.0));
    auto mass_flux = [&](double mach) {
        return mach * std::pow(1.0 + 0.5 * (gamma - 1.0) * mach * mach, exponent);
    };

    const bool axisymmetric = (m_options.flow_type == MocFlowKind::AXISYMMETRIC);
    const double r_throat = (m_throat_radius > 0.0) ? m_throat_radius : 1.0;

    double mdot = 0.0;
    for (size_t i = 1; i < data_line.size(); i++) {
        const CharacteristicPoint& below = data_line[i - 1];
        const CharacteristicPoint& above = data_line[i];
        const double dx = above.x - below.x;
        const double dy = above.y - below.y;
        const double theta_avg = 0.5 * (below.theta + above.theta);
        const double y_avg = 0.5 * (below.y + above.y);
        const double flux = 0.5 * (mass_flux(below.mach) + mass_flux(above.mach));

        // V dot n over the segment, where n is the downstream normal (dy, -dx)/|ds|; the
        // |ds| cancels against the area element, leaving the bare differences.
        const double normal_flux = std::cos(theta_avg) * dy - std::sin(theta_avg) * dx;
        const double area_weight = axisymmetric ? (2.0 * M_PI * y_avg) : 1.0;
        mdot += flux * normal_flux * area_weight;
    }

    // 1-D critical mass flow through the throat, in the same normalization. The
    // axisymmetric case is the full disc; the planar case is the half-height the solver
    // actually meshes.
    const double throat_area = axisymmetric ? (M_PI * r_throat * r_throat) : r_throat;
    const double mdot_reference = mass_flux(1.0) * throat_area;
    if (!(std::abs(mdot_reference) > 0.0)) return std::numeric_limits<double>::quiet_NaN();

    return (mdot - mdot_reference) / mdot_reference;
}

MocInitDiagnostics MocNozzle::record_init_diagnostics(
    const std::vector<CharacteristicPoint>& data_line) const
{
    MocInitDiagnostics diag;
    diag.points = data_line.size();
    if (data_line.size() < 2) return diag;

    // The generators return the line ordered axis to wall (index 0 on the axis for both
    // the Kliegel-Levine line and the marched centered fan).
    const CharacteristicPoint& axis_pt = data_line.front();
    const CharacteristicPoint& wall_pt = data_line.back();
    const double spacing = (m_reference_spacing > 0.0) ? m_reference_spacing : 1.0;

    diag.mach_axis = axis_pt.mach;
    diag.mach_wall = wall_pt.mach;
    diag.mach_ratio = (axis_pt.mach > 0.0) ? wall_pt.mach / axis_pt.mach : 0.0;
    diag.mu_axis = axis_pt.mu;
    diag.mu_wall = wall_pt.mu;
    {
        const double cot_axis = std::tan(axis_pt.mu) > 0.0 ? 1.0 / std::tan(axis_pt.mu) : 0.0;
        const double cot_wall = std::tan(wall_pt.mu) > 0.0 ? 1.0 / std::tan(wall_pt.mu) : 0.0;
        diag.cot_mu_ratio = (cot_wall > 0.0) ? cot_axis / cot_wall : 0.0;
    }

    diag.start_line_used = m_initial_line_family.has_value()
        ? MocStartLine::CENTERED_FAN : MocStartLine::KLIEGEL_LEVINE;

    // K+ of the topmost *interior* point. The Kliegel-Levine line carries a separate wall
    // point (data_line.back()), so its topmost interior point is the one just below it; the
    // centered fan carries no such boundary point within data_line at all (solve() seeds
    // the throat lip anchor outside it), so its own last point already is the topmost
    // interior point.
    diag.kplus_wall_end =
        (diag.start_line_used == MocStartLine::KLIEGEL_LEVINE && data_line.size() >= 2)
            ? data_line[data_line.size() - 2].K_plus
            : wall_pt.K_plus;

    // Fit to the prescribed wall. Two things gate this. A design mode has no prescribed
    // contour to be off by. And a line that lies along one characteristic (the centered
    // fan) does not carry its own wall point at all -- solve() seeds the throat lip
    // separately as the wall anchor, and that lip is on the wall by construction -- so its
    // top point is an interior point and measuring a "gap" from it would be meaningless.
    const NozzleProfile& wall = m_options.nozzle_profile;
    const bool line_owns_wall_point = !m_initial_line_family.has_value();
    if (line_owns_wall_point && wall.size() >= 2 &&
        wall_pt.x >= wall.x_min() && wall_pt.x <= wall.x_max()) {
        diag.wall_gap = std::abs(wall_pt.y - wall.radius_at(wall_pt.x));
        diag.wall_gap_over_spacing = diag.wall_gap / spacing;
        diag.wall_theta_mismatch = std::abs(wall_pt.theta - wall.theta_at(wall_pt.x));

        // How far the raw (pre-correction) series missed this wall angle. The correction
        // has already overwritten wall_pt.theta (which is why wall_theta_mismatch above is
        // ~0), so the value comes from the initializer, which measured it before correcting.
        if (diag.start_line_used == MocStartLine::KLIEGEL_LEVINE) {
            diag.wall_bc_residual = m_init_wall_bc_residual;
        }

        // Distance to the sharpest slope break on the contour. A conical profile's arc runs
        // into a straight cone, and seeding the wall march exactly there means the first
        // wall angle query is a one-sided difference across the break.
        double worst_jump = 0.0;
        double break_x = std::numeric_limits<double>::quiet_NaN();
        for (size_t i = 2; i + 1 < wall.size(); i++) {
            const double jump = std::abs(wall.slope_at_idx(i) - wall.slope_at_idx(i - 1));
            if (jump > worst_jump) {
                worst_jump = jump;
                break_x = wall.x[i - 1];
            }
        }
        // Below this the contour is smooth to within its own discretization and there is no
        // break to be near; report a sentinel rather than a meaningless distance.
        constexpr double slope_break_threshold = 1e-6;
        diag.wall_station_to_tangency = (worst_jump > slope_break_threshold)
            ? (wall_pt.x - break_x) / spacing
            : std::numeric_limits<double>::max();
    }
    else {
        diag.wall_station_to_tangency = std::numeric_limits<double>::max();
    }

    // The shift's natural scale. The transonic region's axial extent goes as
    // sqrt(r_throat * R_curvature), so the same absolute shift is a different fraction of
    // it at every throat curvature.
    if (m_options.geometry.downstream_wall_curvature_radius > 0.0 &&
        m_options.geometry.throat_radius > 0.0)
    {
        const double R = m_options.geometry.downstream_wall_curvature_radius
                       / m_options.geometry.throat_radius;
        diag.shift_over_transonic_length = m_options.initial_line_axial_shift / std::sqrt(R);
    }

    // Where each point's C- reaches the axis, on a straight-ray estimate, and how unevenly
    // graded those arrivals are. A centered fan's points share one location, so its rays
    // still land at distinct stations and the measure stays meaningful.
    std::vector<double> arrivals;
    arrivals.reserve(data_line.size());
    for (const CharacteristicPoint& pt : data_line) {
        const double slope = std::tan(pt.theta - pt.mu);
        arrivals.push_back((std::abs(slope) > 1e-12) ? pt.x + pt.y / std::abs(slope) : pt.x);
    }
    double min_gap = std::numeric_limits<double>::max();
    double max_gap = 0.0;
    for (size_t i = 1; i < arrivals.size(); i++) {
        const double gap = std::abs(arrivals[i] - arrivals[i - 1]);
        min_gap = std::min(min_gap, gap);
        max_gap = std::max(max_gap, gap);
    }
    diag.axis_arrival_grading = (min_gap > 0.0 && min_gap != std::numeric_limits<double>::max())
        ? max_gap / min_gap
        : std::numeric_limits<double>::max();

    // Spacelike margin of the line's own segments, defined exactly as in
    // record_front_diagnostics so the two are directly comparable. A centered fan is
    // degenerate here (all points share a location) and is skipped.
    double min_margin = std::numeric_limits<double>::max();
    for (size_t i = 1; i < data_line.size(); i++) {
        const CharacteristicPoint& below = data_line[i - 1];
        const CharacteristicPoint& above = data_line[i];
        const double dx = above.x - below.x;
        const double dy = above.y - below.y;
        if (dy <= 0.0) continue;
        const double margin_minus = (dy - std::tan(above.theta - above.mu) * dx) / dy;
        const double margin_plus = (dy - std::tan(below.theta + below.mu) * dx) / dy;
        min_margin = std::min({min_margin, margin_minus, margin_plus});
    }
    diag.min_spacelike_margin =
        (min_margin == std::numeric_limits<double>::max()) ? 0.0 : min_margin;

    diag.mass_flow_error = start_line_mass_flow_error(data_line);
    return diag;
}

size_t MocNozzle::control_front_spacing(
    CharacteristicNet& net,
    LeadingEdgeView& plus_edges,
    LeadingEdgeView& minus_edges,
    MocPassDiagnostics& diag)
{
    using Family = ChainMetadata::Family;

    // A minimum-length contour is *defined* by the C- characteristics it absorbs at the
    // wall, so inserting more would change the nozzle being designed rather than the mesh
    // resolving it. That mode also uses fan initialization and does not develop the void.
    if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
        record_front_diagnostics(net, build_front(net, plus_edges, minus_edges), 0.0, diag);
        return 0;
    }
    if (!(m_reference_spacing > 0.0)) return 0;

    std::vector<FrontSegment> front = build_front(net, plus_edges, minus_edges);
    if (front.empty()) {
        record_front_diagnostics(net, front, m_reference_spacing, diag);
        return 0;
    }

    // The front's own outer radius sets how much coarsening the nozzle's expansion
    // justifies. Measuring against a physical length rather than against the front's own
    // median spacing is the whole point: a threshold computed from the spacings it judges
    // is satisfied by any uniformly coarsening mesh, and so never sees a void whose width
    // is set by the flow.
    double front_radius = 0.0;
    for (const FrontSegment& seg : front) {
        front_radius = std::max(front_radius, net.points[seg.upper_pt_idx].y);
    }
    const double radius_ratio = std::max(1.0, front_radius / m_throat_radius);
    const double h_target =
        m_reference_spacing * (1.0 + m_options.front_spacing_growth * (radius_ratio - 1.0));
    const double h_max = m_options.max_front_spacing_factor * h_target;
    const double h_min = m_options.min_front_spacing_factor * h_target;

    const size_t front_cap = m_options.max_front_points != 0
        ? m_options.max_front_points
        : 4 * static_cast<size_t>(m_options.num_characteristics);

    // Bound on how finely one segment may be split in a single pass. A gap that has already
    // grown large is walked down over several passes instead of all at once, which keeps
    // each interpolation spanning a smaller range of the solution.
    constexpr size_t max_subdivisions_per_gap = 8;

    size_t front_points = front.size() + 1;
    size_t inserted = 0;

    for (const FrontSegment& seg : front) {
        if (seg.arc <= h_max) continue;
        if (front_points >= front_cap) {
            log_warning("Marching front reached its size cap ({}); refinement stopped.",
                front_cap);
            break;
        }

        // Copied by value, not bound by reference: insert_rung appends to net.points, which
        // reallocates and would leave a reference into the old buffer dangling.
        const CharacteristicPoint below = net.points[seg.lower_pt_idx];
        const CharacteristicPoint above = net.points[seg.upper_pt_idx];
        const double dx = above.x - below.x;
        const double dy = above.y - below.y;

        // Refuse to refine a front segment that is no longer spacelike -- one whose
        // endpoints the C+ from `below` and the C- from `above` no longer bracket. The
        // interpolation below assumes they do; once that fails the pair is already
        // degenerate and fabricating points on it would convert an honest
        // NON_DOWNSTREAM_POINT failure into silent garbage. Holding the spacing bounded
        // from the first pass should prevent a segment ever reaching this state, so if it
        // does, mesh control has already failed and that is worth saying out loud.
        if (dy - std::tan(below.theta + below.mu) * dx <= 0.0 ||
            dy - std::tan(above.theta - above.mu) * dx <= 0.0) {
            if (!m_warned_non_spacelike) {
                m_warned_non_spacelike = true;
                log_warning("Front segment between plus(x={:.6f},y={:.6f}) and "
                    "minus(x={:.6f},y={:.6f}) is already non-spacelike; mesh control cannot "
                    "refine it. The spacing bound engaged too late.",
                    below.x, below.y, above.x, above.y);
            }
            continue;
        }

        // Subdivide down to the *target*, not down to the trigger. Splitting only until the
        // segment falls under h_max leaves it permanently coarser than the rest of the
        // front, which against an engine that re-stretches it every pass is a losing race.
        size_t subdivisions = static_cast<size_t>(std::ceil(seg.arc / h_target));
        subdivisions = std::min(subdivisions, max_subdivisions_per_gap);
        subdivisions = std::min(subdivisions, front_cap - front_points + 1);
        if (subdivisions < 2) continue;

        const std::vector<size_t> previous_front = previous_front_points(net, seg);

        for (size_t k = 1; k < subdivisions; k++) {
            const double t = static_cast<double>(k) / static_cast<double>(subdivisions);

            // Cross-front interpolation, used only to seed the inverse solve below. On its
            // own it produces a point satisfying neither compatibility relation.
            CharacteristicPoint seed{};
            seed.theta = below.theta + t * (above.theta - below.theta);
            seed.nu = below.nu + t * (above.nu - below.nu);
            seed.mach = below.mach + t * (above.mach - below.mach);

            const double x_new = below.x + t * dx;
            const double y_new = below.y + t * dy;

            // A rung that cannot be given a valid state is simply not inserted: the pass
            // then takes the unrefined (oversized) step exactly as it would otherwise, so a
            // failed solve degrades refinement rather than corrupting the net.
            PointResult solved = solve_inverse_interior_point(
                x_new, y_new, net, previous_front, seed);
            if (solved.error != MocErrorCode::NONE) {
                log_debug("INSERT skipped at (x={:.6f},y={:.6f}): inverse solve failed ({})",
                    x_new, y_new, to_string(solved.error));
                continue;
            }
            CharacteristicPoint pt = solved.point;
            net.insert_rung(pt);
            inserted++;
            front_points++;
            log_debug("INSERT rung (x={:.6f},y={:.6f},th={:.6f},mach={:.6f}) splitting arc "
                "{:.6f} (target {:.6f}) between plus(x={:.6f},y={:.6f}) and "
                "minus(x={:.6f},y={:.6f})",
                pt.x, pt.y, pt.theta, pt.mach, seg.arc, h_target,
                below.x, below.y, above.x, above.y);
        }
    }

    if (inserted > 0) {
        m_inserted_characteristics += inserted;
        update_leading_edges(plus_edges, net, Family::PLUS);
        sort_plus_edges_by_proximity(plus_edges, net);
        update_leading_edges(minus_edges, net, Family::MINUS);
        front = build_front(net, plus_edges, minus_edges);
    }

    // Coarsening. Without it the front's rung count only ever grows -- a wall reflection
    // trades a C+ for a C- and the axis trades back, so nothing retires on its own -- and a
    // compression region (a Rao contour's turn-back above all) drives adjacent rungs
    // together until the cells are ill-conditioned.
    size_t retired = 0;
    if (front.size() > 2) {
        // Never thin the front faster than it can recover, and never retire two neighbours
        // in one pass: both endpoints of a short segment are usually short-spaced, and
        // taking both would leave a hole rather than a merged segment.
        const size_t retire_budget = front.size() / 3;
        bool previous_retired = false;

        for (size_t i = 0; i < front.size() && retired < retire_budget; i++) {
            // Two independent ways for a pairing to become degenerate. Crowding along the
            // front is the obvious one. The other is a cell whose C+ side has collapsed
            // while its diagonal still looks healthy -- the C+ has stalled and now sits
            // essentially on its partner's C- characteristic, so the two points are
            // redundant as C- samples even though they are not close together.
            const std::optional<std::pair<double, double>> sides = cell_sides(net, front[i]);
            const bool crowded = front[i].arc < h_min;
            const bool degenerate = sides.has_value() && sides->first > 0.0 &&
                sides->second / sides->first > m_options.max_cell_aspect_ratio;
            if (!crowded && !degenerate) { previous_retired = false; continue; }
            if (previous_retired) { previous_retired = false; continue; }

            // retire_rung declines any point that does not lead an active chain of each
            // family, which is exactly how the wall and axis anchors protect themselves: a
            // wall point's C+ and an axis point's C- are already terminated.
            const size_t candidate = front[i].upper_pt_idx;
            if (!net.retire_rung(candidate).has_value()) { previous_retired = false; continue; }

            retired++;
            previous_retired = true;
            log_debug("RETIRE rung (x={:.6f},y={:.6f}) {} at arc {:.6f} (min {:.6f}, "
                "aspect {:.3f})",
                net.points[candidate].x, net.points[candidate].y,
                crowded ? "crowding" : "degenerate cell", front[i].arc, h_min,
                sides.has_value() && sides->first > 0.0 ? sides->second / sides->first : 0.0);
        }
    }

    if (retired > 0) {
        m_retired_characteristics += retired;
        update_leading_edges(plus_edges, net, Family::PLUS);
        sort_plus_edges_by_proximity(plus_edges, net);
        update_leading_edges(minus_edges, net, Family::MINUS);
        front = build_front(net, plus_edges, minus_edges);
    }

    // Per-point snapshot of the front the pass will actually march. Aggregate statistics
    // cannot show a front smearing toward horizontal -- the cells stay well-formed while the
    // front itself stops being a sensible spacelike line -- so dump the geometry itself.
    if (m_options.log_level == MocLogLevel::DEBUG) {
        for (size_t i = 0; i < front.size(); i++) {
            const CharacteristicPoint& lo = net.points[front[i].lower_pt_idx];
            const CharacteristicPoint& hi = net.points[front[i].upper_pt_idx];
            log_debug("FRONTPT pass={} i={} lo_idx={} lo_x={:.6f} lo_y={:.6f} lo_th={:.6f} "
                "lo_mu={:.6f} hi_idx={} hi_x={:.6f} hi_y={:.6f} hi_th={:.6f} hi_mu={:.6f} "
                "pchain={} mchain={} arc={:.6f}",
                diag.pass, i, front[i].lower_pt_idx, lo.x, lo.y, lo.theta, lo.mu,
                front[i].upper_pt_idx, hi.x, hi.y, hi.theta, hi.mu,
                front[i].lower_plus_chain, front[i].upper_minus_chain, front[i].arc);
        }
    }

    record_front_diagnostics(net, front, h_target, diag);
    diag.inserted = inserted;
    diag.retired = retired;
    return inserted;
}

} // namespace Goddard