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

bool MocNozzle::is_solved() const {
    return m_is_solved;
}

// -- MocNozzle --
MocResult MocNozzle::solve() {
    // Re-checked here, not just in the constructor: m_options is public and callers do
    // adjust it between construction and solve().
    validate_moc_options(m_options);

    m_log = MocLog(m_options.log_level);
    m_theta_schedule.clear();
    m_initial_line_family.reset();
    m_is_solved = false;
    m_pass_diagnostics.clear();
    m_reference_spacing = 0.0;
    m_throat_radius = 1.0;
    m_init_wall_bc_residual = 0.0;

    CharacteristicNet net;
    std::vector<CharacteristicPoint> data_line;
    m_options.nozzle_profile = setup_nozzle_profile(m_options.geometry);

    ThroatCondition throat{};

    if (m_options.chemistry == GasChemistry::PERFECT_GAS) {
        // Pure algebraic path: no Cantera dependency
        // throat populated with dummy values; they are not used in the perfect gas case.
        throat = ThroatCondition{
            .converged = true,
            .speed_of_sound = 1.0,
            .H_stagnation = 1.0,
            .P_inlet = 1.0, // dimensionless stagnation pressure
            .S_inlet = 1.0,
            .gamma_s = m_options.gamma,
            .dlV_dlP_T = -1.0,
            .dlV_dlT_P = 1.0,
            .state = {}
        };
        m_thermo = MocThermo::perfect_gas(m_options.gamma);
    }
    else {
        // Cantera-backed path for frozen/equilibrium chemistry
        NozzleOptions nozzle_opts;
        nozzle_opts.chemistry = m_options.chemistry;
        Nozzle nozzle(*m_gas, nozzle_opts);
        throat = nozzle.solve_throat_conditions();
        m_thermo = MocThermo::tabulated(*m_gas, m_options.chemistry, throat);
    }

    // Everything a kernel or unit process needs for this solve, besides the points it works
    // on: m_options.nozzle_profile is still where the resolved contour lives at this step
    // (see setup_nozzle_profile above); a later step moves it out of MocOptions.
    m_context.emplace(MocSolveContext{m_options, m_options.nozzle_profile, *m_thermo, m_log});

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

        if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
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
            // CharacteristicNet::add_front); wall_x.front() is F_0's own wall point, not
            // (0, 1), so the throat radius for area_ratio/mesh purposes comes from the
            // geometry directly.
            std::vector<CharacteristicPoint> front0 = build_inverse_initial_front(data_line);
            if (front0.size() < 2) {
                throw ConvergenceError(
                    "Inverse march: initial front has fewer than two points.");
            }
            net.add_front(front0);
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
        result.messages = m_log.messages;
        m_is_solved = true;
        return result;
    }

    // Measured before the march so it is available even when the kernel fails partway --
    // an initialization defect is exactly the case where the march does not finish.
    result.init_diagnostics = record_init_diagnostics(data_line);

    std::optional<MocFailure> kernel_failure = (m_options.mode == MocMode::DESIGN_MIN_LENGTH)
        ? solve_characteristic_kernel(net)
        : solve_inverse_characteristic_kernel(net);

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
            m_log.warning("Flow angle reaches {:.2f} deg at (x={:.3f}, y={:.3f}): a compression is "
                "converging on the axis (a forming shock); the isentropic solution downstream "
                "of it is approximate.", result.min_theta * 180.0 / M_PI,
                result.min_theta_x, result.min_theta_y);
        }
    }

    result.net = net;
    result.messages = m_log.messages;
    if (kernel_failure.has_value()) {
        result.failure = *kernel_failure;
        result.converged = false;
    } else {
        result.converged = true;
    }

    // Populate performance fields
    result.crossings = find_like_characteristic_crossings(net);
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
    else {
        // The inverse march always lands its last front exactly on the exit plane (see
        // solve_inverse_characteristic_kernel's MocStepLimiter::EXIT step) -- there is no
        // ragged outflow staircase to fall short of. A march that failed partway did not
        // reach it; report how far it got instead.
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
    if (!net.wall_x.empty()) {
        result.nozzle_length = net.wall_x.back();
    }
    if (m_options.mode != MocMode::DESIGN_MIN_LENGTH) {
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
        if (m_options.mode != MocMode::DESIGN_MIN_LENGTH) {
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

        // Add the last characteristic for min length nozzle. Guarded by
        // emptiness: a kernel that now fails fast can abort before any axis/wall
        // point has been recorded, where leading_axis_point()/leading_wall_point()
        // would otherwise index an empty vector.
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
            m_log.info("Kliegel-Levine initialization is only valid for axisymmetric flow; "
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
    MocInitialization initializer(geometry, *m_thermo, init_options);

    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // The minimum length nozzle is a special case as the sonic line is straight, and
            // a centered expansion fan is the "exact" solution.
            // Therefore all initialization methods simplify to straight line initialization.

            // A centered expansion fan is collinear along a single C+ characteristic.
            m_initial_line_family = CharacteristicFamily::PLUS;

            auto expansion_line = initializer.initialize_centered_expansion(throat);
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
                        solve_fan_axis_point(*m_context, expansion_point),
                        "Initial data line (min-length axis point)");
                }
                else {
                    upstream_point = unwrap_or_throw(
                        solve_interior_point(*m_context, expansion_point, upstream_point),
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
                    data_line = initializer.initialize_kliegel_levine(throat);
                    m_init_wall_bc_residual = initializer.last_wall_bc_residual;
                }
                catch (const KlWallAngleFallback& e) {
                    // MocOptions::start_line was AUTO (a forced KLIEGEL_LEVINE throws
                    // ConvergenceError instead, which is not this type and is left to
                    // propagate to solve()'s exception boundary as INITIALIZATION_FAILED).
                    // Rebuild as the centered fan below, with theta_max derived exactly as
                    // the eager fan path derives it.
                    m_log.info("Kliegel-Levine series misses the wall angle by {:.6f} rad at "
                             "R = {:.6f}; using centered-fan initialization.", e.miss, e.R);
                    resolved_start_line = MocStartLine::CENTERED_FAN;
                    derive_fan_theta_max();
                    initializer = MocInitialization(geometry, *m_thermo, init_options);
                }
            }
            if (resolved_start_line == MocStartLine::CENTERED_FAN) {
                m_initial_line_family = CharacteristicFamily::PLUS;
                auto expansion_line = initializer.initialize_centered_expansion(throat);
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
                            solve_fan_axis_point(*m_context, expansion_point),
                            "Initial data line (centered-fan axis point)");
                        data_line.push_back(upstream_point);
                    }
                    else {
                        upstream_point = unwrap_or_throw(
                            solve_interior_point(*m_context, expansion_point, upstream_point),
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
        m_log.debug("INIT[{}] x={:.6f} y={:.6f} theta={:.6f} nu={:.6f} mach={:.6f} mu={:.6f}",
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
    using Family = CharacteristicFamily;
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
        m_log.debug("--- kernel pass {}: {} active C+, {} active C- ---",
            iters, plus_edges.chain_indices.size(), minus_edges.chain_indices.size());
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

                PointResult result = solve_interior_point(*m_context, minus_pt, plus_pt);
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
                m_log.debug("PAIR minus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) "
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
                m_log.debug("SKIP plus(x={:.6f},y={:.6f}) waiting for scarce C- partner "
                    "(already claimed this pass)", plus_pt.x, plus_pt.y);
            }
            else [[unlikely]] { // C+ intersects with wall

                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                PointResult result = solve_wall_point(
                    plus_pt,
                    net.leading_wall_point(),
                    static_cast<int>(net.wall_point_indices.size()) - 1
                );

                if (result.error != MocErrorCode::NONE) {
                    failure = MocFailure{
                        result.error,
                        std::format("Wall point solve failed ({}) for plus(x={:.6f},y={:.6f}).",
                            to_string(result.error), plus_pt.x, plus_pt.y),
                        result.point.x, result.point.y,
                        iters
                    };
                    break;
                }
                m_log.debug("WALL plus(x={:.6f},y={:.6f}) -> wall hit at (x={:.6f},y={:.6f})",
                    plus_pt.x, plus_pt.y, result.point.x, result.point.y);
                intersections.push_back({
                    result.point,
                    PointMembership {
                        .c_plus_chain_idx = plus_edges.chain_indices[i],
                        .c_minus_chain_idx = std::nullopt
                    }
                });
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
            PointResult axis_result = solve_axis_point(*m_context, minus_pt);
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
            m_log.debug("AXIS minus(x={:.6f},y={:.6f}) reflects off axis", minus_pt.x, minus_pt.y);
            intersections.push_back({
                axis_result.point,
                PointMembership {
                    .c_plus_chain_idx = std::nullopt,
                    .c_minus_chain_idx = min_y_cminus_idx
                }
            });
        }
        // insert all intersections into net. Minimum-length design only: the wall
        // absorbs each C+ (no reflected wave), which is exactly what produces a
        // minimum-length contour.
        for (auto&& [pt, mem] : intersections) {
            if (!mem.c_minus_chain_idx.has_value()) {
                // C+ reaches the wall: absorb it (no reflected wave). terminate_c_plus_at_wall
                // pushes wall_x/wall_y itself.
                net.terminate_c_plus_at_wall(*mem.c_plus_chain_idx, pt);
            }
            else if (!mem.c_plus_chain_idx.has_value()) {
                // C- reaches the axis: reflect it into a new upward-marching C+.
                net.reflect_c_minus_off_axis(*mem.c_minus_chain_idx, pt);
            }
            else {
                net.add_point(pt, mem);
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

PointResult MocNozzle::solve_wall_point(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    int wall_point_index)
{
    // Minimum length design only: theta is determined by theta schedule.
    // Guard against an empty schedule (the Cantera min-length path does not
    // populate m_theta_schedule) and against indices beyond it (the kernel
    // currently produces more wall hits than scheduled characteristics).
    double theta_wall = m_options.theta_max;
    if (!m_theta_schedule.empty()) {
        size_t k = std::min(static_cast<size_t>(wall_point_index),
                            m_theta_schedule.size() - 1);
        theta_wall = m_options.theta_max - m_theta_schedule[k];
    }

    return solve_wall_point_design(*m_context, interior_parent, previous_wall_point, theta_wall);
}

MocNozzle::LeadingEdgeView MocNozzle::leading_edges(const CharacteristicNet& net, CharacteristicFamily fam) const {
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

void MocNozzle::update_leading_edges(LeadingEdgeView& view, const CharacteristicNet& net, CharacteristicFamily family) const {
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

    // Spacelike margin of the line's own segments: how far a segment sits from being
    // parallel to a characteristic through one of its endpoints, as a fraction of its
    // height. A centered fan is degenerate here (all points share a location) and is
    // skipped.
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

} // namespace Goddard