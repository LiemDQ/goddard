// MoC convergence sweep: one CSV row per configuration.
//
// The companion to instructions/moc_convergence_roadmap.md. That document's elimination
// scoreboard exists because the same measurements kept being re-derived by hand from
// one-off probes; this tool makes the evidence tables reproducible with one command.
//
// Emits every MocInitDiagnostics field, the front-shear summary, and the convergence
// metrics the roadmap says to judge on (exit_coverage and exit_mach -- never area_ratio,
// which mesh control moves by 4% while exit_mach is bit-identical).
//
// Build:  cmake -Dgoddard_BUILD_TOOLS=ON ... && cmake --build build --target moc_sweep
// Run:    pixi run moc-sweep --out sweep.csv

#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/profile.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace Goddard;

namespace {

struct Config {
    std::string mode;      // "analysis" | "rao" | "minlength"
    std::string init;      // "kl" | "fan"
    std::string flow;      // "axi" | "planar"
    double area_ratio = 4.0;
    int n = 15;
    double r_arc = 0.382;  // throat expansion-curve radius AND the KL curvature ratio
    double shift = 0.1;
    double length_frac = 0.8;
    bool mesh_control = false;
    double gamma = 1.4;
    double half_angle = 15.0;
};

MocOptions build_options(const Config& c) {
    MocOptions o;
    o.chemistry = GasChemistry::PERFECT_GAS;
    o.gamma = c.gamma;
    o.num_characteristics = c.n;
    o.geometry.throat_radius = 1.0;
    o.flow_type = (c.flow == "planar") ? MocFlowKind::PLANAR : MocFlowKind::AXISYMMETRIC;
    o.initial_line_axial_shift = c.shift;

    // The KL series' curvature ratio and the contour's throat arc are separate inputs that
    // nothing in the library couples. Setting them together is the only way the start line
    // and the wall it is seeded against describe the same throat.
    o.geometry.downstream_wall_curvature_radius =
        (c.init == "fan") ? -1.0 : c.r_arc;

    if (!c.mesh_control) {
        // There is no disable switch; make every bound non-binding.
        o.max_front_spacing_factor = 1e9;
        o.min_front_spacing_factor = 1e-9;
        o.max_cell_aspect_ratio = 1e9;
    }

    if (c.mode == "minlength") {
        o.mode = MocMode::DESIGN_MIN_LENGTH;
        o.theta_max = c.half_angle * M_PI / 180.0;
        return o;
    }
    if (c.mode == "rao") {
        o.mode = MocMode::DESIGN_RAO;
        o.geometry.expansion_ratio = c.area_ratio;
        o.geometry.length_fraction = c.length_frac;
        return o;
    }

    o.mode = MocMode::ANALYSIS;
    // generate_conical_nozzle sets r_exit = sqrt(arg), so a planar half-height of
    // `area_ratio` needs the argument squared.
    const double profile_arg =
        (c.flow == "planar") ? c.area_ratio * c.area_ratio : c.area_ratio;
    o.nozzle_profile =
        NozzleProfile::generate_conical_nozzle(profile_arg, c.r_arc, 1.0, c.half_angle, 60);
    return o;
}

void write_header(std::ostream& os) {
    os << "mode,init,flow,area_ratio,n,r_arc,shift,length_frac,mesh_control,gamma,"
          "converged,failure_code,fail_x,fail_y,kernel_pass,"
          "exit_mach,exit_coverage,reached_exit_plane,area_ratio_achieved,nozzle_length,"
          "points,chains,inserted,retired,crossings,passes,"
          "init_points,wall_gap,wall_gap_over_spacing,wall_theta_mismatch,"
          "mach_axis,mach_wall,mach_ratio,mu_axis,mu_wall,cot_mu_ratio,"
          "axis_arrival_grading,init_min_spacelike_margin,mass_flow_error,"
          "shift_over_transonic_length,wall_station_to_tangency,"
          "shear_valid,axis_growth,wall_growth,shear_ratio\n";
}

// %.10g keeps the columns readable while staying well inside the precision any of these
// comparisons need. Infinities and NaNs are emitted as bare tokens rather than platform
// spellings so the reporter can parse them uniformly.
std::string num(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%.10g", v);
    return buf;
}

void write_row(std::ostream& os, const Config& c, const MocResult& r) {
    const MocInitDiagnostics& d = r.init_diagnostics;
    const MocFrontShear shear = summarize_front_shear(r.pass_diagnostics);

    os << c.mode << ',' << c.init << ',' << c.flow << ',' << num(c.area_ratio) << ','
       << c.n << ',' << num(c.r_arc) << ',' << num(c.shift) << ',' << num(c.length_frac)
       << ',' << (c.mesh_control ? 1 : 0) << ',' << num(c.gamma) << ','
       << (r.converged ? 1 : 0) << ',' << to_string(r.failure.code) << ','
       << num(r.failure.x) << ',' << num(r.failure.y) << ',' << r.failure.kernel_pass << ','
       << num(r.exit_mach) << ',' << num(r.exit_coverage) << ','
       << (r.reached_exit_plane ? 1 : 0) << ',' << num(r.area_ratio) << ','
       << num(r.nozzle_length) << ','
       << r.net.points.size() << ',' << r.net.c_chains.size() << ','
       << r.inserted_characteristics << ',' << r.retired_characteristics << ','
       << r.crossings.count << ',' << r.pass_diagnostics.size() << ','
       << d.points << ',' << num(d.wall_gap) << ',' << num(d.wall_gap_over_spacing) << ','
       << num(d.wall_theta_mismatch) << ','
       << num(d.mach_axis) << ',' << num(d.mach_wall) << ',' << num(d.mach_ratio) << ','
       << num(d.mu_axis) << ',' << num(d.mu_wall) << ',' << num(d.cot_mu_ratio) << ','
       << num(d.axis_arrival_grading) << ',' << num(d.min_spacelike_margin) << ','
       << num(d.mass_flow_error) << ','
       << num(d.shift_over_transonic_length) << ',' << num(d.wall_station_to_tangency) << ','
       << (shear.valid ? 1 : 0) << ',' << num(shear.axis_growth) << ','
       << num(shear.wall_growth) << ',' << num(shear.shear_ratio) << '\n';
}

std::vector<Config> default_grid() {
    std::vector<Config> grid;

    // Conical analysis: the sweep the roadmap's Sec 5.1 table is built from. r_arc and
    // shift move together because the arc/cone tangency sits at r_arc*sin(theta_n)
    // (0.0989 / 0.259 / 0.518 at theta_n = 15 deg) while the default shift is 0.1 --
    // varying r_arc alone also moves the start line's wall point on and off that corner.
    for (const char* init : {"kl", "fan"}) {
        for (double ar : {2.0, 4.0, 8.0}) {
            for (int n : {8, 15, 31, 61}) {
                for (double r_arc : {0.382, 1.0, 2.0}) {
                    for (double shift : {0.05, 0.1, 0.2}) {
                        Config c;
                        c.mode = "analysis";
                        c.init = init;
                        c.flow = "axi";
                        c.area_ratio = ar;
                        c.n = n;
                        c.r_arc = r_arc;
                        c.shift = shift;
                        grid.push_back(c);
                    }
                }
            }
        }
    }

    // Planar control: it converges, so it is the reference for every axisymmetric row.
    for (double ar : {2.0, 4.0}) {
        for (int n : {8, 15, 31, 61}) {
            Config c;
            c.mode = "analysis";
            c.init = "fan";
            c.flow = "planar";
            c.area_ratio = ar;
            c.n = n;
            grid.push_back(c);
        }
    }

    for (double ar : {5.0, 10.0, 20.0}) {
        for (double lf : {0.6, 0.8, 1.0}) {
            Config c;
            c.mode = "rao";
            c.init = "kl";
            c.flow = "axi";
            c.area_ratio = ar;
            c.n = 15;
            c.length_frac = lf;
            grid.push_back(c);
        }
    }

    // Invariance control: min-length is mode-excluded from mesh control and always uses
    // the fan, so these rows must not move at all across this work.
    for (const char* flow : {"planar", "axi"}) {
        for (int n : {15, 31}) {
            Config c;
            c.mode = "minlength";
            c.init = "fan";
            c.flow = flow;
            c.n = n;
            c.half_angle = 12.0;
            grid.push_back(c);
        }
    }
    return grid;
}

} // namespace

int main(int argc, char** argv) {
    std::string out_path;
    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc) out_path = argv[++i];
        else if (std::strcmp(argv[i], "--help") == 0) {
            std::cout << "usage: moc_sweep [--out FILE]\n"
                         "Writes the MoC convergence sweep as CSV (stdout if --out is absent).\n";
            return 0;
        }
    }

    std::ofstream file;
    if (!out_path.empty()) {
        file.open(out_path);
        if (!file) {
            std::cerr << "moc_sweep: cannot open " << out_path << " for writing\n";
            return 2;
        }
    }
    std::ostream& os = out_path.empty() ? std::cout : file;

    write_header(os);
    const std::vector<Config> grid = default_grid();
    size_t solved = 0, threw = 0;
    for (const Config& c : grid) {
        MocResult r;
        try {
            MocNozzle solver(build_options(c));
            r = solver.solve();
        }
        catch (const std::exception& e) {
            // solve() is contracted not to throw for numerical failures, so anything here
            // is a genuine defect worth seeing rather than a row to silently drop.
            std::cerr << "moc_sweep: " << c.mode << '/' << c.init << " AR=" << c.area_ratio
                      << " N=" << c.n << " r_arc=" << c.r_arc << " threw: " << e.what() << '\n';
            threw++;
            continue;
        }
        write_row(os, c, r);
        solved++;
    }
    os.flush();
    std::cerr << "moc_sweep: " << solved << " rows written, " << threw << " threw\n";
    return 0;
}
