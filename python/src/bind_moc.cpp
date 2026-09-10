#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/optional.h>
#include <nanobind/eigen/dense.h>

#include <algorithm>
#include <cstdint>
#include <format>
#include <string>
#include <vector>

#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/characteristics.hpp"
#include "goddard/characteristic_net.hpp"
#include "goddard/profile.hpp"
#include "goddard/gas.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

namespace {

using IndexArray = Eigen::Matrix<int64_t, Eigen::Dynamic, 1>;

// Bulk numeric data crosses into Python as numpy arrays, copied rather than viewed:
// NozzleProfile::push_back reallocates its backing std::vector, so a zero-copy view would be
// left dangling by an ordinary mutation. Callers that iterate should snapshot once and reuse.
Eigen::VectorXd to_array(const std::vector<double>& values) {
    Eigen::VectorXd out(static_cast<Eigen::Index>(values.size()));
    std::copy(values.begin(), values.end(), out.data());
    return out;
}

IndexArray to_index_array(const std::vector<size_t>& indices) {
    IndexArray out(static_cast<Eigen::Index>(indices.size()));
    for (size_t i = 0; i < indices.size(); i++) {
        out[static_cast<Eigen::Index>(i)] = static_cast<int64_t>(indices[i]);
    }
    return out;
}

// Columnar extraction of one CharacteristicPoint field across a whole net, so the field can be
// handed straight to numpy/matplotlib without a per-point Python loop.
Eigen::VectorXd gather(
    const std::vector<Goddard::CharacteristicPoint>& points,
    double Goddard::CharacteristicPoint::* field)
{
    Eigen::VectorXd out(static_cast<Eigen::Index>(points.size()));
    for (size_t i = 0; i < points.size(); i++) {
        out[static_cast<Eigen::Index>(i)] = points[i].*field;
    }
    return out;
}

std::string_view flow_kind_name(Goddard::MocFlowKind kind) {
    switch (kind) {
        case Goddard::MocFlowKind::PLANAR: return "PLANAR";
        case Goddard::MocFlowKind::AXISYMMETRIC: return "AXISYMMETRIC";
    }
    return "UNKNOWN";
}

std::string_view mode_name(Goddard::MocMode mode) {
    switch (mode) {
        case Goddard::MocMode::DESIGN_MIN_LENGTH: return "DESIGN_MIN_LENGTH";
        case Goddard::MocMode::DESIGN_RAO: return "DESIGN_RAO";
        case Goddard::MocMode::DESIGN_CENTERLINE: return "DESIGN_CENTERLINE";
        case Goddard::MocMode::ANALYSIS: return "ANALYSIS";
    }
    return "UNKNOWN";
}

} // namespace

void bind_moc(nb::module_& m) {
    // ---- Enums ----

    nb::enum_<Goddard::MocFlowKind>(m, "MocFlowKind", DOC(Goddard, MocFlowKind))
        .value("PLANAR", Goddard::MocFlowKind::PLANAR)
        .value("AXISYMMETRIC", Goddard::MocFlowKind::AXISYMMETRIC);

    nb::enum_<Goddard::MocMode>(m, "MocMode", DOC(Goddard, MocMode))
        .value("DESIGN_MIN_LENGTH", Goddard::MocMode::DESIGN_MIN_LENGTH)
        .value("DESIGN_RAO", Goddard::MocMode::DESIGN_RAO)
        .value("DESIGN_CENTERLINE", Goddard::MocMode::DESIGN_CENTERLINE)
        .value("ANALYSIS", Goddard::MocMode::ANALYSIS);

    nb::enum_<Goddard::MocLogLevel>(m, "MocLogLevel", DOC(Goddard, MocLogLevel))
        .value("NORMAL", Goddard::MocLogLevel::NORMAL)
        .value("DEBUG", Goddard::MocLogLevel::DEBUG);

    // Nested in ChainMetadata on the C++ side; flattened to module scope here, since nesting
    // them under the Python class would make every reference a mouthful.
    //
    // Literal docstrings rather than DOC(): both enums are declared inline as part of the
    // member they type (`enum class Family {...} family;`), so it is ambiguous whether
    // pybind11_mkdoc attaches the preceding comment to the type or to the member.
    nb::enum_<Goddard::ChainMetadata::Family>(m, "CharacteristicFamily",
                                              "Which family a characteristic belongs to: "
                                              "C+ (PLUS) or C- (MINUS).")
        .value("UNSPECIFIED", Goddard::ChainMetadata::Family::UNSPECIFIED)
        .value("PLUS", Goddard::ChainMetadata::Family::PLUS)
        .value("MINUS", Goddard::ChainMetadata::Family::MINUS);

    nb::enum_<Goddard::ChainMetadata::TerminationType>(m, "ChainTermination",
                                                       "How a characteristic chain stopped "
                                                       "being marched.")
        .value("NOT_TERMINATED", Goddard::ChainMetadata::TerminationType::NOT_TERMINATED)
        .value("WALL", Goddard::ChainMetadata::TerminationType::WALL)
        .value("AXIS", Goddard::ChainMetadata::TerminationType::AXIS)
        .value("OUTFLOW", Goddard::ChainMetadata::TerminationType::OUTFLOW)
        .value("CORNER_FAN_ORIGIN", Goddard::ChainMetadata::TerminationType::CORNER_FAN_ORIGIN)
        .value("MERGED", Goddard::ChainMetadata::TerminationType::MERGED);

    // MocErrorCode is bound in bind_enums.cpp.

    // ---- NozzleGeometry ----

    // No nb::init<>() overload on purpose: NozzleGeometry has members without default
    // initializers, so a bare value-initializing constructor would hand back throat_radius = 0.
    // Routing every construction through the keyword constructor below keeps the documented
    // defaults as the only defaults.
    nb::class_<Goddard::NozzleGeometry>(m, "NozzleGeometry", DOC(Goddard, NozzleGeometry))
        .def("__init__", [](Goddard::NozzleGeometry* self,
                            double throat_radius,
                            double upstream_wall_curvature_radius,
                            double downstream_wall_curvature_radius,
                            double length_fraction,
                            double expansion_ratio) {
            new (self) Goddard::NozzleGeometry();
            self->throat_radius = throat_radius;
            self->upstream_wall_curvature_radius = upstream_wall_curvature_radius;
            self->downstream_wall_curvature_radius = downstream_wall_curvature_radius;
            self->length_fraction = length_fraction;
            self->expansion_ratio = expansion_ratio;
        },  "throat_radius"_a = 1.0,
            "upstream_wall_curvature_radius"_a = 1.5,
            "downstream_wall_curvature_radius"_a = 0.382,
            "length_fraction"_a = 0.8,
            "expansion_ratio"_a = 5.0)
        .def_rw("throat_radius", &Goddard::NozzleGeometry::throat_radius,
                DOC(Goddard, NozzleGeometry, throat_radius))
        .def_rw("upstream_wall_curvature_radius", &Goddard::NozzleGeometry::upstream_wall_curvature_radius,
                DOC(Goddard, NozzleGeometry, upstream_wall_curvature_radius))
        .def_rw("downstream_wall_curvature_radius", &Goddard::NozzleGeometry::downstream_wall_curvature_radius,
                DOC(Goddard, NozzleGeometry, downstream_wall_curvature_radius))
        .def_rw("length_fraction", &Goddard::NozzleGeometry::length_fraction,
                DOC(Goddard, NozzleGeometry, length_fraction))
        .def_rw("expansion_ratio", &Goddard::NozzleGeometry::expansion_ratio,
                DOC(Goddard, NozzleGeometry, expansion_ratio))
        .def("__repr__", [](const Goddard::NozzleGeometry& self) {
            return std::format(
                "<NozzleGeometry throat_radius={:g} r_up={:g} r_down={:g} "
                "length_fraction={:g} expansion_ratio={:g}>",
                self.throat_radius, self.upstream_wall_curvature_radius,
                self.downstream_wall_curvature_radius, self.length_fraction,
                self.expansion_ratio);
        });

    // ---- NozzleProfile ----

    nb::class_<Goddard::NozzleProfile>(m, "NozzleProfile", DOC(Goddard, NozzleProfile))
        .def(nb::init<>())
        .def("__init__", [](Goddard::NozzleProfile* self,
                            std::vector<double> x,
                            std::vector<double> y,
                            size_t throat_index) {
            if (x.size() != y.size()) {
                throw std::invalid_argument("x and y must have the same length.");
            }
            new (self) Goddard::NozzleProfile();
            self->x = std::move(x);
            self->y = std::move(y);
            self->throat_index = throat_index;
        }, "x"_a, "y"_a, "throat_index"_a = 0)

        // ---- Coordinates ----

        .def_prop_rw("x",
            [](const Goddard::NozzleProfile& self) { return to_array(self.x); },
            [](Goddard::NozzleProfile& self, std::vector<double> values) { self.x = std::move(values); },
            DOC(Goddard, NozzleProfile, x))
        .def_prop_rw("y",
            [](const Goddard::NozzleProfile& self) { return to_array(self.y); },
            [](Goddard::NozzleProfile& self, std::vector<double> values) { self.y = std::move(values); },
            DOC(Goddard, NozzleProfile, y))
        .def_rw("throat_index", &Goddard::NozzleProfile::throat_index,
                DOC(Goddard, NozzleProfile, throat_index))

        // ---- Geometric queries ----

        .def("slope_at", &Goddard::NozzleProfile::slope_at, "x_query"_a,
             DOC(Goddard, NozzleProfile, slope_at))
        .def("theta_at", &Goddard::NozzleProfile::theta_at, "x_query"_a,
             DOC(Goddard, NozzleProfile, theta_at))
        .def("radius_at", &Goddard::NozzleProfile::radius_at, "x_query"_a,
             DOC(Goddard, NozzleProfile, radius_at))
        .def("area_at", &Goddard::NozzleProfile::area_at, "x_query"_a,
             DOC(Goddard, NozzleProfile, area_at))
        .def("slope_at_idx", &Goddard::NozzleProfile::slope_at_idx, "idx"_a,
             DOC(Goddard, NozzleProfile, slope_at_idx))
        .def("theta_at_idx", &Goddard::NozzleProfile::theta_at_idx, "idx"_a,
             DOC(Goddard, NozzleProfile, theta_at_idx))
        .def("max_theta", &Goddard::NozzleProfile::max_theta,
             DOC(Goddard, NozzleProfile, max_theta))
        .def("x_min", &Goddard::NozzleProfile::x_min, DOC(Goddard, NozzleProfile, x_min))
        .def("x_max", &Goddard::NozzleProfile::x_max, DOC(Goddard, NozzleProfile, x_max))
        .def("length", &Goddard::NozzleProfile::length, DOC(Goddard, NozzleProfile, length))
        .def("radius_max", &Goddard::NozzleProfile::radius_max,
             DOC(Goddard, NozzleProfile, radius_max))

        // ---- Container protocol ----

        .def("at", &Goddard::NozzleProfile::at, "idx"_a, DOC(Goddard, NozzleProfile, at))
        .def("size", &Goddard::NozzleProfile::size, DOC(Goddard, NozzleProfile, size))
        .def("push_back", [](Goddard::NozzleProfile& self, double x, double y) {
            self.push_back({x, y});
        }, "x"_a, "y"_a, DOC(Goddard, NozzleProfile, push_back))
        .def("__len__", &Goddard::NozzleProfile::size)
        // Defining __getitem__ is also what makes the profile iterable: Python falls back to
        // the sequence protocol, walking indices until IndexError.
        .def("__getitem__", [](const Goddard::NozzleProfile& self, Py_ssize_t index) {
            Py_ssize_t count = static_cast<Py_ssize_t>(self.size());
            if (index < 0) {
                index += count;
            }
            if (index < 0 || index >= count) {
                throw nb::index_error("NozzleProfile index out of range");
            }
            return self.at(static_cast<size_t>(index));
        }, "index"_a)
        .def("__repr__", [](const Goddard::NozzleProfile& self) {
            if (self.size() == 0) {
                return std::string("<NozzleProfile empty>");
            }
            return std::format(
                "<NozzleProfile {} points, x=[{:.4g}, {:.4g}], r=[{:.4g}, {:.4g}]>",
                self.size(), self.x_min(), self.x_max(),
                self.y.front(), self.radius_max().second);
        })

        // ---- Contour generators ----
        //
        // Shape angles are in degrees here, as the C++ generators define them; every other
        // angle in the MoC API is in radians. The _deg suffix is there to keep that visible at
        // the call site.

        .def_static("conical", &Goddard::NozzleProfile::generate_conical_nozzle,
                    "area_ratio"_a, "r_expansion_curve"_a, "r_throat"_a = 1.0,
                    "angle_deg"_a = 15.0, "n_points"_a = 50,
                    DOC(Goddard, NozzleProfile, generate_conical_nozzle))
        .def_static("rao_top", &Goddard::NozzleProfile::generate_Rao_TOP_nozzle,
                    "area_ratio"_a, "r_throat"_a = 1.0, "length_frac"_a = 0.8,
                    "n_points"_a = 50,
                    DOC(Goddard, NozzleProfile, generate_Rao_TOP_nozzle))
        .def_static("bezier", &Goddard::NozzleProfile::generate_bezier_nozzle,
                    "area_ratio"_a, "theta_n_deg"_a, "theta_e_deg"_a,
                    "r_expansion_curve"_a = 0.382, "r_throat"_a = 1.0,
                    "length_frac"_a = 0.8, "n_points"_a = 50,
                    DOC(Goddard, NozzleProfile, generate_bezier_nozzle))
        .def_static("throat_expansion_curve",
                    &Goddard::NozzleProfile::generate_throat_expansion_curve,
                    "theta_n_deg"_a, "r_expansion_curve"_a, "r_throat"_a, "n_points"_a = 50,
                    DOC(Goddard, NozzleProfile, generate_throat_expansion_curve))

        // ---- I/O ----

        .def_static("load_csv", &Goddard::NozzleProfile::load_profile_csv, "filename"_a,
                    DOC(Goddard, NozzleProfile, load_profile_csv))
        .def("save_csv", &Goddard::NozzleProfile::save_profile_csv, "filename"_a,
             DOC(Goddard, NozzleProfile, save_profile_csv));

    // ---- MocOptions ----

    // Every field is listed three times below: as a constructor parameter, as an assignment,
    // and as a default. That redundancy is deliberate -- nanobind.stubgen turns the named
    // parameters into the signature that _core.pyi and the docs site show, which a **kwargs
    // constructor would erase. Adding a field to MocOptions means editing all three places.
    nb::class_<Goddard::MocOptions>(m, "MocOptions", DOC(Goddard, MocOptions))
        .def("__init__", [](Goddard::MocOptions* self,
                            Goddard::MocFlowKind flow_type,
                            Goddard::GasChemistry chemistry,
                            Goddard::MocMode mode,
                            int num_characteristics,
                            double gamma,
                            Goddard::SolverOptions solver_options,
                            Goddard::NozzleGeometry geometry,
                            double theta_max,
                            double exit_mach,
                            std::vector<double> theta_schedule,
                            Goddard::NozzleProfile nozzle_profile,
                            Goddard::MocLogLevel log_level,
                            double initial_line_axial_shift,
                            double max_front_spacing_factor,
                            double min_front_spacing_factor,
                            double front_spacing_growth,
                            double max_cell_aspect_ratio,
                            size_t max_front_points,
                            Goddard::MocStartLine start_line,
                            double kl_max_wall_angle_error,
                            Goddard::MocMarchScheme march_scheme,
                            double inverse_cfl,
                            double max_wall_turn_per_step,
                            double front_tilt_decay) {
            new (self) Goddard::MocOptions();
            self->flow_type = flow_type;
            self->chemistry = chemistry;
            self->mode = mode;
            self->num_characteristics = num_characteristics;
            self->gamma = gamma;
            self->solver_options = solver_options;
            self->geometry = std::move(geometry);
            self->theta_max = theta_max;
            self->exit_mach = exit_mach;
            self->theta_schedule = std::move(theta_schedule);
            self->nozzle_profile = std::move(nozzle_profile);
            self->log_level = log_level;
            self->initial_line_axial_shift = initial_line_axial_shift;
            self->max_front_spacing_factor = max_front_spacing_factor;
            self->min_front_spacing_factor = min_front_spacing_factor;
            self->front_spacing_growth = front_spacing_growth;
            self->max_cell_aspect_ratio = max_cell_aspect_ratio;
            self->max_front_points = max_front_points;
            self->start_line = start_line;
            self->kl_max_wall_angle_error = kl_max_wall_angle_error;
            self->march_scheme = march_scheme;
            self->inverse_cfl = inverse_cfl;
            self->max_wall_turn_per_step = max_wall_turn_per_step;
            self->front_tilt_decay = front_tilt_decay;
        },  "flow_type"_a = Goddard::MocFlowKind::PLANAR,
            "chemistry"_a = Goddard::GasChemistry::PERFECT_GAS,
            "mode"_a = Goddard::MocMode::DESIGN_MIN_LENGTH,
            "num_characteristics"_a = 10,
            "gamma"_a = 1.4,
            "solver_options"_a = Goddard::SolverOptions{.abstol = 1e-10, .reltol = 1e-5},
            "geometry"_a = Goddard::NozzleGeometry{1.0, 1.5, 0.382, 0.8, 5.0},
            "theta_max"_a = 0.0,
            "exit_mach"_a = 0.0,
            "theta_schedule"_a = std::vector<double>(),
            "nozzle_profile"_a = Goddard::NozzleProfile(),
            "log_level"_a = Goddard::MocLogLevel::NORMAL,
            "initial_line_axial_shift"_a = 0.1,
            "max_front_spacing_factor"_a = 1.5,
            "min_front_spacing_factor"_a = 0.35,
            "front_spacing_growth"_a = 1.0,
            "max_cell_aspect_ratio"_a = 6.0,
            "max_front_points"_a = 0,
            "start_line"_a = Goddard::MocStartLine::AUTO,
            "kl_max_wall_angle_error"_a = 0.25,
            "march_scheme"_a = Goddard::MocMarchScheme::AUTO,
            "inverse_cfl"_a = 0.8,
            "max_wall_turn_per_step"_a = 0.0175,
            "front_tilt_decay"_a = 0.9)
        .def_rw("flow_type", &Goddard::MocOptions::flow_type)
        .def_rw("chemistry", &Goddard::MocOptions::chemistry)
        .def_rw("mode", &Goddard::MocOptions::mode)
        .def_rw("num_characteristics", &Goddard::MocOptions::num_characteristics,
                DOC(Goddard, MocOptions, num_characteristics))
        .def_rw("gamma", &Goddard::MocOptions::gamma, DOC(Goddard, MocOptions, gamma))
        .def_rw("solver_options", &Goddard::MocOptions::solver_options,
                DOC(Goddard, MocOptions, solver_options))
        .def_rw("geometry", &Goddard::MocOptions::geometry, DOC(Goddard, MocOptions, geometry))
        .def_rw("theta_max", &Goddard::MocOptions::theta_max, DOC(Goddard, MocOptions, theta_max))
        .def_rw("exit_mach", &Goddard::MocOptions::exit_mach, DOC(Goddard, MocOptions, exit_mach))
        .def_prop_rw("theta_schedule",
            [](const Goddard::MocOptions& self) { return to_array(self.theta_schedule); },
            [](Goddard::MocOptions& self, std::vector<double> values) {
                self.theta_schedule = std::move(values);
            },
            DOC(Goddard, MocOptions, theta_schedule))
        .def_rw("nozzle_profile", &Goddard::MocOptions::nozzle_profile,
                DOC(Goddard, MocOptions, nozzle_profile))
        .def_rw("log_level", &Goddard::MocOptions::log_level, DOC(Goddard, MocOptions, log_level))
        .def_rw("initial_line_axial_shift", &Goddard::MocOptions::initial_line_axial_shift,
                DOC(Goddard, MocOptions, initial_line_axial_shift))
        .def_rw("max_front_spacing_factor", &Goddard::MocOptions::max_front_spacing_factor,
                DOC(Goddard, MocOptions, max_front_spacing_factor))
        .def_rw("min_front_spacing_factor", &Goddard::MocOptions::min_front_spacing_factor,
                DOC(Goddard, MocOptions, min_front_spacing_factor))
        .def_rw("front_spacing_growth", &Goddard::MocOptions::front_spacing_growth,
                DOC(Goddard, MocOptions, front_spacing_growth))
        .def_rw("max_cell_aspect_ratio", &Goddard::MocOptions::max_cell_aspect_ratio,
                DOC(Goddard, MocOptions, max_cell_aspect_ratio))
        .def_rw("max_front_points", &Goddard::MocOptions::max_front_points,
                DOC(Goddard, MocOptions, max_front_points))
        .def_rw("start_line", &Goddard::MocOptions::start_line,
                DOC(Goddard, MocOptions, start_line))
        .def_rw("kl_max_wall_angle_error", &Goddard::MocOptions::kl_max_wall_angle_error,
                DOC(Goddard, MocOptions, kl_max_wall_angle_error))
        .def_rw("march_scheme", &Goddard::MocOptions::march_scheme,
                DOC(Goddard, MocOptions, march_scheme))
        .def_rw("inverse_cfl", &Goddard::MocOptions::inverse_cfl,
                DOC(Goddard, MocOptions, inverse_cfl))
        .def_rw("max_wall_turn_per_step", &Goddard::MocOptions::max_wall_turn_per_step,
                DOC(Goddard, MocOptions, max_wall_turn_per_step))
        .def_rw("front_tilt_decay", &Goddard::MocOptions::front_tilt_decay,
                DOC(Goddard, MocOptions, front_tilt_decay))
        .def("__repr__", [](const Goddard::MocOptions& self) {
            return std::format("<MocOptions mode={} flow={} N={}>",
                               mode_name(self.mode), flow_kind_name(self.flow_type),
                               self.num_characteristics);
        });

    // ---- CharacteristicPoint ----

    nb::class_<Goddard::CharacteristicPoint>(m, "CharacteristicPoint",
                                             DOC(Goddard, CharacteristicPoint))
        .def(nb::init<>())
        .def_ro("theta", &Goddard::CharacteristicPoint::theta,
                DOC(Goddard, CharacteristicPoint, theta))
        .def_ro("nu", &Goddard::CharacteristicPoint::nu, DOC(Goddard, CharacteristicPoint, nu))
        .def_ro("pressure", &Goddard::CharacteristicPoint::pressure,
                DOC(Goddard, CharacteristicPoint, pressure))
        .def_ro("temperature", &Goddard::CharacteristicPoint::temperature,
                DOC(Goddard, CharacteristicPoint, temperature))
        .def_ro("gamma_s", &Goddard::CharacteristicPoint::gamma_s,
                DOC(Goddard, CharacteristicPoint, gamma_s))
        .def_ro("mach", &Goddard::CharacteristicPoint::mach,
                DOC(Goddard, CharacteristicPoint, mach))
        .def_ro("V", &Goddard::CharacteristicPoint::V, DOC(Goddard, CharacteristicPoint, V))
        .def_ro("mu", &Goddard::CharacteristicPoint::mu, DOC(Goddard, CharacteristicPoint, mu))
        .def_ro("K_plus", &Goddard::CharacteristicPoint::K_plus,
                DOC(Goddard, CharacteristicPoint, K_plus))
        .def_ro("K_minus", &Goddard::CharacteristicPoint::K_minus,
                DOC(Goddard, CharacteristicPoint, K_minus))
        .def_ro("x", &Goddard::CharacteristicPoint::x, DOC(Goddard, CharacteristicPoint, x))
        .def_ro("y", &Goddard::CharacteristicPoint::y, DOC(Goddard, CharacteristicPoint, y))
        // cantera_state intentionally omitted (internal detail)
        .def("__repr__", [](const Goddard::CharacteristicPoint& self) {
            return std::format("<CharacteristicPoint x={:.4g} y={:.4g} M={:.4g} theta={:.4g}>",
                               self.x, self.y, self.mach, self.theta);
        });

    // ---- Net topology ----

    nb::class_<Goddard::PointMembership>(m, "PointMembership", DOC(Goddard, PointMembership))
        .def(nb::init<>())
        .def_rw("c_plus_chain_idx", &Goddard::PointMembership::c_plus_chain_idx,
                DOC(Goddard, PointMembership, c_plus_chain_idx))
        .def_rw("c_minus_chain_idx", &Goddard::PointMembership::c_minus_chain_idx,
                DOC(Goddard, PointMembership, c_minus_chain_idx));

    nb::class_<Goddard::ChainMetadata>(m, "ChainMetadata", DOC(Goddard, ChainMetadata))
        .def(nb::init<>())
        .def_rw("active", &Goddard::ChainMetadata::active, DOC(Goddard, ChainMetadata, active))
        // See the enum bindings above: the inline `enum class ... member;` declarations make
        // the DOC() key ambiguous, so these two carry literal docstrings.
        .def_rw("termination", &Goddard::ChainMetadata::termination,
                "How this chain stopped being marched.")
        .def_rw("family", &Goddard::ChainMetadata::family,
                "Which characteristic family this chain belongs to.")
        .def_rw("origin_point_idx", &Goddard::ChainMetadata::origin_point_idx,
                DOC(Goddard, ChainMetadata, origin_point_idx))
        .def_rw("latest_point_idx", &Goddard::ChainMetadata::latest_point_idx,
                DOC(Goddard, ChainMetadata, latest_point_idx));

    // ---- CharacteristicNet ----

    nb::class_<Goddard::CharacteristicNet>(m, "CharacteristicNet")
        .def(nb::init<>())
        .def_ro("points", &Goddard::CharacteristicNet::points,
                DOC(Goddard, CharacteristicNet, points))
        .def_ro("chain_metadata", &Goddard::CharacteristicNet::chain_metadata,
                DOC(Goddard, CharacteristicNet, chain_metadata))
        .def_ro("membership", &Goddard::CharacteristicNet::membership,
                DOC(Goddard, CharacteristicNet, membership))

        // Ragged, so a list of index arrays rather than a single 2-D array. Combined with the
        // columnar properties below this is what draws the mesh: x[chain], y[chain].
        .def_prop_ro("chains", [](const Goddard::CharacteristicNet& self) {
            std::vector<IndexArray> chains;
            chains.reserve(self.c_chains.size());
            for (const std::vector<size_t>& chain : self.c_chains) {
                chains.push_back(to_index_array(chain));
            }
            return chains;
        }, DOC(Goddard, CharacteristicNet, c_chains))

        .def_prop_ro("wall_x", [](const Goddard::CharacteristicNet& self) {
            return to_array(self.wall_x);
        }, DOC(Goddard, CharacteristicNet, wall_x))
        .def_prop_ro("wall_y", [](const Goddard::CharacteristicNet& self) {
            return to_array(self.wall_y);
        }, DOC(Goddard, CharacteristicNet, wall_y))
        .def_prop_ro("wall_point_indices", [](const Goddard::CharacteristicNet& self) {
            return to_index_array(self.wall_point_indices);
        }, DOC(Goddard, CharacteristicNet, wall_point_indices))
        .def_prop_ro("axis_point_indices", [](const Goddard::CharacteristicNet& self) {
            return to_index_array(self.axis_point_indices);
        }, DOC(Goddard, CharacteristicNet, axis_point_indices))

        // Ragged, like c_chains above; nanobind converts std::vector<std::vector<size_t>>
        // straight to a list of lists of Python ints.
        .def_ro("fronts", &Goddard::CharacteristicNet::fronts,
                DOC(Goddard, CharacteristicNet, fronts))

        // ---- Columnar views over every point, for plotting and post-processing ----

        .def_prop_ro("x", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::x);
        }, "Axial position of every point, in length units.")
        .def_prop_ro("y", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::y);
        }, "Radial (or transverse) position of every point, in length units.")
        .def_prop_ro("theta", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::theta);
        }, "Flow angle at every point, in radians.")
        .def_prop_ro("nu", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::nu);
        }, "Prandtl-Meyer angle at every point, in radians.")
        .def_prop_ro("mach", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::mach);
        }, "Mach number at every point.")
        .def_prop_ro("mu", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::mu);
        }, "Mach angle at every point, in radians.")
        .def_prop_ro("pressure", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::pressure);
        }, "Static pressure at every point.")
        .def_prop_ro("temperature", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::temperature);
        }, "Static temperature at every point.")
        .def_prop_ro("gamma_s", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::gamma_s);
        }, "Local isentropic exponent at every point.")
        .def_prop_ro("V", [](const Goddard::CharacteristicNet& self) {
            return gather(self.points, &Goddard::CharacteristicPoint::V);
        }, "Velocity at every point; m/s, or the Mach number for perfect gas.")

        // ---- Boundary and chain queries ----

        .def("wall_points", &Goddard::CharacteristicNet::wall_points,
             DOC(Goddard, CharacteristicNet, wall_points))
        .def("axis_points", &Goddard::CharacteristicNet::axis_points,
             DOC(Goddard, CharacteristicNet, axis_points))
        .def("outflow_points", &Goddard::CharacteristicNet::outflow_points,
             DOC(Goddard, CharacteristicNet, outflow_points))
        .def("leading_point",
             nb::overload_cast<size_t>(&Goddard::CharacteristicNet::leading_point, nb::const_),
             "chain_idx"_a, DOC(Goddard, CharacteristicNet, leading_point))
        .def("has_active_chains", &Goddard::CharacteristicNet::has_active_chains,
             DOC(Goddard, CharacteristicNet, has_active_chains))
        .def("empty", &Goddard::CharacteristicNet::empty,
             DOC(Goddard, CharacteristicNet, empty))
        .def("__bool__", [](const Goddard::CharacteristicNet& self) { return !self.empty(); })
        .def("__len__", [](const Goddard::CharacteristicNet& self) { return self.points.size(); })
        .def("__repr__", [](const Goddard::CharacteristicNet& self) {
            return std::format("<CharacteristicNet {} points, {} chains, {} wall points>",
                               self.points.size(), self.c_chains.size(),
                               self.wall_point_indices.size());
        });

    // ---- ExitPlane ----

    nb::class_<Goddard::ExitPlane>(m, "ExitPlane", DOC(Goddard, ExitPlane))
        .def(nb::init<>())
        .def_prop_ro("y", [](const Goddard::ExitPlane& self) { return to_array(self.y); },
                     DOC(Goddard, ExitPlane, y))
        .def_prop_ro("mach", [](const Goddard::ExitPlane& self) { return to_array(self.mach); },
                     DOC(Goddard, ExitPlane, mach))
        .def_prop_ro("theta", [](const Goddard::ExitPlane& self) { return to_array(self.theta); },
                     DOC(Goddard, ExitPlane, theta))
        .def_prop_ro("pressure", [](const Goddard::ExitPlane& self) { return to_array(self.pressure); },
                     DOC(Goddard, ExitPlane, pressure))
        .def_prop_ro("temperature", [](const Goddard::ExitPlane& self) { return to_array(self.temperature); },
                     DOC(Goddard, ExitPlane, temperature))
        .def_prop_ro("gamma_s", [](const Goddard::ExitPlane& self) { return to_array(self.gamma_s); },
                     DOC(Goddard, ExitPlane, gamma_s))
        .def_prop_ro("velocity", [](const Goddard::ExitPlane& self) { return to_array(self.velocity); },
                     DOC(Goddard, ExitPlane, velocity))
        .def("__len__", [](const Goddard::ExitPlane& self) { return self.y.size(); });

    // ---- MocFailure ----

    nb::class_<Goddard::MocFailure>(m, "MocFailure", DOC(Goddard, MocFailure))
        .def(nb::init<>())
        .def_ro("code", &Goddard::MocFailure::code, DOC(Goddard, MocFailure, code))
        .def_ro("message", &Goddard::MocFailure::message, DOC(Goddard, MocFailure, message))
        .def_ro("x", &Goddard::MocFailure::x, DOC(Goddard, MocFailure, x))
        .def_ro("y", &Goddard::MocFailure::y, DOC(Goddard, MocFailure, y))
        .def_ro("kernel_pass", &Goddard::MocFailure::kernel_pass,
                DOC(Goddard, MocFailure, kernel_pass))
        .def("__repr__", [](const Goddard::MocFailure& self) {
            if (self.code == Goddard::MocErrorCode::NONE) {
                return std::string("<MocFailure NONE>");
            }
            return std::format("<MocFailure {} at ({:.4g}, {:.4g}) on pass {}: {}>",
                               Goddard::to_string(self.code), self.x, self.y,
                               self.kernel_pass, self.message);
        });

    // ---- MocPassDiagnostics ----

    nb::class_<Goddard::MocPassDiagnostics>(m, "MocPassDiagnostics",
                                            DOC(Goddard, MocPassDiagnostics))
        .def(nb::init<>())
        // `pass` is a Python keyword, so the field is renamed rather than made unreachable.
        .def_ro("pass_index", &Goddard::MocPassDiagnostics::pass,
                DOC(Goddard, MocPassDiagnostics, pass))
        .def_ro("front_points", &Goddard::MocPassDiagnostics::front_points,
                DOC(Goddard, MocPassDiagnostics, front_points))
        .def_ro("min_spacing", &Goddard::MocPassDiagnostics::min_spacing,
                DOC(Goddard, MocPassDiagnostics, min_spacing))
        .def_ro("max_spacing", &Goddard::MocPassDiagnostics::max_spacing,
                DOC(Goddard, MocPassDiagnostics, max_spacing))
        .def_ro("mean_spacing", &Goddard::MocPassDiagnostics::mean_spacing,
                DOC(Goddard, MocPassDiagnostics, mean_spacing))
        .def_ro("target_spacing", &Goddard::MocPassDiagnostics::target_spacing,
                DOC(Goddard, MocPassDiagnostics, target_spacing))
        .def_ro("max_cell_aspect", &Goddard::MocPassDiagnostics::max_cell_aspect,
                DOC(Goddard, MocPassDiagnostics, max_cell_aspect))
        .def_ro("min_spacelike_margin", &Goddard::MocPassDiagnostics::min_spacelike_margin,
                DOC(Goddard, MocPassDiagnostics, min_spacelike_margin))
        .def_ro("inserted", &Goddard::MocPassDiagnostics::inserted,
                DOC(Goddard, MocPassDiagnostics, inserted))
        .def_ro("retired", &Goddard::MocPassDiagnostics::retired,
                DOC(Goddard, MocPassDiagnostics, retired))
        .def_ro("front_axis_x", &Goddard::MocPassDiagnostics::front_axis_x,
                DOC(Goddard, MocPassDiagnostics, front_axis_x))
        .def_ro("front_wall_x", &Goddard::MocPassDiagnostics::front_wall_x,
                DOC(Goddard, MocPassDiagnostics, front_wall_x))
        .def_ro("front_axis_spacing", &Goddard::MocPassDiagnostics::front_axis_spacing,
                DOC(Goddard, MocPassDiagnostics, front_axis_spacing))
        .def_ro("front_wall_spacing", &Goddard::MocPassDiagnostics::front_wall_spacing,
                DOC(Goddard, MocPassDiagnostics, front_wall_spacing))
        .def_ro("step_dx", &Goddard::MocPassDiagnostics::step_dx,
                DOC(Goddard, MocPassDiagnostics, step_dx))
        .def_ro("step_limiter", &Goddard::MocPassDiagnostics::step_limiter,
                DOC(Goddard, MocPassDiagnostics, step_limiter))
        .def("__repr__", [](const Goddard::MocPassDiagnostics& self) {
            return std::format(
                "<MocPassDiagnostics pass={} front_points={} spacing=[{:.3g}, {:.3g}] "
                "target={:.3g} aspect={:.3g}>",
                self.pass, self.front_points, self.min_spacing, self.max_spacing,
                self.target_spacing, self.max_cell_aspect);
        });

    // ---- MocFrontShear ----

    nb::class_<Goddard::MocFrontShear>(m, "MocFrontShear", DOC(Goddard, MocFrontShear))
        .def(nb::init<>())
        .def_ro("axis_growth", &Goddard::MocFrontShear::axis_growth,
                DOC(Goddard, MocFrontShear, axis_growth))
        .def_ro("wall_growth", &Goddard::MocFrontShear::wall_growth,
                DOC(Goddard, MocFrontShear, wall_growth))
        .def_ro("shear_ratio", &Goddard::MocFrontShear::shear_ratio,
                DOC(Goddard, MocFrontShear, shear_ratio))
        .def_ro("valid", &Goddard::MocFrontShear::valid, DOC(Goddard, MocFrontShear, valid))
        .def("__repr__", [](const Goddard::MocFrontShear& self) {
            return std::format(
                "<MocFrontShear axis={:.3g} wall={:.3g} ratio={:.3g} valid={}>",
                self.axis_growth, self.wall_growth, self.shear_ratio, self.valid);
        });

    m.def("summarize_front_shear", &Goddard::summarize_front_shear,
          "pass_diagnostics"_a, "window"_a = 4,
          DOC(Goddard, summarize_front_shear));

    // ---- MocInitDiagnostics ----

    nb::class_<Goddard::MocInitDiagnostics>(m, "MocInitDiagnostics",
                                            DOC(Goddard, MocInitDiagnostics))
        .def(nb::init<>())
        .def_ro("points", &Goddard::MocInitDiagnostics::points,
                DOC(Goddard, MocInitDiagnostics, points))
        .def_ro("wall_gap", &Goddard::MocInitDiagnostics::wall_gap,
                DOC(Goddard, MocInitDiagnostics, wall_gap))
        .def_ro("wall_gap_over_spacing", &Goddard::MocInitDiagnostics::wall_gap_over_spacing,
                DOC(Goddard, MocInitDiagnostics, wall_gap_over_spacing))
        .def_ro("wall_theta_mismatch", &Goddard::MocInitDiagnostics::wall_theta_mismatch,
                DOC(Goddard, MocInitDiagnostics, wall_theta_mismatch))
        .def_ro("mach_axis", &Goddard::MocInitDiagnostics::mach_axis,
                DOC(Goddard, MocInitDiagnostics, mach_axis))
        .def_ro("mach_wall", &Goddard::MocInitDiagnostics::mach_wall,
                DOC(Goddard, MocInitDiagnostics, mach_wall))
        .def_ro("mach_ratio", &Goddard::MocInitDiagnostics::mach_ratio,
                DOC(Goddard, MocInitDiagnostics, mach_ratio))
        .def_ro("mu_axis", &Goddard::MocInitDiagnostics::mu_axis,
                DOC(Goddard, MocInitDiagnostics, mu_axis))
        .def_ro("mu_wall", &Goddard::MocInitDiagnostics::mu_wall,
                DOC(Goddard, MocInitDiagnostics, mu_wall))
        .def_ro("cot_mu_ratio", &Goddard::MocInitDiagnostics::cot_mu_ratio,
                DOC(Goddard, MocInitDiagnostics, cot_mu_ratio))
        .def_ro("axis_arrival_grading", &Goddard::MocInitDiagnostics::axis_arrival_grading,
                DOC(Goddard, MocInitDiagnostics, axis_arrival_grading))
        .def_ro("min_spacelike_margin", &Goddard::MocInitDiagnostics::min_spacelike_margin,
                DOC(Goddard, MocInitDiagnostics, min_spacelike_margin))
        .def_ro("mass_flow_error", &Goddard::MocInitDiagnostics::mass_flow_error,
                DOC(Goddard, MocInitDiagnostics, mass_flow_error))
        .def_ro("shift_over_transonic_length",
                &Goddard::MocInitDiagnostics::shift_over_transonic_length,
                DOC(Goddard, MocInitDiagnostics, shift_over_transonic_length))
        .def_ro("wall_station_to_tangency",
                &Goddard::MocInitDiagnostics::wall_station_to_tangency,
                DOC(Goddard, MocInitDiagnostics, wall_station_to_tangency))
        .def_ro("wall_bc_residual", &Goddard::MocInitDiagnostics::wall_bc_residual,
                DOC(Goddard, MocInitDiagnostics, wall_bc_residual))
        .def_ro("kplus_wall_end", &Goddard::MocInitDiagnostics::kplus_wall_end,
                DOC(Goddard, MocInitDiagnostics, kplus_wall_end))
        .def_ro("start_line_used", &Goddard::MocInitDiagnostics::start_line_used,
                DOC(Goddard, MocInitDiagnostics, start_line_used))
        .def("__repr__", [](const Goddard::MocInitDiagnostics& self) {
            return std::format(
                "<MocInitDiagnostics points={} wall_gap={:.3g} mach=[{:.4g}, {:.4g}] "
                "mass_flow_error={:.3g}>",
                self.points, self.wall_gap, self.mach_axis, self.mach_wall,
                self.mass_flow_error);
        });

    // ---- MocCrossings ----

    nb::class_<Goddard::MocCrossings>(m, "MocCrossings", DOC(Goddard, MocCrossings))
        .def(nb::init<>())
        .def_ro("count", &Goddard::MocCrossings::count, DOC(Goddard, MocCrossings, count))
        .def_ro("first_x", &Goddard::MocCrossings::first_x, DOC(Goddard, MocCrossings, first_x))
        .def_ro("first_y", &Goddard::MocCrossings::first_y, DOC(Goddard, MocCrossings, first_y))
        .def_ro("first_family", &Goddard::MocCrossings::first_family,
                DOC(Goddard, MocCrossings, first_family))
        .def("__repr__", [](const Goddard::MocCrossings& self) {
            if (self.count == 0) {
                return std::string("<MocCrossings none>");
            }
            return std::format("<MocCrossings count={} first at ({:.4g}, {:.4g})>",
                               self.count, self.first_x, self.first_y);
        });

    // ---- MocResult ----

    nb::class_<Goddard::MocResult>(m, "MocResult", DOC(Goddard, MocResult))
        .def(nb::init<>())
        .def_ro("converged", &Goddard::MocResult::converged, DOC(Goddard, MocResult, converged))
        .def_ro("net", &Goddard::MocResult::net, DOC(Goddard, MocResult, net))
        .def_ro("profile", &Goddard::MocResult::profile, DOC(Goddard, MocResult, profile))
        .def_ro("messages", &Goddard::MocResult::messages, DOC(Goddard, MocResult, messages))
        .def_ro("failure", &Goddard::MocResult::failure, DOC(Goddard, MocResult, failure))
        .def_ro("exit_mach", &Goddard::MocResult::exit_mach, DOC(Goddard, MocResult, exit_mach))
        .def_ro("nozzle_length", &Goddard::MocResult::nozzle_length,
                DOC(Goddard, MocResult, nozzle_length))
        .def_ro("area_ratio", &Goddard::MocResult::area_ratio, DOC(Goddard, MocResult, area_ratio))
        .def_ro("inserted_characteristics", &Goddard::MocResult::inserted_characteristics,
                DOC(Goddard, MocResult, inserted_characteristics))
        .def_ro("retired_characteristics", &Goddard::MocResult::retired_characteristics,
                DOC(Goddard, MocResult, retired_characteristics))
        .def_ro("pass_diagnostics", &Goddard::MocResult::pass_diagnostics,
                DOC(Goddard, MocResult, pass_diagnostics))
        .def_ro("init_diagnostics", &Goddard::MocResult::init_diagnostics,
                DOC(Goddard, MocResult, init_diagnostics))
        .def_ro("exit_coverage", &Goddard::MocResult::exit_coverage,
                DOC(Goddard, MocResult, exit_coverage))
        .def_ro("crossings", &Goddard::MocResult::crossings, DOC(Goddard, MocResult, crossings))
        .def_ro("reached_exit_plane", &Goddard::MocResult::reached_exit_plane,
                DOC(Goddard, MocResult, reached_exit_plane))
        .def_ro("min_theta", &Goddard::MocResult::min_theta, DOC(Goddard, MocResult, min_theta))
        .def_ro("min_theta_x", &Goddard::MocResult::min_theta_x,
                DOC(Goddard, MocResult, min_theta_x))
        .def_ro("min_theta_y", &Goddard::MocResult::min_theta_y,
                DOC(Goddard, MocResult, min_theta_y))
        .def_ro("exit_plane", &Goddard::MocResult::exit_plane, DOC(Goddard, MocResult, exit_plane))
        .def("__repr__", [](const Goddard::MocResult& self) {
            return std::format(
                "<MocResult converged={} exit_mach={:.4g} area_ratio={:.4g} length={:.4g}>",
                self.converged, self.exit_mach, self.area_ratio, self.nozzle_length);
        });

    // ---- ThrustCoefficient ----

    nb::class_<Goddard::ThrustCoefficient>(m, "ThrustCoefficient", DOC(Goddard, ThrustCoefficient))
        .def(nb::init<>())
        .def_ro("Cf_vacuum", &Goddard::ThrustCoefficient::Cf_vacuum,
                DOC(Goddard, ThrustCoefficient, Cf_vacuum))
        .def_ro("Cf", &Goddard::ThrustCoefficient::Cf, DOC(Goddard, ThrustCoefficient, Cf))
        .def_ro("momentum_thrust", &Goddard::ThrustCoefficient::momentum_thrust,
                DOC(Goddard, ThrustCoefficient, momentum_thrust))
        .def_ro("pressure_thrust", &Goddard::ThrustCoefficient::pressure_thrust,
                DOC(Goddard, ThrustCoefficient, pressure_thrust))
        .def("__repr__", [](const Goddard::ThrustCoefficient& self) {
            return std::format("<ThrustCoefficient Cf={:.4g} Cf_vacuum={:.4g}>",
                               self.Cf, self.Cf_vacuum);
        });

    // ---- MocNozzle ----

    nb::class_<Goddard::MocNozzle>(m, "MocNozzle", DOC(Goddard, MocNozzle))
        .def("__init__", [](Goddard::MocNozzle* self, Goddard::MocOptions options) {
            new (self) Goddard::MocNozzle(std::move(options));
        }, "options"_a)
        .def("__init__", [](Goddard::MocNozzle* self,
                            std::shared_ptr<Cantera::Solution> sol,
                            Goddard::MocOptions options) {
            Goddard::Gas gas(std::move(sol), options.chemistry);
            new (self) Goddard::MocNozzle(std::move(gas), std::move(options));
        }, "solution"_a, "options"_a)
        .def("__init__", [](Goddard::MocNozzle* self,
                            Goddard::Gas& gas,
                            Goddard::MocOptions options) {
            new (self) Goddard::MocNozzle(gas, std::move(options));
        }, "gas"_a, "options"_a)
        .def_rw("options", &Goddard::MocNozzle::m_options,
                "Solver configuration. Public and re-read by solve(), so it can be adjusted "
                "between solves.")
        .def("solve", &Goddard::MocNozzle::solve, DOC(Goddard, MocNozzle, solve))
        .def("is_solved", &Goddard::MocNozzle::is_solved, DOC(Goddard, MocNozzle, is_solved))
        .def("__repr__", [](const Goddard::MocNozzle& self) {
            return std::format("<MocNozzle mode={} flow={} solved={}>",
                               mode_name(self.m_options.mode),
                               flow_kind_name(self.m_options.flow_type),
                               self.is_solved());
        });

    // ---- Free functions ----

    m.def("compute_thrust_coefficient",
          &Goddard::compute_thrust_coefficient,
          "result"_a, "flow_type"_a, "ambient_pressure_ratio"_a = 0.0,
          DOC(Goddard, compute_thrust_coefficient));

    m.def("find_like_characteristic_crossings",
          &Goddard::find_like_characteristic_crossings,
          "net"_a, DOC(Goddard, find_like_characteristic_crossings));

    m.def("validate_moc_options",
          &Goddard::validate_moc_options,
          "options"_a, DOC(Goddard, validate_moc_options));
}
