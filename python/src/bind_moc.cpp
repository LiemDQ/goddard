#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/pair.h>

#include "goddard/moc.hpp"
#include "goddard/characteristics.hpp"
#include "goddard/gas.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_moc(nb::module_& m) {
    // ---- Enums ----

    nb::enum_<Goddard::MocFlowKind>(m, "MocFlowKind")
        .value("PLANAR", Goddard::MocFlowKind::PLANAR)
        .value("AXISYMMETRIC", Goddard::MocFlowKind::AXISYMMETRIC);

    nb::enum_<Goddard::MocMode>(m, "MocMode")
        .value("DESIGN_MIN_LENGTH", Goddard::MocMode::DESIGN_MIN_LENGTH)
        .value("DESIGN_RAO", Goddard::MocMode::DESIGN_RAO)
        .value("ANALYSIS", Goddard::MocMode::ANALYSIS);

    nb::enum_<Goddard::MocInitialization>(m, "MocInitialization")
        .value("STRAIGHT_SONIC_LINE", Goddard::MocInitialization::STRAIGHT_SONIC_LINE);

    // ---- ThroatGeometry ----

    nb::class_<Goddard::ThroatGeometry>(m, "ThroatGeometry")
        .def(nb::init<>())
        .def("__init__", [](Goddard::ThroatGeometry* self,
                            double throat_radius,
                            double upstream_wall_curvature_radius,
                            double downstream_wall_curvature_radius) {
            new (self) Goddard::ThroatGeometry();
            self->throat_radius = throat_radius;
            self->upstream_wall_curvature_radius = upstream_wall_curvature_radius;
            self->downstream_wall_curvature_radius = downstream_wall_curvature_radius;
        },  "throat_radius"_a = 1.0,
            "upstream_wall_curvature_radius"_a = 1.5,
            "downstream_wall_curvature_radius"_a = 0.382)
        .def_rw("throat_radius", &Goddard::ThroatGeometry::throat_radius)
        .def_rw("upstream_wall_curvature_radius", &Goddard::ThroatGeometry::upstream_wall_curvature_radius)
        .def_rw("downstream_wall_curvature_radius", &Goddard::ThroatGeometry::downstream_wall_curvature_radius);

    // ---- NozzleProfile ----

    nb::class_<Goddard::NozzleProfile>(m, "NozzleProfile")
        .def(nb::init<>())
        .def_rw("x", &Goddard::NozzleProfile::x)
        .def_rw("y", &Goddard::NozzleProfile::y)
        .def("slope_at", &Goddard::NozzleProfile::slope_at, "x_query"_a)
        .def("theta_at", &Goddard::NozzleProfile::theta_at, "x_query"_a)
        .def("max_theta", &Goddard::NozzleProfile::max_theta)
        .def("at", &Goddard::NozzleProfile::at, "idx"_a)
        .def("push_back", [](Goddard::NozzleProfile& self, double x, double y) {
            self.push_back({x, y});
        }, "x"_a, "y"_a)
        .def("__len__", &Goddard::NozzleProfile::size)
        .def_static("load_csv", &Goddard::NozzleProfile::load_profile_csv, "filename"_a)
        .def("save_csv", &Goddard::NozzleProfile::save_profile_csv, "filename"_a);

    // ---- MocOptions ----

    nb::class_<Goddard::MocOptions>(m, "MocOptions")
        .def("__init__", [](Goddard::MocOptions* self,
                            Goddard::MocFlowKind flow_type,
                            Goddard::GasChemistry chemistry,
                            Goddard::MocMode mode,
                            Goddard::MocInitialization initialization,
                            int num_characteristics,
                            double gamma,
                            Goddard::SolverOptions solver_options,
                            Goddard::ThroatGeometry geometry,
                            double theta_max,
                            double exit_mach,
                            std::vector<double> theta_schedule,
                            Goddard::NozzleProfile nozzle_profile) {
            new (self) Goddard::MocOptions();
            self->flow_type = flow_type;
            self->chemistry = chemistry;
            self->mode = mode;
            self->initialization = initialization;
            self->num_characteristics = num_characteristics;
            self->gamma = gamma;
            self->solver_options = solver_options;
            self->geometry = std::move(geometry);
            self->theta_max = theta_max;
            self->exit_mach = exit_mach;
            self->theta_schedule = std::move(theta_schedule);
            self->nozzle_profile = std::move(nozzle_profile);
        },  "flow_type"_a = Goddard::MocFlowKind::PLANAR,
            "chemistry"_a = Goddard::GasChemistry::PERFECT_GAS,
            "mode"_a = Goddard::MocMode::DESIGN_MIN_LENGTH,
            "initialization"_a = Goddard::MocInitialization::STRAIGHT_SONIC_LINE,
            "num_characteristics"_a = 10,
            "gamma"_a = 1.4,
            "solver_options"_a = Goddard::SolverOptions{.abstol = 1e-10, .reltol = 1e-5},
            "geometry"_a = Goddard::ThroatGeometry{1.0, 1.5, 0.382},
            "theta_max"_a = 0.0,
            "exit_mach"_a = 0.0,
            "theta_schedule"_a = std::vector<double>(),
            "nozzle_profile"_a = Goddard::NozzleProfile())
        .def_rw("flow_type", &Goddard::MocOptions::flow_type)
        .def_rw("chemistry", &Goddard::MocOptions::chemistry)
        .def_rw("mode", &Goddard::MocOptions::mode)
        .def_rw("initialization", &Goddard::MocOptions::initialization)
        .def_rw("num_characteristics", &Goddard::MocOptions::num_characteristics)
        .def_rw("gamma", &Goddard::MocOptions::gamma)
        .def_rw("solver_options", &Goddard::MocOptions::solver_options)
        .def_rw("geometry", &Goddard::MocOptions::geometry)
        .def_rw("theta_max", &Goddard::MocOptions::theta_max)
        .def_rw("exit_mach", &Goddard::MocOptions::exit_mach)
        .def_rw("theta_schedule", &Goddard::MocOptions::theta_schedule)
        .def_rw("nozzle_profile", &Goddard::MocOptions::nozzle_profile);

    // ---- CharacteristicPoint ----

    nb::class_<Goddard::CharacteristicPoint>(m, "CharacteristicPoint")
        .def(nb::init<>())
        .def_ro("theta", &Goddard::CharacteristicPoint::theta)
        .def_ro("nu", &Goddard::CharacteristicPoint::nu)
        .def_ro("pressure", &Goddard::CharacteristicPoint::pressure)
        .def_ro("temperature", &Goddard::CharacteristicPoint::temperature)
        .def_ro("gamma_s", &Goddard::CharacteristicPoint::gamma_s)
        .def_ro("mach", &Goddard::CharacteristicPoint::mach)
        .def_ro("V", &Goddard::CharacteristicPoint::V)
        .def_ro("mu", &Goddard::CharacteristicPoint::mu)
        .def_ro("K_plus", &Goddard::CharacteristicPoint::K_plus)
        .def_ro("K_minus", &Goddard::CharacteristicPoint::K_minus)
        .def_ro("x", &Goddard::CharacteristicPoint::x)
        .def_ro("y", &Goddard::CharacteristicPoint::y);
        // cantera_state intentionally omitted (internal detail)

    // ---- CharacteristicNet ----

    nb::class_<Goddard::CharacteristicNet>(m, "CharacteristicNet")
        .def(nb::init<>())
        .def_ro("num_c_plus", &Goddard::CharacteristicNet::num_c_plus)
        .def_ro("num_c_minus", &Goddard::CharacteristicNet::num_c_minus)
        .def_ro("wavefronts", &Goddard::CharacteristicNet::wavefronts)
        .def_ro("wall_x", &Goddard::CharacteristicNet::wall_x)
        .def_ro("wall_y", &Goddard::CharacteristicNet::wall_y)
        .def_ro("wall_points", &Goddard::CharacteristicNet::wall_points);

    // ---- ExitPlane ----

    nb::class_<Goddard::ExitPlane>(m, "ExitPlane")
        .def(nb::init<>())
        .def_ro("y", &Goddard::ExitPlane::y)
        .def_ro("mach", &Goddard::ExitPlane::mach)
        .def_ro("theta", &Goddard::ExitPlane::theta)
        .def_ro("pressure", &Goddard::ExitPlane::pressure)
        .def_ro("temperature", &Goddard::ExitPlane::temperature)
        .def_ro("gamma_s", &Goddard::ExitPlane::gamma_s)
        .def_ro("velocity", &Goddard::ExitPlane::velocity);

    // ---- MocResult ----

    nb::class_<Goddard::MocResult>(m, "MocResult")
        .def(nb::init<>())
        .def_ro("converged", &Goddard::MocResult::converged)
        .def_ro("net", &Goddard::MocResult::net)
        .def_ro("profile", &Goddard::MocResult::profile)
        .def_ro("messages", &Goddard::MocResult::messages)
        .def_ro("exit_mach", &Goddard::MocResult::exit_mach)
        .def_ro("nozzle_length", &Goddard::MocResult::nozzle_length)
        .def_ro("area_ratio", &Goddard::MocResult::area_ratio)
        .def_ro("exit_plane", &Goddard::MocResult::exit_plane);

    // ---- ThrustCoefficient ----

    nb::class_<Goddard::ThrustCoefficient>(m, "ThrustCoefficient")
        .def(nb::init<>())
        .def_ro("Cf_vacuum", &Goddard::ThrustCoefficient::Cf_vacuum)
        .def_ro("Cf", &Goddard::ThrustCoefficient::Cf)
        .def_ro("momentum_thrust", &Goddard::ThrustCoefficient::momentum_thrust)
        .def_ro("pressure_thrust", &Goddard::ThrustCoefficient::pressure_thrust);

    // ---- MocNozzle ----

    nb::class_<Goddard::MocNozzle>(m, "MocNozzle")
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
        .def("solve", &Goddard::MocNozzle::solve);

    // ---- compute_thrust_coefficient ----

    m.def("compute_thrust_coefficient",
          &Goddard::compute_thrust_coefficient,
          "result"_a, "flow_type"_a, "ambient_pressure_ratio"_a = 0.0);
}
