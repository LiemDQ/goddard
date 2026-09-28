#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/pair.h>
#include "goddard/shocks.hpp"
#include "goddard/gas.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_shocks(nb::module_& m) {
    // ShockResult
    nb::class_<Goddard::ShockResult>(m, "ShockResult", DOC(Goddard, ShockResult))
        .def(nb::init<>())
        .def_ro("valid", &Goddard::ShockResult::valid, DOC(Goddard, ShockResult, valid))
        .def_ro("mach_in", &Goddard::ShockResult::mach_in, DOC(Goddard, ShockResult, mach_in))
        .def_ro("mach_out", &Goddard::ShockResult::mach_out, DOC(Goddard, ShockResult, mach_out))
        .def_ro("static_pressure_ratio", &Goddard::ShockResult::static_pressure_ratio,
                DOC(Goddard, ShockResult, static_pressure_ratio))
        .def_ro("static_temperature_ratio", &Goddard::ShockResult::static_temperature_ratio,
                DOC(Goddard, ShockResult, static_temperature_ratio))
        .def_ro("density_ratio", &Goddard::ShockResult::density_ratio,
                DOC(Goddard, ShockResult, density_ratio))
        .def_ro("total_pressure_ratio", &Goddard::ShockResult::total_pressure_ratio,
                DOC(Goddard, ShockResult, total_pressure_ratio));

    // ReflectedShockResult
    nb::class_<Goddard::ReflectedShockResult>(m, "ReflectedShockResult", DOC(Goddard, ReflectedShockResult))
        .def(nb::init<>())
        .def_ro("valid", &Goddard::ReflectedShockResult::valid, DOC(Goddard, ReflectedShockResult, valid))
        .def_ro("incident", &Goddard::ReflectedShockResult::incident,
                DOC(Goddard, ReflectedShockResult, incident))
        .def_ro("reflected", &Goddard::ReflectedShockResult::reflected,
                DOC(Goddard, ReflectedShockResult, reflected));

    // ObliqueShockResult
    nb::class_<Goddard::ObliqueShockResult>(m, "ObliqueShockResult", DOC(Goddard, ObliqueShockResult))
        .def(nb::init<>())
        .def_ro("valid", &Goddard::ObliqueShockResult::valid, DOC(Goddard, ObliqueShockResult, valid))
        .def_ro("mach_in", &Goddard::ObliqueShockResult::mach_in, DOC(Goddard, ObliqueShockResult, mach_in))
        .def_ro("mach_out", &Goddard::ObliqueShockResult::mach_out, DOC(Goddard, ObliqueShockResult, mach_out))
        .def_ro("shock", &Goddard::ObliqueShockResult::shock, DOC(Goddard, ObliqueShockResult, shock))
        .def_ro("beta", &Goddard::ObliqueShockResult::beta, DOC(Goddard, ObliqueShockResult, beta))
        .def_ro("theta", &Goddard::ObliqueShockResult::theta, DOC(Goddard, ObliqueShockResult, theta));

    // Perfect-gas free functions
    m.def("normal_shock",
          nb::overload_cast<double, double>(&Goddard::normal_shock),
          "mach"_a, "gamma"_a, DOC(Goddard, normal_shock));
    m.def("reflected_shock",
          nb::overload_cast<double, double>(&Goddard::reflected_shock),
          "mach"_a, "gamma"_a, DOC(Goddard, reflected_shock));
    m.def("oblique_shock_wave_angle",
          &Goddard::oblique_shock_wave_angle,
          "mach"_a, "deflection_angle"_a, "gamma"_a, DOC(Goddard, oblique_shock_wave_angle));
    m.def("oblique_shock_deflection_angle",
          &Goddard::oblique_shock_deflection_angle,
          "mach"_a, "wave_angle"_a, "gamma"_a, DOC(Goddard, oblique_shock_deflection_angle));
    m.def("oblique_shock_max_deflection_wave_angle",
          &Goddard::oblique_shock_max_deflection_wave_angle,
          "mach"_a, "gamma"_a, DOC(Goddard, oblique_shock_max_deflection_wave_angle));
    m.def("oblique_shock_max_deflection",
          &Goddard::oblique_shock_max_deflection,
          "mach"_a, "gamma"_a, DOC(Goddard, oblique_shock_max_deflection));
    m.def("oblique_shock_from_wave_angle",
          nb::overload_cast<double, double, double>(&Goddard::oblique_shock_from_wave_angle),
          "mach"_a, "wave_angle"_a, "gamma"_a, DOC(Goddard, oblique_shock_from_wave_angle));
    m.def("oblique_shock_from_deflection",
          nb::overload_cast<double, double, double, bool>(&Goddard::oblique_shock_from_deflection),
          "mach"_a, "deflection_angle"_a, "gamma"_a, "weak"_a = true,
          DOC(Goddard, oblique_shock_from_deflection));

    // ShockSolver — chemistry-dispatched solver with state management
    nb::class_<Goddard::ShockSolver>(m, "ShockSolver", DOC(Goddard, ShockSolver))
        .def(nb::init<Goddard::Gas, Goddard::SolverOptions>(),
             "gas"_a, "options"_a = Goddard::SolverOptions(), DOC(Goddard, ShockSolver, ShockSolver))
        .def("normal_shock", &Goddard::ShockSolver::normal_shock,
             "mach"_a, DOC(Goddard, ShockSolver, normal_shock))
        .def("normal_shock_from_velocity", &Goddard::ShockSolver::normal_shock_from_velocity,
             "velocity"_a, DOC(Goddard, ShockSolver, normal_shock_from_velocity))
        .def("reflected_shock",
             nb::overload_cast<double>(&Goddard::ShockSolver::reflected_shock),
             "mach"_a, DOC(Goddard, ShockSolver, reflected_shock))
        .def("reflected_shock",
             nb::overload_cast<double, Goddard::GasChemistry, Goddard::GasChemistry>(
                 &Goddard::ShockSolver::reflected_shock),
             "mach"_a, "incident_chemistry"_a, "reflected_chemistry"_a,
             DOC(Goddard, ShockSolver, reflected_shock, 2))
        .def("reflected_shock_from_velocity",
             nb::overload_cast<double>(&Goddard::ShockSolver::reflected_shock_from_velocity),
             "velocity"_a, DOC(Goddard, ShockSolver, reflected_shock_from_velocity))
        .def("reflected_shock_from_velocity",
             nb::overload_cast<double, Goddard::GasChemistry, Goddard::GasChemistry>(
                 &Goddard::ShockSolver::reflected_shock_from_velocity),
             "velocity"_a, "incident_chemistry"_a, "reflected_chemistry"_a,
             DOC(Goddard, ShockSolver, reflected_shock_from_velocity, 2))
        .def("oblique_shock_from_wave_angle",
             &Goddard::ShockSolver::oblique_shock_from_wave_angle,
             "mach"_a, "wave_angle"_a, DOC(Goddard, ShockSolver, oblique_shock_from_wave_angle))
        .def("oblique_shock_from_deflection",
             &Goddard::ShockSolver::oblique_shock_from_deflection,
             "mach"_a, "deflection"_a, "weak"_a = true,
             DOC(Goddard, ShockSolver, oblique_shock_from_deflection))
        .def("max_deflection", &Goddard::ShockSolver::max_deflection,
             "mach"_a, DOC(Goddard, ShockSolver, max_deflection))
        .def("pre_shock_state", &Goddard::ShockSolver::pre_shock_state,
             nb::rv_policy::reference_internal, DOC(Goddard, ShockSolver, pre_shock_state))
        .def("post_shock_state", &Goddard::ShockSolver::post_shock_state,
             nb::rv_policy::reference_internal, DOC(Goddard, ShockSolver, post_shock_state));
}
