#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/pair.h>
#include "goddard/shocks.hpp"
#include "goddard/gas.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_shocks(nb::module_& m) {
    // Gas — thin wrapper around Cantera::Solution with chemistry-aware properties
    nb::class_<Goddard::Gas>(m, "Gas")
        .def("__init__",
             [](Goddard::Gas* self,
                std::shared_ptr<Cantera::Solution> sol,
                Goddard::GasChemistry chemistry) {
                 new (self) Goddard::Gas(sol, chemistry);
             },
             "solution"_a, "chemistry"_a = Goddard::GasChemistry::FROZEN)
        .def("temperature", &Goddard::Gas::temperature)
        .def("pressure", &Goddard::Gas::pressure)
        .def("density", &Goddard::Gas::density)
        .def("enthalpy_mass", &Goddard::Gas::enthalpy_mass)
        .def("entropy_mass", &Goddard::Gas::entropy_mass)
        .def("cp_mass", &Goddard::Gas::cp_mass)
        .def("cv_mass", &Goddard::Gas::cv_mass)
        .def("molecular_weight", &Goddard::Gas::molecular_weight)
        .def("gamma_s", &Goddard::Gas::gamma_s)
        .def("speed_of_sound", &Goddard::Gas::speed_of_sound)
        .def("mach", &Goddard::Gas::mach, "velocity"_a)
        .def("snapshot", &Goddard::Gas::snapshot)
        .def_rw("chemistry", &Goddard::Gas::chemistry);

    // ShockResult
    nb::class_<Goddard::ShockResult>(m, "ShockResult")
        .def(nb::init<>())
        .def_ro("valid", &Goddard::ShockResult::valid)
        .def_ro("mach_in", &Goddard::ShockResult::mach_in)
        .def_ro("mach_out", &Goddard::ShockResult::mach_out)
        .def_ro("static_pressure_ratio", &Goddard::ShockResult::static_pressure_ratio)
        .def_ro("static_temperature_ratio", &Goddard::ShockResult::static_temperature_ratio)
        .def_ro("total_pressure_ratio", &Goddard::ShockResult::total_pressure_ratio);

    // ObliqueShockResult
    nb::class_<Goddard::ObliqueShockResult>(m, "ObliqueShockResult")
        .def(nb::init<>())
        .def_ro("valid", &Goddard::ObliqueShockResult::valid)
        .def_ro("mach_in", &Goddard::ObliqueShockResult::mach_in)
        .def_ro("mach_out", &Goddard::ObliqueShockResult::mach_out)
        .def_ro("shock", &Goddard::ObliqueShockResult::shock)
        .def_ro("beta", &Goddard::ObliqueShockResult::beta)
        .def_ro("theta", &Goddard::ObliqueShockResult::theta);

    // Perfect-gas free functions
    m.def("normal_shock",
          nb::overload_cast<double, double>(&Goddard::normal_shock),
          "mach"_a, "gamma"_a);
    m.def("reflected_shock",
          nb::overload_cast<double, double>(&Goddard::reflected_shock),
          "mach"_a, "gamma"_a);
    m.def("oblique_shock_wave_angle",
          &Goddard::oblique_shock_wave_angle,
          "mach"_a, "deflection_angle"_a, "gamma"_a);
    m.def("oblique_shock_deflection_angle",
          &Goddard::oblique_shock_deflection_angle,
          "mach"_a, "wave_angle"_a, "gamma"_a);
    m.def("oblique_shock_from_wave_angle",
          nb::overload_cast<double, double, double>(&Goddard::oblique_shock_from_wave_angle),
          "mach"_a, "wave_angle"_a, "gamma"_a);
    m.def("oblique_shock_from_deflection",
          nb::overload_cast<double, double, double, bool>(&Goddard::oblique_shock_from_deflection),
          "mach"_a, "deflection_angle"_a, "gamma"_a, "weak"_a = true);

    // ShockSolver — chemistry-dispatched solver with state management
    nb::class_<Goddard::ShockSolver>(m, "ShockSolver")
        .def(nb::init<Goddard::Gas, Goddard::SolverOptions>(),
             "gas"_a, "options"_a = Goddard::SolverOptions())
        .def("normal_shock", &Goddard::ShockSolver::normal_shock,
             "mach"_a)
        .def("reflected_shock", &Goddard::ShockSolver::reflected_shock,
             "mach"_a)
        .def("oblique_shock_from_wave_angle",
             &Goddard::ShockSolver::oblique_shock_from_wave_angle,
             "mach"_a, "wave_angle"_a)
        .def("oblique_shock_from_deflection",
             &Goddard::ShockSolver::oblique_shock_from_deflection,
             "mach"_a, "deflection"_a, "weak"_a = true)
        .def("pre_shock_state", &Goddard::ShockSolver::pre_shock_state,
             nb::rv_policy::reference_internal)
        .def("post_shock_state", &Goddard::ShockSolver::post_shock_state,
             nb::rv_policy::reference_internal);
}
