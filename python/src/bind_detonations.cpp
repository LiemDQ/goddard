#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/detonations.hpp"
#include "goddard/gas.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_detonations(nb::module_& m) {
    nb::enum_<Goddard::DetonationBranch>(m, "DetonationBranch", DOC(Goddard, DetonationBranch))
        .value("OVERDRIVEN", Goddard::DetonationBranch::OVERDRIVEN)
        .value("UNDERDRIVEN", Goddard::DetonationBranch::UNDERDRIVEN);

    // DetonationResult
    nb::class_<Goddard::DetonationResult>(m, "DetonationResult", DOC(Goddard, DetonationResult))
        .def(nb::init<>())
        .def_ro("valid", &Goddard::DetonationResult::valid, DOC(Goddard, DetonationResult, valid))
        .def_ro("velocity", &Goddard::DetonationResult::velocity, DOC(Goddard, DetonationResult, velocity))
        .def_ro("drive_factor", &Goddard::DetonationResult::drive_factor,
                DOC(Goddard, DetonationResult, drive_factor))
        .def_ro("mach_in", &Goddard::DetonationResult::mach_in, DOC(Goddard, DetonationResult, mach_in))
        .def_ro("mach_out", &Goddard::DetonationResult::mach_out, DOC(Goddard, DetonationResult, mach_out))
        .def_ro("static_pressure_ratio", &Goddard::DetonationResult::static_pressure_ratio,
                DOC(Goddard, DetonationResult, static_pressure_ratio))
        .def_ro("static_temperature_ratio", &Goddard::DetonationResult::static_temperature_ratio,
                DOC(Goddard, DetonationResult, static_temperature_ratio))
        .def_ro("density_ratio", &Goddard::DetonationResult::density_ratio,
                DOC(Goddard, DetonationResult, density_ratio))
        .def_ro("molecular_weight_ratio", &Goddard::DetonationResult::molecular_weight_ratio,
                DOC(Goddard, DetonationResult, molecular_weight_ratio))
        .def_ro("total_pressure_ratio", &Goddard::DetonationResult::total_pressure_ratio,
                DOC(Goddard, DetonationResult, total_pressure_ratio))
        .def_ro("von_neumann", &Goddard::DetonationResult::von_neumann,
                DOC(Goddard, DetonationResult, von_neumann));

    // ReflectedDetonationResult
    nb::class_<Goddard::ReflectedDetonationResult>(m, "ReflectedDetonationResult",
                                                   DOC(Goddard, ReflectedDetonationResult))
        .def(nb::init<>())
        .def_ro("valid", &Goddard::ReflectedDetonationResult::valid,
                DOC(Goddard, ReflectedDetonationResult, valid))
        .def_ro("incident", &Goddard::ReflectedDetonationResult::incident,
                DOC(Goddard, ReflectedDetonationResult, incident))
        .def_ro("reflected", &Goddard::ReflectedDetonationResult::reflected,
                DOC(Goddard, ReflectedDetonationResult, reflected));

    // Perfect-gas free functions
    m.def("chapman_jouguet_detonation", &Goddard::chapman_jouguet_detonation,
          "gamma"_a, "heat_release"_a, DOC(Goddard, chapman_jouguet_detonation));
    m.def("detonation",
          nb::overload_cast<double, double, double, Goddard::DetonationBranch>(&Goddard::detonation),
          "drive_factor"_a, "gamma"_a, "heat_release"_a, "branch"_a = Goddard::DetonationBranch::OVERDRIVEN,
          DOC(Goddard, detonation));
    m.def("reflected_detonation",
          nb::overload_cast<double, double, double, Goddard::DetonationBranch>(&Goddard::reflected_detonation),
          "drive_factor"_a, "gamma"_a, "heat_release"_a, "branch"_a = Goddard::DetonationBranch::OVERDRIVEN,
          DOC(Goddard, reflected_detonation));

    // DetonationSolver — real-gas detonations with state management
    nb::class_<Goddard::DetonationSolver>(m, "DetonationSolver", DOC(Goddard, DetonationSolver))
        .def(nb::init<const Goddard::Gas&, Goddard::SolverOptions>(),
             "gas"_a, "options"_a = Goddard::SolverOptions(), DOC(Goddard, DetonationSolver, DetonationSolver))
        .def("chapman_jouguet", &Goddard::DetonationSolver::chapman_jouguet,
             DOC(Goddard, DetonationSolver, chapman_jouguet))
        .def("detonation", &Goddard::DetonationSolver::detonation,
             "drive_factor"_a, "branch"_a = Goddard::DetonationBranch::OVERDRIVEN,
             DOC(Goddard, DetonationSolver, detonation))
        .def("detonation_from_velocity", &Goddard::DetonationSolver::detonation_from_velocity,
             "velocity"_a, "branch"_a = Goddard::DetonationBranch::OVERDRIVEN,
             DOC(Goddard, DetonationSolver, detonation_from_velocity))
        .def("reflected_detonation", &Goddard::DetonationSolver::reflected_detonation,
             "drive_factor"_a = 1.0, "branch"_a = Goddard::DetonationBranch::OVERDRIVEN,
             DOC(Goddard, DetonationSolver, reflected_detonation))
        .def("pre_detonation_state", &Goddard::DetonationSolver::pre_detonation_state, DOC(Goddard, DetonationSolver, pre_detonation_state))
        .def("post_detonation_state", &Goddard::DetonationSolver::post_detonation_state, DOC(Goddard, DetonationSolver, post_detonation_state))
        .def("von_neumann_state", &Goddard::DetonationSolver::von_neumann_state, DOC(Goddard, DetonationSolver, von_neumann_state));
}
