#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/nozzle.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_nozzle(nb::module_& m) {
    // NozzleBase (abstract -- no direct construction from Python)
    nb::class_<Goddard::NozzleBase>(m, "NozzleBase")
        .def("solve",
             nb::overload_cast<Goddard::ExpansionType, double>(
                 &Goddard::NozzleBase::solve),
             "expansion_type"_a, "ratio"_a = 1.0)
        .def("solve_ratios",
             nb::overload_cast<Goddard::ExpansionType, const std::vector<double>&>(
                 &Goddard::NozzleBase::solve),
             "expansion_type"_a, "ratios"_a)
        .def("solve_throat_conditions",
             &Goddard::NozzleBase::solve_throat_conditions,
             "abstol"_a = 4e-4)
        .def("reset_state", &Goddard::NozzleBase::reset_state)
        .def("set_inlet_state", &Goddard::NozzleBase::set_inlet_state, "state"_a)
        .def("get_inlet_state", &Goddard::NozzleBase::get_inlet_state)
        .def_rw("inlet_state", &Goddard::NozzleBase::inlet_state);

    // EquilibriumNozzle
    // Constructor takes Solution& but calls shared_from_this(), so we wrap with a lambda
    nb::class_<Goddard::EquilibriumNozzle, Goddard::NozzleBase>(m, "EquilibriumNozzle")
        .def("__init__", [](Goddard::EquilibriumNozzle* self,
                            std::shared_ptr<Cantera::Solution> sol) {
            new (self) Goddard::EquilibriumNozzle(*sol);
        }, "solution"_a)
        .def("__init__", [](Goddard::EquilibriumNozzle* self,
                            std::shared_ptr<Cantera::Solution> sol,
                            std::vector<double> state) {
            new (self) Goddard::EquilibriumNozzle(*sol, std::move(state));
        }, "solution"_a, "state"_a);

    // FrozenNozzle
    nb::class_<Goddard::FrozenNozzle, Goddard::NozzleBase>(m, "FrozenNozzle")
        .def("__init__", [](Goddard::FrozenNozzle* self,
                            std::shared_ptr<Cantera::Solution> sol) {
            new (self) Goddard::FrozenNozzle(*sol);
        }, "solution"_a)
        .def("__init__", [](Goddard::FrozenNozzle* self,
                            std::shared_ptr<Cantera::Solution> sol,
                            std::vector<double> state) {
            new (self) Goddard::FrozenNozzle(*sol, std::move(state));
        }, "solution"_a, "state"_a);
}
