#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/nozzle.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_nozzle(nb::module_& m) {
    // Constructor takes Solution& but calls shared_from_this(), so we wrap with a lambda
    nb::class_<Goddard::Nozzle>(m, "Nozzle")
        .def("__init__", [](Goddard::Nozzle* self,
                            std::shared_ptr<Cantera::Solution> sol,
                            Goddard::NozzleChemistryType chemistry) {
            new (self) Goddard::Nozzle(*sol, chemistry);
        }, "solution"_a, "chemistry"_a)
        .def("__init__", [](Goddard::Nozzle* self,
                            std::shared_ptr<Cantera::Solution> sol,
                            Goddard::NozzleChemistryType chemistry,
                            std::vector<double> state) {
            new (self) Goddard::Nozzle(*sol, chemistry, std::move(state));
        }, "solution"_a, "chemistry"_a, "state"_a)
        .def("solve",
             nb::overload_cast<Goddard::ExpansionType, double>(
                 &Goddard::Nozzle::solve),
             "expansion_type"_a, "ratio"_a = 1.0)
        .def("solve_ratios",
             nb::overload_cast<Goddard::ExpansionType, const std::vector<double>&>(
                 &Goddard::Nozzle::solve),
             "expansion_type"_a, "ratios"_a)
        .def("solve_throat_conditions",
             &Goddard::Nozzle::solve_throat_conditions,
             "abstol"_a = 4e-4)
        .def("reset_state", &Goddard::Nozzle::reset_state)
        .def("set_inlet_state", &Goddard::Nozzle::set_inlet_state, "state"_a)
        .def("get_inlet_state", &Goddard::Nozzle::get_inlet_state)
        .def_rw("inlet_state", &Goddard::Nozzle::inlet_state);
}
