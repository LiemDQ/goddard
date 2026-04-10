#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/nozzle.hpp"
#include "goddard/profile.hpp"
#include "goddard/gas.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_nozzle(nb::module_& m) {
    // Constructor takes Solution& but calls shared_from_this(), so we wrap with a lambda
    nb::class_<Goddard::Nozzle>(m, "Nozzle")
   
        // Gas-based constructors
        .def("__init__", [](Goddard::Nozzle* self, 
                            const Goddard::Gas& gas, 
                            Goddard::NozzleOptions options) {
            new (self) Goddard::Nozzle(gas, options);
        }, "gas"_a, "options"_a=Goddard::NozzleOptions{})
        .def("__init__", [](Goddard::Nozzle* self, const Goddard::Gas& gas,
                            std::vector<double> state, Goddard::NozzleOptions options) {
            new (self) Goddard::Nozzle(gas, std::move(state), options);
        }, "gas"_a, "state"_a, "options"_a=Goddard::NozzleOptions{})
        .def("solve",
             nb::overload_cast<Goddard::ExpansionType, double>(
                 &Goddard::Nozzle::solve),
             "expansion_type"_a, "ratio"_a = 1.0)
        .def("solve_ratios",
             nb::overload_cast<Goddard::ExpansionType, const std::vector<double>&>(
                 &Goddard::Nozzle::solve),
             "expansion_type"_a, "ratios"_a)
        .def("solve_profile",
             nb::overload_cast<const Goddard::NozzleProfile&, int>(
                 &Goddard::Nozzle::solve),
             "profile"_a, "num_stations"_a = 50)
        .def("solve_throat_conditions",
             &Goddard::Nozzle::solve_throat_conditions,
             "abstol"_a = 4e-4)
        .def("reset_state", &Goddard::Nozzle::reset_state)
        .def("set_inlet_state", &Goddard::Nozzle::set_inlet_state, "state"_a)
        .def("get_inlet_state", &Goddard::Nozzle::get_inlet_state)
        .def_rw("inlet_state", &Goddard::Nozzle::inlet_state);
}
