#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/nozzle.hpp"
#include "goddard/profile.hpp"
#include "goddard/gas.hpp"
#include "goddard/combustor.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_nozzle(nb::module_& m) {
    nb::class_<Goddard::Nozzle>(m, "Nozzle", DOC(Goddard, Nozzle))
   
        // Gas-based constructors
        .def("__init__", [](Goddard::Nozzle* self, 
                            const Goddard::Gas& gas, 
                            Goddard::NozzleOptions options) {
            new (self) Goddard::Nozzle(gas, options);
        }, "gas"_a, "options"_a=Goddard::NozzleOptions{}, DOC(Goddard, Nozzle, Nozzle))
        .def("__init__", [](Goddard::Nozzle* self, const Goddard::Gas& gas,
                            std::vector<double> state, Goddard::NozzleOptions options) {
            new (self) Goddard::Nozzle(gas, std::move(state), options);
        }, "gas"_a, "state"_a, "options"_a=Goddard::NozzleOptions{},
           DOC(Goddard, Nozzle, Nozzle, 2))
        .def("solve",
             nb::overload_cast<Goddard::ExpansionType, double>(
                 &Goddard::Nozzle::solve),
             "expansion_type"_a, "ratio"_a, DOC(Goddard, Nozzle, solve))
        .def("solve_ratios",
             nb::overload_cast<Goddard::ExpansionType, const std::vector<double>&>(
                 &Goddard::Nozzle::solve),
             "expansion_type"_a, "ratios"_a, DOC(Goddard, Nozzle, solve, 2))
        .def("solve_profile",
             nb::overload_cast<const Goddard::NozzleProfile&, int>(
                 &Goddard::Nozzle::solve),
             "profile"_a, "num_stations"_a = 50, DOC(Goddard, Nozzle, solve, 3))
        .def("solve_throat_conditions",
             &Goddard::Nozzle::solve_throat_conditions,
             "abstol"_a = 4e-4, DOC(Goddard, Nozzle, solve_throat_conditions))
        .def("solve_stations",
             nb::overload_cast<const Goddard::ThroatCondition&, Goddard::ExpansionType,
                 const std::vector<double>&>(&Goddard::Nozzle::solve_stations),
             "throat_condition"_a, "expansion_type"_a, "ratios"_a,
             DOC(Goddard, Nozzle, solve_stations))
        .def("solve_stations",
             nb::overload_cast<const Goddard::FiniteAreaChamber&, Goddard::ExpansionType,
                 const std::vector<double>&>(&Goddard::Nozzle::solve_stations),
             "chamber"_a, "expansion_type"_a, "ratios"_a,
             DOC(Goddard, Nozzle, solve_stations, 2))
        .def("solve_finite_area_chamber",
             &Goddard::Nozzle::solve_finite_area_chamber,
             "injector_state"_a, "type"_a, "value"_a, "reltol"_a = 1e-6,
             DOC(Goddard, Nozzle, solve_finite_area_chamber))
        .def("reset_state", &Goddard::Nozzle::reset_state, DOC(Goddard, Nozzle, reset_state))
        .def_prop_rw("inlet_state", &Goddard::Nozzle::get_inlet_state,
             &Goddard::Nozzle::set_inlet_state, DOC(Goddard, Nozzle, set_inlet_state));
}
