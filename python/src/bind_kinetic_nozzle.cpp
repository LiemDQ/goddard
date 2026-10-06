#include <nanobind/nanobind.h>
#include <format>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/kinetic_nozzle.hpp"
#include "goddard/gas.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_kinetic_nozzle(nb::module_& m) {
    // KineticNozzleStation — spatially-resolved output at one axial position
    nb::class_<Goddard::KineticNozzleStation>(m, "KineticNozzleStation",
            DOC(Goddard, KineticNozzleStation))
        .def_ro("x",              &Goddard::KineticNozzleStation::x,
                DOC(Goddard, KineticNozzleStation, x))
        .def_ro("velocity",       &Goddard::KineticNozzleStation::velocity,
                DOC(Goddard, KineticNozzleStation, velocity))
        .def_ro("mach",           &Goddard::KineticNozzleStation::mach,
                DOC(Goddard, KineticNozzleStation, mach))
        .def_ro("area_ratio",     &Goddard::KineticNozzleStation::area_ratio,
                DOC(Goddard, KineticNozzleStation, area_ratio))
        .def_ro("thermo",         &Goddard::KineticNozzleStation::thermo,
                DOC(Goddard, KineticNozzleStation, thermo))
        .def_ro("state",          &Goddard::KineticNozzleStation::state,
                DOC(Goddard, KineticNozzleStation, state))
        .def_ro("damkohler",      &Goddard::KineticNozzleStation::damkohler,
                DOC(Goddard, KineticNozzleStation, damkohler))
        .def_ro("Da_min",         &Goddard::KineticNozzleStation::Da_min,
                DOC(Goddard, KineticNozzleStation, Da_min))
        .def_ro("min_Da_species", &Goddard::KineticNozzleStation::min_Da_species,
                DOC(Goddard, KineticNozzleStation, min_Da_species))
        .def("__repr__", [](const Goddard::KineticNozzleStation& self) {
            return std::format("<KineticNozzleStation x={:.6g} A/At={:.6g} T={:.6g} K P={:.6g} Pa M={:.6g} Da_min={:.3g} ({})>", self.x, self.area_ratio, self.thermo.temperature, self.thermo.pressure, self.mach, self.Da_min, self.min_Da_species);
        });

    // KineticNozzleResults — full solver output
    nb::class_<Goddard::KineticNozzleResults>(m, "KineticNozzleResults",
            DOC(Goddard, KineticNozzleResults))
        .def_ro("throat",   &Goddard::KineticNozzleResults::throat,
                DOC(Goddard, KineticNozzleResults, throat))
        .def_ro("stations", &Goddard::KineticNozzleResults::stations,
                DOC(Goddard, KineticNozzleResults, stations))
        // False if solve() stopped at max_steps before the exit of the profile.
        .def_ro("reached_exit", &Goddard::KineticNozzleResults::reached_exit,
                DOC(Goddard, KineticNozzleResults, reached_exit));

    // KineticNozzle — 1D kinetic nozzle solver using Cantera IdealGasMoleReactor
    // NozzleProfile is passed by value (it is copyable).
    // `options.chemistry` (NozzleOptions) selects the throat model: EQUILIBRIUM or FROZEN.
    // KINETIC is not a valid throat model and raises ValueError. The expansion itself is
    // always finite-rate.
    nb::class_<Goddard::KineticNozzle>(m, "KineticNozzle", DOC(Goddard, KineticNozzle))
    
        // Gas-based constructors
        .def("__init__",
             [](Goddard::KineticNozzle* self,
                Goddard::Gas gas,
                Goddard::NozzleProfile profile,
                double mdot,
                Goddard::NozzleOptions options) {
                 new (self) Goddard::KineticNozzle(gas, profile, mdot, std::move(options));
             },
             "gas"_a, "profile"_a, "mdot"_a, "options"_a=Goddard::NozzleOptions{},
             DOC(Goddard, KineticNozzle, KineticNozzle, 3))
        .def("__init__",
             [](Goddard::KineticNozzle* self,
                Goddard::Gas gas,
                Goddard::NozzleProfile profile,
                double mdot,
                std::vector<double> state,
                Goddard::NozzleOptions options) {
                 new (self) Goddard::KineticNozzle(gas, profile, mdot, std::move(state),
                                                    std::move(options));
             },
             "gas"_a, "profile"_a, "mdot"_a, "state"_a, "options"_a=Goddard::NozzleOptions{},
             DOC(Goddard, KineticNozzle, KineticNozzle, 4))
        .def("solve", &Goddard::KineticNozzle::solve,
             "dt_max"_a = 1e-6, "dx_max"_a = 1e-3, "max_steps"_a = 100000,
             DOC(Goddard, KineticNozzle, solve))
        .def_rw("profile", &Goddard::KineticNozzle::profile,
                DOC(Goddard, KineticNozzle, profile))
        .def_rw("mdot",    &Goddard::KineticNozzle::mdot,
                DOC(Goddard, KineticNozzle, mdot));
}
