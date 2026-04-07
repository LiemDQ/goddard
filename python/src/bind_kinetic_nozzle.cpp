#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/kinetic_nozzle.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_kinetic_nozzle(nb::module_& m) {
    // KineticNozzleStation — spatially-resolved output at one axial position
    nb::class_<Goddard::KineticNozzleStation>(m, "KineticNozzleStation")
        .def(nb::init<>())
        .def_ro("x",              &Goddard::KineticNozzleStation::x)
        .def_ro("velocity",       &Goddard::KineticNozzleStation::velocity)
        .def_ro("mach",           &Goddard::KineticNozzleStation::mach)
        .def_ro("area_ratio",     &Goddard::KineticNozzleStation::area_ratio)
        .def_ro("state",          &Goddard::KineticNozzleStation::state)
        .def_ro("damkohler",      &Goddard::KineticNozzleStation::damkohler)
        .def_ro("Da_min",         &Goddard::KineticNozzleStation::Da_min)
        .def_ro("min_Da_species", &Goddard::KineticNozzleStation::min_Da_species);

    // KineticNozzleResults — full solver output
    nb::class_<Goddard::KineticNozzleResults>(m, "KineticNozzleResults")
        .def(nb::init<>())
        .def_ro("throat",   &Goddard::KineticNozzleResults::throat)
        .def_ro("stations", &Goddard::KineticNozzleResults::stations);

    // KineticNozzle — 1D kinetic nozzle solver using Cantera IdealGasMoleReactor
    // NozzleProfile is passed by value (it is copyable).
    // The GasChemistry parameter selects the throat model: EQUILIBRIUM or FROZEN.
    // KINETIC is not a valid throat model and will throw.
    nb::class_<Goddard::KineticNozzle>(m, "KineticNozzle")
        .def("__init__",
             [](Goddard::KineticNozzle* self,
                std::shared_ptr<Cantera::Solution> sol,
                Goddard::NozzleProfile profile,
                double mdot,
                Goddard::GasChemistry chemistry) {
                 new (self) Goddard::KineticNozzle(*sol, profile, mdot, chemistry);
             },
             "solution"_a, "profile"_a, "mdot"_a,
             "chemistry"_a = Goddard::GasChemistry::EQUILIBRIUM)
        .def("__init__",
             [](Goddard::KineticNozzle* self,
                std::shared_ptr<Cantera::Solution> sol,
                Goddard::NozzleProfile profile,
                double mdot,
                std::vector<double> state,
                Goddard::GasChemistry chemistry) {
                 new (self) Goddard::KineticNozzle(*sol, profile, mdot,
                                                    std::move(state), chemistry);
             },
             "solution"_a, "profile"_a, "mdot"_a, "state"_a,
             "chemistry"_a = Goddard::GasChemistry::EQUILIBRIUM)
        .def("solve", &Goddard::KineticNozzle::solve,
             "dt_max"_a = 1e-6, "dx_max"_a = 1e-3, "max_steps"_a = 100000)
        .def_rw("m_profile", &Goddard::KineticNozzle::m_profile)
        .def_rw("m_mdot",    &Goddard::KineticNozzle::m_mdot);
}
