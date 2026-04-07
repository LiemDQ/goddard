#include <nanobind/nanobind.h>
#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/chemistry.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_enums(nb::module_& m) {
    nb::enum_<Goddard::CombustorType>(m, "CombustorType")
        .value("INFINITE_AREA", Goddard::CombustorType::INFINITE_AREA)
        .value("FINITE_MASS_FLUX", Goddard::CombustorType::FINITE_MASS_FLUX)
        .value("FINITE_CONTRACTION_RATIO", Goddard::CombustorType::FINITE_CONTRACTION_RATIO)
        .value("NONE", Goddard::CombustorType::NONE);

    nb::enum_<Goddard::GasChemistry>(m, "GasChemistry")
        .value("PERFECT_GAS", Goddard::GasChemistry::PERFECT_GAS)
        .value("FROZEN", Goddard::GasChemistry::FROZEN)
        .value("EQUILIBRIUM", Goddard::GasChemistry::EQUILIBRIUM)
        .value("KINETIC", Goddard::GasChemistry::KINETIC);

    nb::enum_<Goddard::ExpansionType>(m, "ExpansionType")
        .value("SUPERSONIC_AREA_RATIO", Goddard::ExpansionType::SUPERSONIC_AREA_RATIO)
        .value("SUBSONIC_AREA_RATIO", Goddard::ExpansionType::SUBSONIC_AREA_RATIO)
        .value("PRESSURE_RATIO", Goddard::ExpansionType::PRESSURE_RATIO);

    nb::class_<Goddard::SolverOptions>(m, "SolverOptions")
        .def("__init__", [](Goddard::SolverOptions* self,
                            double abstol,
                            double reltol,
                            int max_iterations) {
            new (self) Goddard::SolverOptions();
            self->abstol = abstol;
            self->reltol = reltol;
            self->max_iterations = max_iterations;
        },  "abstol"_a = 1e-6,
            "reltol"_a = 1e-5,
            "max_iterations"_a = 100)
        .def_rw("abstol", &Goddard::SolverOptions::abstol)
        .def_rw("reltol", &Goddard::SolverOptions::reltol)
        .def_rw("max_iterations", &Goddard::SolverOptions::max_iterations);
}
