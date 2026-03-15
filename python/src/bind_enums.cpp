#include <nanobind/nanobind.h>
#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"

namespace nb = nanobind;

void bind_enums(nb::module_& m) {
    nb::enum_<Goddard::CombustorType>(m, "CombustorType")
        .value("INFINITE_AREA", Goddard::CombustorType::INFINITE_AREA)
        .value("FINITE_MASS_FLUX", Goddard::CombustorType::FINITE_MASS_FLUX)
        .value("FINITE_CONTRACTION_RATIO", Goddard::CombustorType::FINITE_CONTRACTION_RATIO)
        .value("NONE", Goddard::CombustorType::NONE);

    nb::enum_<Goddard::NozzleChemistryType>(m, "NozzleChemistryType")
        .value("FROZEN", Goddard::NozzleChemistryType::FROZEN)
        .value("EQUILIBRIUM", Goddard::NozzleChemistryType::EQUILIBRIUM)
        .value("KINETIC", Goddard::NozzleChemistryType::KINETIC);

    nb::enum_<Goddard::ExpansionType>(m, "ExpansionType")
        .value("SUPERSONIC_AREA_RATIO", Goddard::ExpansionType::SUPERSONIC_AREA_RATIO)
        .value("SUBSONIC_AREA_RATIO", Goddard::ExpansionType::SUBSONIC_AREA_RATIO)
        .value("PRESSURE_RATIO", Goddard::ExpansionType::PRESSURE_RATIO);
}
