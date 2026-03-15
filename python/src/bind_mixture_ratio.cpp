#include <nanobind/nanobind.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/eigen/dense.h>
#include "goddard/mixture_ratio.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_mixture_ratio(nb::module_& m) {
    // MixtureRatio (single O/F ratio)
    nb::class_<Goddard::MixtureRatio>(m, "MixtureRatio")
        .def(nb::init<double, double, double>(),
             "OF"_a, "fuel_molar_mass"_a, "oxidizer_molar_mass"_a)
        .def("fuel_mole_frac", &Goddard::MixtureRatio::fuel_mole_frac)
        .def("oxidizer_mole_frac", &Goddard::MixtureRatio::oxidizer_mole_frac)
        .def("molar_ratio", &Goddard::MixtureRatio::molar_ratio)
        .def_rw("OF_ratio", &Goddard::MixtureRatio::OF_ratio)
        .def_rw("M_fuel", &Goddard::MixtureRatio::M_fuel)
        .def_rw("M_ox", &Goddard::MixtureRatio::M_ox);

    // MixtureRatios (array of O/F ratios)
    nb::class_<Goddard::MixtureRatios>(m, "MixtureRatios")
        .def(nb::init<double, double, double>(),
             "OF"_a, "fuel_molar_mass"_a, "oxidizer_molar_mass"_a)
        .def(nb::init<const Eigen::ArrayXd&, double, double>(),
             "OF"_a, "fuel_molar_mass"_a, "oxidizer_molar_mass"_a)
        .def("OF_ratio", &Goddard::MixtureRatios::OF_ratio)
        .def("fuel_mass_frac", &Goddard::MixtureRatios::fuel_mass_frac)
        .def("oxidizer_mass_frac", &Goddard::MixtureRatios::oxidizer_mass_frac)
        .def("fuel_mole_frac", &Goddard::MixtureRatios::fuel_mole_frac)
        .def("oxidizer_mole_frac", &Goddard::MixtureRatios::oxidizer_mole_frac)
        .def("molar_ratio", &Goddard::MixtureRatios::molar_ratio)
        .def("size", &Goddard::MixtureRatios::size)
        .def_rw("M_fuel", &Goddard::MixtureRatios::M_fuel)
        .def_rw("M_ox", &Goddard::MixtureRatios::M_ox);
}
