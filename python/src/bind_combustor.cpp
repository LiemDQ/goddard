#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/map.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/eigen/dense.h>
#include "goddard/combustor.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_combustor(nb::module_& m) {
    nb::class_<Goddard::BaseCombustor>(m, "BaseCombustor")
        .def("get_combustion_species", &Goddard::BaseCombustor::get_combustion_species);

    nb::class_<Goddard::Combustor, Goddard::BaseCombustor>(m, "Combustor")
        .def(nb::init<Goddard::Gas, const Goddard::Composition&, const Goddard::Composition&>(),
             "gas"_a, "fuel"_a, "oxidizer"_a)
        .def("solve",
             nb::overload_cast<const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               const Goddard::CombustorOptions&>(
                 &Goddard::Combustor::solve),
             "temperatures"_a, "pressures"_a, "mixture_ratios"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("solve_adiabatic",
             nb::overload_cast<double, double,
                               const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               const Goddard::CombustorOptions&>(
                 &Goddard::Combustor::solve),
             "fuel_temperature"_a, "oxidizer_temperature"_a,
             "pressures"_a, "mixture_ratios"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("generate_mole_fraction_matrix",
             &Goddard::Combustor::generate_mole_fraction_matrix,
             "mixture_ratios"_a, "type"_a = Goddard::MixtureRatioType::OF_RATIO)
        .def("generate_mass_fraction_matrix",
             &Goddard::Combustor::generate_mass_fraction_matrix,
             "mixture_ratios"_a, "type"_a = Goddard::MixtureRatioType::OF_RATIO);

    nb::class_<Goddard::DilutedCombustor, Goddard::BaseCombustor>(m, "DilutedCombustor")
        .def(nb::init<Goddard::Gas, const Goddard::Composition&, const Goddard::Composition&,
                       const Goddard::Composition&>(),
             "gas"_a, "fuel"_a, "oxidizer"_a, "flue"_a)
        .def("solve",
             nb::overload_cast<const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               double,
                               const Goddard::CombustorOptions&>(
                 &Goddard::DilutedCombustor::solve),
             "temperatures"_a, "pressures"_a, "mixture_ratios"_a,
             "recirculation_ratio"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("solve_adiabatic",
             nb::overload_cast<double, double, double,
                               const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               double,
                               const Goddard::CombustorOptions&>(
                 &Goddard::DilutedCombustor::solve),
             "fuel_temperature"_a, "oxidizer_temperature"_a, "flue_temperature"_a,
             "pressures"_a, "mixture_ratios"_a,
             "recirculation_ratio"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("generate_mole_fraction_matrix",
             &Goddard::DilutedCombustor::generate_mole_fraction_matrix,
             "mixture_ratios"_a, "type"_a = Goddard::MixtureRatioType::OF_RATIO,
             "recirculation_ratio"_a = 0.0)
        .def("generate_mass_fraction_matrix",
             &Goddard::DilutedCombustor::generate_mass_fraction_matrix,
             "mixture_ratios"_a, "type"_a = Goddard::MixtureRatioType::OF_RATIO,
             "recirculation_ratio"_a = 0.0);
}
