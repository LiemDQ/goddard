#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/eigen/dense.h>
#include "goddard/combustor.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_combustor(nb::module_& m) {
    nb::class_<Goddard::BaseCombustor>(m, "BaseCombustor")
        .def("get_combustion_species", &Goddard::BaseCombustor::get_combustion_species);

    nb::class_<Goddard::Combustor, Goddard::BaseCombustor>(m, "Combustor")
        .def(nb::init<std::shared_ptr<Cantera::Solution>,
                       std::vector<double>&,
                       std::vector<double>&>(),
             "solution"_a, "fuel_state"_a, "oxidizer_state"_a)
        .def("solve",
             nb::overload_cast<const Eigen::ArrayXd&,
                               const Goddard::MixtureRatios&,
                               const Goddard::CombustorOptions&>(
                 &Goddard::Combustor::solve),
             "pressures"_a, "mixture_ratios"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("solve_with_temperatures",
             nb::overload_cast<const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               const Goddard::MixtureRatios&,
                               const Goddard::CombustorOptions&>(
                 &Goddard::Combustor::solve),
             "temperatures"_a, "pressures"_a, "mixture_ratios"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("generate_mole_fraction_matrix",
             &Goddard::Combustor::generate_mole_fraction_matrix,
             "mixture_ratios"_a)
        .def("generate_mass_fraction_matrix",
             &Goddard::Combustor::generate_mass_fraction_matrix,
             "mixture_ratios"_a)
        .def_rw("fuel_state", &Goddard::Combustor::fuel_state)
        .def_rw("oxidizer_state", &Goddard::Combustor::oxidizer_state);

    nb::class_<Goddard::DilutedCombustor, Goddard::BaseCombustor>(m, "DilutedCombustor")
        .def(nb::init<std::shared_ptr<Cantera::Solution>,
                       std::vector<double>&,
                       std::vector<double>&,
                       std::vector<double>&>(),
             "solution"_a, "fuel_state"_a, "oxidizer_state"_a,
             "flue_state"_a)
        .def("solve",
             nb::overload_cast<const Eigen::ArrayXd&,
                               const Goddard::MixtureRatios&,
                               double,
                               const Goddard::CombustorOptions&>(
                 &Goddard::DilutedCombustor::solve),
             "pressures"_a, "mixture_ratios"_a, "recirculation_ratio"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("solve_with_temperatures",
             nb::overload_cast<const Eigen::ArrayXd&,
                               const Eigen::ArrayXd&,
                               const Goddard::MixtureRatios&,
                               double,
                               const Goddard::CombustorOptions&>(
                 &Goddard::DilutedCombustor::solve),
             "temperatures"_a, "pressures"_a, "mixture_ratios"_a,
             "recirculation_ratio"_a,
             "options"_a = Goddard::CombustorOptions{})
        .def("generate_mole_fraction_matrix",
             &Goddard::DilutedCombustor::generate_mole_fraction_matrix,
             "mixture_ratios"_a, "recirculation_ratio"_a)
        .def("generate_mass_fraction_matrix",
             &Goddard::DilutedCombustor::generate_mass_fraction_matrix,
             "mixture_ratios"_a, "recirculation_ratio"_a)
        .def_rw("fuel_state", &Goddard::DilutedCombustor::fuel_state)
        .def_rw("oxidizer_state", &Goddard::DilutedCombustor::oxidizer_state)
        .def_rw("flue_state", &Goddard::DilutedCombustor::flue_state);
}
