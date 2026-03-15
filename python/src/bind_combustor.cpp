#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/eigen/dense.h>
#include "goddard/combustor.hpp"

namespace nb = nanobind;
using namespace nb::literals;

void bind_combustor(nb::module_& m) {
    nb::class_<Goddard::Combustor>(m, "Combustor")
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
        .def("get_combustion_species", &Goddard::Combustor::get_combustion_species)
        .def_rw("fuel_state", &Goddard::Combustor::fuel_state)
        .def_rw("oxidizer_state", &Goddard::Combustor::oxidizer_state);
}
