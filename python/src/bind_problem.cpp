#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/optional.h>
#include "goddard/problem.hpp"
#include "goddard/rocket_results.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_problem(nb::module_& m) {
    // RocketProblem
    nb::class_<Goddard::RocketProblem>(m, "RocketProblem", DOC(Goddard, RocketProblem))
        .def(nb::init<const Goddard::ChemicalParameters&,
                       const std::vector<Goddard::RocketCaseParameters>&,
                       const std::string&, bool, bool, double>(),
             "chem_params"_a, "cases"_a,
             "phase_name"_a = "", "transport"_a = false,
             "ionized_species"_a = false, "trace"_a = 0.0,
             DOC(Goddard, RocketProblem, RocketProblem))
        .def("solve", &Goddard::RocketProblem::solve, DOC(Goddard, RocketProblem, solve))
        .def_rw("include_transport", &Goddard::RocketProblem::include_transport,
             DOC(Goddard, RocketProblem, include_transport))
        .def_rw("include_ionized_species", &Goddard::RocketProblem::include_ionized_species,
             DOC(Goddard, RocketProblem, include_ionized_species))
        .def_rw("trace_cutoff", &Goddard::RocketProblem::trace_cutoff,
             DOC(Goddard, RocketProblem, trace_cutoff))
        .def_rw("problem_cases", &Goddard::RocketProblem::problem_cases,
             DOC(Goddard, RocketProblem, problem_cases))
        .def_rw("chemical_params", &Goddard::RocketProblem::chemical_params,
             DOC(Goddard, RocketProblem, chemical_params));

    // RocketProblemResults
    nb::class_<Goddard::RocketProblemResults>(m, "RocketProblemResults",
            DOC(Goddard, RocketProblemResults))
        .def_prop_ro("stations", &Goddard::RocketProblemResults::stations,
             DOC(Goddard, RocketProblemResults, stations))
        .def("stations_of_type", &Goddard::RocketProblemResults::stations_of_type,
             "type"_a, "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, stations_of_type))
        // Everything after of_index is keyword-only: pressure_index was inserted before
        // exit_index and case_name, so an old positional call must fail instead of
        // silently reading another operating point.
        .def("chamber", &Goddard::RocketProblemResults::chamber,
             "of_index"_a = 0, nb::kw_only(), "pressure_index"_a = 0, "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, chamber),
             nb::rv_policy::reference_internal)
        .def("throat", &Goddard::RocketProblemResults::throat,
             "of_index"_a = 0, nb::kw_only(), "pressure_index"_a = 0, "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, throat),
             nb::rv_policy::reference_internal)
        .def("stagnation", &Goddard::RocketProblemResults::stagnation,
             "of_index"_a = 0, nb::kw_only(), "pressure_index"_a = 0, "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, stagnation),
             nb::rv_policy::reference_internal)
        .def("combustion_end", &Goddard::RocketProblemResults::combustion_end,
             "of_index"_a = 0, nb::kw_only(), "pressure_index"_a = 0, "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, combustion_end),
             nb::rv_policy::reference_internal)
        .def("exits", &Goddard::RocketProblemResults::exits,
             "of_index"_a = 0, nb::kw_only(), "pressure_index"_a = 0, "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, exits))
        .def("performance", &Goddard::RocketProblemResults::performance,
             "of_index"_a = 0, nb::kw_only(), "pressure_index"_a = 0, "exit_index"_a = 0,
             "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, performance))
        .def("case_names", &Goddard::RocketProblemResults::case_names,
             DOC(Goddard, RocketProblemResults, case_names))
        .def("of_ratios", &Goddard::RocketProblemResults::of_ratios,
             "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, of_ratios),
             nb::rv_policy::reference_internal)
        .def("report", &Goddard::RocketProblemResults::report,
             "case_name"_a = "",
             DOC(Goddard, RocketProblemResults, report))
        .def_static("calculate_performance",
             &Goddard::RocketProblemResults::calculate_performance,
             "chamber"_a, "throat"_a, "exit"_a,
             DOC(Goddard, RocketProblemResults, calculate_performance));
}
