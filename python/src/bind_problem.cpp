#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/unordered_set.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/problem.hpp"
#include "goddard/rocket_results.hpp"
#include "goddard/speciate.hpp"
#include "cantera/core.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_problem(nb::module_& m) {
    // Opaque handle to Cantera::Solution
    nb::class_<Cantera::Solution>(m, "SolutionHandle",
        "Opaque handle to an internal Cantera Solution object. "
        "Obtain via RocketProblem.solution() or create_solution().");

    // Factory function to create a Cantera Solution from a YAML file
    m.def("create_solution", [](const std::string& yaml_file,
                                const std::string& phase_name,
                                const std::unordered_set<std::string>& species) {
        auto root_node = Goddard::select_species(yaml_file, species);
        const Cantera::AnyMap& phase_node = root_node.at("phases").getMapWhere("name", phase_name);
        auto sln = Cantera::newSolution(phase_node, root_node);
        sln->setSource(yaml_file);
        return sln;
    }, "yaml_file"_a, "phase_name"_a = "", "species"_a = std::unordered_set<std::string>{},
       "Create a Cantera Solution handle from a YAML thermodynamic data file.");

    // RocketProblem
    nb::class_<Goddard::RocketProblem>(m, "RocketProblem")
        .def(nb::init<const Goddard::ChemicalParameters&,
                       const std::vector<Goddard::RocketCaseParameters>&,
                       const std::string&, bool, bool, double>(),
             "chem_params"_a, "cases"_a,
             "name"_a = "", "transport"_a = false,
             "ionized_species"_a = false, "trace"_a = 1e-6)
        .def("solve", &Goddard::RocketProblem::solve)
        .def("solution", &Goddard::RocketProblem::solution)
        .def_rw("include_transport", &Goddard::RocketProblem::include_transport)
        .def_rw("include_ionized_species", &Goddard::RocketProblem::include_ionized_species)
        .def_rw("trace_cutoff", &Goddard::RocketProblem::trace_cutoff)
        .def_rw("problem_cases", &Goddard::RocketProblem::problem_cases)
        .def_rw("chemical_params", &Goddard::RocketProblem::chemical_params);

    // RocketProblemResults
    nb::class_<Goddard::RocketProblemResults>(m, "RocketProblemResults")
        .def_prop_ro("stations", &Goddard::RocketProblemResults::stations,
             "Flat list of all RocketStation objects across all cases.")
        .def("stations_of_type", &Goddard::RocketProblemResults::stations_of_type,
             "type"_a, "case_name"_a = "",
             "Return all stations of a given StationType, optionally filtered to a case.")
        .def("chamber", &Goddard::RocketProblemResults::chamber,
             "of_index"_a = 0, "case_name"_a = "",
             "Chamber station for the given O/F index. "
             "case_name may be omitted when there is only one case.",
             nb::rv_policy::reference_internal)
        .def("throat", &Goddard::RocketProblemResults::throat,
             "of_index"_a = 0, "case_name"_a = "",
             "Throat station for the given O/F index. "
             "case_name may be omitted when there is only one case.",
             nb::rv_policy::reference_internal)
        .def("exits", &Goddard::RocketProblemResults::exits,
             "of_index"_a = 0, "case_name"_a = "",
             "Exit stations for the given O/F index, ordered by expansion ratio. "
             "case_name may be omitted when there is only one case.")
        .def("performance", &Goddard::RocketProblemResults::performance,
             "of_index"_a = 0, "exit_index"_a = 0, "case_name"_a = "",
             "Compute rocket performance for the given O/F index and exit station.")
        .def("case_names", &Goddard::RocketProblemResults::case_names,
             "List all case names.")
        .def("of_ratios", &Goddard::RocketProblemResults::of_ratios,
             "case_name"_a = "",
             "O/F ratios used for a case.",
             nb::rv_policy::reference_internal)
        .def("report", &Goddard::RocketProblemResults::report,
             "case_name"_a = "",
             "Generate a CEA-style formatted text report. "
             "If case_name is empty, all cases are reported.")
        .def_static("calculate_performance",
             &Goddard::RocketProblemResults::calculate_performance,
             "chamber"_a, "throat"_a, "exit"_a,
             "Compute performance metrics from individual ThermoStateInfo objects.");
}
