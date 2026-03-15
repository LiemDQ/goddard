#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/unordered_set.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include "goddard/problem.hpp"
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

    // RocketProblemCaseResult
    nb::class_<Goddard::RocketProblemCaseResult>(m, "RocketProblemCaseResult")
        .def_ro("problem_type", &Goddard::RocketProblemCaseResult::problem_type)
        .def_ro("inlet_states", &Goddard::RocketProblemCaseResult::inlet_states)
        .def_ro("chemistry", &Goddard::RocketProblemCaseResult::chemistry)
        .def_ro("nozzle_states", &Goddard::RocketProblemCaseResult::nozzle_states);

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
        .def("extract_thermo_info", &Goddard::RocketProblemResults::extract_thermo_info,
             "case_name"_a, "index"_a)
        .def("get_chamber_state", &Goddard::RocketProblemResults::get_chamber_state,
             "case_name"_a, "of_index"_a)
        .def("get_throat_state", &Goddard::RocketProblemResults::get_throat_state,
             "case_name"_a, "of_index"_a)
        .def("get_exit_states", &Goddard::RocketProblemResults::get_exit_states,
             "case_name"_a, "of_index"_a)
        .def_static("calculate_performance",
             &Goddard::RocketProblemResults::calculate_performance,
             "chamber"_a, "throat"_a, "exit"_a)
        .def("get_case", [](Goddard::RocketProblemResults& self, const std::string& name)
                -> Goddard::RocketProblemCaseResult& {
            return self.cases.at(name);
        }, nb::rv_policy::reference_internal, "name"_a,
           "Get case result by name. Raises KeyError if not found.")
        .def("case_names", [](Goddard::RocketProblemResults& self) {
            std::vector<std::string> names;
            names.reserve(self.cases.size());
            for (auto& [k, v] : self.cases) {
                names.push_back(k);
            }
            return names;
        }, "Get list of case names.");
}
