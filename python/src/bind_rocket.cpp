#include <nanobind/nanobind.h>
#include <format>
#include <nanobind/stl/map.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <optional>
#include <vector>

#include "goddard/rocket.hpp"
#include "bind_status.hpp"
#include "goddard_docstrings.h"

namespace nb = nanobind;
using namespace nb::literals;

void bind_rocket(nb::module_& m) {
    // RocketPerformance
    nb::class_<Goddard::RocketPerformance>(m, "RocketPerformance", DOC(Goddard, RocketPerformance))
        .def_ro("status", &Goddard::RocketPerformance::status,
             DOC(Goddard, RocketPerformance, status))
        .def_ro("pressure_ratio", &Goddard::RocketPerformance::pressure_ratio,
             DOC(Goddard, RocketPerformance, pressure_ratio))
        .def_ro("area_ratio", &Goddard::RocketPerformance::area_ratio,
             DOC(Goddard, RocketPerformance, area_ratio))
        .def_ro("mach_number", &Goddard::RocketPerformance::mach_number,
             DOC(Goddard, RocketPerformance, mach_number))
        .def_ro("cstar", &Goddard::RocketPerformance::cstar,
             DOC(Goddard, RocketPerformance, cstar))
        .def_ro("CF", &Goddard::RocketPerformance::CF, DOC(Goddard, RocketPerformance, CF))
        .def_ro("isp", &Goddard::RocketPerformance::isp, DOC(Goddard, RocketPerformance, isp))
        .def_ro("ivac", &Goddard::RocketPerformance::ivac, DOC(Goddard, RocketPerformance, ivac))
        .def_prop_ro("isp_s", &Goddard::RocketPerformance::isp_s,
             DOC(Goddard, RocketPerformance, isp_s))
        .def_prop_ro("ivac_s", &Goddard::RocketPerformance::ivac_s,
             DOC(Goddard, RocketPerformance, ivac_s))
        .def("__repr__", [](const Goddard::RocketPerformance& self) {
            return std::format("<RocketPerformance {} Pc/Pe={:.6g} Ae/At={:.6g} M={:.6g} "
                               "cstar={:.6g} m/s CF={:.6g} Isp={:.6g} m/s Ivac={:.6g} m/s>",
                               status_tag(self.status), self.pressure_ratio, self.area_ratio,
                               self.mach_number, self.cstar, self.CF, self.isp, self.ivac);
        });

    m.def("rocket_performance", &Goddard::rocket_performance,
          "stagnation"_a, "throat"_a, "exit"_a, DOC(Goddard, rocket_performance));

    // RocketPoint
    nb::class_<Goddard::RocketPoint>(m, "RocketPoint", DOC(Goddard, RocketPoint))
        .def("__init__", [](Goddard::RocketPoint* self,
                            const std::optional<Goddard::PhaseSpecification>& fuel,
                            const std::optional<Goddard::PhaseSpecification>& oxidizer,
                            double mixture_ratio,
                            double chamber_pressure,
                            const std::optional<Goddard::CombustorOptions>& combustor,
                            const std::optional<Goddard::NozzleOptions>& nozzle,
                            Goddard::ExpansionType expansion_type,
                            std::vector<double> expansion_ratios,
                            std::string label) {
            new (self) Goddard::RocketPoint();
            self->fuel = fuel.value_or(Goddard::PhaseSpecification());
            self->oxidizer = oxidizer.value_or(Goddard::PhaseSpecification());
            self->mixture_ratio = mixture_ratio;
            self->chamber_pressure = chamber_pressure;
            self->combustor = combustor.value_or(Goddard::CombustorOptions());
            self->nozzle = nozzle.value_or(Goddard::NozzleOptions());
            self->expansion_type = expansion_type;
            self->expansion_ratios = std::move(expansion_ratios);
            self->label = std::move(label);
        },  // Defaults of bound types are passed as None: a bound object held as a default
            // argument is reported as leaked when the interpreter exits.
            "fuel"_a = nb::none(),
            "oxidizer"_a = nb::none(),
            "mixture_ratio"_a = Goddard::NaN,
            "chamber_pressure"_a = Goddard::NaN,
            "combustor"_a = nb::none(),
            "nozzle"_a = nb::none(),
            "expansion_type"_a = Goddard::ExpansionType::SUPERSONIC_AREA_RATIO,
            "expansion_ratios"_a = std::vector<double>(),
            "label"_a = "",
            "Create a point; each keyword sets the field of the same name.")
        .def_rw("fuel", &Goddard::RocketPoint::fuel, DOC(Goddard, RocketPoint, fuel))
        .def_rw("oxidizer", &Goddard::RocketPoint::oxidizer, DOC(Goddard, RocketPoint, oxidizer))
        .def_rw("mixture_ratio", &Goddard::RocketPoint::mixture_ratio,
             DOC(Goddard, RocketPoint, mixture_ratio))
        .def_rw("chamber_pressure", &Goddard::RocketPoint::chamber_pressure,
             DOC(Goddard, RocketPoint, chamber_pressure))
        .def_rw("combustor", &Goddard::RocketPoint::combustor,
             DOC(Goddard, RocketPoint, combustor))
        .def_rw("nozzle", &Goddard::RocketPoint::nozzle, DOC(Goddard, RocketPoint, nozzle))
        .def_rw("expansion_type", &Goddard::RocketPoint::expansion_type,
             DOC(Goddard, RocketPoint, expansion_type))
        .def_rw("expansion_ratios", &Goddard::RocketPoint::expansion_ratios,
             DOC(Goddard, RocketPoint, expansion_ratios))
        .def_rw("label", &Goddard::RocketPoint::label, DOC(Goddard, RocketPoint, label))
        .def("__repr__", [](const Goddard::RocketPoint& self) {
            return std::format("<RocketPoint '{}' mixture_ratio={:.6g} chamber_pressure={:.6g} Pa "
                               "exits={}>", self.label, self.mixture_ratio, self.chamber_pressure,
                               self.expansion_ratios.size());
        });

    // RocketResult
    nb::class_<Goddard::RocketResult>(m, "RocketResult", DOC(Goddard, RocketResult))
        .def_ro("point", &Goddard::RocketResult::point, DOC(Goddard, RocketResult, point))
        .def_ro("chamber", &Goddard::RocketResult::chamber, DOC(Goddard, RocketResult, chamber))
        .def_ro("stagnation", &Goddard::RocketResult::stagnation,
             DOC(Goddard, RocketResult, stagnation))
        .def_ro("combustion_end", &Goddard::RocketResult::combustion_end,
             DOC(Goddard, RocketResult, combustion_end))
        .def_ro("throat", &Goddard::RocketResult::throat, DOC(Goddard, RocketResult, throat))
        .def_ro("exits", &Goddard::RocketResult::exits, DOC(Goddard, RocketResult, exits))
        .def_ro("performance", &Goddard::RocketResult::performance,
             DOC(Goddard, RocketResult, performance))
        .def_ro("contraction_ratio", &Goddard::RocketResult::contraction_ratio,
             DOC(Goddard, RocketResult, contraction_ratio))
        .def_ro("mass_flux", &Goddard::RocketResult::mass_flux,
             DOC(Goddard, RocketResult, mass_flux))
        .def_ro("status", &Goddard::RocketResult::status, DOC(Goddard, RocketResult, status))
        .def("report", [](const Goddard::RocketResult& self) { return Goddard::rocket_report(self); },
             DOC(Goddard, rocket_report))
        .def("__repr__", [](const Goddard::RocketResult& self) {
            return std::format("<RocketResult {} '{}' mixture_ratio={:.6g} Tc={:.6g} K exits={}>",
                               status_tag(self.status), self.point.label, self.point.mixture_ratio,
                               self.chamber.thermo.temperature, self.exits.size());
        });

    // RocketSolver
    nb::class_<Goddard::RocketSolver>(m, "RocketSolver", DOC(Goddard, RocketSolver))
        .def(nb::init<const Goddard::Gas&, std::string>(),
             "products"_a, "reactant_file"_a = "", DOC(Goddard, RocketSolver, RocketSolver))
        .def("solve", nb::overload_cast<const Goddard::RocketPoint&>(&Goddard::RocketSolver::solve),
             "point"_a, DOC(Goddard, RocketSolver, solve))
        // Any iterable of points, a generator for example, not only a list.
        .def("solve", [](Goddard::RocketSolver& self, const nb::iterable& points) {
                std::vector<Goddard::RocketPoint> list;
                for (nb::handle point : points) {
                    list.push_back(nb::cast<Goddard::RocketPoint>(point));
                }
                return self.solve(list);
            }, "points"_a, DOC(Goddard, RocketSolver, solve, 2))
        .def("validate", &Goddard::RocketSolver::validate, "point"_a,
             DOC(Goddard, RocketSolver, validate))
        .def_rw("include_transport", &Goddard::RocketSolver::include_transport,
             DOC(Goddard, RocketSolver, include_transport))
        .def_rw("include_ionized_species", &Goddard::RocketSolver::include_ionized_species,
             DOC(Goddard, RocketSolver, include_ionized_species))
        .def_rw("trace_cutoff", &Goddard::RocketSolver::trace_cutoff,
             DOC(Goddard, RocketSolver, trace_cutoff));

    m.def("rocket_report",
          nb::overload_cast<const std::vector<Goddard::RocketResult>&>(&Goddard::rocket_report),
          "results"_a, DOC(Goddard, rocket_report, 2));
    m.def("rocket_report", nb::overload_cast<const Goddard::RocketResult&>(&Goddard::rocket_report),
          "result"_a, DOC(Goddard, rocket_report));
}
