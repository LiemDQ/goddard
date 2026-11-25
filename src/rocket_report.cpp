#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "goddard/format.hpp"
#include "goddard/rocket.hpp"

namespace Goddard {

namespace {

constexpr double PA_TO_BAR  = 1e-5;
constexpr double PA_TO_PSIA = 1.0 / 6894.757;
constexpr double J_TO_KJ    = 1e-3;

/** Extra header and column values needed only for a finite-area combustor page. */
struct FacReportInfo {
    double injector_pressure_pa;   ///< Pinj [Pa].
    double stagnation_pressure_pa; ///< Pinf [Pa].
    double mass_flux;              ///< mdot/Ac [kg/(m^2 s)].
    double contraction_ratio;      ///< Ac/At [-].
};

std::string mixture_label(MixtureRatioType type) {
    switch (type) {
        case MixtureRatioType::OF_RATIO: return "O/F=";
        case MixtureRatioType::FUEL_FRAC: return "FUEL FRACTION=";
        case MixtureRatioType::PHI_RATIO: return "PHI=";
    }
    return "O/F=";
}

std::string build_report_page(
    GasChemistry chemistry,
    CombustionProcess process,
    const std::vector<ThermodynamicState>& states,
    std::size_t throat_col,
    const ThermodynamicState& stagnation,
    MixtureRatioType mixture_type,
    double mixture_ratio,
    double display_pressure_pa,
    double initial_pressure_pa,
    int frozen_NFZ,
    const FacReportInfo* fac = nullptr)
{
    std::string page;

    std::string combustor_description;
    switch (process) {
        case CombustionProcess::ISOBARIC:
            combustor_description = (fac != nullptr) ? "FINITE AREA COMBUSTOR" : "INFINITE AREA COMBUSTOR";
            break;
        case CombustionProcess::ISOCHORIC:
            combustor_description = "CONSTANT-VOLUME COMBUSTOR";
            break;
    }

    switch (chemistry) {
        case GasChemistry::EQUILIBRIUM:
            page += "         THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM\n\n";
            page += "      COMPOSITION DURING EXPANSION FROM " + combustor_description + "\n\n";
            break;
        case GasChemistry::FROZEN: {
            // Nozzle station numbering: 0 is the chamber (the combustion end for a finite-area
            // combustor), 1 the throat, k >= 2 the (k-1)th exit.
            std::string freezing_station;
            if (frozen_NFZ <= 0) {
                freezing_station = (fac != nullptr) ? "THE COMBUSTION END" : "THE CHAMBER";
            } else if (frozen_NFZ == 1) {
                freezing_station = "THE THROAT";
            } else {
                freezing_station = "EXIT " + std::to_string(frozen_NFZ - 1);
            }
            page += "         THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION\n";
            page += "                      FROZEN AT " + freezing_station + "\n\n";
            if (process == CombustionProcess::ISOCHORIC) {
                page += "      EXPANSION FROM " + combustor_description + "\n\n";
            }
            break;
        }
        default: break;
    }

    page += "Pin = " + format_fixed(display_pressure_pa * PA_TO_PSIA, 7, 1) + " PSIA\n";
    if (process == CombustionProcess::ISOCHORIC) {
        page += "Pinitial = " + format_fixed(initial_pressure_pa * PA_TO_PSIA, 7, 1) + " PSIA\n";
    }
    if (fac != nullptr) {
        page += "MDOT/Ac = " + format_fixed(fac->mass_flux, 9, 3) + " (KG/S)/M**2      Pinj/Pinf = "
              + format_fixed(fac->injector_pressure_pa / fac->stagnation_pressure_pa, 9, 6) + "\n";
        page += "Ac/At = " + format_fixed(fac->contraction_ratio, 9, 4) + "\n";
    }
    page += mixture_label(mixture_type) + format_fixed(mixture_ratio, 11, 5) + "\n\n";

    TextTable table(18, 10);
    std::vector<std::string> headers;
    headers.reserve(states.size());
    if (fac != nullptr) {
        headers.push_back("INJECTOR");
        headers.push_back("COMB END");
    } else {
        headers.push_back("CHAMBER");
    }
    headers.push_back("THROAT");
    for (std::size_t i = throat_col + 1; i < states.size(); i++) {
        headers.push_back("EXIT");
    }
    table.set_headers(headers);

    auto row_vals = [&](auto fn) {
        std::vector<std::string> vals;
        vals.reserve(states.size());
        for (const auto& s : states) {
            vals.push_back(fn(s));
        }
        return vals;
    };

    if (fac != nullptr) {
        table.add_row("Pinj/P", row_vals([&](const ThermodynamicState& s){ return format_fixed(fac->injector_pressure_pa / s.pressure, 10, 4); }));
        table.add_row("Pinf/P", row_vals([&](const ThermodynamicState& s){ return format_fixed(stagnation.pressure / s.pressure, 10, 4); }));
    } else {
        table.add_row("Pinf/P", row_vals([&](const ThermodynamicState& s){ return format_fixed(states[0].pressure / s.pressure, 10, 4); }));
    }
    table.add_row("P, BAR",       row_vals([](const ThermodynamicState& s){ return format_fixed(s.pressure * PA_TO_BAR, 10, 4); }));
    table.add_row("T, K",         row_vals([](const ThermodynamicState& s){ return format_fixed(s.temperature, 10, 2); }));
    table.add_row("RHO, KG/CU M", row_vals([](const ThermodynamicState& s){ return format_cea_engineering(s.density, 10); }));
    table.add_row("H, KJ/KG",     row_vals([](const ThermodynamicState& s){ return format_fixed(s.enthalpy * J_TO_KJ, 10, 2); }));
    table.add_row("U, KJ/KG",     row_vals([](const ThermodynamicState& s){ return format_fixed(s.internal_energy * J_TO_KJ, 10, 2); }));
    table.add_row("G, KJ/KG",     row_vals([](const ThermodynamicState& s){ return format_fixed(s.gibbs * J_TO_KJ, 10, 1); }));
    table.add_row("S, KJ/(KG)(K)",row_vals([](const ThermodynamicState& s){ return format_fixed(s.entropy * J_TO_KJ, 10, 4); }));

    table.add_blank_line();
    table.add_row("M, (1/n)", row_vals([](const ThermodynamicState& s){ return format_fixed(s.molecular_weight, 10, 3); }));

    // CEA reports the mixture molecular weight, which counts the condensed moles as well, only
    // when a condensed phase is present. "M, (1/n)" above counts the gas alone.
    const bool has_condensed = std::any_of(states.begin(), states.end(),
        [](const ThermodynamicState& s){ return s.gas_mass_fraction < 1.0; });
    if (has_condensed) {
        table.add_row("MW, MOL WT", row_vals([](const ThermodynamicState& s){
            return format_fixed(s.mixture_molecular_weight, 10, 3); }));
    }

    if (chemistry == GasChemistry::EQUILIBRIUM) {
        table.add_row("(dLV/dLP)t", row_vals([](const ThermodynamicState& s){ return format_fixed(s.dlV_dlP_T, 10, 5); }));
        table.add_row("(dLV/dLT)p", row_vals([](const ThermodynamicState& s){ return format_fixed(s.dlV_dlT_P, 10, 4); }));
    }

    // The equilibrium specific heat is infinite at a pinned phase transition; CEA prints zero.
    table.add_row("Cp, KJ/(KG)(K)", row_vals([](const ThermodynamicState& s){
        return format_fixed(s.pinned_transition ? 0.0 : s.cp * J_TO_KJ, 10, 4); }));
    table.add_row("GAMMAs",          row_vals([](const ThermodynamicState& s){ return format_fixed(s.gamma_s, 10, 4); }));
    table.add_row("SON VEL,M/SEC",   row_vals([](const ThermodynamicState& s){ return format_fixed(s.speed_of_sound, 10, 1); }));

    // Mach number: 0 at the first column, 1 at the throat, computed for every other column.
    {
        std::vector<std::string> vals;
        for (std::size_t i = 0; i < states.size(); i++) {
            double mach = 0.0;
            if (i == throat_col) {
                mach = 1.0;
            } else if (i != 0) {
                mach = rocket_performance(stagnation, states[throat_col], states[i]).mach_number;
            }
            vals.push_back(format_fixed(mach, 10, 3));
        }
        table.add_row("MACH NUMBER", vals);
    }

    table.add_blank_line();
    table.add_section_header("PERFORMANCE PARAMETERS");
    table.add_blank_line();

    // Every performance row is blank in the first column.
    auto performance_row = [&](const std::string& label, auto value) {
        std::vector<std::string> vals;
        vals.push_back("");
        for (std::size_t i = 1; i < states.size(); i++) {
            vals.push_back(value(i, rocket_performance(stagnation, states[throat_col], states[i])));
        }
        table.add_row(label, vals);
    };
    performance_row("Ae/At", [&](std::size_t i, const RocketPerformance& perf) {
        return i == throat_col ? format_fixed(1.0, 10, 4) : format_fixed(perf.area_ratio, 10, 3);
    });
    // c* is a property of the throat, the same in every column.
    const double cstar = rocket_performance(stagnation, states[throat_col], states[throat_col]).cstar;
    performance_row("CSTAR, M/SEC", [&](std::size_t, const RocketPerformance&) {
        return format_fixed(cstar, 10, 1);
    });
    performance_row("CF", [](std::size_t, const RocketPerformance& perf) {
        return format_fixed(perf.CF, 10, 4);
    });
    performance_row("Ivac, M/SEC", [](std::size_t, const RocketPerformance& perf) {
        return format_fixed(perf.ivac, 10, 1);
    });
    performance_row("Isp, M/SEC", [](std::size_t, const RocketPerformance& perf) {
        return format_fixed(perf.isp, 10, 1);
    });

    table.add_blank_line();
    table.add_blank_line();
    table.add_section_header("MASS FRACTIONS");
    table.add_blank_line();

    if (chemistry == GasChemistry::EQUILIBRIUM) {
        std::set<std::string> all_species;
        for (const auto& s : states) {
            for (const auto& [name, frac] : s.composition) {
                if (frac > 1e-10) all_species.insert(name);
            }
        }
        for (const auto& species : all_species) {
            std::vector<std::string> vals;
            for (const auto& s : states) {
                auto it = s.composition.find(species);
                double val = (it != s.composition.end()) ? it->second : 0.0;
                vals.push_back(format_mass_fraction(val, 10));
            }
            table.add_row(species, vals);
        }
    } else {
        // The composition the flow keeps downstream of the freezing station, which the last
        // station always has.
        const auto& frozen_comp = states.back().composition;
        std::vector<std::pair<std::string, double>> species_list;
        for (const auto& [name, frac] : frozen_comp) {
            if (frac > 1e-10) species_list.push_back({name, frac});
        }
        std::sort(species_list.begin(), species_list.end());

        constexpr int entries_per_row = 3;
        constexpr int name_width = 16;
        constexpr int frac_width = 10;
        for (std::size_t i = 0; i < species_list.size(); i += entries_per_row) {
            std::string line;
            for (std::size_t j = i; j < i + entries_per_row && j < species_list.size(); j++) {
                std::string name = species_list[j].first;
                if (static_cast<int>(name.size()) < name_width)
                    name += std::string(name_width - name.size(), ' ');
                line += name + format_mass_fraction(species_list[j].second, frac_width);
                if (j + 1 < i + entries_per_row && j + 1 < species_list.size()) {
                    line += "   ";
                }
            }
            table.add_section_header(line);
        }
    }

    page += table.render();
    page += "\n\n";
    return page;
}

} // namespace

std::string rocket_report(const RocketResult& result) {
    const RocketPoint& point = result.point;
    const bool is_fac = point.combustor.type != CombustorType::INFINITE_AREA;

    std::string page;
    if (!point.label.empty()) {
        page += "CASE = " + point.label + "\n\n";
    }

    // Stations that were not solved, and whether one up to the throat is among them: the table
    // then has nothing to show.
    std::string failures;
    bool upstream_failed = false;
    auto note_failure = [&](const SolveStatus& status, const std::string& label, bool upstream) {
        if (status.ok()) {
            return;
        }
        upstream_failed = upstream_failed || upstream;
        failures += "   " + label + ": " + status.describe() + "\n";
    };
    note_failure(result.chamber.status, is_fac ? "INJECTOR" : "CHAMBER", true);
    if (is_fac) {
        note_failure(result.stagnation.status, "STAGNATION", true);
        note_failure(result.combustion_end.status, "COMB END", true);
    }
    note_failure(result.throat.status, "THROAT", true);
    for (std::size_t i = 0; i < result.exits.size(); i++) {
        note_failure(result.exits[i].status, "EXIT " + std::to_string(i + 1), false);
    }

    if (upstream_failed) {
        page += "         OPERATING POINT NOT SOLVED\n\n";
        page += "Pin = " + format_fixed(point.chamber_pressure * PA_TO_PSIA, 7, 1) + " PSIA\n";
        page += mixture_label(point.combustor.mixture_type)
              + format_fixed(point.mixture_ratio, 11, 5) + "\n\n";
        page += failures + "\n\n";
        return page;
    }

    // Columns: chamber or injector, [combustion end, for a finite-area combustor], throat, exits.
    std::vector<ThermodynamicState> states;
    states.push_back(result.chamber.thermo);
    if (is_fac) {
        states.push_back(result.combustion_end.thermo);
    }
    states.push_back(result.throat.thermo);
    const std::size_t throat_col = states.size() - 1;
    for (const NozzleStation& exit : result.exits) {
        states.push_back(exit.thermo);
    }

    FacReportInfo fac_info{};
    const FacReportInfo* fac = nullptr;
    if (is_fac) {
        fac_info.injector_pressure_pa = result.chamber.thermo.pressure;
        fac_info.stagnation_pressure_pa = result.stagnation.thermo.pressure;
        fac_info.contraction_ratio = result.contraction_ratio;
        fac_info.mass_flux = result.mass_flux;
        fac = &fac_info;
    }

    page += build_report_page(
        point.nozzle.chemistry,
        point.combustor.process,
        states,
        throat_col,
        result.stagnation.thermo,
        point.combustor.mixture_type,
        point.mixture_ratio,
        states[0].pressure,
        point.chamber_pressure,
        point.nozzle.frozen_NFZ,
        fac);
    if (!failures.empty()) {
        page += " STATIONS NOT SOLVED (NaN ABOVE)\n" + failures + "\n\n";
    }
    return page;
}

std::string rocket_report(const std::vector<RocketResult>& results) {
    std::string report;
    for (const RocketResult& result : results) {
        report += rocket_report(result);
    }
    return report;
}

} // namespace Goddard
