

#include <algorithm>
#include <set>

#include "goddard/rocket_results.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/format.hpp"

namespace Goddard {


RocketProblemResults::RocketProblemResults(
    std::unordered_map<std::string, RocketProblemCaseResult>&& case_results,
    std::shared_ptr<Cantera::Solution> sln)
    : m_sln(std::move(sln))
{
    auto tmo = m_sln->thermo();
    const std::vector<std::string> species_names = tmo->speciesNames();
    const std::size_t n_species = tmo->nSpecies();

    // Extract species mass fractions from a raw Cantera state vector.
    auto make_composition = [&](const std::vector<double>& state) {
        std::unordered_map<std::string, double> comp;
        comp.reserve(n_species);
        for (std::size_t i = 0; i < n_species; i++) {
            comp[species_names[i]] = state[i + 2];
        }
        return comp;
    };

    // Restore Cantera state, then read all thermodynamic quantities.
    auto read_thermo = [&](const std::vector<double>& state,
                           double gamma, double dlV_dlP_T, double dlV_dlT_P) {
        tmo->restoreState(state);
        return ThermoStateInfo{
            tmo->pressure(),
            tmo->temperature(),
            tmo->density(),
            tmo->enthalpy_mass(),
            tmo->intEnergy_mass(),
            tmo->gibbs_mass(),
            tmo->entropy_mass(),
            tmo->meanMolecularWeight(),
            tmo->cp_mass(),
            gamma,
            dlV_dlP_T,
            dlV_dlT_P,
            gas_sonic_velocity(*tmo, gamma),
            0.0,
            make_composition(state)
        };
    };

    for (auto& [name, case_result] : case_results) {
        m_case_meta[name] = CaseMeta{
            case_result.OF_ratios,
            case_result.pressures,
            case_result.expansion_ratios,
            case_result.chemistry,
            case_result.expansion_type
        };

        const std::size_t N_of = case_result.OF_ratios.size();
        const std::size_t N_p  = case_result.pressures.size();

        for (std::size_t of_idx = 0; of_idx < N_of; of_idx++) {
            for (std::size_t p_idx = 0; p_idx < N_p; p_idx++) {
                const std::size_t state_idx = of_idx * N_p + p_idx;
                if (state_idx >= static_cast<std::size_t>(case_result.inlet_states.size())) {
                    continue;
                }

                std::vector<double> inlet_state = case_result.inlet_states.get_state(state_idx);
                const NozzleResults& nozzle = case_result.nozzle_states[state_idx];

                // Chamber gamma must be computed (not stored like throat/exit).
                double inlet_gamma, inlet_dlP, inlet_dlT;
                switch (case_result.chemistry) {
                    case GasChemistry::EQUILIBRIUM: {
                        tmo->restoreState(inlet_state);
                        auto props = get_thermo_equilibrium_properties(*tmo);
                        inlet_gamma = props.gamma_s;
                        inlet_dlP   = props.dlogV_dlogP_T;
                        inlet_dlT   = props.dlogV_dlogT_P;
                        break;
                    }
                    case GasChemistry::FROZEN: {
                        tmo->restoreState(inlet_state);
                        inlet_gamma = tmo->cp_mass() / tmo->cv_mass();
                        inlet_dlP   = -1.0;
                        inlet_dlT   =  1.0;
                        break;
                    }
                    default:
                        throw NotImplementedError("Chemistry type not implemented.");
                }

                // Chamber station
                {
                    RocketStation s;
                    s.case_name       = name;
                    s.type            = StationType::CHAMBER;
                    s.of_index        = of_idx;
                    s.pressure_index  = p_idx;
                    s.expansion_index = 0;
                    s.area_ratio      = 0.0;
                    s.converged       = true;
                    s.thermo          = read_thermo(inlet_state, inlet_gamma, inlet_dlP, inlet_dlT);
                    m_stations.push_back(std::move(s));
                }

                // Throat station — gamma and derivatives are stored in ThroatCondition.
                {
                    const ThroatCondition& tc = nozzle.throat;
                    RocketStation s;
                    s.case_name       = name;
                    s.type            = StationType::THROAT;
                    s.of_index        = of_idx;
                    s.pressure_index  = p_idx;
                    s.expansion_index = 0;
                    s.area_ratio      = 1.0;
                    s.converged       = tc.converged;
                    s.thermo          = read_thermo(tc.state, tc.gamma_s, tc.dlV_dlP_T, tc.dlV_dlT_P);
                    m_stations.push_back(std::move(s));
                }

                // Exit stations — gamma and derivatives are stored in NozzleResult.
                for (std::size_t exp_idx = 0; exp_idx < nozzle.expansions.size(); exp_idx++) {
                    const NozzleStation& exp = nozzle.expansions[exp_idx];
                    RocketStation s;
                    s.case_name       = name;
                    s.type            = StationType::EXIT;
                    s.of_index        = of_idx;
                    s.pressure_index  = p_idx;
                    s.expansion_index = exp_idx;
                    s.area_ratio      = (exp_idx < case_result.expansion_ratios.size())
                                        ? case_result.expansion_ratios[exp_idx] : 0.0;
                    s.converged       = exp.converged;
                    s.thermo          = read_thermo(exp.state, exp.gamma_s, exp.dlV_dlP_T, exp.dlV_dlT_P);
                    m_stations.push_back(std::move(s));
                }
            }
        }
    }
}

std::string RocketProblemResults::resolve_case(const std::string& case_name) const {
    if (!case_name.empty()) {
        if (m_case_meta.count(case_name) == 0) {
            throw std::runtime_error("Case '" + case_name + "' not found in results.");
        }
        return case_name;
    }
    if (m_case_meta.size() == 1) {
        return m_case_meta.begin()->first;
    }
    if (m_case_meta.empty()) {
        throw std::runtime_error("No cases in results.");
    }
    throw std::runtime_error("Multiple cases present; specify case_name explicitly.");
}

std::vector<RocketStation> RocketProblemResults::stations_of_type(
    StationType type, const std::string& case_name) const
{
    const std::string resolved = resolve_case(case_name);
    std::vector<RocketStation> result;
    for (const auto& s : m_stations) {
        if (s.case_name == resolved && s.type == type) {
            result.push_back(s);
        }
    }
    return result;
}

const RocketStation& RocketProblemResults::chamber(
    std::size_t of_index, const std::string& case_name) const
{
    const std::string resolved = resolve_case(case_name);
    for (const auto& s : m_stations) {
        if (s.case_name == resolved && s.type == StationType::CHAMBER && s.of_index == of_index) {
            return s;
        }
    }
    throw std::runtime_error("Chamber state not found for case '" + resolved +
                              "', of_index=" + std::to_string(of_index));
}

const RocketStation& RocketProblemResults::throat(
    std::size_t of_index, const std::string& case_name) const
{
    const std::string resolved = resolve_case(case_name);
    for (const auto& s : m_stations) {
        if (s.case_name == resolved && s.type == StationType::THROAT && s.of_index == of_index) {
            return s;
        }
    }
    throw std::runtime_error("Throat state not found for case '" + resolved +
                              "', of_index=" + std::to_string(of_index));
}

std::vector<RocketStation> RocketProblemResults::exits(
    std::size_t of_index, const std::string& case_name) const
{
    const std::string resolved = resolve_case(case_name);
    std::vector<RocketStation> result;
    for (const auto& s : m_stations) {
        if (s.case_name == resolved && s.type == StationType::EXIT && s.of_index == of_index) {
            result.push_back(s);
        }
    }
    std::sort(result.begin(), result.end(),
              [](const RocketStation& a, const RocketStation& b) {
                  return a.expansion_index < b.expansion_index;
              });
    return result;
}

RocketPerformance RocketProblemResults::performance(
    std::size_t of_index, std::size_t exit_index, const std::string& case_name) const
{
    const RocketStation& c = chamber(of_index, case_name);
    const RocketStation& t = throat(of_index, case_name);
    auto exit_stations = exits(of_index, case_name);
    if (exit_index >= exit_stations.size()) {
        throw std::runtime_error(
            "exit_index " + std::to_string(exit_index) +
            " out of range (" + std::to_string(exit_stations.size()) + " exits available)");
    }
    return calculate_performance(c.thermo, t.thermo, exit_stations[exit_index].thermo);
}

std::vector<std::string> RocketProblemResults::case_names() const {
    std::vector<std::string> names;
    names.reserve(m_case_meta.size());
    for (const auto& [name, _] : m_case_meta) {
        names.push_back(name);
    }
    return names;
}

const std::vector<double>& RocketProblemResults::of_ratios(const std::string& case_name) const {
    return m_case_meta.at(resolve_case(case_name)).of_ratios;
}

RocketPerformance RocketProblemResults::calculate_performance(
    const ThermoStateInfo& chamber,
    const ThermoStateInfo& throat,
    const ThermoStateInfo& exit) {

    // Pressure ratios
    double pressure_ratio = chamber.pressure / exit.pressure;
    double throat_pressure_ratio = chamber.pressure / throat.pressure;

    // Area ratio from mass flux continuity: (rho*v) is constant at throat and exit
    double area_ratio = (throat.density * throat.speed_of_sound) /
                        (exit.density * exit.speed_of_sound) *
                        (throat_pressure_ratio / pressure_ratio);

    // Characteristic velocity c* = P_c * A_t / m_dot
    double gamma = throat.gamma_s;
    double c_star = throat.speed_of_sound * std::sqrt(
        std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (gamma - 1.0)) / gamma
    );

    // Exit velocity from energy conservation: v_e = sqrt(2*(h_c - h_e))
    double exit_velocity = std::sqrt(2.0 * (chamber.enthalpy - exit.enthalpy));

    // Thrust coefficient: CF = v_e / c*  (matched nozzle approximation)
    double CF = exit_velocity / c_star;

    constexpr double g0 = 9.80665;
    double isp  = exit_velocity / g0;
    double ivac = isp + (exit.pressure / chamber.pressure) * area_ratio * c_star / g0;

    double mach_number = exit_velocity / exit.speed_of_sound;

    return RocketPerformance{pressure_ratio, area_ratio, mach_number, c_star, CF, isp, ivac};
}

namespace {

constexpr double PA_TO_BAR  = 1e-5;
constexpr double PA_TO_PSIA = 1.0 / 6894.757;
constexpr double J_TO_KJ    = 1e-3;

std::string build_report_page(
    GasChemistry chemistry,
    const std::vector<ThermoStateInfo>& states,
    double of_ratio,
    double chamber_pressure_pa)
{
    std::string page;

    switch (chemistry) {
        case GasChemistry::EQUILIBRIUM:
            page += "         THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM\n\n";
            page += "      COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR\n\n";
            break;
        case GasChemistry::FROZEN:
            page += "         THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION\n\n";
            break;
        default: break;
    }

    page += "Pin = " + format_fixed(chamber_pressure_pa * PA_TO_PSIA, 7, 1) + " PSIA\n";
    page += "O/F=" + format_fixed(of_ratio, 11, 5) + "\n\n";

    TextTable table(18, 10);
    std::vector<std::string> headers;
    headers.push_back("CHAMBER");
    headers.push_back("THROAT");
    for (std::size_t i = 2; i < states.size(); i++) {
        headers.push_back("EXIT");
    }
    table.set_headers(headers);

    const auto& ch = states[0];

    auto row_vals = [&](auto fn) {
        std::vector<std::string> vals;
        for (const auto& s : states) {
            vals.push_back(fn(s));
        }
        return vals;
    };

    table.add_row("Pinf/P",       row_vals([&](const ThermoStateInfo& s){ return format_fixed(ch.pressure / s.pressure, 10, 4); }));
    table.add_row("P, BAR",       row_vals([](const ThermoStateInfo& s){ return format_fixed(s.pressure * PA_TO_BAR, 10, 4); }));
    table.add_row("T, K",         row_vals([](const ThermoStateInfo& s){ return format_fixed(s.temperature, 10, 2); }));
    table.add_row("RHO, KG/CU M", row_vals([](const ThermoStateInfo& s){ return format_cea_engineering(s.density, 10); }));
    table.add_row("H, KJ/KG",     row_vals([](const ThermoStateInfo& s){ return format_fixed(s.enthalpy * J_TO_KJ, 10, 2); }));
    table.add_row("U, KJ/KG",     row_vals([](const ThermoStateInfo& s){ return format_fixed(s.internal_energy * J_TO_KJ, 10, 2); }));
    table.add_row("G, KJ/KG",     row_vals([](const ThermoStateInfo& s){ return format_fixed(s.gibbs * J_TO_KJ, 10, 1); }));
    table.add_row("S, KJ/(KG)(K)",row_vals([](const ThermoStateInfo& s){ return format_fixed(s.entropy * J_TO_KJ, 10, 4); }));

    table.add_blank_line();
    table.add_row("M, (1/n)", row_vals([](const ThermoStateInfo& s){ return format_fixed(s.molecular_weight, 10, 3); }));

    if (chemistry == GasChemistry::EQUILIBRIUM) {
        table.add_row("(dLV/dLP)t", row_vals([](const ThermoStateInfo& s){ return format_fixed(s.dlV_dlP_T, 10, 5); }));
        table.add_row("(dLV/dLT)p", row_vals([](const ThermoStateInfo& s){ return format_fixed(s.dlV_dlT_P, 10, 4); }));
    }

    table.add_row("Cp, KJ/(KG)(K)", row_vals([](const ThermoStateInfo& s){ return format_fixed(s.cp * J_TO_KJ, 10, 4); }));
    table.add_row("GAMMAs",          row_vals([](const ThermoStateInfo& s){ return format_fixed(s.gamma_s, 10, 4); }));
    table.add_row("SON VEL,M/SEC",   row_vals([](const ThermoStateInfo& s){ return format_fixed(s.speed_of_sound, 10, 1); }));

    // Mach number: 0 at chamber, 1 at throat, computed for exits
    {
        std::vector<std::string> vals;
        vals.push_back(format_fixed(0.0, 10, 3));
        vals.push_back(format_fixed(1.0, 10, 3));
        for (std::size_t i = 2; i < states.size(); i++) {
            auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
            vals.push_back(format_fixed(perf.mach_number, 10, 3));
        }
        table.add_row("MACH NUMBER", vals);
    }

    table.add_blank_line();
    table.add_section_header("PERFORMANCE PARAMETERS");
    table.add_blank_line();

    // Ae/At
    {
        std::vector<std::string> vals;
        vals.push_back("");
        vals.push_back(format_fixed(1.0, 10, 4));
        for (std::size_t i = 2; i < states.size(); i++) {
            auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
            vals.push_back(format_fixed(perf.area_ratio, 10, 3));
        }
        table.add_row("Ae/At", vals);
    }

    // CSTAR — same for all supersonic stations
    double cstar_val = 0.0;
    if (states.size() >= 2) {
        cstar_val = RocketProblemResults::calculate_performance(states[0], states[1], states[1]).cstar;
    }
    {
        std::vector<std::string> vals;
        vals.push_back("");
        for (std::size_t i = 1; i < states.size(); i++) {
            vals.push_back(format_fixed(cstar_val, 10, 1));
        }
        table.add_row("CSTAR, M/SEC", vals);
    }

    // CF
    {
        std::vector<std::string> vals;
        vals.push_back("");
        auto throat_perf = RocketProblemResults::calculate_performance(states[0], states[1], states[1]);
        vals.push_back(format_fixed(throat_perf.CF, 10, 4));
        for (std::size_t i = 2; i < states.size(); i++) {
            auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
            vals.push_back(format_fixed(perf.CF, 10, 4));
        }
        table.add_row("CF", vals);
    }

    // Ivac
    {
        std::vector<std::string> vals;
        vals.push_back("");
        auto throat_perf = RocketProblemResults::calculate_performance(states[0], states[1], states[1]);
        vals.push_back(format_fixed(throat_perf.ivac, 10, 1));
        for (std::size_t i = 2; i < states.size(); i++) {
            auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
            vals.push_back(format_fixed(perf.ivac, 10, 1));
        }
        table.add_row("Ivac, M/SEC", vals);
    }

    // Isp
    {
        std::vector<std::string> vals;
        vals.push_back("");
        auto throat_perf = RocketProblemResults::calculate_performance(states[0], states[1], states[1]);
        vals.push_back(format_fixed(throat_perf.isp, 10, 1));
        for (std::size_t i = 2; i < states.size(); i++) {
            auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
            vals.push_back(format_fixed(perf.isp, 10, 1));
        }
        table.add_row("Isp, M/SEC", vals);
    }

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
        const auto& chamber_comp = states[0].composition;
        std::vector<std::pair<std::string, double>> species_list;
        for (const auto& [name, frac] : chamber_comp) {
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

} // anonymous namespace

std::string RocketProblemResults::report(const std::string& case_name_arg) const {
    std::vector<std::string> cases_to_report;
    if (case_name_arg.empty()) {
        for (const auto& [name, _] : m_case_meta) {
            cases_to_report.push_back(name);
        }
        std::sort(cases_to_report.begin(), cases_to_report.end());
    } else {
        if (m_case_meta.find(case_name_arg) == m_case_meta.end()) {
            throw std::runtime_error("Case '" + case_name_arg + "' not found.");
        }
        cases_to_report.push_back(case_name_arg);
    }

    std::string result;
    for (const auto& name : cases_to_report) {
        const CaseMeta& meta = m_case_meta.at(name);
        double prev_pressure = -1.0;

        for (std::size_t of_idx = 0; of_idx < meta.of_ratios.size(); of_idx++) {
            for (std::size_t p_idx = 0; p_idx < meta.pressures.size(); p_idx++) {
                // Skip duplicate pressures from the Cantera SolutionArray workaround.
                if (p_idx > 0 && meta.pressures[p_idx] == prev_pressure) continue;
                prev_pressure = meta.pressures[p_idx];

                // Collect stations for this (of_idx, p_idx) in order: chamber, throat, exits.
                std::vector<ThermoStateInfo> thermo_states;
                for (const auto& s : m_stations) {
                    if (s.case_name == name && s.type == StationType::CHAMBER &&
                        s.of_index == of_idx && s.pressure_index == p_idx) {
                        thermo_states.push_back(s.thermo);
                        break;
                    }
                }
                for (const auto& s : m_stations) {
                    if (s.case_name == name && s.type == StationType::THROAT &&
                        s.of_index == of_idx && s.pressure_index == p_idx) {
                        thermo_states.push_back(s.thermo);
                        break;
                    }
                }
                for (std::size_t exp_idx = 0; ; exp_idx++) {
                    bool found = false;
                    for (const auto& s : m_stations) {
                        if (s.case_name == name && s.type == StationType::EXIT &&
                            s.of_index == of_idx && s.pressure_index == p_idx &&
                            s.expansion_index == exp_idx) {
                            thermo_states.push_back(s.thermo);
                            found = true;
                            break;
                        }
                    }
                    if (!found) break;
                }

                if (thermo_states.size() < 2) continue;

                result += build_report_page(
                    meta.chemistry,
                    thermo_states,
                    meta.of_ratios[of_idx],
                    thermo_states[0].pressure);
            }
        }
    }
    return result;
}

} // namespace Goddard
