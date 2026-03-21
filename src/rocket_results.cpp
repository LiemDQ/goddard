

#include <optional>

#include "goddard/rocket_results.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/format.hpp"

namespace Goddard {


RocketProblemResults::RocketProblemResults(
    std::unordered_map<std::string, RocketProblemCaseResult>&& case_results, 
    std::shared_ptr<Cantera::Solution> sln) 
    : cases(case_results), m_sln(std::move(sln)) {
}

std::vector<ThermoStateInfo> RocketProblemResults::extract_thermo_info(const std::string& case_name, std::size_t index){
    RocketProblemCaseResult& case_result = cases.at(case_name);
    auto tmo = thermo();

    if (index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        throw std::runtime_error("Provided index exceeds length of ThermoArray.");
    }

    std::vector<double> state = case_result.inlet_states.get_state(static_cast<int>(index));
    NozzleResults nozzle_results = case_result.nozzle_states[index];
    

    std::unordered_map<std::string, double> inlet_compositions; 
    std::vector<std::string> names = tmo->speciesNames();
    

    for (std::size_t i = 0; i < tmo->nSpecies(); i++) {
        inlet_compositions[names[i]] = state[i+2];
    }

    tmo->restoreState(state);
    double inlet_gamma = 0;
    double inlet_dlogV_dlogP_T = 0;
    double inlet_dlogV_dlogT_P = 0;

    switch (case_result.chemistry) {
        case NozzleChemistryType::EQUILIBRIUM: {
            auto props = get_thermo_equilibrium_properties(*tmo);
            inlet_gamma = props.gamma_s;
            inlet_dlogV_dlogP_T = props.dlogV_dlogP_T;
            inlet_dlogV_dlogT_P = props.dlogV_dlogT_P;
            break;
        }
        case NozzleChemistryType::FROZEN: {
            inlet_gamma = tmo->cp_mass()/tmo->cv_mass();
            inlet_dlogV_dlogP_T = -1.0; //by definition -- see NASA reference publication 1311 eq 6.38
            inlet_dlogV_dlogT_P = 1.0; //by definition -- see NASA reference publication 1311 eq 6.37
            break;
        }
        default: throw NotImplementedError("This chemistry type has not been implemented.");
    }
        
    ThermoStateInfo inlet_state{
        tmo->pressure(),
        tmo->temperature(),
        tmo->density(),
        tmo->enthalpy_mass(),
        tmo->intEnergy_mass(),
        tmo->gibbs_mass(),
        tmo->entropy_mass(),
        tmo->meanMolecularWeight(),
        tmo->cp_mass(),
        inlet_gamma,
        inlet_dlogV_dlogP_T,
        inlet_dlogV_dlogT_P,
        gas_sonic_velocity(*tmo, inlet_gamma),
        inlet_compositions
    };

    //throat
    ThroatCondition throat = nozzle_results.throat;
    std::unordered_map<std::string, double> throat_compositions; 
    tmo->restoreState(throat.state);

    for (std::size_t i = 0; i < tmo->nSpecies(); i++) {
        throat_compositions[names[i]] = state[i+2];
    }

    
    ThermoStateInfo throat_state {
        tmo->pressure(),
        tmo->temperature(),
        tmo->density(),
        tmo->enthalpy_mass(),
        tmo->intEnergy_mass(),
        tmo->gibbs_mass(),
        tmo->entropy_mass(),
        tmo->meanMolecularWeight(),
        tmo->cp_mass(),
        throat.gamma_s,
        throat.dlV_dlP_T,
        throat.dlV_dlT_P,
        gas_sonic_velocity(*tmo, throat.gamma_s),
        throat_compositions
    };

    std::vector<ThermoStateInfo> state_info {inlet_state, throat_state}; 

    for (auto&& exp: nozzle_results.expansions) {
        tmo->restoreState(exp.state);
        std::unordered_map<std::string, double> compositions;
        for (std::size_t i = 0; i < tmo->nSpecies(); i++) {
           compositions[names[i]] = state[i+2];
        }

        ThermoStateInfo exp_state {
            tmo->pressure(),
            tmo->temperature(),
            tmo->density(),
            tmo->enthalpy_mass(),
            tmo->intEnergy_mass(),
            tmo->gibbs_mass(),
            tmo->entropy_mass(),
            tmo->meanMolecularWeight(),
            tmo->cp_mass(),
            exp.gamma_s,
            exp.dlV_dlP_T,
            exp.dlV_dlT_P,
            gas_sonic_velocity(*tmo, exp.gamma_s),
            compositions
        };

        state_info.push_back(exp_state);

    }

    return state_info;
}

std::optional<ThermoStateInfo> RocketProblemResults::get_chamber_state(
    const std::string& case_name, std::size_t of_index) {

    auto it = cases.find(case_name);
    if (it == cases.end()) {
        return std::nullopt;
    }

    RocketProblemCaseResult& case_result = it->second;
    if (of_index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        return std::nullopt;
    }

    auto all_states = extract_thermo_info(case_name, of_index);
    if (all_states.empty()) {
        return std::nullopt;
    }

    return all_states[0]; // Chamber is first state
}

std::optional<ThermoStateInfo> RocketProblemResults::get_throat_state(
    const std::string& case_name, std::size_t of_index) {

    auto it = cases.find(case_name);
    if (it == cases.end()) {
        return std::nullopt;
    }

    RocketProblemCaseResult& case_result = it->second;
    if (of_index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        return std::nullopt;
    }

    auto all_states = extract_thermo_info(case_name, of_index);
    if (all_states.size() < 2) {
        return std::nullopt;
    }

    return all_states[1]; // Throat is second state
}

std::vector<ThermoStateInfo> RocketProblemResults::get_exit_states(
    const std::string& case_name, std::size_t of_index) {

    auto it = cases.find(case_name);
    if (it == cases.end()) {
        return {};
    }

    RocketProblemCaseResult& case_result = it->second;
    if (of_index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        return {};
    }

    auto all_states = extract_thermo_info(case_name, of_index);
    if (all_states.size() <= 2) {
        return {}; // No exit states
    }

    // Exit states start at index 2
    return std::vector<ThermoStateInfo>(all_states.begin() + 2, all_states.end());
}

RocketPerformance RocketProblemResults::calculate_performance(
    const ThermoStateInfo& chamber,
    const ThermoStateInfo& throat,
    const ThermoStateInfo& exit) {

    // Pressure ratios
    double pressure_ratio = chamber.pressure / exit.pressure;
    double throat_pressure_ratio = chamber.pressure / throat.pressure;

    // Area ratio from isentropic flow relations
    // A/A* = (1/M) * [(2/(gamma+1)) * (1 + (gamma-1)/2 * M^2)]^((gamma+1)/(2*(gamma-1)))
    // For now, estimate from density ratio (approximate)
    double area_ratio = (throat.density * throat.speed_of_sound) /
                        (exit.density * exit.speed_of_sound) *
                        (throat_pressure_ratio / pressure_ratio);

    // Characteristic velocity (c*)
    // c* = P_c * A_t / m_dot = sqrt(gamma * R * T_c) / gamma * sqrt((2/(gamma+1))^((gamma+1)/(gamma-1)))
    // Simplified: c* = throat.speed_of_sound / sqrt(throat.gamma_s) * factor
    double gamma = throat.gamma_s;
    double cstar = throat.speed_of_sound * std::sqrt(
        std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (gamma - 1.0)) / gamma
    );

    // Exit velocity from enthalpy difference
    // v_e = sqrt(2 * (h_chamber - h_exit))
    double exit_velocity = std::sqrt(2.0 * (chamber.enthalpy - exit.enthalpy));

    // Thrust coefficient CF = v_e / c* + (P_e - P_amb) * A_e / (P_c * A_t)
    // For vacuum: CF_vac = v_e / c* + P_e * A_e / (P_c * A_t)
    // Simplified (assuming matched nozzle): CF ≈ v_e / c*
    double CF = exit_velocity / cstar;

    // Specific impulse Isp = v_e / g0
    constexpr double g0 = 9.80665; // m/s^2
    double isp = exit_velocity / g0;

    // Vacuum specific impulse (includes pressure thrust term)
    // Ivac = Isp + P_e * A_e / (m_dot * g0)
    // Approximation: Ivac ≈ Isp * (1 + P_e/(P_c) * area_ratio * some_factor)
    double ivac = isp + (exit.pressure / chamber.pressure) * area_ratio * cstar / g0;

    // Mach number at exit (approximate from speed of sound)
    double mach_number = exit_velocity / exit.speed_of_sound;

    return RocketPerformance{
        pressure_ratio,
        area_ratio,
        mach_number,
        cstar,
        CF,
        isp,
        ivac
    };
}

namespace {

constexpr double PA_TO_BAR = 1e-5;
constexpr double PA_TO_PSIA = 1.0 / 6894.757;
constexpr double J_TO_KJ = 1e-3;

std::string build_report_page(
    const RocketProblemCaseResult& case_result,
    const std::vector<ThermoStateInfo>& states,
    double of_ratio,
    double chamber_pressure_pa)
{
    std::string page;

    // Title
    switch (case_result.chemistry) {
        case NozzleChemistryType::EQUILIBRIUM:
            page += "         THEORETICAL ROCKET PERFORMANCE ASSUMING EQUILIBRIUM\n\n";
            page += "      COMPOSITION DURING EXPANSION FROM INFINITE AREA COMBUSTOR\n\n";
            break;
        case NozzleChemistryType::FROZEN:
            page += "         THEORETICAL ROCKET PERFORMANCE ASSUMING FROZEN COMPOSITION\n\n";
            break;
        default: break;
    }

    page += "Pin = " + format_fixed(chamber_pressure_pa * PA_TO_PSIA, 7, 1) + " PSIA\n";
    page += "O/F=" + format_fixed(of_ratio, 11, 5) + "\n\n";

    // Column headers
    TextTable table(18, 10);
    std::vector<std::string> headers;
    headers.push_back("CHAMBER");
    headers.push_back("THROAT");
    for (std::size_t i = 2; i < states.size(); i++) {
        headers.push_back("EXIT");
    }
    table.set_headers(headers);

    const auto& chamber = states[0];

    // Pinf/P
    std::vector<std::string> vals;
    for (auto& s : states) {
        vals.push_back(format_fixed(chamber.pressure / s.pressure, 10, 4));
    }
    table.add_row("Pinf/P", vals);

    // P, BAR
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.pressure * PA_TO_BAR, 10, 4));
    }
    table.add_row("P, BAR", vals);

    // T, K
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.temperature, 10, 2));
    }
    table.add_row("T, K", vals);

    // RHO, KG/CU M (CEA engineering notation)
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_cea_engineering(s.density, 10));
    }
    table.add_row("RHO, KG/CU M", vals);

    // H, KJ/KG
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.enthalpy * J_TO_KJ, 10, 2));
    }
    table.add_row("H, KJ/KG", vals);

    // U, KJ/KG
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.internal_energy * J_TO_KJ, 10, 2));
    }
    table.add_row("U, KJ/KG", vals);

    // G, KJ/KG
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.gibbs * J_TO_KJ, 10, 1));
    }
    table.add_row("G, KJ/KG", vals);

    // S, KJ/(KG)(K)
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.entropy * J_TO_KJ, 10, 4));
    }
    table.add_row("S, KJ/(KG)(K)", vals);

    table.add_blank_line();

    // M, (1/n)
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.molecular_weight, 10, 3));
    }
    table.add_row("M, (1/n)", vals);

    // (dLV/dLP)t and (dLV/dLT)p — equilibrium only
    if (case_result.chemistry == NozzleChemistryType::EQUILIBRIUM) {
        vals.clear();
        for (auto& s : states) {
            vals.push_back(format_fixed(s.dlV_dlP_T, 10, 5));
        }
        table.add_row("(dLV/dLP)t", vals);

        vals.clear();
        for (auto& s : states) {
            vals.push_back(format_fixed(s.dlV_dlT_P, 10, 4));
        }
        table.add_row("(dLV/dLT)p", vals);
    }

    // Cp, KJ/(KG)(K)
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.cp * J_TO_KJ, 10, 4));
    }
    table.add_row("Cp, KJ/(KG)(K)", vals);

    // GAMMAs
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.gamma_s, 10, 4));
    }
    table.add_row("GAMMAs", vals);

    // SON VEL, M/SEC
    vals.clear();
    for (auto& s : states) {
        vals.push_back(format_fixed(s.speed_of_sound, 10, 1));
    }
    table.add_row("SON VEL,M/SEC", vals);

    // MACH NUMBER: chamber=0, throat=1, exits computed from performance
    vals.clear();
    vals.push_back(format_fixed(0.0, 10, 3));
    vals.push_back(format_fixed(1.0, 10, 3));
    for (std::size_t i = 2; i < states.size(); i++) {
        auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
        vals.push_back(format_fixed(perf.mach_number, 10, 3));
    }
    table.add_row("MACH NUMBER", vals);

    // Performance parameters section
    table.add_blank_line();
    table.add_section_header("PERFORMANCE PARAMETERS");
    table.add_blank_line();

    // Ae/At
    vals.clear();
    vals.push_back(""); // chamber — no Ae/At
    vals.push_back(format_fixed(1.0, 10, 4)); // throat
    for (std::size_t i = 2; i < states.size(); i++) {
        auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
        vals.push_back(format_fixed(perf.area_ratio, 10, 3));
    }
    table.add_row("Ae/At", vals);

    // CSTAR, M/SEC — same for all stations
    double cstar_val = 0.0;
    if (states.size() >= 2) {
        // Use throat as exit to compute cstar (cstar depends only on throat conditions)
        auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[1]);
        cstar_val = perf.cstar;
    }
    vals.clear();
    vals.push_back(""); // chamber
    for (std::size_t i = 1; i < states.size(); i++) {
        vals.push_back(format_fixed(cstar_val, 10, 1));
    }
    table.add_row("CSTAR, M/SEC", vals);

    // CF
    vals.clear();
    vals.push_back(""); // chamber
    {
        auto throat_perf = RocketProblemResults::calculate_performance(states[0], states[1], states[1]);
        vals.push_back(format_fixed(throat_perf.CF, 10, 4));
    }
    for (std::size_t i = 2; i < states.size(); i++) {
        auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
        vals.push_back(format_fixed(perf.CF, 10, 4));
    }
    table.add_row("CF", vals);

    // Ivac, M/SEC
    vals.clear();
    vals.push_back(""); // chamber
    {
        auto throat_perf = RocketProblemResults::calculate_performance(states[0], states[1], states[1]);
        vals.push_back(format_fixed(throat_perf.ivac, 10, 1));
    }
    for (std::size_t i = 2; i < states.size(); i++) {
        auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
        vals.push_back(format_fixed(perf.ivac, 10, 1));
    }
    table.add_row("Ivac, M/SEC", vals);

    // Isp, M/SEC
    vals.clear();
    vals.push_back(""); // chamber
    {
        auto throat_perf = RocketProblemResults::calculate_performance(states[0], states[1], states[1]);
        vals.push_back(format_fixed(throat_perf.isp, 10, 1));
    }
    for (std::size_t i = 2; i < states.size(); i++) {
        auto perf = RocketProblemResults::calculate_performance(states[0], states[1], states[i]);
        vals.push_back(format_fixed(perf.isp, 10, 1));
    }
    table.add_row("Isp, M/SEC", vals);

    // Mass fractions section
    table.add_blank_line();
    table.add_blank_line();
    table.add_section_header("MASS FRACTIONS");
    table.add_blank_line();

    if (case_result.chemistry == NozzleChemistryType::EQUILIBRIUM) {
        // Collect all species that appear in any station
        std::set<std::string> all_species;
        for (auto& s : states) {
            for (auto& [name, frac] : s.composition) {
                if (frac > 1e-10) all_species.insert(name);
            }
        }

        // One row per species, values at each station
        for (const auto& species : all_species) {
            vals.clear();
            for (auto& s : states) {
                auto it = s.composition.find(species);
                double val = (it != s.composition.end()) ? it->second : 0.0;
                vals.push_back(format_mass_fraction(val, 10));
            }
            table.add_row(species, vals);
        }
    } else {
        // Frozen: compact 3-per-row format (composition constant across stations)
        const auto& chamber_comp = states[0].composition;

        // Collect and sort species with non-trivial fractions
        std::vector<std::pair<std::string, double>> species_list;
        for (auto& [name, frac] : chamber_comp) {
            if (frac > 1e-10) {
                species_list.push_back({name, frac});
            }
        }
        std::sort(species_list.begin(), species_list.end());

        // Format 3 per row: "species  value   species  value   species  value"
        constexpr int entries_per_row = 3;
        constexpr int name_width = 16;
        constexpr int frac_width = 10;

        for (std::size_t i = 0; i < species_list.size(); i += entries_per_row) {
            std::string line;
            for (std::size_t j = i; j < i + entries_per_row && j < species_list.size(); j++) {
                std::string name = species_list[j].first;
                // Pad name to name_width
                if (static_cast<int>(name.size()) < name_width)
                    name += std::string(name_width - name.size(), ' ');
                line += name + format_mass_fraction(species_list[j].second, frac_width);
                if (j + 1 < i + entries_per_row && j + 1 < species_list.size()) {
                    line += "   "; // separator between entries
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

std::string RocketProblemResults::report(const std::string& case_name) {
    // Determine which cases to report
    std::vector<std::string> case_names;
    if (case_name.empty()) {
        for (auto& [name, _] : cases) {
            case_names.push_back(name);
        }
        std::sort(case_names.begin(), case_names.end());
    } else {
        if (cases.find(case_name) == cases.end()) {
            throw std::runtime_error("Case '" + case_name + "' not found.");
        }
        case_names.push_back(case_name);
    }

    std::string result;

    for (const auto& name : case_names) {
        RocketProblemCaseResult& case_result = cases.at(name);

        std::size_t num_of = case_result.OF_ratios.size();
        std::size_t num_pressures = case_result.pressures.size();

        // Track previous pressure to skip duplicates from the Cantera SolutionArray workaround
        double prev_pressure = -1.0;

        for (std::size_t of_idx = 0; of_idx < num_of; of_idx++) {
            for (std::size_t p_idx = 0; p_idx < num_pressures; p_idx++) {
                // Skip duplicate pressures from the Cantera workaround
                if (p_idx > 0 && case_result.pressures[p_idx] == prev_pressure) {
                    continue;
                }
                prev_pressure = case_result.pressures[p_idx];

                // ThermoArray shape is {N_of, N_p}: each O/F composition is
                // paired with its corresponding enthalpy, not broadcast.
                std::size_t state_idx = of_idx * num_pressures + p_idx;

                if (state_idx >= static_cast<std::size_t>(case_result.inlet_states.size())) {
                    continue;
                }

                auto thermo_states = extract_thermo_info(name, state_idx);
                if (thermo_states.size() < 2) {
                    continue;
                }

                result += build_report_page(case_result, thermo_states,
                    case_result.OF_ratios[of_idx],
                    case_result.pressures[p_idx]);
            }
        }
    }

    return result;
}

} //namespace Goddard