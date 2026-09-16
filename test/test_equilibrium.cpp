#include "goddard/equilibrium.hpp"
#include "goddard/config.h"
#include "goddard/gas.hpp"
#include "goddard/utils.hpp"
#include "goddard/numerics.hpp"

#include "condensed_test_helpers.hpp"

#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "gtest/gtest.h"
#include <memory>
#include <vector>
#include <cmath>
#include <map>
#include <string>
#include <algorithm>

class DerivativeTests: public ::testing::Test {
    protected:
    DerivativeTests() {
        this->sln = Cantera::newSolution("h2o2.yaml", "ohmech");
        double temp = 2400.0; //K
        double pressure = 50.0*Cantera::OneAtm;
        auto gas = sln->thermo();
        
        n_species = sln->thermo()->nSpecies();
        n_elements = sln->thermo()->nElements();

        gas->setState_TPX(temp, pressure, "H2O:1, N2:1, O2:1. AR:0.1"); //completely random composition lol
    }
    std::shared_ptr<Cantera::Solution> sln;
    std::vector<std::string> elements = {"O", "H", "Ar", "N"};
    std::vector<std::string> species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"};
    size_t n_species;
    size_t n_elements;
    

};

using namespace Goddard;


TEST_F(DerivativeTests, stoichiometricCoeffsAreCorrect){
    auto coeffs = Goddard::get_stoichiometric_coeffs(*sln->thermo());
    Eigen::ArrayXXd expected_coeffs{
        {0, 0, 1, 2, 1, 1, 2, 2, 0, 0},
        {2, 1, 0, 0, 1, 2, 1, 2, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
    };

    long n_elements_signed = static_cast<long>(n_elements); //suppress compiler warnings
    long n_species_signed = static_cast<long>(n_species);

    EXPECT_EQ(coeffs.size(), expected_coeffs.size());

    expected_coeffs.transposeInPlace();

    for (int i = 0; i < n_elements_signed; i++){
        for (int j = 0; j < n_species_signed; j++){
            std::string stoich_id = "Species: ";
            stoich_id += species[j];
            stoich_id += ", element: ";
            stoich_id += elements[i];
            
            EXPECT_EQ(coeffs(j,i), expected_coeffs(j,i)) << stoich_id;
        }
    }
}

TEST_F(DerivativeTests, moleVectorIsCorrect){
    auto moles = Goddard::get_mole_vector(*sln->thermo());
    
    std::vector<double> expected_moles{0., 0., 0.,  0.01219185, 0., 0.01219185, 0. , 0. , 0.001219185, 0.01219185};
    ASSERT_EQ(moles.size(), n_species);

    for (size_t i = 0; i < n_species; i++){
        double tol = max_fp_error(expected_moles[i], 1e-5, 1e-7);
        EXPECT_NEAR(moles(i), expected_moles[i], tol);
    }
}

TEST_F(DerivativeTests, moleVectorIsCorrectAfterEquilibration){
    sln->thermo()->equilibrate("HP");
    auto moles = Goddard::get_mole_vector(*sln->thermo());
    Eigen::ArrayXd expected_moles(n_species);
    expected_moles << 9.11244396e-06, 1.13250699e-06, 2.20601379e-05, 1.21177599e-02,
        2.63718243e-04, 1.20488303e-02, 2.45871008e-06, 2.53459510e-07,
        1.21918510e-03, 1.21918510e-02;

    ASSERT_EQ(moles.size(), n_species);

    for (size_t i = 0; i < n_species; i++){
        double tol = max_fp_error(expected_moles(i));
        EXPECT_NEAR(moles(i), expected_moles(i), tol);
    }
}

TEST_F(DerivativeTests, stdEnthalpiesRTAreCorrect){
    auto enthalpies = Goddard::get_enthalpyRT_vector(*sln->thermo());
    Eigen::ArrayXd exp_enthalpies(n_species);
    exp_enthalpies << 3.35335209, 13.11402496, 14.69431338,  3.73354955,  5.37563442,
       -7.39428256,  5.8582109 ,  0.05367515,  2.18942708,  3.54038084;
    
    ASSERT_EQ(enthalpies.size(), n_species);
    for (size_t i = 0; i < n_species; i++){
        double tol = max_fp_error(exp_enthalpies[i]);
        EXPECT_NEAR(enthalpies(i), exp_enthalpies(i), tol);
    }
}

TEST_F(DerivativeTests, thermoDerivativesAreCorrect){
    EXPECT_NEAR(sln->thermo()->temperature(), 2400.0, 1e-1);
    EXPECT_NEAR(sln->thermo()->enthalpy_mass(), 23985.583414274723, 1e-1);
    
    Goddard::EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    ASSERT_EQ(derivs.dpi_dlogP_T.size(),n_elements);
    ASSERT_EQ(derivs.dpi_dlogT_P.size(), n_elements);

    Eigen::ArrayXd exp_dpi_dlogT_P(n_elements);
    exp_dpi_dlogT_P << -1.86677477,  4.63052867, -2.18942708, -1.77019042;
    
    double exp_dlogn_dlogT_P = -1.8359417928589634e-16;

    Eigen::ArrayXd exp_dpi_dlogP_T(n_elements);
    exp_dpi_dlogP_T << 0.5, 0.25, 1., 0.5;

    double exp_dlogn_dlogP_T = 9.179708964294817e-17;

    for (size_t i = 0; i < n_elements; i++){
        double tol = max_fp_error(exp_dpi_dlogT_P(i));
        EXPECT_NEAR(derivs.dpi_dlogT_P(i), exp_dpi_dlogT_P(i), tol);
    }

    EXPECT_NEAR(derivs.dlogn_dlogT_P, exp_dlogn_dlogT_P, max_fp_error(exp_dlogn_dlogT_P));

    for (size_t i = 0; i < n_elements; i++){
        double tol = max_fp_error(exp_dpi_dlogP_T(i));
        EXPECT_NEAR(derivs.dpi_dlogP_T(i), exp_dpi_dlogP_T(i), tol);
    }

    EXPECT_NEAR(derivs.dlogn_dlogP_T, exp_dlogn_dlogP_T, max_fp_error(exp_dlogn_dlogP_T));

}

TEST_F(DerivativeTests, thermoDerivativesAreCorrectAfterEquilibration){
    sln->thermo()->equilibrate("HP");

    EXPECT_NEAR(sln->thermo()->temperature(), 2367.8424567893862, 1e-1);
    EXPECT_NEAR(sln->thermo()->enthalpy_mass(), 23985.583414274723, 1e-1);
    

    Goddard::EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    ASSERT_EQ(derivs.dpi_dlogP_T.size(),n_elements);
    ASSERT_EQ(derivs.dpi_dlogT_P.size(), n_elements);

    double exp_dlogn_dlogT_P = 0.0198466495473248;
    double exp_dlogn_dlogP_T = -0.0006600595545953462;

    EXPECT_NEAR(derivs.dlogn_dlogT_P, exp_dlogn_dlogT_P, max_fp_error(exp_dlogn_dlogT_P));
    EXPECT_NEAR(derivs.dlogn_dlogP_T, exp_dlogn_dlogP_T, max_fp_error(exp_dlogn_dlogP_T));
}

class PropertyTests: public ::testing::Test {
    protected:
    PropertyTests() {
        this->sln = Cantera::newSolution("h2o2.yaml", "ohmech");
        double temp = 2400.0; //K
        double pressure = 50.0*Cantera::OneAtm;
        auto gas = sln->thermo();
        
        n_species = sln->thermo()->nSpecies();
        n_elements = sln->thermo()->nElements();

        gas->setState_TPX(temp, pressure, "H2O:1, N2:1, O2:1. AR:0.1"); //completely random composition lol
    }
    std::shared_ptr<Cantera::Solution> sln;
    Goddard::ExpansionProperties expected_props{0.9999999999999998,-0.9999999999999999,1604.459611106932,1.2435582661124545};
    Goddard::ExpansionProperties expected_eq_props{1.0198466495473248,-1.0006600595545954,1796.3940426543525,1.222008621037549};
    size_t n_species;
    size_t n_elements;
};

TEST_F(PropertyTests, equilibriumPropertiesAreCorrect) {
    auto derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    Goddard::ExpansionProperties props = Goddard::get_thermo_equilibrium_properties(*sln->thermo(), derivs);
    
    EXPECT_NEAR(props.dlogV_dlogT_P, expected_props.dlogV_dlogT_P, max_fp_error(expected_props.dlogV_dlogT_P));
    EXPECT_NEAR(props.dlogV_dlogP_T, expected_props.dlogV_dlogP_T, max_fp_error(expected_props.dlogV_dlogP_T));
    EXPECT_NEAR(props.spec_heat_p, expected_props.spec_heat_p, max_fp_error(expected_props.spec_heat_p));
    EXPECT_NEAR(props.gamma_s, expected_props.gamma_s, max_fp_error(expected_props.gamma_s));
}

TEST_F(PropertyTests, equilibriumPropertiesAreCorrectAfterEquilibrium) {
    sln->thermo()->equilibrate("HP", "gibbs");
    auto derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    Goddard::ExpansionProperties props = Goddard::get_thermo_equilibrium_properties(*sln->thermo(), derivs);
    
    EXPECT_NEAR(props.dlogV_dlogT_P, expected_eq_props.dlogV_dlogT_P, max_fp_error(expected_eq_props.dlogV_dlogT_P));
    EXPECT_NEAR(props.dlogV_dlogP_T, expected_eq_props.dlogV_dlogP_T, max_fp_error(expected_eq_props.dlogV_dlogP_T));
    EXPECT_NEAR(props.spec_heat_p, expected_eq_props.spec_heat_p, max_fp_error(expected_eq_props.spec_heat_p));
    EXPECT_NEAR(props.gamma_s, expected_eq_props.gamma_s, max_fp_error(expected_eq_props.gamma_s));
}
// ---- Condensed-phase derivatives and properties ----

namespace {

using Goddard::test_helpers::mixture_enthalpy;
using Goddard::test_helpers::mixture_entropy;
using Goddard::test_helpers::mixture_specific_volume;
using Goddard::test_helpers::solve_multiphase_TP;

constexpr const char* NASA9_GAS = DATA_DIR "/nasa9_gas.yaml";
constexpr const char* NASA9_CONDENSED = DATA_DIR "/nasa9_condensed.yaml";

/** Relative step used by every finite-difference check below. */
constexpr double FD_STEP = 1e-4;

/** Central-difference reference values for the equilibrium expansion properties. */
struct FiniteDifferences {
    double spec_heat_p = 0.0;
    double dlogV_dlogT_P = 0.0;
    double dlogV_dlogP_T = 0.0;
    double gamma_s = 0.0;
};

/** Re-equilibrate a copy of `base` at temperature `T` [K] and pressure `P` [Pa]. */
Gas solved_at(const Gas& base, double T, double P) {
    Gas copy = base.clone();
    solve_multiphase_TP(copy, T, P);
    return copy;
}

/**
 * State on the isentrope through `base` at pressure `P` [Pa]: a secant iteration on temperature
 * that holds the mixture entropy (gas plus condensed) constant.
 */
Gas isentropic_state(const Gas& base, double P) {
    const double target_entropy = mixture_entropy(base);

    double T_previous = base.temperature();
    double residual_previous = mixture_entropy(solved_at(base, T_previous, P)) - target_entropy;

    double T_current = T_previous * (1.0 + 1e-4);
    Gas current = solved_at(base, T_current, P);
    double residual_current = mixture_entropy(current) - target_entropy;

    for (int iteration = 0; iteration < 15; iteration++) {
        const double slope = residual_current - residual_previous;
        if (slope == 0.0) break;
        const double T_next = T_current - residual_current * (T_current - T_previous) / slope;

        T_previous = T_current;
        residual_previous = residual_current;
        T_current = T_next;
        current = solved_at(base, T_current, P);
        residual_current = mixture_entropy(current) - target_entropy;

        if (std::abs(T_current - T_previous) < 1e-11 * T_current) break;
    }
    return current;
}

/** Central differences of the equilibrium expansion properties around the state of `base`. */
FiniteDifferences finite_differences(const Gas& base) {
    const double T = base.temperature();
    const double P = base.pressure();

    const Gas hot = solved_at(base, T * (1.0 + FD_STEP), P);
    const Gas cold = solved_at(base, T * (1.0 - FD_STEP), P);
    const Gas compressed = solved_at(base, T, P * (1.0 + FD_STEP));
    const Gas expanded = solved_at(base, T, P * (1.0 - FD_STEP));
    const Gas isentropic_high = isentropic_state(base, P * (1.0 + FD_STEP));
    const Gas isentropic_low = isentropic_state(base, P * (1.0 - FD_STEP));

    FiniteDifferences differences;
    differences.spec_heat_p = (mixture_enthalpy(hot) - mixture_enthalpy(cold))
        / (hot.temperature() - cold.temperature());
    differences.dlogV_dlogT_P =
        std::log(mixture_specific_volume(hot) / mixture_specific_volume(cold))
        / std::log(hot.temperature() / cold.temperature());
    differences.dlogV_dlogP_T =
        std::log(mixture_specific_volume(compressed) / mixture_specific_volume(expanded))
        / std::log(compressed.pressure() / expanded.pressure());
    differences.gamma_s =
        -std::log(isentropic_high.pressure() / isentropic_low.pressure())
        / std::log(mixture_specific_volume(isentropic_high) / mixture_specific_volume(isentropic_low));
    return differences;
}

/** Compare the analytic equilibrium properties of `gas` with central differences. */
void expect_matches_finite_differences(const Gas& gas, double reltol = 1e-4) {
    const ExpansionProperties props = Goddard::get_thermo_equilibrium_properties(gas);
    const FiniteDifferences differences = finite_differences(gas);

    EXPECT_NEAR(props.spec_heat_p, differences.spec_heat_p,
                reltol * std::abs(differences.spec_heat_p)) << "cp_e [J/(kg.K)]";
    EXPECT_NEAR(props.dlogV_dlogT_P, differences.dlogV_dlogT_P,
                reltol * std::abs(differences.dlogV_dlogT_P)) << "(dlnV/dlnT)_P";
    EXPECT_NEAR(props.dlogV_dlogP_T, differences.dlogV_dlogP_T,
                reltol * std::abs(differences.dlogV_dlogP_T)) << "(dlnV/dlnP)_T";
    EXPECT_NEAR(props.gamma_s, differences.gamma_s,
                reltol * std::abs(differences.gamma_s)) << "gamma_s";
}

/** CH4/O2 at O/F = 0.5 by mass, which leaves graphite in the products. */
Gas make_methane_oxygen_gas(double T, double P) {
    // A CEA-style product list rather than every C/H/O species of the database: with all of them
    // offered, `MultiPhaseEquil` does not converge at this rich, cool state.
    Gas gas = Gas::create_from_species(
        NASA9_GAS, "gas", {"CH4", "CO", "CO2", "H", "H2", "H2O", "O", "O2", "OH", "C"},
        GasChemistry::EQUILIBRIUM);
    gas.set_state_TPX(T, P, Composition{{"CH4", (1.0 / 1.5) / 16.04246},
                                        {"O2", (0.5 / 1.5) / 31.9988}});
    gas.add_condensed_species(NASA9_CONDENSED, {"C(gr)", "H2O(L)", "H2O(cr)"});
    return gas;
}

/** H2/O2 in a 100:60 molar ratio, cold enough for the product water to condense. */
Gas make_hydrogen_oxygen_gas(double T, double P) {
    Gas gas = Gas::create_from_elements(NASA9_GAS, "gas", {"H", "O"}, GasChemistry::EQUILIBRIUM);
    gas.set_state_TPX(T, P, Composition{{"H2", 100.0}, {"O2", 60.0}});
    gas.add_condensed_species(NASA9_CONDENSED, {"H2O(L)", "H2O(cr)"});
    return gas;
}

/**
 * AP/Al/binder propellant per kg: 72 % NH4ClO4, 10 % C1H1.86955O0.031256 binder, 18 % Al.
 * Molar masses [kg/kmol] are those of the elemental formulas above.
 */
Gas make_ap_aluminium_gas(double T, double P, const std::vector<std::string>& condensed) {
    Gas gas = Gas::create_from_elements(NASA9_GAS, "gas", {"N", "H", "Cl", "O", "C", "Al"},
                                        GasChemistry::EQUILIBRIUM);

    const double perchlorate = 0.72 / 117.48906;
    const double binder = 0.10 / 14.39543;
    const double aluminium = 0.18 / 26.9815;

    const std::map<std::string, double> amounts{
        {"N", perchlorate},
        {"H", 4.0 * perchlorate + 1.86955 * binder},
        {"Cl", perchlorate},
        {"O", 4.0 * perchlorate + 0.031256 * binder},
        {"C", binder},
        {"Al", aluminium},
    };

    const std::vector<std::string> element_names = gas.element_names();
    Eigen::ArrayXd element_moles(static_cast<long>(element_names.size()));
    for (size_t m = 0; m < element_names.size(); m++) {
        element_moles(static_cast<long>(m)) = amounts.at(element_names[m]);
    }
    gas.set_element_moles(element_moles, T, P);

    gas.add_condensed_species(NASA9_CONDENSED, condensed);
    return gas;
}

} // namespace

TEST(CondensedEquilibriumTests, MethaneOxygenGraphiteMatchesFiniteDifferences) {
    Gas gas = make_methane_oxygen_gas(1500.0, 1.0e5);
    solve_multiphase_TP(gas, 1500.0, 1.0e5);

    ASSERT_TRUE(gas.has_condensed_phases());
    EXPECT_GT(gas.condensed_moles()[0], 0.0) << "C(gr) should be present";

    const EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(gas);
    EXPECT_EQ(derivs.dn_condensed_dlogT_P.size(), 1);
    EXPECT_FALSE(derivs.pinned_transition);

    expect_matches_finite_differences(gas);
}

TEST(CondensedEquilibriumTests, TemperatureDerivativesConserveElements) {
    Gas gas = make_methane_oxygen_gas(1500.0, 1.0e5);
    solve_multiphase_TP(gas, 1500.0, 1.0e5);
    ASSERT_TRUE(gas.has_condensed_phases());

    const EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(gas);
    const Eigen::ArrayXXd gas_stoich = Goddard::get_stoichiometric_coeffs(*gas.thermo());
    const Eigen::ArrayXd gas_moles =
        Goddard::get_mole_vector(*gas.thermo()) * gas.gas_mass_fraction();
    const Eigen::ArrayXd enthalpy_RT = Goddard::get_enthalpyRT_vector(*gas.thermo());
    const Eigen::ArrayXXd condensed_stoich = gas.condensed_stoich_coeffs();
    const Eigen::ArrayXd element_moles = gas.element_moles();

    const long n_species = gas_moles.size();
    const long n_elements = gas_stoich.cols();

    // d ln n_j / d ln T = H_j/RT + d ln n / d ln T + sum_i a_ij dpi_i/d ln T
    Eigen::ArrayXd dlogn_species(n_species);
    for (long j = 0; j < n_species; j++) {
        double value = enthalpy_RT(j) + derivs.dlogn_dlogT_P;
        for (long i = 0; i < n_elements; i++) {
            value += gas_stoich(j, i) * derivs.dpi_dlogT_P(i);
        }
        dlogn_species(j) = value;
    }

    for (long k = 0; k < n_elements; k++) {
        double residual = (gas_stoich.col(k) * gas_moles * dlogn_species).sum();
        for (long c = 0; c < condensed_stoich.rows(); c++) {
            residual += condensed_stoich(c, k) * derivs.dn_condensed_dlogT_P(c);
        }
        EXPECT_NEAR(residual, 0.0, 1e-8 * element_moles(k))
            << "element " << gas.element_names()[static_cast<size_t>(k)];
    }
}

TEST(CondensedEquilibriumTests, HydrogenOxygenLiquidWaterMatchesFiniteDifferences) {
    const double pressure = 0.05 * Cantera::OneAtm;
    Gas gas = make_hydrogen_oxygen_gas(300.0, pressure);
    solve_multiphase_TP(gas, 300.0, pressure);

    ASSERT_TRUE(gas.has_condensed_phases());
    EXPECT_GT(gas.condensed_moles()[0], 0.0) << "H2O(L) should be present";
    EXPECT_DOUBLE_EQ(gas.condensed_moles()[1], 0.0) << "H2O(cr) is out of range at 300 K";

    expect_matches_finite_differences(gas);
}

TEST(CondensedEquilibriumTests, HydrogenOxygenIceMatchesFiniteDifferences) {
    const double pressure = 0.05 * Cantera::OneAtm;
    Gas gas = make_hydrogen_oxygen_gas(265.0, pressure);
    solve_multiphase_TP(gas, 265.0, pressure);

    ASSERT_TRUE(gas.has_condensed_phases());
    EXPECT_GT(gas.condensed_moles()[1], 0.0) << "H2O(cr) should be present";
    EXPECT_DOUBLE_EQ(gas.condensed_moles()[0], 0.0) << "H2O(L) is out of range at 265 K";

    expect_matches_finite_differences(gas);
}

TEST(CondensedEquilibriumTests, AluminisedPropellantMatchesFiniteDifferences) {
    const double pressure = 3.45e5;
    Gas gas = make_ap_aluminium_gas(3000.0, pressure, {"AL2O3(a)", "AL2O3(L)"});
    solve_multiphase_TP(gas, 3000.0, pressure);

    ASSERT_TRUE(gas.has_condensed_phases());
    EXPECT_GT(gas.condensed_moles()[1], 0.0) << "AL2O3(L) should be present";

    expect_matches_finite_differences(gas);
}

TEST(CondensedEquilibriumTests, PinnedAluminaTransitionMergesThePolymorphPair) {
    const double transition_temperature = 2327.0;
    const double pressure = 3.45e5;

    Gas gas = make_ap_aluminium_gas(transition_temperature, pressure, {"AL2O3(a)", "AL2O3(L)"});
    solve_multiphase_TP(gas, transition_temperature, pressure, {"AL2O3(a)"});

    std::vector<double> moles = gas.condensed_moles();
    ASSERT_EQ(moles.size(), 2u);
    ASSERT_GT(moles[0], 0.0) << "AL2O3(a) should be present";

    // Reference state: the same gas with the whole group as the low-temperature polymorph.
    const Gas unpinned = gas.clone();
    const EquilibriumDerivatives reference_derivs =
        Goddard::get_thermo_equilibrium_derivatives(unpinned);
    const ExpansionProperties reference =
        Goddard::get_thermo_equilibrium_properties(unpinned, reference_derivs);
    EXPECT_FALSE(reference.pinned_transition);

    // Split the group between the two polymorphs and pin it, as the multiphase solver will.
    moles[0] *= 0.5;
    moles[1] = moles[0];
    gas.set_condensed_moles(moles);
    gas.set_phase_transition("AL2O3(a)", "AL2O3(L)");
    ASSERT_TRUE(gas.at_phase_transition());
    ASSERT_EQ(gas.pinned_polymorphs(), std::make_pair(0L, 1L));

    const EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(gas);
    EXPECT_TRUE(derivs.pinned_transition);
    EXPECT_TRUE(std::isnan(derivs.dlogn_dlogT_P));
    ASSERT_EQ(derivs.dn_condensed_dlogT_P.size(), 2);
    EXPECT_TRUE(std::isnan(derivs.dn_condensed_dlogT_P(0)));
    ASSERT_EQ(derivs.dn_condensed_dlogP_T.size(), 2);
    // The merged pressure derivative is reported on the lower-temperature polymorph.
    EXPECT_DOUBLE_EQ(derivs.dn_condensed_dlogP_T(1), 0.0);
    EXPECT_NEAR(derivs.dn_condensed_dlogP_T(0), reference_derivs.dn_condensed_dlogP_T(0),
                1e-10 * std::abs(reference_derivs.dn_condensed_dlogP_T(0)));

    const ExpansionProperties props = Goddard::get_thermo_equilibrium_properties(gas, derivs);
    EXPECT_TRUE(props.pinned_transition);
    EXPECT_TRUE(std::isinf(props.spec_heat_p));
    EXPECT_TRUE(std::isinf(props.spec_heat_v));
    EXPECT_TRUE(std::isinf(props.dlogV_dlogT_P));
    EXPECT_TRUE(std::isfinite(props.gamma_s));
    EXPECT_GT(props.gamma_s, 0.0);
    EXPECT_DOUBLE_EQ(props.gamma_s, -1.0 / props.dlogV_dlogP_T);
    EXPECT_NEAR(props.gamma_s, -1.0 / reference.dlogV_dlogP_T,
                1e-10 * std::abs(reference.dlogV_dlogP_T));
    EXPECT_GT(props.speed_of_sound, 0.0);
}

TEST(CondensedEquilibriumTests, MixtureOverloadsReduceToTheGasOnlyOnes) {
    Gas gas = make_methane_oxygen_gas(1500.0, 1.0e5);
    // Candidates are attached but nothing has condensed yet.
    ASSERT_FALSE(gas.has_condensed_phases());

    const EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(gas);
    const EquilibriumDerivatives expected_derivs =
        Goddard::get_thermo_equilibrium_derivatives(*gas.thermo());
    EXPECT_DOUBLE_EQ(derivs.dlogn_dlogT_P, expected_derivs.dlogn_dlogT_P);
    EXPECT_DOUBLE_EQ(derivs.dlogn_dlogP_T, expected_derivs.dlogn_dlogP_T);
    EXPECT_EQ(derivs.dn_condensed_dlogT_P.size(), 0);
    EXPECT_EQ(derivs.dn_condensed_dlogP_T.size(), 0);

    const ExpansionProperties props = Goddard::get_thermo_equilibrium_properties(gas);
    const ExpansionProperties expected =
        Goddard::get_thermo_equilibrium_properties(*gas.thermo());
    EXPECT_DOUBLE_EQ(props.spec_heat_p, expected.spec_heat_p);
    EXPECT_DOUBLE_EQ(props.gamma_s, expected.gamma_s);
    EXPECT_DOUBLE_EQ(props.dlogV_dlogT_P, expected.dlogV_dlogT_P);
    EXPECT_DOUBLE_EQ(props.dlogV_dlogP_T, expected.dlogV_dlogP_T);
    EXPECT_DOUBLE_EQ(props.total_moles, props.gas_moles);
    EXPECT_DOUBLE_EQ(gas.gas_mass_fraction(), 1.0);
}

TEST_F(PropertyTests, GasOverloadReproducesTheGoldenValues) {
    Gas gas(sln, GasChemistry::EQUILIBRIUM);
    const ExpansionProperties props = Goddard::get_thermo_equilibrium_properties(gas);

    EXPECT_NEAR(props.dlogV_dlogT_P, expected_props.dlogV_dlogT_P,
                max_fp_error(expected_props.dlogV_dlogT_P));
    EXPECT_NEAR(props.dlogV_dlogP_T, expected_props.dlogV_dlogP_T,
                max_fp_error(expected_props.dlogV_dlogP_T));
    EXPECT_NEAR(props.spec_heat_p, expected_props.spec_heat_p,
                max_fp_error(expected_props.spec_heat_p));
    EXPECT_NEAR(props.gamma_s, expected_props.gamma_s, max_fp_error(expected_props.gamma_s));

    sln->thermo()->equilibrate("HP", "gibbs");
    const ExpansionProperties eq_props = Goddard::get_thermo_equilibrium_properties(gas);
    EXPECT_NEAR(eq_props.dlogV_dlogT_P, expected_eq_props.dlogV_dlogT_P,
                max_fp_error(expected_eq_props.dlogV_dlogT_P));
    EXPECT_NEAR(eq_props.dlogV_dlogP_T, expected_eq_props.dlogV_dlogP_T,
                max_fp_error(expected_eq_props.dlogV_dlogP_T));
    EXPECT_NEAR(eq_props.spec_heat_p, expected_eq_props.spec_heat_p,
                max_fp_error(expected_eq_props.spec_heat_p));
    EXPECT_NEAR(eq_props.gamma_s, expected_eq_props.gamma_s,
                max_fp_error(expected_eq_props.gamma_s));
}

TEST_F(DerivativeTests, GasOverloadReproducesTheGoldenDerivatives) {
    Gas gas(sln, GasChemistry::EQUILIBRIUM);
    const EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(gas);

    Eigen::ArrayXd exp_dpi_dlogT_P(n_elements);
    exp_dpi_dlogT_P << -1.86677477, 4.63052867, -2.18942708, -1.77019042;
    Eigen::ArrayXd exp_dpi_dlogP_T(n_elements);
    exp_dpi_dlogP_T << 0.5, 0.25, 1., 0.5;

    for (size_t i = 0; i < n_elements; i++) {
        EXPECT_NEAR(derivs.dpi_dlogT_P(i), exp_dpi_dlogT_P(i), max_fp_error(exp_dpi_dlogT_P(i)));
        EXPECT_NEAR(derivs.dpi_dlogP_T(i), exp_dpi_dlogP_T(i), max_fp_error(exp_dpi_dlogP_T(i)));
    }
    EXPECT_NEAR(derivs.dlogn_dlogT_P, -1.8359417928589634e-16, max_fp_error(0.0));
    EXPECT_NEAR(derivs.dlogn_dlogP_T, 9.179708964294817e-17, max_fp_error(0.0));
    EXPECT_EQ(derivs.dn_condensed_dlogT_P.size(), 0);
    EXPECT_EQ(derivs.dn_condensed_dlogP_T.size(), 0);
    EXPECT_FALSE(derivs.pinned_transition);
}
