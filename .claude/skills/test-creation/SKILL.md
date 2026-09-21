---
name: test-creation
description: Writes tests for Goddard using the recommended best practices.
---
# Skill: Writing Tests for Goddard

## Overview

Goddard is a C++/Python rocket engine simulation toolkit. Tests must verify that the physics is right, not just that the code runs. Every test you write should answer: "What physical or mathematical property would break if this code were wrong?"

## Ground Rules

### Choose the right language for the job

- **C++ (GTest)**: Testing fundamental implementation correctness — invariants, convergence, state management, error handling, thermodynamic identities, and algebraic relationships that hold regardless of specific numerical values. Located in `test/` with pattern `test_*.cpp`.
- **Python (pytest)**: Numerical validation against reference data (CEA, literature values, analytical solutions). Python tests live in `python/tests/` and use the shared infrastructure in `conftest.py`. The `cea` package may not be installed; use `pytest.importorskip("cea")` to skip gracefully.

If you're unsure: if the test checks "does this quantity match a known number from CEA or a textbook?", write it in Python. If it checks "does this operation preserve an invariant of the physics?", write it in C++.

### Test physics, not plumbing

Every test must verify something meaningful about the simulation's correctness. Ask yourself:

1. **Would a physically wrong implementation fail this test?** If a function returned random numbers in the right range and the test could still pass, the test is worthless.
2. **Is this testing the claim in the test name?** A test named `EnergyConservation` that only checks `converged == true` is lying about what it tests.
3. **Could this test catch a real regression?** A test that only checks array sizes or that outputs are nonzero is not a regression test — it's a smoke test at best.

Bad tests are actively harmful because they create false confidence. A small number of well-chosen tests is better than many weak ones.

### What makes a good test for numerical simulation code

Prefer tests that check **invariants** — properties that must hold regardless of the specific numerical values:

- **Conservation laws**: Stagnation enthalpy is conserved through an adiabatic nozzle. Mass fractions sum to 1. Element mass is conserved across combustion.
- **Thermodynamic identities**: Entropy is constant in isentropic flow. Gibbs free energy is minimized at equilibrium.
- **Ordering and bounding relationships**: Frozen flow Isp ≤ equilibrium flow Isp for the same conditions (equilibrium allows further energy extraction via recombination). Throat pressure < chamber pressure. Exit temperature < chamber temperature for an expanding nozzle. Mach number at throat ≈ 1.
- **Consistency checks**: The same problem solved via pressure ratio and area ratio should give consistent exit states. A `DilutedCombustor` with zero recirculation should match a plain `Combustor`.
- **Limiting cases**: At very high temperature, composition approaches fully dissociated. At very low area ratio (just past throat), exit conditions approach throat conditions.
- **Known solutions**: In some cases, reference solutions to algorithms or techniques are available in textbooks or papers. For compressible gas dynamics in particular, Zucrow, M. J., Hoffman, J. D. (1976). Gas dynamics : volume 2 has a wealth of worked numerical examples, particularly for the Method of Characteristics.

### What to avoid

- **Tests that only check convergence or that outputs exist.** `EXPECT_TRUE(result.converged)` alone proves nothing about correctness.
- **Tests that hardcode implementation details.** Don't check the number of solver iterations, internal array sizes that aren't part of the API contract, or intermediate variables.
- **Regression tests against magic numbers with no provenance.** If you hardcode an expected value, document where it comes from (CEA output, analytical formula, Cantera baseline). Include units.
- **Trivial identity tests.** Setting a temperature and reading it back tests Cantera, not Goddard. The same goes for simply verifying that `setState_TPX` followed by `temperature()` returns what you put in — that's Cantera's job.

## C++ Test Conventions

### Framework and structure

Tests use GTest and GMock. The test binary is `goddardTests`, built by `test/CMakeLists.txt`. New test files must be added to the `add_executable` list in that file.

```cpp
#include "goddard/nozzle.hpp"   // the component under test
#include "goddard/global.hpp"   // for setup_defaults()
#include "gtest/gtest.h"

// Use test fixtures for shared setup
class NozzleTests : public ::testing::Test {
protected:
    NozzleTests() {
        Goddard::setup_defaults();  // required for data file paths
        gas = Cantera::newSolution("h2o2.yaml", "ohmech");
        // Set up a realistic combustion product state
        gas->thermo()->setState_TPX(3000.0, 70.0 * Cantera::OneAtm,
                                     "H2O:0.7, H2:0.1, O2:0.05, OH:0.1, H:0.05");
        gas->thermo()->equilibrate("HP");
        gas->thermo()->saveState(inlet_state);
    }

    std::shared_ptr<Cantera::Solution> gas;
    std::vector<double> inlet_state;
};
```

### Floating point comparison

Use the project's `max_fp_error` helper for tolerance calculation:

```cpp
#include "goddard/utils.hpp"  // provides max_fp_error

// max_fp_error(expected, reltol, abstol) returns max(|expected| * reltol, abstol)
double tol = max_fp_error(expected_entropy, 1e-5, 1e-4);
EXPECT_NEAR(actual_entropy, expected_entropy, tol);
```

Choose tolerances thoughtfully. Solver tolerances are typically 1e-6 to 1e-8; test tolerances should be somewhat looser (1e-4 to 1e-5) to account for floating point accumulation, but not so loose that real errors slip through.

### Invariant test examples

**Energy conservation through a nozzle:**

```cpp
TEST_F(NozzleTests, StagnationEnthalpyConserved) {
    EquilibriumNozzle nozzle(*gas);
    double H_inlet = gas->thermo()->enthalpy_mass();

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 20.0);
    ASSERT_TRUE(results.throat.converged);

    // Stagnation enthalpy = static enthalpy + kinetic energy
    // At the inlet (zero velocity), H_stagnation = H_static
    // This must be conserved at every station in adiabatic flow
    EXPECT_NEAR(results.throat.H_stagnation, H_inlet,
                max_fp_error(H_inlet, 1e-5, 1e-4))
        << "Stagnation enthalpy must be conserved (first law)";
}
```

**Frozen flow bounds equilibrium flow:**

```cpp
TEST_F(NozzleTests, FrozenIspBoundsEquilibriumIsp) {
    gas->thermo()->restoreState(inlet_state);
    EquilibriumNozzle eq_nozzle(*gas);

    gas->thermo()->restoreState(inlet_state);
    FrozenNozzle frozen_nozzle(*gas);

    double area_ratio = 20.0;
    NozzleResults eq_results = eq_nozzle.solve(
        ExpansionType::SUPERSONIC_AREA_RATIO, area_ratio);
    NozzleResults frz_results = frozen_nozzle.solve(
        ExpansionType::SUPERSONIC_AREA_RATIO, area_ratio);

    ASSERT_TRUE(eq_results.expansions.front().converged);
    ASSERT_TRUE(frz_results.expansions.front().converged);

    // Extract exit velocities from enthalpy difference: v_e = sqrt(2 * (h0 - h_e))
    gas->thermo()->restoreState(eq_results.expansions.front().state);
    double h_exit_eq = gas->thermo()->enthalpy_mass();

    gas->thermo()->restoreState(frz_results.expansions.front().state);
    double h_exit_frz = gas->thermo()->enthalpy_mass();

    double H0 = eq_results.throat.H_stagnation;
    double v_exit_eq = std::sqrt(2.0 * (H0 - h_exit_eq));
    double v_exit_frz = std::sqrt(2.0 * (H0 - h_exit_frz));

    EXPECT_GE(v_exit_eq, v_exit_frz * 0.999)
        << "Equilibrium exit velocity should be >= frozen "
           "(recombination recovers energy in equilibrium flow)";
}
```

**Isentropic expansion:**

```cpp
TEST_F(NozzleTests, NozzleExpansionIsIsentropic) {
    EquilibriumNozzle nozzle(*gas);
    double inlet_entropy = gas->thermo()->entropy_mass();

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 15.0);
    ASSERT_TRUE(results.expansions.front().converged);

    gas->thermo()->restoreState(results.expansions.front().state);
    double exit_entropy = gas->thermo()->entropy_mass();

    EXPECT_NEAR(exit_entropy, inlet_entropy,
                max_fp_error(inlet_entropy, 1e-4, 1e-3))
        << "Nozzle expansion must be isentropic";
}
```

### Marking incomplete tests

If a test is a placeholder, a smoke test pending real assertions, or is temporarily set to always pass, mark it with a `TODO` comment explaining what's missing:

```cpp
TEST_F(CombustorTests, EquilibriumTemperaturesArePhysical) {
    auto results = combustor->solve(temperatures, pressures, *MRs, options);

    // TODO: expand this test to compare combustion temperatures against
    // adiabatic flame temperature bounds or CEA reference values.
    // Currently only checks that combustion occurred (T > T_initial).
    for (int i = 0; i < results.size(); i++) {
        sln->thermo()->restoreState(results.get_state(i));
        EXPECT_GT(sln->thermo()->temperature(), 300.0);
    }
}
```

Do the same if a test uses `GTEST_SKIP()` or `SUCCEED()` as a workaround:

```cpp
// TODO: This test is a workflow demonstration that always passes.
// Replace with actual numerical comparisons once the CEA loader
// fully populates thermodynamic states.
SUCCEED();
```

## Python Test Conventions

### Use the shared infrastructure in `conftest.py`

`conftest.py` provides `RocketTestCase`, `ComparisonTolerances`, `build_goddard_problem`, `solve_cea_problem`, `compare_thermo_states`, and helper assertions. Use these rather than writing ad hoc comparison code.

Define new test cases declaratively:

```python
RP1_LOX_EQUILIBRIUM = RocketTestCase(
    name="rp1_lox_eq",
    fuel=ReactantSpec(cea_name="RP-1", cantera_composition="C12H24:1", temperature=300.0),
    oxidizer=ReactantSpec(cea_name="O2(L)", cantera_composition="O2:1", temperature=90.0),
    of_ratio=2.6,
    chamber_pressure_bar=68.9,
    area_ratios=[10.0, 25.0],
    thermo_file="rp1_lox.yaml",
    phase_name="gas",
    species={"CO2", "CO", "H2O", "H2", "OH", "O2", "H", "O"},
    nozzle_chemistry="equilibrium",
)
```

### Comparing against CEA

Before writing or debugging a comparison, read [pycea.md](pycea.md): call signatures, station order, `n_frz` numbering, and CEA behaviors that look like bugs.

Every numerical comparison needs clear tolerances and labels:

```python
@pytest.mark.parametrize("case", EQUILIBRIUM_CASES, ids=lambda c: c.name)
def test_chamber_temperature(case: RocketTestCase):
    """Chamber temperature should match CEA within 1%."""
    problem = build_goddard_problem(case)
    results = problem.solve()
    chamber = results.get_chamber_state(case.name, 0)

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)

    tol = ComparisonTolerances()
    assert_close_rel(chamber.temperature, cea_sol.T[stations["chamber"]],
                     tol.temperature_rel, "chamber temperature")
```

### Performance metrics

Compute performance from first principles, not from convenience methods that may embed their own approximations:

```python
# c* from definition: c* = P_c / (rho_t * a_t)
goddard_cstar = chamber.pressure / (throat.density * throat.speed_of_sound)

# Isp from energy: v_e = sqrt(2*(h_c - h_e)), Isp = v_e / g0
exit_velocity = math.sqrt(2.0 * (chamber.enthalpy - exit_state.enthalpy))
```

### Marking incomplete Python tests

Same principle as C++. If a test is a stub or known to be weak:

```python
def test_frozen_gamma_at_throat(case):
    """Frozen gamma at throat should match CEA."""
    # TODO: enable once frozen nozzle gamma_s calculation is validated.
    # Currently skipped because of a known discrepancy in the frozen
    # gamma_s definition between Goddard and CEA.
    pytest.skip("pending gamma_s fix")
```

## Checklist Before Committing a Test

1. **Does the test name say what physical property it checks?** Good: `StagnationEnthalpyConserved`, `FrozenIspBoundsEquilibriumIsp`. Bad: `NozzleTest1`, `BasicCheck`.
2. **Would the test fail if the computation were wrong?** Mentally substitute a plausible bug (e.g., sign error in enthalpy, missing factor of 2) and check that the test would catch it.
3. **Are expected values documented?** If from CEA, say so. If from an analytical formula, write the formula in a comment. If from a previous Cantera run, say "Cantera baseline" and note the version.
4. **Are tolerances justified?** Don't use `1e-2` because it's round. The solver converges to `reltol=1e-6`, so test at `1e-4` to `1e-5` with a rationale for the margin.
5. **Is the test complete?** If it's a stub or smoke test, add a `TODO` comment describing what assertions should be added later.
6. **Right language?** Numerical comparisons against CEA or literature → Python. Invariants, convergence, state management → C++.

## Build and Run

New tests should be run and confirmed to pass, or explicitly annotated with `TODO` comments as described previously.

```bash
pixi shell -e test

# C++ tests
pixi run compile
pixi run ctest
# or directly:
cmake --build build --target goddardTests
./build/test/goddardTests

# Python tests (requires editable install)
pixi run pip-install
pixi run test-python
# or directly:
cd python && pytest tests/ -v
```
