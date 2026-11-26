# CEA as a test oracle

Goddard's comparison tests use NASA's `cea` Python package (3.3.4, in the `test` pixi environment) as the reference. CEA is not correct everywhere, and some of its failures are silent. This page lists the defects found so far, how the tests detect them, and how they work around them. `python/tests/test_propellant_cea.py` is the main user.

## Known defects of the oracle

| Defect | Where it was found | What CEA returns | Workaround in the tests |
|---|---|---|---|
| A frozen expansion stops "more than 50 K outside the range of a condensed species" although the species has data at that temperature | Fuel-rich RP-1/LOX at 5 bar (O/F 0.25–0.6), graphite in the products | The stations before the stop, then nothing; `converged` is False | Compare only the stations CEA produced; the missing ones are listed in `KNOWN_CEA_ROCKET_FAILURES` |
| An equilibrium expansion with graphite at low pressure fails with "Too many singular matrices" | Fuel-rich RP-1/LOX at 5 bar, O/F 0.25 and 0.3 | NaN or non-physical stations (3800 K, or a temperature that rises along the nozzle) while `converged` is True | Listed in `KNOWN_CEA_ROCKET_FAILURES`; the exits are checked against CEA's TP equilibrium instead |
| A TP equilibrium solved afresh does not converge where the rocket solver, continuing from the previous station, does | LF2/LH2 at O/F 9 and 70 bar, area ratio 100 (804 K, 0.0286 bar), even at CEA's own station state | `converged` False | That exit is compared with CEA's rocket station only |

The frozen stops that are a defect fit one rule. Freezing above 600 K, CEA stops every expansion that reaches 550 K or below: all failing exits are at 464–509 K, and every exit CEA solves is at 584 K or above. 600 K is a break between the temperature intervals of the graphite data. That is consistent with CEA treating each interval of a condensed species as its own species and applying its 50 K rule to the interval it froze in. This is inferred from the pattern, not read from CEA's source. Goddard treats C(gr) as one species with data from 200 to 6000 K and continues.

Where CEA's rocket solver can be used, Goddard agrees with it to about 1e-7 at the chamber and the throat and about 1e-6 at the exits in the cases examined (fuel-rich RP-1/LOX, LF2/LH2).

## Frozen expansions stop at CEA's data limits, and drop the stations after the stop

CEA stops a frozen expansion at a temperature below 0.8 T_min of a gas species ("temperature below lower-limit (0.8*Tmin=160.000)" for data starting at 200 K) or 50 K below the range of a condensed one, and does not compute any later station. That is the same kind of rule as Goddard's range check (`Gas::check_temperature_range`), and not a defect, but it has a trap. pycea requires at least one pressure-ratio station (`pi_p`), and CEA solves those before the area-ratio stations. A placeholder such as Pc/P = 1000, used only to satisfy the requirement, can take the frozen flow below the limit, and every area-ratio station is then silently missing. That happened for N2H4 and fuel-rich RP-1/LOX. The tests use Pc/P = 3, which stays warm.

## How the tests detect a failed oracle

- `converged` is not enough. `cea_problem` checks each exit station CEA returns: a station must exist at that area ratio (within 1e-3, not merely the nearest one), be finite, and be colder than the throat.
- An exit station that fails the check fails the test, unless it is listed in `KNOWN_CEA_ROCKET_FAILURES` with the reason. A listed exit that CEA then solves also fails the test, so the list cannot go stale when CEA is fixed.
- `test_propellant_exit_entropy` checks every equilibrium exit independently of CEA's rocket solver. CEA's TP equilibrium at Goddard's exit temperature and pressure must have the chamber's entropy, to 1e-6.

## Differences that are not defects

- **Atomic weights.** CEA uses H = 1.00794, O = 15.9994, N = 14.0067; Cantera's default table has the newer 1.008, 15.999, 14.007. Near a stoichiometric mixture this moves the leftover minor species by about 1%: it moved exit temperatures by 2e-4 for LF2/LH2 at O/F 19. The NASA-9 data files carry CEA's weights, so Goddard uses them too (see `scripts/parse_nasa_data.py`).
- **Gas constant.** CEA uses 8.314510 J/(mol K) and Cantera 8.314462618, a 5.7e-6 offset in mass-basis properties. It is built into Cantera and remains.
- **Reference pressure.** CEA's NASA-9 data are referenced to 1 bar, Cantera assumes 1 atm unless told otherwise. The data files set 1 bar.

## Calling the oracle

Practical details of the `cea` package (crashes that need a forked process, call signatures, station order, `n_frz`) are in `.claude/skills/test-creation/pycea.md`. A running record of everything found while extending the propellant coverage, in Goddard and its data as well as in CEA, is kept in `instructions/propellant_coverage_report.md` (a local working document). Two more found here: `RocketSolver.solve` rejects an empty `supar` list (omit the argument instead), and `EqSolution.entropy` and `RocketSolution.entropy` are in kJ/(kg K).
