
import goddard


bar = 1e5

fuel_state = goddard.ThermodynamicState(
    T = 273.15,
    P = 101325.0,
    composition="H2:1"
)

oxidizer_state = goddard.ThermodynamicState(
    T = 273.15,
    P = 101325.0,
    composition="O2:1"
)

params = goddard.ChemicalParameters(
    thermo_file = "h2o2.yaml",
    species =  {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"},
    OF_ratios = [3.0, 6.0, 9.0],
    cantera_fuel_state = fuel_state,
    cantera_oxidizer_state = oxidizer_state
)


combust_opts = goddard.CombustorOptions(
    type=goddard.CombustorType.INFINITE_AREA,
    pressures=[10*bar, 50*bar],
)

nozzle_opts = goddard.NozzleOptions(
    chemistry=goddard.GasChemistry.EQUILIBRIUM,
    expansion_type=goddard.ExpansionType.SUPERSONIC_AREA_RATIO,
    expansion_ratios=[2.0, 5.0, 10.0, 50.0],
)

case1 = goddard.RocketCaseParameters(
    name = "h2o2eq",
    problem_type = "",
    combustor_options = combust_opts,
    nozzle_options = nozzle_opts
)

cases = [
    case1
]


problem = goddard.RocketProblem(
    params,
    cases,
    name = "ohmech"
)

result = problem.solve()

print(result.report())