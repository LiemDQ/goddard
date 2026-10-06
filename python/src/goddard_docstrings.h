/*
  This file contains docstrings for use in the Python bindings.
  Do not edit! They were automatically extracted by pybind11_mkdoc.
 */

#define __EXPAND(x)                                      x
#define __COUNT(_1, _2, _3, _4, _5, _6, _7, COUNT, ...)  COUNT
#define __VA_SIZE(...)                                   __EXPAND(__COUNT(__VA_ARGS__, 7, 6, 5, 4, 3, 2, 1))
#define __CAT1(a, b)                                     a ## b
#define __CAT2(a, b)                                     __CAT1(a, b)
#define __DOC1(n1)                                       __doc_##n1
#define __DOC2(n1, n2)                                   __doc_##n1##_##n2
#define __DOC3(n1, n2, n3)                               __doc_##n1##_##n2##_##n3
#define __DOC4(n1, n2, n3, n4)                           __doc_##n1##_##n2##_##n3##_##n4
#define __DOC5(n1, n2, n3, n4, n5)                       __doc_##n1##_##n2##_##n3##_##n4##_##n5
#define __DOC6(n1, n2, n3, n4, n5, n6)                   __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6
#define __DOC7(n1, n2, n3, n4, n5, n6, n7)               __doc_##n1##_##n2##_##n3##_##n4##_##n5##_##n6##_##n7
#define DOC(...)                                         __EXPAND(__EXPAND(__CAT2(__DOC, __VA_SIZE(__VA_ARGS__)))(__VA_ARGS__))

#if defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#endif


static const char *__doc_Goddard_BaseCombustor =
R"doc(Base class for isobaric (HP) and isochoric (UV) combustion. Provides
shared utilities for stream mixing and equilibration.)doc";

static const char *__doc_Goddard_BaseCombustor_BaseCombustor =
R"doc(Parameter ``gas``:
    product gas. The solver works on its own copy of `gas`; the
    caller's `Gas` is not modified. */)doc";

static const char *__doc_Goddard_BaseCombustor_combust =
R"doc(Equilibrate reactant states in place and return them.

Throws:
    NotImplementedError for `CombustionProcess::ISOCHORIC` when
    `states` carries candidate condensed species.)doc";

static const char *__doc_Goddard_BaseCombustor_get_combustion_species = R"doc()doc";

static const char *__doc_Goddard_BaseCombustor_m_gas = R"doc()doc";

static const char *__doc_Goddard_BaseCombustor_reactant_states =
R"doc(Array of reactant states with the given shape, carrying the candidate
condensed species of `m_gas` with every amount set to zero (the
reactants are all gas phase). Sets the condensed amounts of `m_gas` to
zero as well.)doc";

static const char *__doc_Goddard_BaseCombustor_require_nonempty =
R"doc(Reject an empty input grid before any entry of it is read.

Parameter ``values``:
    Grid of input values, e.g. pressures [Pa] or mixture ratios [-].

Parameter ``name``:
    Name of the grid, for the error message.

Throws:
    std::invalid_argument if `values` is empty.)doc";

static const char *__doc_Goddard_BaseCombustor_set_mixture_composition = R"doc()doc";

static const char *__doc_Goddard_BaseCombustor_validate_options =
R"doc(Check that the combustor options describe a supported combustor.

The combustor itself is the same for every `CombustorType`: for
`FINITE_MASS_FLUX` and `FINITE_CONTRACTION_RATIO` it produces the
injector-face state, and the finite-area chamber is solved afterwards
by `Nozzle::solve_finite_area_chamber`.

Throws:
    std::invalid_argument for a finite-area type with
    `contraction_ratio` <= 1 or `mass_flux` <= 0, or combined with
    `CombustionProcess::ISOCHORIC`.)doc";

static const char *__doc_Goddard_ChainMetadata =
R"doc(Bookkeeping for one characteristic chain in a CharacteristicNet.

A chain is a single characteristic line, stored as an ordered list of
point indices. It is born at the initial data line or at a wall/axis
reflection, and stops when it terminates.)doc";

static const char *__doc_Goddard_ChainMetadata_active =
R"doc(False once the chain has terminated; only active chains are still
marched.)doc";

static const char *__doc_Goddard_ChainMetadata_family = R"doc(Which characteristic family this chain belongs to.)doc";

static const char *__doc_Goddard_ChainMetadata_latest_point_idx =
R"doc(Index into CharacteristicNet::points of the chain's leading (most
downstream) point.)doc";

static const char *__doc_Goddard_ChainMetadata_origin_point_idx = R"doc(Index into CharacteristicNet::points of the chain's first point.)doc";

static const char *__doc_Goddard_ChainMetadata_termination = R"doc(How this chain stopped being marched.)doc";

static const char *__doc_Goddard_ChainTermination = R"doc(How a characteristic chain stopped being marched. */)doc";

static const char *__doc_Goddard_ChainTermination_AXIS = R"doc()doc";

static const char *__doc_Goddard_ChainTermination_NOT_TERMINATED = R"doc()doc";

static const char *__doc_Goddard_ChainTermination_WALL = R"doc()doc";

static const char *__doc_Goddard_CharacteristicFamily = R"doc(Which family a characteristic belongs to: C+ (PLUS) or C- (MINUS). */)doc";

static const char *__doc_Goddard_CharacteristicFamily_MINUS = R"doc()doc";

static const char *__doc_Goddard_CharacteristicFamily_PLUS = R"doc()doc";

static const char *__doc_Goddard_CharacteristicFamily_UNSPECIFIED = R"doc()doc";

static const char *__doc_Goddard_CharacteristicNet =
R"doc(The method-of-characteristics flow field mesh: every solved point, the
characteristic chains and marching fronts that connect them, and the
wall/axis boundary bookkeeping. Built up incrementally by the two
kernels (the chain-pairing ladder for minimum-length design, the
reference-plane march for analysis and Rao design) via the mutation
methods below; queried afterward via the columnar point data and the
boundary/chain accessors.)doc";

static const char *__doc_Goddard_CharacteristicNet_active_chains =
R"doc(View of the metadata of active chains that belong to a given family
(or all active chains, if `family` is UNSPECIFIED).)doc";

static const char *__doc_Goddard_CharacteristicNet_add_front =
R"doc(Add a marching front's points (axis to wall) to the net: no chain
membership (the front-based kernels do not use c_chains), but
registers the front's axis and wall points and records it in `fronts`.

Returns:
    Indices into `points` of the added front, axis to wall.)doc";

static const char *__doc_Goddard_CharacteristicNet_add_initial_characteristic =
R"doc(Add points from an initial data line that happens to be along a
characteristic. This is primarily used when the initial data line
originates from a centered expansion.)doc";

static const char *__doc_Goddard_CharacteristicNet_add_initial_data_line =
R"doc(Seed the net from an initial data line.

Handles both data-line topologies. When `on_characteristic` is empty
the points are assumed to lie on distinct characteristics (a non-
collinear transonic start line, e.g. Kliegel-Levine): every interior
point seeds both a C+ and a C- chain, mirroring the fan-init topology,
and the last point (already at the wall) immediately reflects into the
first C- chain, exactly as a wall reflection would during marching.
The first point (the axis bootstrap) seeds only a C+: giving it a C-
would immediately re-reflect it off the axis as a degenerate point.
When `on_characteristic` is set the points are collinear along a
single characteristic of that family (e.g. a centered expansion fan),
and are seeded as one shared chain via add_initial_characteristic.

The caller must supply the topology; it cannot be inferred from the
points, since a Riemann invariant is constant along a characteristic
only for planar flow.)doc";

static const char *__doc_Goddard_CharacteristicNet_add_initialization_point =
R"doc(Add a point along an initial data line.

Returns:
    index of the added point and corresponding characteristics)doc";

static const char *__doc_Goddard_CharacteristicNet_add_point =
R"doc(Add a point and its membership to list of points, and appends it to
its member characteristic chains.

Returns:
    index of the added point)doc";

static const char *__doc_Goddard_CharacteristicNet_axis_point_indices =
R"doc(Indices into `points` of the points lying on the centerline, in march
order. */)doc";

static const char *__doc_Goddard_CharacteristicNet_axis_points =
R"doc(Every point at `axis_point_indices`, in march order. Kernel-
independent: populated for both the chain-pairing ladder and the
front-based (analysis/Rao) kernels.)doc";

static const char *__doc_Goddard_CharacteristicNet_c_chains =
R"doc(Each entry is the chain of point indices along one characteristic.
Note that wall and axis points terminate a chain, and also start the
next chain.)doc";

static const char *__doc_Goddard_CharacteristicNet_chain_metadata = R"doc(Metadata for each chain in `c_chains`, in the same order. */)doc";

static const char *__doc_Goddard_CharacteristicNet_create_chain =
R"doc(Create a new characteristic with `pt_idx` as the starting point, and
add it to the tracking lists.

Returns:
    Index of the chain)doc";

static const char *__doc_Goddard_CharacteristicNet_empty = R"doc(True when the net holds no points at all. */)doc";

static const char *__doc_Goddard_CharacteristicNet_fronts =
R"doc(Point indices of every marching front, axis to wall, in order.
Populated by the front-based nets (analysis and Rao design):
`fronts.front()` is the initial front F_0 and `fronts.back()` is the
exit plane once the march has completed. Empty for the chain-pairing
ladder (minimum-length design), which has no synchronized front -- its
topology lives entirely in `c_chains`.)doc";

static const char *__doc_Goddard_CharacteristicNet_has_active_chains = R"doc(True while at least one chain is still being marched. */)doc";

static const char *__doc_Goddard_CharacteristicNet_leading_axis_point = R"doc(Get the leading axis point.)doc";

static const char *__doc_Goddard_CharacteristicNet_leading_axis_point_2 = R"doc()doc";

static const char *__doc_Goddard_CharacteristicNet_leading_point = R"doc(Get the leading point of the chain at the given index.)doc";

static const char *__doc_Goddard_CharacteristicNet_leading_point_2 = R"doc()doc";

static const char *__doc_Goddard_CharacteristicNet_leading_wall_point = R"doc(Get the leading wall point.)doc";

static const char *__doc_Goddard_CharacteristicNet_leading_wall_point_2 = R"doc()doc";

static const char *__doc_Goddard_CharacteristicNet_membership =
R"doc(Given a point index, which chains does it belong to? `membership[i]`
describes `points[i]`.)doc";

static const char *__doc_Goddard_CharacteristicNet_points =
R"doc(Every point in the net, in creation order. The index is opaque; the
flow-field topology lives in `c_chains`, and the boundaries in
`wall_point_indices`/`axis_point_indices`.)doc";

static const char *__doc_Goddard_CharacteristicNet_push_chain =
R"doc(Add a characteristic chain and its metadata to be tracked.

Returns:
    Index of the chain)doc";

static const char *__doc_Goddard_CharacteristicNet_reflect_c_minus_off_axis =
R"doc(Terminate a C- characteristic off the central axis and generate a new
reflected C+ characteristic.

Returns:
    Index of axis point, index of new C+ chain)doc";

static const char *__doc_Goddard_CharacteristicNet_reflect_c_plus_off_wall =
R"doc(Terminate a C+ characteristic off a wall and generate a new reflected
C- characteristic.

Returns:
    Index of wall point, index of new C- chain)doc";

static const char *__doc_Goddard_CharacteristicNet_seed_wall_point =
R"doc(Seed an initial wall point (e.g. the throat lip) that anchors the wall
march. The point owns no characteristic chain; it only bootstraps
leading_wall_point() and the wall coordinate lists.

Returns:
    index of the seeded point)doc";

static const char *__doc_Goddard_CharacteristicNet_terminate_c_plus_at_wall =
R"doc(Add a point where the C+ characteristic hits the wall without emitting
a reflected C- characteristic. This is used when designing minimum
length nozzles.)doc";

static const char *__doc_Goddard_CharacteristicNet_terminate_chain = R"doc(Mark a chain as inactive. */)doc";

static const char *__doc_Goddard_CharacteristicNet_update_and_terminate_chain = R"doc(Mark a chain as inactive while updating the last point. */)doc";

static const char *__doc_Goddard_CharacteristicNet_wall_point_indices =
R"doc(Indices into `points` of the points lying on the nozzle wall, in march
order. */)doc";

static const char *__doc_Goddard_CharacteristicNet_wall_points =
R"doc(Every point at `wall_point_indices`, in march order. Kernel-
independent: populated for both the chain-pairing ladder and the
front-based (analysis/Rao) kernels.)doc";

static const char *__doc_Goddard_CharacteristicNet_wall_x =
R"doc(Axial coordinates of the wall points, in length units. The net owns
the invariant that `wall_x`/`wall_y` stay parallel to
`wall_point_indices` (same length, same order); every method that
appends to `wall_point_indices` appends to these too.)doc";

static const char *__doc_Goddard_CharacteristicNet_wall_y =
R"doc(Radial coordinates of the wall points, in length units. See `wall_x`.
*/)doc";

static const char *__doc_Goddard_CharacteristicPoint =
R"doc(One node of the characteristic mesh: the flow state at a point, plus
the two Riemann invariants carried through it.

For perfect-gas solves the thermodynamic quantities are normalized by
their stagnation values and `V` holds the Mach number; for frozen and
equilibrium chemistry they are dimensional SI.

Plain data: the thermodynamic state is filled in by MocThermo (see
moc_thermo.hpp), not by a method on this type.)doc";

static const char *__doc_Goddard_CharacteristicPoint_K_minus = R"doc(< Riemann invariant carried along the C- characteristic.)doc";

static const char *__doc_Goddard_CharacteristicPoint_K_plus = R"doc(< Riemann invariant carried along the C+ characteristic.)doc";

static const char *__doc_Goddard_CharacteristicPoint_V = R"doc(< Velocity in m/s; equal to the Mach number for perfect gas.)doc";

static const char *__doc_Goddard_CharacteristicPoint_cantera_state =
R"doc(Serialized Cantera state; empty for perfect-gas solves. Internal
detail.)doc";

static const char *__doc_Goddard_CharacteristicPoint_gamma_s = R"doc(< Local isentropic exponent.)doc";

static const char *__doc_Goddard_CharacteristicPoint_mach = R"doc(< Local Mach number.)doc";

static const char *__doc_Goddard_CharacteristicPoint_mu =
R"doc(Mach angle asin(1/M), in radians. Sets the characteristic slopes theta
+/- mu.)doc";

static const char *__doc_Goddard_CharacteristicPoint_nu = R"doc(< Prandtl-Meyer angle (or generalized PM function), in radians.)doc";

static const char *__doc_Goddard_CharacteristicPoint_pressure =
R"doc(< Static pressure; Pa, or normalized by the stagnation pressure for
perfect gas.)doc";

static const char *__doc_Goddard_CharacteristicPoint_temperature =
R"doc(< Static temperature; K, or normalized by the stagnation temperature
for perfect gas.)doc";

static const char *__doc_Goddard_CharacteristicPoint_theta = R"doc(< Flow angle, in radians.)doc";

static const char *__doc_Goddard_CharacteristicPoint_update_Ks = R"doc(Set K_plus/K_minus from the current theta and nu. */)doc";

static const char *__doc_Goddard_CharacteristicPoint_x = R"doc(< Axial position, in length units.)doc";

static const char *__doc_Goddard_CharacteristicPoint_y = R"doc(< Radial (or transverse) position, in length units.)doc";

static const char *__doc_Goddard_ChemicalParameters = R"doc()doc";

static const char *__doc_Goddard_ChemicalParameters_OF_ratios =
R"doc(Mixture ratios of the problem, interpreted through
`CombustorOptions::mixture_type` of each case (O/F ratios by default).)doc";

static const char *__doc_Goddard_ChemicalParameters_all_condensed_species =
R"doc(Offer every species of `condensed_file` whose elements the product
phase has. */)doc";

static const char *__doc_Goddard_ChemicalParameters_condensed_file =
R"doc(Cantera YAML file holding candidate condensed product species, e.g.
`data/nasa9_condensed.yaml`. */)doc";

static const char *__doc_Goddard_ChemicalParameters_condensed_species =
R"doc(Condensed species of `condensed_file` to offer as candidates. Ignored
if `all_condensed_species`. */)doc";

static const char *__doc_Goddard_ChemicalParameters_fuel_state =
R"doc(Fuel stream state. `composition` is always a **mole-fraction** map.
Its keys are product species names when `reactant_file` is empty, and
species of `reactant_file` otherwise (e.g. `H2(L)`, `RP-1`). CEA-style
weight percentages must be converted to mole fractions by the caller.)doc";

static const char *__doc_Goddard_ChemicalParameters_fuel_weight_percentages =
R"doc(Not implemented: `RocketProblem::solve` raises `NotImplementedError`
if non-empty. */)doc";

static const char *__doc_Goddard_ChemicalParameters_mixtures =
R"doc(Not implemented: `RocketProblem::solve` raises `NotImplementedError`
if non-empty. */)doc";

static const char *__doc_Goddard_ChemicalParameters_oxidizer_state =
R"doc(Oxidizer stream state; see `fuel_state` for how `composition` is
interpreted. */)doc";

static const char *__doc_Goddard_ChemicalParameters_phi_ratios =
R"doc(Not implemented: `RocketProblem::solve` raises `NotImplementedError`
if non-empty. */)doc";

static const char *__doc_Goddard_ChemicalParameters_reactant_file =
R"doc(Cantera YAML file holding the reactant species, e.g.
`data/nasa9_reactants.yaml`.

When set, the fuel and oxidizer are built as separate `Gas` streams
from this file and their element amounts and enthalpies are
transferred to the products, so reactants need not be product species.
When empty, the reactant compositions refer to product species and the
combustor blends them in the product phase, as before.)doc";

static const char *__doc_Goddard_ChemicalParameters_species = R"doc()doc";

static const char *__doc_Goddard_ChemicalParameters_thermo_file = R"doc()doc";

static const char *__doc_Goddard_CombustionProcess =
R"doc(Thermodynamic constraint held fixed while the reactants burn to
equilibrium.)doc";

static const char *__doc_Goddard_CombustionProcess_ISOBARIC = R"doc(< Constant enthalpy and pressure (HP).)doc";

static const char *__doc_Goddard_CombustionProcess_ISOCHORIC = R"doc(< Constant internal energy and specific volume (UV).)doc";

static const char *__doc_Goddard_Combustor = R"doc(Handles isobaric combustion reactions with fuel and oxidizer streams.)doc";

static const char *__doc_Goddard_CombustorOptions = R"doc()doc";

static const char *__doc_Goddard_CombustorOptions_contraction_ratio =
R"doc(Contraction ratio A_c/A_t [-] of a `FINITE_CONTRACTION_RATIO` chamber.
*/)doc";

static const char *__doc_Goddard_CombustorOptions_mass_flux =
R"doc(Mass flux through a `FINITE_MASS_FLUX` chamber, mdot/A_c [kg/(m^2 s)].
*/)doc";

static const char *__doc_Goddard_CombustorOptions_mixture_type = R"doc()doc";

static const char *__doc_Goddard_CombustorOptions_pressures =
R"doc(Pressures [Pa]. For `CombustionProcess::ISOBARIC` these are the
chamber pressures. For `CombustionProcess::ISOCHORIC` they are the
initial pressures of the unburnt reactants; the chamber pressure is
the result of the constant-volume combustion.)doc";

static const char *__doc_Goddard_CombustorOptions_process =
R"doc(Combustion constraint. With `ISOCHORIC` in a `RocketProblem`, the
constant-volume equilibrium state is used as the stagnation state of a
steady isentropic nozzle expansion, so the reported performance is
that idealization. It is not a Chapman-Jouguet detonation model.)doc";

static const char *__doc_Goddard_CombustorOptions_type = R"doc()doc";

static const char *__doc_Goddard_CombustorType = R"doc()doc";

static const char *__doc_Goddard_CombustorType_FINITE_CONTRACTION_RATIO = R"doc()doc";

static const char *__doc_Goddard_CombustorType_FINITE_MASS_FLUX = R"doc()doc";

static const char *__doc_Goddard_CombustorType_INFINITE_AREA = R"doc()doc";

static const char *__doc_Goddard_Combustor_Combustor = R"doc()doc";

static const char *__doc_Goddard_Combustor_Combustor_2 = R"doc()doc";

static const char *__doc_Goddard_Combustor_Combustor_3 =
R"doc(Construct from reactant streams given as `Gas` objects.

The fuel and oxidizer streams carry their own temperature, pressure
and composition, set beforehand with
`Gas::set_state_TPX`/`set_state_TPY`. Their species need not be
product species: element amounts [kmol per kg of stream] and specific
enthalpies [J/kg] are transferred to `products` by element name. A
stream built from a condensed reactant such as `H2(L)` or `RP-1` is an
ideal-gas phase whose NASA9 polynomials give the correct enthalpy and
element amounts; its density and entropy are meaningless and are only
read on the isochoric path, which is therefore restricted to gaseous
reactants.

Parameter ``products``:
    Product gas, which also defines the element set of the problem.

Parameter ``fuel``:
    Fuel stream at its own state.

Parameter ``oxidizer``:
    Oxidizer stream at its own state.)doc";

static const char *__doc_Goddard_Combustor_generate_mass_fraction_matrix = R"doc()doc";

static const char *__doc_Goddard_Combustor_generate_mole_fraction_matrix = R"doc()doc";

static const char *__doc_Goddard_Combustor_m_fuel_composition = R"doc()doc";

static const char *__doc_Goddard_Combustor_m_fuel_gas = R"doc()doc";

static const char *__doc_Goddard_Combustor_m_oxidizer_composition = R"doc()doc";

static const char *__doc_Goddard_Combustor_m_oxidizer_gas = R"doc()doc";

static const char *__doc_Goddard_Combustor_solve = R"doc()doc";

static const char *__doc_Goddard_Combustor_solve_2 = R"doc()doc";

static const char *__doc_Goddard_Combustor_solve_3 =
R"doc(Burn the reactant streams given to the `Gas`-stream constructor over a
grid of pressures and mixture ratios.

For every mixture ratio the fuel mass fraction f follows from
`options.mixture_type` (`OF_RATIO`: f = 1/(1+O/F); `FUEL_FRAC`: f is
the value itself). The element amounts b [kmol/kg of mixture] and the
specific enthalpy h [J/kg] of the mixture are the mass blends f*x_fuel
+ (1-f)*x_oxidizer of the two stream states, with the element amounts
mapped by element name onto the element order of the product gas. Each
combination is then equilibrated at constant enthalpy and pressure
(`ISOBARIC`) or at constant internal energy and specific volume
(`ISOCHORIC`), warm-started from the previous pressure of the same
mixture ratio.

Parameter ``pressures``:
    Pressures [Pa]: chamber pressures for `ISOBARIC`, initial reactant
    pressures for `ISOCHORIC`.

Parameter ``mixture_ratios``:
    Mixture ratios, interpreted according to `options.mixture_type`.

Parameter ``options``:
    Combustor settings.

Returns:
    Array of combustion states with shape (mixture ratios, pressures).

Throws:
    NotImplementedError if the combustor was not built from reactant
    `Gas` streams, if `options.mixture_type` is `PHI_RATIO` (CEA's
    valence rule is not implemented on this path), if `options.type`
    is `NONE`, or if the process is `ISOCHORIC` while the product gas
    carries candidate condensed species.

Throws:
    FmtError if a reactant contains an element the product gas does
    not have.

Throws:
    std::invalid_argument if `pressures` or `mixture_ratios` is empty.)doc";

static const char *__doc_Goddard_Combustor_stream_element_moles =
R"doc(Element amounts of a reactant stream [kmol per kg of stream] re-
indexed onto the element order of the product gas. Product elements
the stream does not have are zero.)doc";

static const char *__doc_Goddard_CondensedPhaseSet =
R"doc(Candidate condensed species of a `Gas`. Implementation detail; see
`condensed.hpp`. */)doc";

static const char *__doc_Goddard_CondensedPhaseSet_2 = R"doc()doc";

static const char *__doc_Goddard_DetonationBranch =
R"doc(Root of the Rayleigh line and the equilibrium Hugoniot for a wave
faster than the Chapman-Jouguet (CJ) speed.

Above the CJ speed the Rayleigh line cuts the Hugoniot twice; at the
CJ speed the two roots merge, and below it there is no steady
solution. Both branches therefore take a drive factor u1/u_CJ > 1.)doc";

static const char *__doc_Goddard_DetonationBranch_OVERDRIVEN =
R"doc(< Strong root: higher pressure, products subsonic relative to the
wave.)doc";

static const char *__doc_Goddard_DetonationBranch_UNDERDRIVEN = R"doc(< Weak root: lower pressure, products supersonic relative to the wave.)doc";

static const char *__doc_Goddard_DetonationResult =
R"doc(Jump conditions across a planar detonation, in the frame of the wave.

State 1 is the unburned gas ahead of the wave and state 2 the products
at chemical equilibrium behind it. Ratios are state 2 over state 1.

An invalid result (a wave slower than the CJ speed, or a mixture that
releases no heat) has `valid = false` and every other field set to -1.)doc";

static const char *__doc_Goddard_DetonationResult_density_ratio =
R"doc(Density ratio rho2/rho1 [-]. It is also the ratio u1/u2 of the
velocities relative to the wave. */)doc";

static const char *__doc_Goddard_DetonationResult_drive_factor =
R"doc(Wave speed over the CJ speed of the mixture [-]. 1 for a CJ
detonation. */)doc";

static const char *__doc_Goddard_DetonationResult_mach_in =
R"doc(Wave Mach number u1/a1, with the frozen sound speed of the unburned
gas [-]. */)doc";

static const char *__doc_Goddard_DetonationResult_mach_out =
R"doc(Mach number of the products relative to the wave, u2/a2, with the
equilibrium sound speed [-]. 1 for a CJ detonation, below 1 when
overdriven and above 1 when under-driven.)doc";

static const char *__doc_Goddard_DetonationResult_molecular_weight_ratio = R"doc(Molecular weight ratio M2/M1 [-], with M = 1/n the gas molar mass. */)doc";

static const char *__doc_Goddard_DetonationResult_static_pressure_ratio = R"doc(Static pressure ratio P2/P1 [-]. */)doc";

static const char *__doc_Goddard_DetonationResult_static_temperature_ratio = R"doc(Static temperature ratio T2/T1 [-]. */)doc";

static const char *__doc_Goddard_DetonationResult_total_pressure_ratio =
R"doc(Stagnation pressure ratio P02/P01 [-], in the frame of the wave. P01
is the frozen stagnation pressure of the unburned gas and P02 the
equilibrium stagnation pressure of the products.)doc";

static const char *__doc_Goddard_DetonationResult_valid = R"doc(False if the inputs admit no detonation. */)doc";

static const char *__doc_Goddard_DetonationResult_velocity =
R"doc(Wave speed relative to the unburned gas [m/s]. NaN from the perfect-
gas functions, which are dimensionless.)doc";

static const char *__doc_Goddard_DetonationResult_von_neumann =
R"doc(Frozen (non-reacting) shock at the same wave speed: the von Neumann
spike at the head of the ZND reaction zone. For an under-driven
detonation it is not on a steady ZND path.)doc";

static const char *__doc_Goddard_DetonationSolver =
R"doc(Planar detonations in a real gas: Chapman-Jouguet (CJ) detonations
following Gordon & McBride, NASA RP-1311 Part I, chapter 8, and
overdriven and under-driven detonations at a given wave speed.

The unburned gas (state 1) is the state of the `Gas` at construction.
It is taken as given: it need not be at equilibrium, and its sound
speed is the frozen one. The products are always at chemical
equilibrium, whatever the chemistry of the `Gas`. The CJ detonation is
solved once, at construction.

Every result also reports the von Neumann state, the frozen shock at
the same wave speed.

The solver works on its own copy of the `Gas` passed to the
constructor, and the state accessors return independent copies, so
neither the caller's `Gas` nor an earlier accessor result changes when
the solver runs.)doc";

static const char *__doc_Goddard_DetonationSolver_DetonationSolver =
R"doc(Parameter ``gas``:
    unburned gas: FROZEN or EQUILIBRIUM chemistry, with no condensed
    species

Parameter ``options``:
    Newton tolerance on the log pressure and temperature ratios, and
    iteration limit

Throws:
    std::invalid_argument for PERFECT_GAS or KINETIC chemistry, or if
    condensed species are present

Throws:
    ConvergenceError if the CJ detonation cannot be solved)doc";

static const char *__doc_Goddard_DetonationSolver_chapman_jouguet =
R"doc(Chapman-Jouguet detonation: the slowest steady detonation, whose
products leave at their equilibrium sound speed.

The result is the one solved at construction.
`post_detonation_state()` and `von_neumann_state()` are then the CJ
states.

Returns:
    jump conditions; invalid if the mixture releases no heat)doc";

static const char *__doc_Goddard_DetonationSolver_detonation =
R"doc(Overdriven or under-driven detonation at a multiple of the CJ speed.

The branch is chosen by the initial guess, and checked against the
Mach number of the products after the solve.

Parameter ``drive_factor``:
    wave speed over the CJ speed [-]

Parameter ``branch``:
    root of the Rayleigh line and the Hugoniot

Returns:
    jump conditions; invalid if drive_factor < 1. A drive factor of 1
    returns the CJ detonation.

Throws:
    std::invalid_argument if condensed species form behind the wave

Throws:
    ConvergenceError if the jump conditions cannot be solved on the
    requested branch, e.g. for a drive factor so close to 1 that the
    two branches merge)doc";

static const char *__doc_Goddard_DetonationSolver_detonation_from_velocity =
R"doc(Overdriven or under-driven detonation at a given wave speed.

Parameter ``velocity``:
    wave speed relative to the unburned gas [m/s]

Parameter ``branch``:
    root of the Rayleigh line and the Hugoniot

Returns:
    jump conditions; invalid if the velocity is below the CJ speed)doc";

static const char *__doc_Goddard_DetonationSolver_m_chapman_jouguet = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_chapman_jouguet_state = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_chapman_jouguet_von_neumann_state = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_chemistry = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_gas = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_options = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_post_detonation_state = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_pre_detonation_state = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_m_von_neumann_state = R"doc()doc";

static const char *__doc_Goddard_DetonationSolver_post_detonation_state =
R"doc(Gas at the products of the most recent solve, with EQUILIBRIUM
chemistry. After an invalid result it is the unburned state. An
independent copy.)doc";

static const char *__doc_Goddard_DetonationSolver_pre_detonation_state =
R"doc(Gas at the unburned state, with the chemistry it had at construction.
An independent copy. */)doc";

static const char *__doc_Goddard_DetonationSolver_reflected_detonation =
R"doc(Detonation reflected as a shock from the closed end of a tube, with
the products at equilibrium on both sides of the reflected shock.

`post_detonation_state()` is then state 5, the products brought to
rest.

Parameter ``drive_factor``:
    wave speed of the incident detonation over the CJ speed [-]

Parameter ``branch``:
    root of the Rayleigh line and the Hugoniot for the incident
    detonation

Returns:
    incident detonation and reflected shock; invalid if the incident
    detonation is invalid)doc";

static const char *__doc_Goddard_DetonationSolver_reset = R"doc(Restore the unburned state and the construction chemistry. */)doc";

static const char *__doc_Goddard_DetonationSolver_solve_chapman_jouguet =
R"doc(Solve the CJ detonation from the unburned state and store it with its
states. */)doc";

static const char *__doc_Goddard_DetonationSolver_solve_driven =
R"doc(Overdriven or under-driven detonation at `velocity`, whose drive
factor is `drive_factor`. Leaves the products in `m_gas` and records
the post and von Neumann states.)doc";

static const char *__doc_Goddard_DetonationSolver_store_unburned_states =
R"doc(Record the unburned state as both the post and the von Neumann state,
for an invalid result. */)doc";

static const char *__doc_Goddard_DetonationSolver_von_neumann_state =
R"doc(Gas behind the leading frozen shock of the most recent solve, with
FROZEN chemistry. After an invalid result it is the unburned state. An
independent copy.)doc";

static const char *__doc_Goddard_DilutedCombustor =
R"doc(Handles isobaric combustion with fuel, oxidizer, and recirculated flue
gas streams. The recirculation ratio is defined as r = m_flue /
(m_fuel + m_ox). The flue gas state is user-supplied (fixed), not
iteratively solved.)doc";

static const char *__doc_Goddard_DilutedCombustor_DilutedCombustor = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_DilutedCombustor_2 = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_generate_mass_fraction_matrix = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_generate_mole_fraction_matrix = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_m_flue_composition = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_m_fuel_composition = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_m_oxidizer_composition = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_solve = R"doc()doc";

static const char *__doc_Goddard_DilutedCombustor_solve_2 = R"doc()doc";

static const char *__doc_Goddard_EquilibriumDerivatives = R"doc()doc";

static const char *__doc_Goddard_EquilibriumDerivatives_dlogn_dlogP_T = R"doc((d log n / d log P)_T [-]. */)doc";

static const char *__doc_Goddard_EquilibriumDerivatives_dlogn_dlogT_P =
R"doc((d log n / d log T)_P [-], n being the moles of gas per kg of mixture.
*/)doc";

static const char *__doc_Goddard_EquilibriumDerivatives_dn_condensed_dlogP_T =
R"doc((d n_k / d log P)_T [kmol per kg of mixture], same ordering as
`dn_condensed_dlogT_P`.

At a pinned phase transition only the sum over the two coexisting
polymorphs is determined: it is reported on the lower-temperature
polymorph and the higher-temperature one is 0.)doc";

static const char *__doc_Goddard_EquilibriumDerivatives_dn_condensed_dlogT_P =
R"doc((d n_k / d log T)_P [kmol per kg of mixture], one entry per condensed
species that is currently present, in `Gas::condensed_species_names()`
order restricted to those species. Empty for a gas-only mixture.)doc";

static const char *__doc_Goddard_EquilibriumDerivatives_dpi_dlogP_T = R"doc((d pi_i / d log P)_T [-], one entry per element of the gas phase. */)doc";

static const char *__doc_Goddard_EquilibriumDerivatives_dpi_dlogT_P = R"doc((d pi_i / d log T)_P [-], one entry per element of the gas phase. */)doc";

static const char *__doc_Goddard_EquilibriumDerivatives_pinned_transition =
R"doc(True when the mixture sits exactly at a condensed phase transition
with both polymorphs present (see `Gas::at_phase_transition()`). The
temperature derivatives `dpi_dlogT_P`, `dlogn_dlogT_P` and
`dn_condensed_dlogT_P` are then NaN.)doc";

static const char *__doc_Goddard_EquilibriumOptions =
R"doc(Settings for the multiphase equilibrium solver used by `Gas`.

`rtol` and `max_steps` are passed to Cantera's Gibbs
(`MultiPhaseEquil`) solver for each constant-temperature-and-pressure
solve. The remaining entries control the outer one-dimensional
temperature root find that implements the HP and SP problems, which
Cantera cannot solve directly once the set of condensed phases depends
on temperature.)doc";

static const char *__doc_Goddard_EquilibriumOptions_T_default =
R"doc(Temperature [K] used to start the HP/SP root find when no better guess
is available. */)doc";

static const char *__doc_Goddard_EquilibriumOptions_T_max = R"doc(Highest temperature [K] the HP/SP root find will consider. */)doc";

static const char *__doc_Goddard_EquilibriumOptions_T_min = R"doc(Lowest temperature [K] the HP/SP root find will consider. */)doc";

static const char *__doc_Goddard_EquilibriumOptions_T_rel_tol =
R"doc(Relative convergence tolerance on temperature for the HP/SP outer root
find [-]. */)doc";

static const char *__doc_Goddard_EquilibriumOptions_max_bracket_steps =
R"doc(Maximum number of doubling steps used to bracket the HP/SP temperature
root. */)doc";

static const char *__doc_Goddard_EquilibriumOptions_max_steps =
R"doc(Maximum number of steps of the inner Gibbs solve. The Cantera default
(1000) is not enough near a dew point, where a barely-present
condensed phase converges slowly.)doc";

static const char *__doc_Goddard_EquilibriumOptions_rtol =
R"doc(Relative tolerance of the inner constant-T, constant-P Gibbs solve
[-]. */)doc";

static const char *__doc_Goddard_EquilibriumProperty =
R"doc(Property held constant together with pressure by the multiphase
equilibrium solver.

Cantera cannot solve the HP and SP problems once the set of condensed
phases depends on temperature, so `Gas` solves them as a one-
dimensional root find on temperature over constant-temperature,
constant-pressure solves. This enum selects which of the two.)doc";

static const char *__doc_Goddard_EquilibriumProperty_ENTHALPY =
R"doc(//!< Specific enthalpy [J/kg of mixture] is held fixed (the HP
problem).)doc";

static const char *__doc_Goddard_EquilibriumProperty_ENTROPY =
R"doc(//!< Specific entropy [J/(kg.K) of mixture] is held fixed (the SP
problem).)doc";

static const char *__doc_Goddard_ExitPlane =
R"doc(Flow state sampled across the nozzle exit plane, as parallel arrays
ordered from the axis outward. This is what
compute_thrust_coefficient() integrates over.)doc";

static const char *__doc_Goddard_ExitPlane_gamma_s = R"doc(< Local isentropic exponent.)doc";

static const char *__doc_Goddard_ExitPlane_mach = R"doc(< Local Mach number.)doc";

static const char *__doc_Goddard_ExitPlane_pressure =
R"doc(< Static pressure; Pa for frozen/equilibrium, normalized by the
stagnation pressure for perfect gas.)doc";

static const char *__doc_Goddard_ExitPlane_temperature =
R"doc(< Static temperature; K for frozen/equilibrium, normalized by the
stagnation temperature for perfect gas.)doc";

static const char *__doc_Goddard_ExitPlane_theta = R"doc(< Local flow angle, in radians.)doc";

static const char *__doc_Goddard_ExitPlane_velocity =
R"doc(< Dimensional V (m/s) for frozen/equilibrium; equal to the Mach number
for perfect gas.)doc";

static const char *__doc_Goddard_ExitPlane_y =
R"doc(< Radial (or transverse) station, in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_ExpansionProperties = R"doc()doc";

static const char *__doc_Goddard_ExpansionProperties_density = R"doc(Mixture density [kg/m^3], P / (gas_moles * R * T). */)doc";

static const char *__doc_Goddard_ExpansionProperties_dlV_dlP_T = R"doc((d log V / d log P)_T [-]. */)doc";

static const char *__doc_Goddard_ExpansionProperties_dlV_dlT_P = R"doc((d log V / d log T)_P [-]. */)doc";

static const char *__doc_Goddard_ExpansionProperties_frozen_gamma = R"doc(Frozen-composition ratio of specific heats cp/cv [-]. */)doc";

static const char *__doc_Goddard_ExpansionProperties_frozen_spec_heat_p = R"doc(Frozen-composition constant-pressure specific heat [J/(kg.K)]. */)doc";

static const char *__doc_Goddard_ExpansionProperties_gamma_s = R"doc(Isentropic exponent -(d log P / d log V)_s [-]. */)doc";

static const char *__doc_Goddard_ExpansionProperties_gas_moles = R"doc(Moles of gas per kg of mixture [kmol/kg] (CEA's 1/M). */)doc";

static const char *__doc_Goddard_ExpansionProperties_pinned_transition = R"doc(True when the mixture sits exactly at a condensed phase transition. */)doc";

static const char *__doc_Goddard_ExpansionProperties_spec_heat_p = R"doc(Equilibrium constant-pressure specific heat [J/(kg.K)]. */)doc";

static const char *__doc_Goddard_ExpansionProperties_spec_heat_v = R"doc(Equilibrium constant-volume specific heat [J/(kg.K)]. */)doc";

static const char *__doc_Goddard_ExpansionProperties_speed_of_sound =
R"doc(Equilibrium speed of sound [m/s], sqrt(gas_moles * R * T * gamma_s).
*/)doc";

static const char *__doc_Goddard_ExpansionProperties_total_moles =
R"doc(Moles of gas plus condensed species per kg of mixture [kmol/kg] (CEA's
1/MW). */)doc";

static const char *__doc_Goddard_ExpansionType = R"doc()doc";

static const char *__doc_Goddard_ExpansionType_PRESSURE_RATIO = R"doc()doc";

static const char *__doc_Goddard_ExpansionType_SUBSONIC_AREA_RATIO = R"doc()doc";

static const char *__doc_Goddard_ExpansionType_SUPERSONIC_AREA_RATIO = R"doc()doc";

static const char *__doc_Goddard_FiniteAreaChamber =
R"doc(Converged chamber of a finite-area combustor.

The chamber has constant area A_c, and propellant enters it with no
axial momentum. The flow from the combustion end onward is the
isentropic expansion from a hypothetical stagnation state "inf" with
the injector enthalpy and pressure P_inf < P_inj.)doc";

static const char *__doc_Goddard_FiniteAreaChamber_combustion_end =
R"doc(Combustion-end station, the subsonic point at area ratio
`contraction_ratio`. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_contraction_ratio = R"doc(Contraction ratio A_c/A_t [-]; solved for in mass-flux mode. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_injector = R"doc(Injector-face station; always in equilibrium. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_injector_pressure = R"doc(Injector-face pressure P_inj [Pa]. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_iterations = R"doc(Number of momentum-balance iterations [-]. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_mass_flux =
R"doc(Mass flux through the chamber, mdot/A_c [kg/(m^2 s)]; derived in
contraction-ratio mode. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_stagnation =
R"doc(Stagnation station "inf": equilibrium at the injector enthalpy and
P_inf. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_stagnation_pressure = R"doc(Stagnation pressure P_inf [Pa]. */)doc";

static const char *__doc_Goddard_FiniteAreaChamber_throat =
R"doc(Throat solved from the stagnation state; `throat.P_inlet` is P_inf
[Pa]. */)doc";

static const char *__doc_Goddard_Gas =
R"doc(Forward declaration. `equilibrium.hpp` is included by `gas.hpp` for
the result structs below, so it must not include `gas.hpp`. The `const
Gas&` overloads are defined in `equilibrium.cpp`.)doc";

static const char *__doc_Goddard_Gas_2 =
R"doc(Main Goddard class for querying thermodynamic information. `Gas` wraps
Cantera's `Solution` object and adds several convenience methods for
processes involving reactive flow.

`Gas` is itself a wrapper class that carries very little data, and is
cheap to copy.)doc";

static const char *__doc_Goddard_Gas_3 = R"doc()doc";

static const char *__doc_Goddard_GasChemistry = R"doc()doc";

static const char *__doc_Goddard_GasChemistry_EQUILIBRIUM = R"doc()doc";

static const char *__doc_Goddard_GasChemistry_FROZEN = R"doc()doc";

static const char *__doc_Goddard_GasChemistry_KINETIC = R"doc()doc";

static const char *__doc_Goddard_GasChemistry_PERFECT_GAS = R"doc()doc";

static const char *__doc_Goddard_Gas_Gas =
R"doc(@note The created `Gas` will reference the same `Solution` object
pointed to by the `shared_ptr`.)doc";

static const char *__doc_Goddard_Gas_Gas_2 =
R"doc(Create a `Gas` that carries a set of candidate condensed species.

@note Both the `Solution` and the condensed set are referenced, not
copied. Copies of the resulting `Gas` share both; `clone()` deep-
copies both.)doc";

static const char *__doc_Goddard_Gas_Gas_3 = R"doc(@note This constructor copies the `Solution` object.)doc";

static const char *__doc_Goddard_Gas_Gas_4 = R"doc(@note This constructor will reference the same `Solution` object.)doc";

static const char *__doc_Goddard_Gas_Gas_5 = R"doc(Create `Gas` from input data file.)doc";

static const char *__doc_Goddard_Gas_add_all_condensed_species =
R"doc(Add every species of `infile` whose elements are a subset of this
phase's elements. */)doc";

static const char *__doc_Goddard_Gas_add_condensed_species =
R"doc(@name Condensed species

A `Gas` optionally carries a set of candidate condensed (solid or
liquid) product species and the current amount of each. The set is
shared by copies of the `Gas` exactly as the underlying `Solution` is,
and deep-copied by `clone()`.

With a set attached the scalar property getters above become mixture
properties (per kg of gas plus condensed phases). `mole_fractions()`
and `mass_fractions()` stay gas-phase; `mixture_mass_fractions()`
covers the whole mixture. @{

Add named condensed species from a data file as candidates.

Parameter ``infile``:
    Cantera YAML data file containing the species.

Parameter ``names``:
    Species names, matching the data file exactly.

Throws:
    FmtError if a name is absent from the file, or if a species
    contains an element the gas phase does not have.)doc";

static const char *__doc_Goddard_Gas_area_per_mdot =
R"doc(Nozzle area per unit mass flow rate at the current state, RT/(P M u)
(RP-1311 eq. 6.12).

Parameter ``velocity``:
    flow velocity [m/s]

Returns:
    area per mass flow rate [m^2 s/kg])doc";

static const char *__doc_Goddard_Gas_at_phase_transition =
R"doc(True if the mixture sits exactly at a condensed phase transition (a
pinned state). */)doc";

static const char *__doc_Goddard_Gas_check_for_valid_cantera = R"doc()doc";

static const char *__doc_Goddard_Gas_chemistry = R"doc()doc";

static const char *__doc_Goddard_Gas_clear_phase_transition = R"doc(Clear the pinned phase transition, if any. */)doc";

static const char *__doc_Goddard_Gas_clone =
R"doc(Create a new Gas with a deep copy of the underlying Solution object
and of the condensed species set, if any. The chemistry mode and
stored reference stagnation enthalpy and entropy are preserved.

@note Copying a `Gas` normally is shallow: the copy references the
same `Solution` and the same condensed species set.)doc";

static const char *__doc_Goddard_Gas_condensed_cp_R =
R"doc(Reference-state molar heat capacity divided by R [-] of each *present*
condensed species. */)doc";

static const char *__doc_Goddard_Gas_condensed_enthalpy_RT =
R"doc(Reference-state molar enthalpy divided by R*T [-] of each *present*
condensed species. */)doc";

static const char *__doc_Goddard_Gas_condensed_molar_masses = R"doc(Molar mass [kg/kmol] of each *present* condensed species. */)doc";

static const char *__doc_Goddard_Gas_condensed_moles =
R"doc(Amount of each candidate condensed species [kmol per kg of mixture],
in candidate order. */)doc";

static const char *__doc_Goddard_Gas_condensed_species_names = R"doc(Names of all candidate condensed species, in candidate order. */)doc";

static const char *__doc_Goddard_Gas_condensed_stoich_coeffs =
R"doc(Stoichiometric coefficients of the *present* condensed species: a |C|
x l array whose rows are species and columns are elements, in the
element order of `element_names()`.)doc";

static const char *__doc_Goddard_Gas_copy_state = R"doc()doc";

static const char *__doc_Goddard_Gas_cp_mass =
R"doc(Constant-pressure specific heat [J/(kg.K) of mixture]: w_g cp_gas +
sum_k n_k cp_k(T). */)doc";

static const char *__doc_Goddard_Gas_create = R"doc()doc";

static const char *__doc_Goddard_Gas_create_from_elements = R"doc()doc";

static const char *__doc_Goddard_Gas_create_from_species =
R"doc(Create a `Gas` whose phase holds the named species of an input file.

The phase lists the species in the order of the file, whatever the
order of `species`, so species indices are reproducible.

Parameter ``infile``:
    Cantera YAML input file.

Parameter ``name``:
    Name of the new phase.

Parameter ``species``:
    Species names, matched exactly.

Parameter ``chemistry``:
    Chemistry model of the `Gas`.

Throws:
    std::invalid_argument if a species is not in the file.)doc";

static const char *__doc_Goddard_Gas_cstar =
R"doc(Characteristic velocity from the perfect-gas formula, using the
current temperature, molecular weight and `gamma_s()`. This is an
estimate: the c* of a rocket case is the one in `RocketPerformance`,
computed from the throat mass flux.

Returns:
    c* [m/s])doc";

static const char *__doc_Goddard_Gas_cv_mass =
R"doc(Constant-volume specific heat [J/(kg.K) of mixture], frozen
composition: `cp_mass() - n R`, with n the moles of gas per kg of
mixture.)doc";

static const char *__doc_Goddard_Gas_density =
R"doc(Current mass density [kg/m^3] of the mixture: rho_gas /
gas_mass_fraction(). */)doc";

static const char *__doc_Goddard_Gas_element_moles =
R"doc(Amount of each element [kmol per kg of mixture], gas plus condensed,
in the element order of `element_names()`.)doc";

static const char *__doc_Goddard_Gas_element_names = R"doc(Element names of the gas phase, in the phase's element order. */)doc";

static const char *__doc_Goddard_Gas_enthalpy_mass = R"doc(Specific enthalpy [J/kg of mixture]: w_g h_gas + sum_k n_k H_k(T). */)doc";

static const char *__doc_Goddard_Gas_entropy_mass =
R"doc(Specific entropy [J/(kg.K) of mixture]: w_g s_gas + sum_k n_k S_k(T).
*/)doc";

static const char *__doc_Goddard_Gas_equilibrate =
R"doc(Equilibrate the gas mixture, holding two thermodynamic properties
constant.

Parameter ``XY``:
    two-letter property pair to hold constant, e.g. "HP", "TP", "SP".

Parameter ``solver``:
    Cantera equilibrium solver name (default "gibbs").

With candidate condensed species attached the call is routed to
`equilibrate_TP()`, `equilibrate_HP()` or `equilibrate_SP()` with the
current value of the held properties.

Throws:
    std::invalid_argument if `solver` is "vcs", which gives wrong
    answers with condensed phases, or if `XY` is not one of "TP",
    "HP", "SP", "UV".

Throws:
    NotImplementedError for "UV" with condensed candidates attached.)doc";

static const char *__doc_Goddard_Gas_equilibrate_HP =
R"doc(Equilibrate at fixed specific enthalpy [J/kg of mixture] and pressure
[Pa].

With condensed candidates this is a bracketed root find on temperature
over `equilibrate_TP()` solves, warm-started from the current
temperature. When the requested enthalpy falls inside the latent heat
of a phase transition the temperature is pinned at the transition and
the two polymorphs are split to match; see `at_phase_transition()`.

Throws:
    ConvergenceError if the root cannot be bracketed in
    `[EquilibriumOptions::T_min, EquilibriumOptions::T_max]` or the
    iteration does not converge.)doc";

static const char *__doc_Goddard_Gas_equilibrate_SP =
R"doc(Equilibrate at fixed specific entropy [J/(kg.K) of mixture] and
pressure [Pa].

See also:
    equilibrate_HP() for the condensed-phase algorithm and the failure
    modes.)doc";

static const char *__doc_Goddard_Gas_equilibrate_TP =
R"doc(Equilibrate at fixed temperature [K] and pressure [Pa].

With candidate condensed species attached this is one multiphase Gibbs
minimization over the gas phase and the candidates whose thermodynamic
data covers `T`. Amounts of candidates that have left their data range
are first moved to the in-range polymorph of their group.

Throws:
    FmtError if a condensed species is present and no polymorph of its
    group has data at `T`.

Throws:
    ConvergenceError if the solve fails after all restart strategies.)doc";

static const char *__doc_Goddard_Gas_equilibrium_options =
R"doc(Settings for the equilibrium solvers used by `equilibrate_TP/HP/SP`.
*/)doc";

static const char *__doc_Goddard_Gas_equivalence_ratio = R"doc()doc";

static const char *__doc_Goddard_Gas_equivalence_ratio_2 = R"doc()doc";

static const char *__doc_Goddard_Gas_equivalence_ratio_3 = R"doc()doc";

static const char *__doc_Goddard_Gas_expansion_properties =
R"doc(Compute expansion-related thermodynamic derivatives at the current
state: gamma_s, dlnV/dlnT|_P, dlnV/dlnP|_T, and cp.)doc";

static const char *__doc_Goddard_Gas_fuel_fraction = R"doc()doc";

static const char *__doc_Goddard_Gas_fuel_fraction_2 = R"doc()doc";

static const char *__doc_Goddard_Gas_gamma = R"doc(Convenience for the standard cp/cv adiabatic index.)doc";

static const char *__doc_Goddard_Gas_gamma_s =
R"doc(Effective ratio of specific heats used for compressible-flow
calculations.

The result depends on the active GasChemistry mode: PERFECT_GAS uses
cp/cv, FROZEN holds composition fixed, EQUILIBRIUM accounts for
shifting equilibrium.)doc";

static const char *__doc_Goddard_Gas_gas_mass_fraction =
R"doc(Gas mass fraction w_g = 1 - sum_k n_k M_k [-]. 1 for a gas-only
mixture. */)doc";

static const char *__doc_Goddard_Gas_get_reference_entropy = R"doc(Get the stored reference specific entropy [J/(kg.K)]. */)doc";

static const char *__doc_Goddard_Gas_get_stagnation_enthalpy =
R"doc(Get the stored reference stagnation enthalpy [J/kg].

See set_stagnation_enthalpy() for the semantics of the reference
value.)doc";

static const char *__doc_Goddard_Gas_has_cantera_sln = R"doc(True if underlying Cantera Solution object has been initialized.)doc";

static const char *__doc_Goddard_Gas_has_condensed_candidates = R"doc(True if any candidate condensed species is attached. */)doc";

static const char *__doc_Goddard_Gas_has_condensed_phases =
R"doc(True if any candidate condensed species is present in a nonzero
amount. */)doc";

static const char *__doc_Goddard_Gas_isenthalpic_velocity =
R"doc(Velocity that produces a given stagnation enthalpy from the current
static state.

Parameter ``H_stagnation``:
    target stagnation enthalpy [J/kg]

Returns:
    velocity [m/s])doc";

static const char *__doc_Goddard_Gas_isenthalpic_velocity_2 =
R"doc(Velocity computed from the stored reference stagnation enthalpy (see
set_stagnation_enthalpy()).

Returns:
    velocity [m/s])doc";

static const char *__doc_Goddard_Gas_isp =
R"doc(Isenthalpic velocity from the stored stagnation enthalpy, the specific
impulse of an exit at the current state.

Returns:
    specific impulse [m/s])doc";

static const char *__doc_Goddard_Gas_ivac =
R"doc(Vacuum specific impulse of an exit at the current state, u + RT/(M u).

Returns:
    vacuum specific impulse [m/s])doc";

static const char *__doc_Goddard_Gas_kinetics = R"doc()doc";

static const char *__doc_Goddard_Gas_last_equilibrium_solve_count =
R"doc(Number of constant-temperature, constant-pressure multiphase solves
consumed by the most recent `equilibrate_TP/HP/SP` call [-]. Zero
before the first call; a diagnostic only.)doc";

static const char *__doc_Goddard_Gas_m_H_stagnation = R"doc()doc";

static const char *__doc_Goddard_Gas_m_S0 = R"doc()doc";

static const char *__doc_Goddard_Gas_m_condensed =
R"doc(Candidate condensed species and their amounts. Null means the `Gas` is
gas-only. */)doc";

static const char *__doc_Goddard_Gas_m_equilibrium_solve_count =
R"doc(Constant-T, constant-P solves used by the most recent `equilibrate_*`
call. */)doc";

static const char *__doc_Goddard_Gas_m_sol = R"doc()doc";

static const char *__doc_Goddard_Gas_mach =
R"doc(Mach number at the current state for a given velocity.

Parameter ``velocity``:
    flow velocity [m/s])doc";

static const char *__doc_Goddard_Gas_mass_fractions =
R"doc(Mass fractions of the *gas* phase [-], summing to 1 over the gas
alone. See `mixture_mass_fractions()` for fractions of the whole
mixture.)doc";

static const char *__doc_Goddard_Gas_mixture_mass_fractions =
R"doc(Mass fractions of the whole mixture [-]: the gas-phase mass fractions
scaled by `gas_mass_fraction()`, followed by n_k M_k for every
candidate condensed species in candidate order. Sums to 1.)doc";

static const char *__doc_Goddard_Gas_mixture_molecular_weight =
R"doc(Mixture molecular weight [kg/kmol], CEA's "MW": 1 / (n + sum_k n_k)
with n the moles of gas per kg of mixture. `molecular_weight()`
remains CEA's "M" = 1/n, the gas-phase value.)doc";

static const char *__doc_Goddard_Gas_mole_fractions =
R"doc(Mole fractions of the *gas* phase [-]. Condensed amounts are reported
separately. */)doc";

static const char *__doc_Goddard_Gas_molecular_weight =
R"doc(Molecular weight [kg/kmol], CEA's "M" = 1/n with n the moles of *gas*
per kg of *mixture*. Equal to the gas-phase mean molecular weight
divided by `gas_mass_fraction()`, and to the gas-phase value when no
condensed phase is present. See `mixture_molecular_weight()` for CEA's
"MW", which counts the condensed moles as well.)doc";

static const char *__doc_Goddard_Gas_name = R"doc()doc";

static const char *__doc_Goddard_Gas_num_species = R"doc()doc";

static const char *__doc_Goddard_Gas_pinned_polymorphs =
R"doc(Indices of the two coexisting polymorphs at a pinned phase transition,
lower-temperature polymorph first. The indices refer to the condensed
species that are currently *present* (nonzero moles), in candidate
order. `{-1, -1}` when no group is pinned.)doc";

static const char *__doc_Goddard_Gas_pressure = R"doc(Current pressure [Pa]. */)doc";

static const char *__doc_Goddard_Gas_report = R"doc(Generate report string of underlying Cantera `ThermoPhase` object.)doc";

static const char *__doc_Goddard_Gas_restore_state =
R"doc(Restore a thermodynamic state previously captured by save_state().

Accepts either the extended length described in `save_state()` or the
bare Cantera state length, in which case all condensed amounts are set
to zero.

Throws:
    std::invalid_argument if the vector has neither length.)doc";

static const char *__doc_Goddard_Gas_save_state =
R"doc(Snapshot the current thermodynamic state into an opaque vector that
can later be passed to restore_state().

The vector is the Cantera state of the gas phase followed by one entry
per candidate condensed species (its amount in kmol per kg of
mixture). With no candidates attached it is the Cantera state alone,
exactly as before.

@note The pinned-transition flag is not stored: it is re-derived on
restore from the fact that two polymorphs of one group are present
simultaneously.)doc";

static const char *__doc_Goddard_Gas_set_OF_ratio = R"doc()doc";

static const char *__doc_Goddard_Gas_set_OF_ratio_2 = R"doc()doc";

static const char *__doc_Goddard_Gas_set_condensed_moles =
R"doc(Set the amount of each candidate condensed species [kmol per kg of
mixture].

Throws:
    std::invalid_argument if `moles` does not have one entry per
    candidate.)doc";

static const char *__doc_Goddard_Gas_set_current_state_as_reference = R"doc()doc";

static const char *__doc_Goddard_Gas_set_element_moles =
R"doc(Set the composition from elemental amounts, without equilibrating.

The state becomes a "basis" composition holding the requested
elements: each element is assigned to a single single-element species
(the homonuclear diatomic if the phase has one, otherwise the
monatomic, otherwise any species made of that element alone). This is
the starting point a subsequent `equilibrate_*` call refines.

Parameter ``element_moles``:
    Amount of each element [kmol per kg of mixture], in the element
    order of `element_names()`.

Parameter ``T``:
    Temperature [K] of the resulting state.

Parameter ``P``:
    Pressure [Pa] of the resulting state.

Throws:
    FmtError if the phase contains no single-element species for one
    of the requested elements. @note All condensed amounts are set to
    zero and any pinned transition is cleared, so the requested
    element amounts are carried entirely by the gas phase.)doc";

static const char *__doc_Goddard_Gas_set_equivalence_ratio = R"doc()doc";

static const char *__doc_Goddard_Gas_set_equivalence_ratio_2 = R"doc()doc";

static const char *__doc_Goddard_Gas_set_fuel_fraction =
R"doc(Set the composition from a fuel fraction, an equivalence ratio or an
O/F ratio, keeping temperature and pressure.

These setters also make the new state the reference: the stored
stagnation enthalpy and reference entropy are reset to those of the
new mixture.)doc";

static const char *__doc_Goddard_Gas_set_fuel_fraction_2 = R"doc()doc";

static const char *__doc_Goddard_Gas_set_name = R"doc()doc";

static const char *__doc_Goddard_Gas_set_phase_transition =
R"doc(Mark the polymorph group shared by two *present* condensed species as
pinned at its phase transition, the state `pinned_polymorphs()`
reports.

The solver sets this itself; the setter exists so a state built by
hand (a test fixture, or a caller that restores a split it computed
elsewhere) can declare the pin.

Parameter ``low``:
    Index of the lower-temperature polymorph among the condensed
    species that are currently present (nonzero moles), in candidate
    order.

Parameter ``high``:
    Index of the higher-temperature polymorph, same convention.

Throws:
    std::invalid_argument if either index is out of range or the two
    species do not belong to the same polymorph group.)doc";

static const char *__doc_Goddard_Gas_set_phase_transition_2 =
R"doc(Mark the polymorph group of two named candidate condensed species as
pinned at its phase transition.

Parameter ``low``:
    Name of the lower-temperature polymorph, as it appears in the data
    file.

Parameter ``high``:
    Name of the higher-temperature polymorph.

Throws:
    std::invalid_argument if either name is not a candidate of this
    `Gas` or the two are not polymorphs of one another.)doc";

static const char *__doc_Goddard_Gas_set_reference_entropy = R"doc(Set the stored reference specific entropy [J/(kg.K)]. */)doc";

static const char *__doc_Goddard_Gas_set_stagnation_enthalpy =
R"doc(Set the stored reference stagnation enthalpy [J/kg].

The reference value is used by isenthalpic_velocity() and other
stagnation helpers and must be set explicitly. It does not change when
the thermodynamic state of `Gas` changes.)doc";

static const char *__doc_Goddard_Gas_set_state_HP = R"doc(Set state from specific enthalpy [J/kg] and pressure [Pa]. */)doc";

static const char *__doc_Goddard_Gas_set_state_SP = R"doc(Set state from specific entropy [J/(kg.K)] and pressure [Pa]. */)doc";

static const char *__doc_Goddard_Gas_set_state_TD = R"doc(Set state from temperature [K] and density [kg/m^3]. */)doc";

static const char *__doc_Goddard_Gas_set_state_TP = R"doc(Set state from temperature [K] and pressure [Pa]. */)doc";

static const char *__doc_Goddard_Gas_set_state_TPX =
R"doc(Set state from temperature [K], pressure [Pa], and mole-fraction
composition string. */)doc";

static const char *__doc_Goddard_Gas_set_state_TPX_2 =
R"doc(Set state from temperature [K], pressure [Pa], and mole-fraction map.
*/)doc";

static const char *__doc_Goddard_Gas_set_state_TPX_3 = R"doc()doc";

static const char *__doc_Goddard_Gas_set_state_TPY =
R"doc(Set state from temperature [K], pressure [Pa], and mass-fraction
composition string. */)doc";

static const char *__doc_Goddard_Gas_set_state_TPY_2 =
R"doc(Set state from temperature [K], pressure [Pa], and mass-fraction map.
*/)doc";

static const char *__doc_Goddard_Gas_set_state_TPY_3 = R"doc()doc";

static const char *__doc_Goddard_Gas_set_state_UV =
R"doc(Set state from specific internal energy [J/kg] and specific volume
[m^3/kg]. */)doc";

static const char *__doc_Goddard_Gas_snapshot =
R"doc(Return a ThermodynamicState struct containing all current properties.

All properties are mixture properties (per kg of gas plus condensed
phases) and the composition is the mixture mass fractions of
`mixture_mass_fractions()`, with the condensed species appended after
the gas species.)doc";

static const char *__doc_Goddard_Gas_solution =
R"doc(Underlying Cantera Solution handle. Pass to solver constructors that
accept one. */)doc";

static const char *__doc_Goddard_Gas_solve_frozen_XP =
R"doc(Newton iteration on temperature for the frozen
`set_state_HP`/`set_state_SP` problems with condensed phases present:
composition and condensed amounts are held fixed.

Throws:
    FmtError if the solution lies outside the temperature range of a
    present condensed species, as CEA reports for a frozen expansion.)doc";

static const char *__doc_Goddard_Gas_solve_multiphase_TP =
R"doc(One constant-temperature, constant-pressure multiphase Gibbs
minimization.

Moves the amounts of every group onto the candidate that is offered to
the solver, builds a `Cantera::MultiPhase` from the gas phase and the
offered candidates on a one kilogram basis, equilibrates it with the
"gibbs" solver and reads the phase amounts back. Mass and element
balances are verified afterwards; a `Cantera::CanteraError` or a
failed balance restarts the solve from a different initial guess.

Parameter ``T``:
    Temperature [K].

Parameter ``P``:
    Pressure [Pa].

Parameter ``offered_override``:
    Candidate indices to offer, one per polymorph group at most. Empty
    means the automatic choice of `CondensedPhaseSet::offered_at()`;
    the outer HP/SP loop uses it to offer the upper polymorph of a
    group at its transition temperature.)doc";

static const char *__doc_Goddard_Gas_solve_multiphase_XP =
R"doc(Outer temperature root find implementing the multiphase HP and SP
problems.

Parameter ``property``:
    Which property `target` refers to.

Parameter ``target``:
    Specific enthalpy [J/kg] or entropy [J/(kg.K)] of the mixture.

Parameter ``P``:
    Pressure [Pa].

Returns:
    The number of constant-T, constant-P solves used.)doc";

static const char *__doc_Goddard_Gas_species_names = R"doc()doc";

static const char *__doc_Goddard_Gas_speed_of_sound =
R"doc(Speed of sound [m/s] at the current state, consistent with gamma_s().
*/)doc";

static const char *__doc_Goddard_Gas_stagnation_enthalpy =
R"doc(Stagnation (total) enthalpy h0 = h + v^2/2.

Parameter ``velocity``:
    flow velocity [m/s]

Returns:
    stagnation enthalpy [J/kg])doc";

static const char *__doc_Goddard_Gas_stagnation_pressure =
R"doc(Stagnation pressure obtained by isentropic deceleration from the
current state.

The deceleration follows the active GasChemistry mode: EQUILIBRIUM re-
equilibrates the mixture at constant entropy, and every other mode
holds the composition fixed. The current state is not modified.

Parameter ``velocity``:
    flow velocity [m/s]

Returns:
    stagnation pressure [Pa])doc";

static const char *__doc_Goddard_Gas_state_size = R"doc(Size of the state vector.)doc";

static const char *__doc_Goddard_Gas_stoich_OF_ratio = R"doc()doc";

static const char *__doc_Goddard_Gas_stoich_OF_ratio_2 = R"doc()doc";

static const char *__doc_Goddard_Gas_temperature = R"doc(Current temperature [K]. */)doc";

static const char *__doc_Goddard_Gas_thermo = R"doc()doc";

static const char *__doc_Goddard_Gas_transport = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle = R"doc()doc";

static const char *__doc_Goddard_KineticNozzleResults = R"doc()doc";

static const char *__doc_Goddard_KineticNozzleResults_reached_exit =
R"doc(True if the integration reached the exit of the profile. False if it
stopped after `max_steps` steps, in which case `stations` ends
upstream of the exit.)doc";

static const char *__doc_Goddard_KineticNozzleResults_stations = R"doc()doc";

static const char *__doc_Goddard_KineticNozzleResults_throat = R"doc()doc";

static const char *__doc_Goddard_KineticNozzleStation = R"doc(One station of a kinetic nozzle integration. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_Da_min = R"doc(Smallest non-zero Damkohler number [-]. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_area_ratio = R"doc(Area ratio A/A_t [-]. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_damkohler = R"doc(Damkohler number of each species [-]; 0 for trace species. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_mach = R"doc(Mach number [-], with the frozen speed of sound. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_min_Da_species =
R"doc(Name of the species with the smallest non-zero Damkohler number; empty
if there is none. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_state = R"doc(Raw state vector, for `Gas::restore_state`. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_thermo = R"doc(Mixture state at the station, with frozen derivatives. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_velocity = R"doc(Flow velocity [m/s]. */)doc";

static const char *__doc_Goddard_KineticNozzleStation_x = R"doc(Axial position, in the length unit of the profile. */)doc";

static const char *__doc_Goddard_KineticNozzle_KineticNozzle = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_KineticNozzle_2 = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_KineticNozzle_3 = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_KineticNozzle_4 = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_m_gas = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_m_inlet_state = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_m_throat_solver = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_mdot = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_opts = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_profile = R"doc()doc";

static const char *__doc_Goddard_KineticNozzle_solve = R"doc()doc";

static const char *__doc_Goddard_MixtureRatioType = R"doc()doc";

static const char *__doc_Goddard_MixtureRatioType_FUEL_FRAC = R"doc()doc";

static const char *__doc_Goddard_MixtureRatioType_OF_RATIO = R"doc()doc";

static const char *__doc_Goddard_MixtureRatioType_PHI_RATIO = R"doc()doc";

static const char *__doc_Goddard_MocCrossings =
R"doc(Result of scanning a characteristic net for intersections between
characteristics of the same family.

Two characteristics of one family meeting is the discrete signature of
an oblique shock: the wave family is coalescing rather than diverging.
An isentropic MoC solution is only valid where this does not happen,
so a converged solve reporting a nonzero count has produced a field
the method cannot represent -- regardless of whether every unit
process succeeded.)doc";

static const char *__doc_Goddard_MocCrossings_count = R"doc(< Number of same-family characteristic crossings found.)doc";

static const char *__doc_Goddard_MocCrossings_first_family = R"doc(< Family that crossed there.)doc";

static const char *__doc_Goddard_MocCrossings_first_x =
R"doc(< Position of the most upstream crossing, in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_MocCrossings_first_y = R"doc(< See `first_x`.)doc";

static const char *__doc_Goddard_MocErrorCode =
R"doc(Classification of the numerical/physical failure modes that can occur
while resolving a single unit process (interior/wall/axis point solve)
or marching the method-of-characteristics net as a whole.)doc";

static const char *__doc_Goddard_MocErrorCode_INITIALIZATION_FAILED =
R"doc(< Construction of the initial data line (transonic start line) failed
to converge.)doc";

static const char *__doc_Goddard_MocErrorCode_MAX_ITERATIONS_REACHED =
R"doc(< The characteristic kernel reached its iteration safety cap before
all chains terminated.)doc";

static const char *__doc_Goddard_MocErrorCode_NEGATIVE_NU =
R"doc(< Prandtl-Meyer angle (or generalized PM function) nu is negative
beyond tolerance.)doc";

static const char *__doc_Goddard_MocErrorCode_NEGATIVE_THETA = R"doc(< Flow angle theta is negative beyond tolerance.)doc";

static const char *__doc_Goddard_MocErrorCode_NONE = R"doc(< No error; the point/solve is valid.)doc";

static const char *__doc_Goddard_MocErrorCode_NONFINITE_VALUE =
R"doc(< One of the point's numeric fields (x, y, theta, nu, mach, mu) is NaN
or infinite.)doc";

static const char *__doc_Goddard_MocErrorCode_NON_DOWNSTREAM_POINT =
R"doc(< The computed intersection lies at or behind (upstream of) one of its
parent points.)doc";

static const char *__doc_Goddard_MocErrorCode_PM_INVERSION_FAILED =
R"doc(< The Prandtl-Meyer inversion (nu -> Mach), or an equivalent Mach
rootfind, failed to converge.)doc";

static const char *__doc_Goddard_MocErrorCode_SUBSONIC_MACH =
R"doc(< Mach number is below 1.0; the supersonic compatibility relations no
longer apply.)doc";

static const char *__doc_Goddard_MocErrorCode_TABLE_RANGE_EXCEEDED =
R"doc(< A PrandtlMeyerTable lookup (by nu, Mach, or velocity) fell outside
the tabulated range.)doc";

static const char *__doc_Goddard_MocErrorCode_WALL_QUERY_OUT_OF_BOUNDS =
R"doc(< A wall-profile query (e.g. theta_at) fell outside the profile's
domain.)doc";

static const char *__doc_Goddard_MocFailure =
R"doc(Describes a numerical or physical failure encountered while solving a
MocNozzle problem. When MocResult::converged is false, this identifies
what went wrong and where, instead of the failure being silently
swallowed or thrown as an exception.)doc";

static const char *__doc_Goddard_MocFailure_code = R"doc(< Failure classification; NONE if no failure occurred.)doc";

static const char *__doc_Goddard_MocFailure_kernel_pass =
R"doc(< Kernel marching pass on which the failure occurred; -1 for pre-
kernel (initialization) failures.)doc";

static const char *__doc_Goddard_MocFailure_message = R"doc(< Human-readable description of the failure.)doc";

static const char *__doc_Goddard_MocFailure_x =
R"doc(< x-coordinate of the failing point (or a parent's, if the point
itself could not be computed), in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_MocFailure_y =
R"doc(< y-coordinate of the failing point (or a parent's, if the point
itself could not be computed), in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_MocFlowKind =
R"doc(Dimensionality of the flow field the characteristics are marched
through. */)doc";

static const char *__doc_Goddard_MocFlowKind_AXISYMMETRIC = R"doc(< Flow about a centerline; the source term is singular on the axis.)doc";

static const char *__doc_Goddard_MocFlowKind_PLANAR =
R"doc(< Two-dimensional flow; the compatibility relations have no source
term.)doc";

static const char *__doc_Goddard_MocFrontShear =
R"doc(How differently the two ends of the marching front advanced over a
solve.

Both ends of an axisymmetric front accelerate; the failure is that
they do so at very different rates, so the front shears. Comparing the
mean axial step over the first few passes with the mean over the last
few reduces that history to one number per end.

Measured on conical AR=8, N=15: the axis end's step grew x20.6 while
the wall end's *shrank* to x0.84 -- a shear ratio near 24. A
converging planar solve on the same geometry gives x7.9 and x44.6, a
ratio near 5.6. The sign differs because planar's wall end outruns its
axis end; what matters is the magnitude of the disparity.)doc";

static const char *__doc_Goddard_MocFrontShear_axis_growth = R"doc(< Late mean axial step at the axis end / early mean.)doc";

static const char *__doc_Goddard_MocFrontShear_shear_ratio =
R"doc(< max(axis_growth, wall_growth) / min(...); 1 means the ends kept
pace.)doc";

static const char *__doc_Goddard_MocFrontShear_valid = R"doc(< False when there were too few passes to form both windows.)doc";

static const char *__doc_Goddard_MocFrontShear_wall_growth = R"doc(< Same at the wall end.)doc";

static const char *__doc_Goddard_MocInitDiagnostics =
R"doc(Properties of the initial data line, measured once before the march
begins.

The initialization decides whether an axisymmetric analysis converges
at all: with a centered-fan start line the solver grid-converges, and
with the Kliegel-Levine start line it gets *worse* under refinement.
Anti-convergence means an error that is fixed in absolute terms while
the cell size shrinks, so the quantities here are reported both raw
and normalized by the characteristic spacing -- a raw value that stays
put while its normalized twin grows with N is the signature.

All of these are cheap and are filled in for every mode and both
initializers, so a healthy solve's numbers are available as the
reference for a sick one's.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_axis_arrival_grading =
R"doc(Largest / smallest gap between the points where the data line's C-
characteristics reach the axis, on a straight-ray estimate. The start
line is uniform in y but its characteristics need not be uniform in
where they land; a large value means the axis is fed by a burst of
arrivals and then starved.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_cot_mu_ratio =
R"doc(cot(mu_axis) / cot(mu_wall). Marching step length scales with cot(mu),
so this is how much faster one end of the front is licensed to advance
than the other at pass 0.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_kplus_wall_end =
R"doc(K+ = theta - nu of the topmost interior point (the data line point
just below the wall end), after the wall-consistency correction. A
large negative value here over-expands the first wall solve.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_mach_axis = R"doc(< Mach at the data line's axis end.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_mach_ratio = R"doc(< mach_wall / mach_axis; 1 for a constant-Mach line.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_mach_wall = R"doc(< Mach at its wall end.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_mass_flow_error =
R"doc(Relative error of the mass flow integrated across the data line
against the 1-D critical mass flow through the throat. This is the one
measurement that checks the start line against physics rather than
against its own construction.

NaN for frozen and equilibrium chemistry, where the density needed for
the integrand is not recoverable from the stored point state.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_min_spacelike_margin =
R"doc(Smallest spacelike margin over the data line's own segments: how far a
segment sits from being parallel to a characteristic through one of
its endpoints, as a fraction of its height. A value <= 0 means the
line is crossed by its own characteristics, i.e. it is not a valid
Cauchy surface for the march.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_mu_axis = R"doc(< Mach angle at the axis end, radians.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_mu_wall = R"doc(< Mach angle at the wall end, radians.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_points = R"doc(< Points on the initial data line.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_shift_over_transonic_length =
R"doc(MocOptions::initial_line_axial_shift divided by sqrt(R), the transonic
length scale. The shift is specified as an absolute offset in throat
radii, but the axial extent of the transonic region scales as
sqrt(r_throat * R_curvature), so a fixed shift means different things
at different throat curvatures. Zero when the shift is unused.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_start_line_used =
R"doc(Which start line the solve actually used -- meaningful when
MocOptions::start_line is AUTO and the Kliegel-Levine series was
rejected in favor of the fan.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_wall_bc_residual =
R"doc(How far the raw (uncorrected) Kliegel-Levine series missed its own
wall boundary condition at the data line's wall end: 1 -
theta_series/theta_wall, measured before initialize_kliegel_levine's
wall-consistency correction is applied. Zero for the centered fan (no
series to miss).)doc";

static const char *__doc_Goddard_MocInitDiagnostics_wall_gap =
R"doc(Distance from the data line's wall end to the prescribed wall contour,
in the unit of NozzleGeometry::throat_radius. The line's top point is
registered as *the* wall point, so a nonzero value means the wall
march begins from a point that is not on the wall. Zero in design
modes, which have no prescribed contour.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_wall_gap_over_spacing =
R"doc(wall_gap divided by the characteristic spacing. Growth with N is the
anti-convergence signature.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_wall_station_to_tangency =
R"doc(Signed distance, in characteristic spacings, from the data line's wall
end to the nearest curvature discontinuity on the contour (negative =
upstream of it). Near zero means the wall march is seeded exactly on a
corner, where the contour's angle query is a one-sided difference
across a slope jump. Large sentinel value for a smooth contour.)doc";

static const char *__doc_Goddard_MocInitDiagnostics_wall_theta_mismatch =
R"doc(|theta at the line's wall end - the contour's own angle there|, in
radians.)doc";

static const char *__doc_Goddard_MocLogLevel =
R"doc(How much diagnostic output a solve collects into MocResult::messages.
*/)doc";

static const char *__doc_Goddard_MocLogLevel_DEBUG = R"doc(< also collect (and echo live to stderr) a verbose kernel trace --)doc";

static const char *__doc_Goddard_MocLogLevel_NORMAL = R"doc(< Default: only MocLog::warning()/MocLog::info() messages collected.)doc";

static const char *__doc_Goddard_MocMode =
R"doc(What the solver is being asked to do: derive a wall contour, or
resolve the flow through a contour that is already known.)doc";

static const char *__doc_Goddard_MocMode_ANALYSIS = R"doc(< Wall contour is input; the flow field is solved through it.)doc";

static const char *__doc_Goddard_MocMode_DESIGN_CENTERLINE =
R"doc(< Optimal nozzle from prescribed centerline values. @warning Currently
unimplemented; throws NotImplementedError.)doc";

static const char *__doc_Goddard_MocMode_DESIGN_MIN_LENGTH = R"doc(< Minimum length nozzle with uniform exit flow.)doc";

static const char *__doc_Goddard_MocMode_DESIGN_RAO = R"doc(< Rao-type length-optimized thrust nozzle.)doc";

static const char *__doc_Goddard_MocNozzle =
R"doc(Main class for performing 2D nozzle supersonic flow simulations using
the method of characteristics (MoC). Simulation specification is done
with `MocOptions`.

There are two modes of usage: - Design mode, where a specific mach
number or other end goal is specified and the nozzle profile is
determined by the solver. - Analysis mode, where an arbitrary nozzle
profile is provided and the solver determines the flow field and
performance metrics.

solve() is the orchestrator: it validates `options`, resolves the
throat and thermo model (MocThermo), builds and measures the start
line (moc_initialization.hpp), then dispatches to one of the two
kernels by MocMode -- DirectMarch for DESIGN_MIN_LENGTH, InverseMarch
for ANALYSIS and DESIGN_RAO (moc_direct_march.hpp,
moc_inverse_march.hpp) -- and assembles the result.)doc";

static const char *__doc_Goddard_MocNozzle_MocNozzle = R"doc()doc";

static const char *__doc_Goddard_MocNozzle_MocNozzle_2 =
R"doc(Solver for real-gas chemistry. The solver works on its own copy of
`gas`. */)doc";

static const char *__doc_Goddard_MocNozzle_is_solved = R"doc(True once solve() has run to completion on this instance. */)doc";

static const char *__doc_Goddard_MocNozzle_m_gas = R"doc()doc";

static const char *__doc_Goddard_MocNozzle_m_inverse_front_override =
R"doc(Test-only hook (kept protected; set only by a test subclass): when
present, solve() seeds the inverse kernel's first front directly from
these points instead of running build_start_line, so a hand-built
front with a known exact solution (e.g. uniform flow, or a
manufactured source-flow field) can be marched without requiring a
throat/KL/fan construction consistent with it. See
InverseMarchUniformFlow.StaysUniform and
InverseMarchSourceFlow.SecondOrderConvergence
(test/test_moc_inverse_march.cpp). The points' coordinates are in the
unit of NozzleGeometry::throat_radius, like
MocOptions::nozzle_profile.)doc";

static const char *__doc_Goddard_MocNozzle_m_is_solved = R"doc()doc";

static const char *__doc_Goddard_MocNozzle_options =
R"doc(Solver configuration. Public and re-read by solve(), so it can be
adjusted between solves.)doc";

static const char *__doc_Goddard_MocNozzle_solve =
R"doc(March the characteristic net and return the solved flow field.

Never throws for numerical failures: a solve that breaks down returns
a MocResult with `converged == false` and `failure` identifying what
went wrong and where. It does throw std::invalid_argument if `options`
is self-inconsistent (checked here as well as in the constructor,
since `options` is public and callers do adjust it in between), and
NotImplementedError for MocMode::DESIGN_CENTERLINE.)doc";

static const char *__doc_Goddard_MocOptions = R"doc(Options for method of characteristics simulations.)doc";

static const char *__doc_Goddard_MocOptions_chemistry = R"doc()doc";

static const char *__doc_Goddard_MocOptions_exit_mach =
R"doc(Target exit Mach number, for the design modes. Not implemented: the
design target is set through theta_max, and `validate_moc_options`
raises NotImplementedError if this is non-zero.)doc";

static const char *__doc_Goddard_MocOptions_flow_type = R"doc()doc";

static const char *__doc_Goddard_MocOptions_front_tilt_decay =
R"doc(Inverse march only: fraction of the marching front's shape retained
per pass, in [0, 1]. Each new front is the previous one translated
downstream, with every point's axial offset from the wall point
multiplied by this factor, so a Kliegel-Levine start line (axis end
downstream of its wall end) relaxes toward a vertical plane as the
march proceeds. The relaxation is additionally capped so it moves no
point by more than half a step, which keeps the step bounds valid; on
a front with large offsets that cap, not this factor, is what limits
the relaxation.)doc";

static const char *__doc_Goddard_MocOptions_gamma =
R"doc(Ratio of specific heats. Used only when `chemistry` is
GasChemistry::PERFECT_GAS.)doc";

static const char *__doc_Goddard_MocOptions_geometry = R"doc(Throat and contour geometry.)doc";

static const char *__doc_Goddard_MocOptions_initial_line_axial_shift =
R"doc(Downstream shift (dimensionless, in throat radii) applied to every
station of the Kliegel-Levine transonic start line, used only for
axisymmetric ANALYSIS/DESIGN_RAO (the KL-init path; see
NozzleGeometry::downstream_wall_curvature_radius). The raw sonic
(zero-radial-velocity) locus that Kliegel-Levine solves for is not
usable as a dual-family (C+ and C-) seeding line: near the axis its
Mach angle mu approaches 90 deg, and rigidly translating every station
downstream by this amount raises the Mach number (lowering mu)
everywhere while preserving the locus's own near-axis curvature -- a
per-station constant-Mach lift does not, and empirically produces
invalid, behind-parent seeding. Default 0.1 throat radii is a moderate
lift validated against the default geometry; a larger shift trades
numerical margin for accuracy, since it extrapolates the KL series
further from the throat plane it is expanded about.)doc";

static const char *__doc_Goddard_MocOptions_inverse_cfl =
R"doc(Inverse march only: fraction of the domain-of-dependence step taken
each pass, in (0, 1]. The full domain-of-dependence step (cfl = 1) is
the largest step for which every new front point's characteristics
still trace back to a point strictly inside the previous front; a
fraction below 1 leaves margin against the linearization error in that
estimate.)doc";

static const char *__doc_Goddard_MocOptions_kl_max_wall_angle_error =
R"doc(Largest |theta_series - theta_wall| (rad) at the Kliegel-Levine line's
wall end for which the line is still corrected to the contour (a
multiplicative rescaling of the flow angle, see
MocInitialization::initialize_kliegel_levine). Above it AUTO falls
back to the centered fan and a forced KLIEGEL_LEVINE fails with
INITIALIZATION_FAILED.

Default 0.25 rad (14 deg). The default throat (r_arc = 0.382) misses
by 0.145 rad, and the corrected line is measurably better there than
the fan: monotone wall Mach, 1% start-line mass-flow error against the
fan's 46%, and higher coverage on the Rao contour. The threshold
exists to catch a series that is not describing the throat at all, not
to reject the correction where it works.)doc";

static const char *__doc_Goddard_MocOptions_log_level = R"doc(Set to MocLogLevel::DEBUG for a verbose kernel trace.)doc";

static const char *__doc_Goddard_MocOptions_max_wall_turn_per_step =
R"doc(Inverse march only: largest change of wall angle tolerated per step,
in radians. Caps the step length on a curving contour (e.g. a throat
expansion arc) so the wall point sampling stays fine enough to resolve
the turn.)doc";

static const char *__doc_Goddard_MocOptions_mode = R"doc()doc";

static const char *__doc_Goddard_MocOptions_nozzle_profile =
R"doc(Wall contour to march against in MocMode::ANALYSIS, in the length unit
of NozzleGeometry::throat_radius, with the throat at x = 0: its
smallest radius must equal `geometry.throat_radius` (within 1e-3
relative). Ignored by the design modes, which generate their own
contour and report it in MocResult::profile instead.)doc";

static const char *__doc_Goddard_MocOptions_num_characteristics = R"doc(Number of C+ lines seeded from the initial expansion fan.)doc";

static const char *__doc_Goddard_MocOptions_solver_options = R"doc(Tolerances for the iterative unit processes.)doc";

static const char *__doc_Goddard_MocOptions_start_line = R"doc(Which initial data line to build; see MocStartLine. */)doc";

static const char *__doc_Goddard_MocOptions_theta_max = R"doc(Maximum wall angle (radians), for MocMode::DESIGN_MIN_LENGTH.)doc";

static const char *__doc_Goddard_MocOptions_theta_schedule =
R"doc(Optional user-supplied theta schedule for the expansion fan.

When provided, overrides the auto-generated schedule and
num_characteristics is inferred from the schedule size. Each entry is
the flow angle (radians) of a C- characteristic from the expansion.
Must be monotonically increasing, with the last entry equal to
theta_max.)doc";

static const char *__doc_Goddard_MocPassDiagnostics =
R"doc(Per-pass record of the marching front's geometry, recorded every
kernel pass.

The front is where every known axisymmetric failure mode shows up
first, and none of it is visible in the finished net: a sampling void
contains no cells, so per-cell statistics look healthy while it grows.
These are the quantities that expose it.

Every length here (spacings, positions, step_dx) is in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_front_axis_spacing = R"doc(< Arc length of the bottom-most front segment.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_front_axis_x = R"doc(< x of the front's lowest (nearest-axis) point.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_front_points = R"doc(< Number of rungs on the front at the start of the pass.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_front_wall_spacing = R"doc(< Arc length of the top-most front segment.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_front_wall_x = R"doc(< x of the front's highest (nearest-wall) point.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_max_spacing = R"doc(< Longest front-segment arc length. The void shows up here.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_mean_spacing = R"doc(< Mean front-segment arc length.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_min_spacing = R"doc(< Shortest front-segment arc length.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_pass = R"doc(< Kernel pass this record describes.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_step_dx =
R"doc(Inverse march only: axial step length taken this pass. 0 for minimum-
length design.)doc";

static const char *__doc_Goddard_MocPassDiagnostics_step_limiter =
R"doc(Inverse march only: which limiter bound step_dx this pass; NONE for
minimum-length design.)doc";

static const char *__doc_Goddard_MocResult =
R"doc(Everything a MocNozzle::solve() produces: the flow field, the contour,
and the diagnostics.

Every length (net point and wall coordinates, `profile`,
`nozzle_length`, `exit_plane.y`, `min_theta_x`/`min_theta_y`, crossing
and failure locations, and the length fields of the diagnostics) is in
the unit of NozzleGeometry::throat_radius, recorded here as
`throat_radius`. Coordinates quoted inside `messages` are in throat
radii.)doc";

static const char *__doc_Goddard_MocResult_area_ratio = R"doc(< Exit area divided by throat area.)doc";

static const char *__doc_Goddard_MocResult_converged =
R"doc(< True iff failure.code == MocErrorCode::NONE and the kernel completed
without hitting the iteration cap.)doc";

static const char *__doc_Goddard_MocResult_crossings =
R"doc(Same-family characteristic crossings in the finished net; nonzero
means coalescence. Meaningful for chain nets, i.e. minimum-length
design; always zero for the front-based nets of analysis/Rao, which
prescribe every front directly and build no characteristic chains
(CharacteristicNet::c_chains) for this scan to see; see
CharacteristicNet::fronts for their own mesh record.)doc";

static const char *__doc_Goddard_MocResult_exit_coverage =
R"doc(Fraction of the target exit radius the solved net reached.

Always 1.0 for MocMode::DESIGN_MIN_LENGTH: the minimum-length march is
defined to stop exactly at theta_max, so there is no partial-coverage
case to report. For the front-based modes (ANALYSIS, DESIGN_RAO) the
reference-plane march (InverseMarch) places its last front exactly on
the exit plane when it converges (1.0); on a solve that fails partway,
this is instead the fraction of the target exit radius the net's last
wall point actually reached.)doc";

static const char *__doc_Goddard_MocResult_exit_mach = R"doc(< Mach number at the exit plane.)doc";

static const char *__doc_Goddard_MocResult_exit_plane = R"doc()doc";

static const char *__doc_Goddard_MocResult_failure = R"doc(< Populated when converged is false; MocErrorCode::NONE otherwise.)doc";

static const char *__doc_Goddard_MocResult_init_diagnostics = R"doc(Properties of the initial data line, measured before the march begins.)doc";

static const char *__doc_Goddard_MocResult_messages = R"doc(< Warnings and error information collected during the solve.)doc";

static const char *__doc_Goddard_MocResult_min_theta =
R"doc(Smallest flow angle anywhere in the net (radians) and where it occurs.
In a diverging nozzle a markedly negative value marks a compression
converging on the axis, i.e. a forming shock, which an isentropic
march can only pass through approximately; the solution downstream of
it is not to be trusted to better than the size of the dip. Conical
nozzles with a circular-arc throat are known to form such a shock
(Darwell & Badham 1963; Migdal & Kosson 1965). Reported for every
scheme.)doc";

static const char *__doc_Goddard_MocResult_min_theta_x =
R"doc(< Axial position of min_theta, in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_MocResult_min_theta_y =
R"doc(< Radial position of min_theta, in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_MocResult_net =
R"doc(< The solved characteristic mesh. Populated even when the march fails
partway.)doc";

static const char *__doc_Goddard_MocResult_nozzle_length =
R"doc(< Distance from throat to exit plane, in the unit of
NozzleGeometry::throat_radius.)doc";

static const char *__doc_Goddard_MocResult_pass_diagnostics =
R"doc(Per-pass marching-front geometry; one entry per kernel pass, in order.
Recorded by the inverse march only; empty for minimum-length design.)doc";

static const char *__doc_Goddard_MocResult_profile = R"doc(< Wall contour: computed in design mode, echoed back in analysis mode.)doc";

static const char *__doc_Goddard_MocResult_reached_exit_plane =
R"doc(True when the net actually reached the exit plane; false only for a
front-based (ANALYSIS/DESIGN_RAO) solve that failed partway -- always
true for minimum-length design.)doc";

static const char *__doc_Goddard_MocResult_throat_radius =
R"doc(The NozzleGeometry::throat_radius of the solve: the length unit of
every length in this result, and the reference
compute_thrust_coefficient() normalizes by.)doc";

static const char *__doc_Goddard_MocStartLine =
R"doc(Which initial data line to build. AUTO keeps today's rule (negative
NozzleGeometry::downstream_wall_curvature_radius or planar flow
selects the fan) and additionally falls back to the fan when the
Kliegel-Levine series misses the wall angle by more than
kl_max_wall_angle_error.)doc";

static const char *__doc_Goddard_MocStartLine_AUTO = R"doc()doc";

static const char *__doc_Goddard_MocStartLine_CENTERED_FAN = R"doc()doc";

static const char *__doc_Goddard_MocStartLine_KLIEGEL_LEVINE = R"doc()doc";

static const char *__doc_Goddard_MocStepLimiter =
R"doc(What limited the length of the last inverse-march step
(MocPassDiagnostics::step_limiter). */)doc";

static const char *__doc_Goddard_MocStepLimiter_CFL =
R"doc(< Bounded by the domain-of-dependence CFL condition
(MocOptions::inverse_cfl).)doc";

static const char *__doc_Goddard_MocStepLimiter_EXIT =
R"doc(< Bounded by the distance remaining to the exit plane; the front this
step builds is the last one.)doc";

static const char *__doc_Goddard_MocStepLimiter_NONE = R"doc(< No step has been taken yet.)doc";

static const char *__doc_Goddard_MocStepLimiter_WALL_FOOT =
R"doc(< Bounded so the top interior point's C- foot stays below the previous
wall point.)doc";

static const char *__doc_Goddard_MocStepLimiter_WALL_TURN = R"doc(< Bounded by MocOptions::max_wall_turn_per_step on a curving contour.)doc";

static const char *__doc_Goddard_Nozzle = R"doc()doc";

static const char *__doc_Goddard_NozzleGeometry =
R"doc(Throat and contour geometry the solver is anchored to.

`throat_radius` fixes the length unit of the whole solve:
MocOptions::nozzle_profile is read in the unit `throat_radius` is
given in, and every length in MocResult is reported in it. The two
curvature radii are dimensionless multiples of `throat_radius`.
`length_fraction` and `expansion_ratio` are consumed only by
MocMode::DESIGN_RAO, which uses them to generate the Rao contour it
then marches.)doc";

static const char *__doc_Goddard_NozzleGeometry_downstream_wall_curvature_radius =
R"doc(Wall radius of curvature downstream of the throat, as a multiple of
`throat_radius`. Set to a negative value to request centered-fan
initialization instead of a Kliegel-Levine transonic start line.)doc";

static const char *__doc_Goddard_NozzleGeometry_expansion_ratio = R"doc(Rao design only: exit-to-throat area ratio of the contour to generate.)doc";

static const char *__doc_Goddard_NozzleGeometry_length_fraction =
R"doc(Rao design only: length as a fraction of a comparable 15-degree
conical nozzle.)doc";

static const char *__doc_Goddard_NozzleGeometry_throat_radius =
R"doc(Throat radius (axisymmetric) or throat half-height (planar), in any
length unit; must be positive. The solver works in throat radii
internally and converts the input contour and every output length with
it, so its value only sets the unit: the dimensionless results do not
depend on it.)doc";

static const char *__doc_Goddard_NozzleGeometry_upstream_wall_curvature_radius =
R"doc(Wall radius of curvature upstream of the throat, as a multiple of
`throat_radius`. Not implemented: `validate_moc_options` raises
NotImplementedError for any other value.)doc";

static const char *__doc_Goddard_NozzleOptions = R"doc()doc";

static const char *__doc_Goddard_NozzleOptions_chemistry = R"doc()doc";

static const char *__doc_Goddard_NozzleOptions_expansion_ratios =
R"doc(Station ratios, interpreted by `expansion_type`: area ratios A/A_t
[-], or pressure ratios P_inlet/P [-] (chamber over station pressure,
as CEA's `pi/p`).)doc";

static const char *__doc_Goddard_NozzleOptions_expansion_type = R"doc()doc";

static const char *__doc_Goddard_NozzleOptions_frozen_NFZ =
R"doc(Freezing station for FROZEN chemistry [-], 0-based: 0 is the chamber
(the default, and CEA's default freezing point), 1 the throat, and 2
onward the expansion stations in the order they are solved. The flow
is in equilibrium up to and including this station and keeps its
composition downstream of it. With 0 the chamber derivatives are
frozen as well.

For an infinite-area combustor CEA's `nfz` is `frozen_NFZ + 1`. For a
finite-area combustor the count starts at the combustion end, so `nfz`
is `frozen_NFZ + 3`; 1 is the throat in both cases.)doc";

static const char *__doc_Goddard_NozzleProfile =
R"doc(Geometric representation of a nozzle wall profile. Used as both input
(analysis) or output (design).

The contour is stored as two parallel coordinate arrays in a
meridional plane: `x` runs downstream along the nozzle axis and `y` is
the wall radius (axisymmetric) or half-height (planar). Both carry
whatever length unit the caller supplied for the throat radius, so a
profile built with the default `r_throat = 1.0` is expressed in throat
radii.

Angles returned by the query methods are in **radians**; the static
generators take their shape angles in **degrees**, matching the
convention used to tabulate them.)doc";

static const char *__doc_Goddard_NozzleProfile_area_at =
R"doc(Calculate the cross-sectional area of nozzle at x-position assuming
the nozzle is axisymmetric.

Parameter ``x_query``:
    Axial position (length units); must lie within the profile's
    domain.

Returns:
    Cross-sectional area (length units squared).)doc";

static const char *__doc_Goddard_NozzleProfile_at =
R"doc(Coordinates of one wall point.

Parameter ``idx``:
    Index into `x`/`y`.

Returns:
    The (x, y) pair at that index, in length units.)doc";

static const char *__doc_Goddard_NozzleProfile_find_index = R"doc()doc";

static const char *__doc_Goddard_NozzleProfile_generate_Rao_TOP_nozzle =
R"doc(Generate a thrust-optimized parabolic (TOP) nozzle based on the
approximations by Rao.

@note The Rao approximation generates parameters for a Bézier curve
construction. For more control over the geometry, use
`generate_bezier_nozzle` instead.

## References

1. G. V. R. Rao, “Exhaust Nozzle Contour for Optimum Thrust,” Journal
of Jet Propulsion, vol. 28, no. 6, pp. 377–382, Jun. 1958, doi:
10.2514/8.7324.

2. G. V. R. Rao, “Approximation of optimum thrust nozzle contour,” Ars
Journal, vol. 30, no. 6, p. 561, 1960.

Parameter ``area_ratio``:
    Ratio of the exit area to the throat area. Must be at least 3.55,
    the lower bound of Rao's tabulated data.

Parameter ``r_throat``:
    Throat radius, in length units.

Parameter ``length_frac``:
    Fraction of length of comparable 15-degree conical nozzle. Must be
    between 0.6 and 1.0.

Parameter ``n_points``:
    Number of points in the profile.

Throws:
    std::invalid_argument if any argument is outside its usable range.)doc";

static const char *__doc_Goddard_NozzleProfile_generate_bezier_nozzle =
R"doc(Generate a parabolic nozzle parametrized by Bezier curves.

Parameter ``area_ratio``:
    Ratio of the exit area to the throat area.

Parameter ``theta_n``:
    Maximum expansion angle from centerline, in degrees.

Parameter ``theta_e``:
    Expansion angle from centerline at the nozzle exit, in degrees.

Parameter ``r_expansion_curve``:
    Radius of curvature of the expansion region as a fraction of the
    throat radius.

Parameter ``r_throat``:
    Throat radius, in length units.

Parameter ``length_frac``:
    Fraction of length of comparable 15-degree conical nozzle.

Parameter ``n_points``:
    Number of points in the profile.

Throws:
    std::invalid_argument if any argument is outside its usable range.)doc";

static const char *__doc_Goddard_NozzleProfile_generate_conical_nozzle =
R"doc(Generate a conical nozzle profile.

The contour is a circular expansion arc off the throat, followed by a
straight cone tangent to it. When `r_expansion_curve` is zero the arc
is omitted and the cone starts at a sharp throat corner.

Parameter ``area_ratio``:
    Ratio of the exit area to the throat area.

Parameter ``r_expansion_curve``:
    Radius of curvature of the throat expansion arc, as a multiple of
    `r_throat`. Zero gives a sharp throat corner.

Parameter ``r_throat``:
    Throat radius, in length units.

Parameter ``theta_n``:
    Conical expansion angle from the centerline, in degrees.

Parameter ``n_points``:
    Number of points in the profile.

Throws:
    std::invalid_argument if any argument is outside its usable range,
    or if the expansion arc alone already exceeds the exit radius.)doc";

static const char *__doc_Goddard_NozzleProfile_generate_throat_expansion_curve =
R"doc(Generate just the circular expansion arc that turns the flow off the
throat.

This is the leading segment shared by the conical and Bézier contours,
useful on its own as the wall geometry a method-of-characteristics
solve is seeded against.

Parameter ``theta_n``:
    Wall angle the arc turns to, measured from the centerline, in
    degrees.

Parameter ``r_expansion_curve``:
    Radius of curvature of the arc, as a fraction of `r_throat`.

Parameter ``r_throat``:
    Throat radius, in length units.

Parameter ``n_points``:
    Number of points in the arc.)doc";

static const char *__doc_Goddard_NozzleProfile_length =
R"doc(Get the length of the nozzle, measured from the throat, in length
units. */)doc";

static const char *__doc_Goddard_NozzleProfile_load_profile_csv =
R"doc(Read a profile from a two-column (x, y) CSV file.

Parameter ``filename``:
    Path to the CSV file.)doc";

static const char *__doc_Goddard_NozzleProfile_max_theta = R"doc(Find the maximum wall angle across all segments, in radians. */)doc";

static const char *__doc_Goddard_NozzleProfile_populate_throat_expansion_curve =
R"doc(Populate profile with expansion curve from the throat, taking into
account curvature radius at the throat.

.. warning::
    Appends to the existing profile!)doc";

static const char *__doc_Goddard_NozzleProfile_push_back = R"doc(Append a wall point, given as an (x, y) pair in length units. */)doc";

static const char *__doc_Goddard_NozzleProfile_radius_at =
R"doc(Interpolate the wall radius at a given axial position.

Parameter ``x_query``:
    Axial position (length units); must lie within the profile's
    domain.

Returns:
    Wall radius (length units).)doc";

static const char *__doc_Goddard_NozzleProfile_radius_max =
R"doc(Largest wall radius on the profile.

Returns:
    Index of the widest point and its radius (length units).)doc";

static const char *__doc_Goddard_NozzleProfile_save_profile_csv =
R"doc(Write the profile to a two-column (x, y) CSV file.

Parameter ``filename``:
    Path to write to.)doc";

static const char *__doc_Goddard_NozzleProfile_scaled =
R"doc(Copy of this profile in another length unit: every `x` and `y`
multiplied by `factor`. Angles and slopes are unchanged;
`throat_index` is kept.

Parameter ``factor``:
    Length scale factor, e.g. 1 / r_throat to express the profile in
    throat radii. Must be positive.

Throws:
    std::invalid_argument if `factor` is not positive.)doc";

static const char *__doc_Goddard_NozzleProfile_size = R"doc(Number of points on the profile. */)doc";

static const char *__doc_Goddard_NozzleProfile_slope_at =
R"doc(Interpolate wall slope dy/dx at a given axial position. Only used in
analysis mode.

Parameter ``x_query``:
    Axial position (length units); must lie within the profile's
    domain.

Returns:
    Dimensionless slope dy/dx.)doc";

static const char *__doc_Goddard_NozzleProfile_slope_at_idx =
R"doc(Wall slope dy/dx at a given point index, by finite difference.

Second-order (weighted two-sided) wherever a downstream neighbour
exists, first-order (backward) at the last point.

Parameter ``idx``:
    Index into `x`/`y`; must be at least 1.

Returns:
    Dimensionless slope dy/dx.)doc";

static const char *__doc_Goddard_NozzleProfile_theta_at =
R"doc(Interpolate the wall angle at a given axial position.

Parameter ``x_query``:
    Axial position (length units); must lie within the profile's
    domain.

Returns:
    Wall angle atan(dy/dx), in radians.)doc";

static const char *__doc_Goddard_NozzleProfile_theta_at_idx =
R"doc(Wall angle at a given point index.

Parameter ``idx``:
    Index into `x`/`y`; must be at least 1.

Returns:
    Wall angle atan(dy/dx), in radians.)doc";

static const char *__doc_Goddard_NozzleProfile_theta_at_interpolated =
R"doc(Wall angle at a given axial position, linearly interpolated between
the vertex angles theta_at_idx() returns.

theta_at() looks up the bracketing facet [x[idx-1], x[idx]] and
returns its single finite-difference angle theta_at_idx(idx) unchanged
for every query inside that facet, i.e. it is piecewise constant in x:
it reports the angle of the facet a query point belongs to. A caller
that samples the wall at stations finer than the profile's own facets
(the inverse march's per-pass wall point) instead needs an angle that
is continuous in x, so this interpolates linearly between the same
theta_at_idx() vertex values across each facet rather than snapping to
one of them.

Parameter ``x_query``:
    Axial position (length units); must lie within the profile's
    domain.

Returns:
    Wall angle atan(dy/dx), continuous in x, in radians.)doc";

static const char *__doc_Goddard_NozzleProfile_throat_index =
R"doc(Index into `x`/`y` of the throat, i.e. the origin `length()` is
measured from. */)doc";

static const char *__doc_Goddard_NozzleProfile_x = R"doc(Axial coordinates of the wall points, ascending (length units). */)doc";

static const char *__doc_Goddard_NozzleProfile_x_max =
R"doc(Axial coordinate of the last (most downstream) point, in length units.
*/)doc";

static const char *__doc_Goddard_NozzleProfile_x_min =
R"doc(Axial coordinate of the first (most upstream) point, in length units.
*/)doc";

static const char *__doc_Goddard_NozzleProfile_y =
R"doc(Wall radius (axisymmetric) or half-height (planar) at each point in
`x` (length units). */)doc";

static const char *__doc_Goddard_NozzleResults = R"doc()doc";

static const char *__doc_Goddard_NozzleResults_expansions = R"doc()doc";

static const char *__doc_Goddard_NozzleResults_inlet =
R"doc(Nozzle inlet (chamber) station. Its derivatives are in the chemistry
of station 0: frozen only for FROZEN chemistry with `frozen_NFZ` == 0,
otherwise equilibrium.)doc";

static const char *__doc_Goddard_NozzleResults_throat = R"doc()doc";

static const char *__doc_Goddard_NozzleStation = R"doc(One solved nozzle station.)doc";

static const char *__doc_Goddard_NozzleStation_area_ratio =
R"doc(Area ratio A/A_t [-], from the mass flux relative to the throat. 0 for
a station where the flow is at rest (the chamber of an infinite-area
combustor, the injector face and the stagnation state of a finite-area
one).)doc";

static const char *__doc_Goddard_NozzleStation_mach = R"doc(Mach number [-], from `velocity` and `thermo.speed_of_sound`. */)doc";

static const char *__doc_Goddard_NozzleStation_state = R"doc(Raw state vector, for `Gas::restore_state`. */)doc";

static const char *__doc_Goddard_NozzleStation_thermo =
R"doc(Mixture state at the station. `gamma_s`, `dlV_dlP_T`, `dlV_dlT_P`,
`speed_of_sound` and `pinned_transition` are in the chemistry of the
station: frozen downstream of the freezing station, equilibrium
otherwise. `pinned_transition` is true when the station sits exactly
at a condensed phase transition (see
`ThroatCondition::pinned_transition`).)doc";

static const char *__doc_Goddard_NozzleStation_velocity =
R"doc(Flow velocity [m/s], from the enthalpy drop below the stagnation
enthalpy. */)doc";

static const char *__doc_Goddard_Nozzle_Nozzle =
R"doc(Parameter ``gas``:
    gas at the inlet (chamber) state. The solver works on its own copy
    of `gas`; the caller's `Gas` is not modified.

Parameter ``options``:
    chemistry, freezing station and stations to solve.

Throws:
    std::invalid_argument for KINETIC chemistry.

Throws:
    NotImplementedError for PERFECT_GAS chemistry.)doc";

static const char *__doc_Goddard_Nozzle_Nozzle_2 = R"doc()doc";

static const char *__doc_Goddard_Nozzle_get_gamma_s = R"doc()doc";

static const char *__doc_Goddard_Nozzle_get_inlet_state = R"doc()doc";

static const char *__doc_Goddard_Nozzle_inlet_state = R"doc()doc";

static const char *__doc_Goddard_Nozzle_is_equilibrium_station =
R"doc(True when station `station` is in equilibrium, from
`NozzleOptions::frozen_NFZ`. */)doc";

static const char *__doc_Goddard_Nozzle_iterate_area_expansion = R"doc()doc";

static const char *__doc_Goddard_Nozzle_iterate_temperature =
R"doc(Temperature [K] of a frozen station: the isentrope at fixed
composition.

With condensed products the amounts of the condensed species are
frozen too, so the temperature cannot leave the data range of any
species that is present.

Parameter ``throat_condition``:
    Throat state the expansion starts from.

Parameter ``pressure_ratio``:
    Chamber pressure divided by the station pressure [-].

Parameter ``T_guess``:
    Starting temperature [K].

Parameter ``composition``:
    Gas-phase mole fractions held fixed [-].

Parameter ``station``:
    Index of the station being solved, for error messages [-].

Parameter ``abstol``:
    Convergence tolerance on d(log T) [-].

Returns:
    Station temperature [K].

Throws:
    FmtError if a present condensed species leaves its temperature
    range, as CEA reports for a frozen expansion carried too far.)doc";

static const char *__doc_Goddard_Nozzle_m_gas = R"doc()doc";

static const char *__doc_Goddard_Nozzle_m_opts = R"doc()doc";

static const char *__doc_Goddard_Nozzle_reset_state = R"doc()doc";

static const char *__doc_Goddard_Nozzle_set_inlet_state = R"doc()doc";

static const char *__doc_Goddard_Nozzle_set_station_chemistry =
R"doc(Set the gas chemistry for station `station`.

Returns:
    True for an equilibrium station.)doc";

static const char *__doc_Goddard_Nozzle_solve = R"doc()doc";

static const char *__doc_Goddard_Nozzle_solve_2 = R"doc()doc";

static const char *__doc_Goddard_Nozzle_solve_3 = R"doc()doc";

static const char *__doc_Goddard_Nozzle_solve_finite_area_chamber =
R"doc(Solve the chamber of a finite-area combustor.

Iterates on the stagnation pressure P_inf until the chamber momentum
balance P_inj = P_c + rho_c u_c^2 holds, where c is the combustion
end. The chamber is always in equilibrium. On return the nozzle inlet
state is the stagnation state, so `solve_stations(result.throat, ...)`
continues the expansion, and the gas holds the stagnation state.

The combustion-end velocity follows from the enthalpy drop h_inj -
h_c. At very large contraction ratios (Mach numbers of order 1e-3 and
below) that drop is within the station solve's error, and a warning
says the combustion-end velocity and Mach number are not resolved.
Pressures and the momentum balance are unaffected.

Parameter ``injector_state``:
    Equilibrium state at the injector face (HP at P_inj).

Parameter ``type``:
    `FINITE_CONTRACTION_RATIO` or `FINITE_MASS_FLUX`.

Parameter ``value``:
    Contraction ratio A_c/A_t [-] or mass flux mdot/A_c [kg/(m^2 s)].

Parameter ``reltol``:
    Tolerance on |1 - P_inj,calc / P_inj| [-].

Returns:
    Converged chamber.

Throws:
    std::invalid_argument if `type` is not a finite-area type, if the
    contraction ratio is not above 1, or if the mass flux thermally
    chokes the chamber: it exceeds the perfect-gas limit by more than
    2 %, or A_c/A_t reaches 1 in the real-gas iteration.

Throws:
    NotImplementedError for FROZEN chemistry with `frozen_NFZ` == 0.

Throws:
    ConvergenceError if the momentum balance does not converge.)doc";

static const char *__doc_Goddard_Nozzle_solve_pressure_ratio = R"doc()doc";

static const char *__doc_Goddard_Nozzle_solve_stations =
R"doc(Solve expansion stations from an already solved throat.

Parameter ``throat_condition``:
    Throat of the current inlet state.

Parameter ``expansion_type``:
    How `ratios` are interpreted.

Parameter ``ratios``:
    Area ratios A/A_t [-] or pressure ratios P_inlet/P [-].

Returns:
    One station per ratio, in order.)doc";

static const char *__doc_Goddard_Nozzle_solve_stations_2 =
R"doc(Solve expansion stations downstream of a finite-area combustor.

Parameter ``chamber``:
    Chamber returned by `solve_finite_area_chamber`.

Parameter ``expansion_type``:
    How `ratios` are interpreted.

Parameter ``ratios``:
    Area ratios A/A_t [-] or pressure ratios P_inj/P [-]. Pressure
    ratios are taken from the injector pressure, as CEA does, not from
    the stagnation pressure.

Returns:
    One station per ratio, in order.)doc";

static const char *__doc_Goddard_Nozzle_solve_subsonic_area_expansion = R"doc()doc";

static const char *__doc_Goddard_Nozzle_solve_supersonic_area_expansion = R"doc()doc";

static const char *__doc_Goddard_Nozzle_solve_throat_conditions = R"doc()doc";

static const char *__doc_Goddard_Nozzle_station_at_current_state =
R"doc(Station at the current gas state, with derivatives `props` already
solved for it.

Parameter ``H_stagnation``:
    Stagnation enthalpy the velocity is measured from [J/kg].

Parameter ``area_per_mdot_throat``:
    Throat area per mass flow rate [m^2 s/kg]; 0 when there is no
    throat yet, which leaves `area_ratio` at 0.)doc";

static const char *__doc_Goddard_Nozzle_station_at_state =
R"doc(Station at `state` with derivatives in the chemistry of station
`station`. Leaves the gas at `state`.)doc";

static const char *__doc_Goddard_Nozzle_throat_area_per_mdot =
R"doc(Throat area per mass flow rate [m^2 s/kg], from the throat state and
velocity. */)doc";

static const char *__doc_Goddard_Nozzle_throw_invalid_expansion_ratio = R"doc()doc";

static const char *__doc_Goddard_ObliqueShockResult =
R"doc(Oblique shock: wave and deflection angles, and the jump across the
normal component. */)doc";

static const char *__doc_Goddard_ObliqueShockResult_beta = R"doc(Wave angle between the shock and the upstream flow [rad]. */)doc";

static const char *__doc_Goddard_ObliqueShockResult_mach_in = R"doc(Upstream Mach number [-]. */)doc";

static const char *__doc_Goddard_ObliqueShockResult_mach_out = R"doc(Downstream Mach number, including the tangential velocity [-]. */)doc";

static const char *__doc_Goddard_ObliqueShockResult_shock = R"doc(Jump across the normal component of the flow. */)doc";

static const char *__doc_Goddard_ObliqueShockResult_theta = R"doc(Flow deflection angle [rad]. */)doc";

static const char *__doc_Goddard_ObliqueShockResult_valid = R"doc(False if the inputs admit no attached oblique shock. */)doc";

static const char *__doc_Goddard_PhaseSpecification = R"doc(Convenience class for initializing a thermodynamic state.)doc";

static const char *__doc_Goddard_PhaseSpecification_P = R"doc()doc";

static const char *__doc_Goddard_PhaseSpecification_PhaseSpecification = R"doc()doc";

static const char *__doc_Goddard_PhaseSpecification_PhaseSpecification_2 = R"doc()doc";

static const char *__doc_Goddard_PhaseSpecification_PhaseSpecification_3 = R"doc()doc";

static const char *__doc_Goddard_PhaseSpecification_T = R"doc()doc";

static const char *__doc_Goddard_PhaseSpecification_composition = R"doc()doc";

static const char *__doc_Goddard_PointMembership =
R"doc(The chains one point belongs to.

An interior point lies on exactly one characteristic of each family,
so it leads both. A wall, axis, or initial-line point may belong to
only one.)doc";

static const char *__doc_Goddard_PointMembership_c_minus_chain_idx =
R"doc(Index into CharacteristicNet::c_chains of the C- chain through this
point, if any.)doc";

static const char *__doc_Goddard_PointMembership_c_plus_chain_idx =
R"doc(Index into CharacteristicNet::c_chains of the C+ chain through this
point, if any.)doc";

static const char *__doc_Goddard_PrandtlMeyerTable = R"doc(Table of thermodynamic data for Prandtl-Meyer expansion fans)doc";

static const char *__doc_Goddard_PrandtlMeyerTable_build_table =
R"doc(Build table of (velocity, nu) pairs along an isentropic expansion.
Steps in pressure from throat conditions to a specified pressure
ratio.

Parameter ``thermo``:
    Cantera ThermoPhase, set to throat conditions on entry

Parameter ``equilibrium``:
    If true, equilibrate at each point; otherwise frozen composition

Parameter ``s0``:
    Stagnation entropy (J/kg/K)

Parameter ``h0``:
    Stagnation enthalpy (J/kg)

Parameter ``a_throat``:
    Sound speed at the throat (m/s)

Parameter ``pressure_ratio``:
    Ratio P_exit/P_throat to expand to (default: 1e-4)

Parameter ``num_points``:
    Number of table points)doc";

static const char *__doc_Goddard_PrandtlMeyerTable_built = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_dnu_dV = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_enthalpies = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_equil = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_find_V_index_and_weight = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_find_closest_nMv_index =
R"doc(Find the index closest to a queried value using binary search. This
works for nu, Mach number, and velocity as they monotonically increase
with expansion.)doc";

static const char *__doc_Goddard_PrandtlMeyerTable_find_mach_index_and_weight = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_find_nu_index_and_weight = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_gamma_s = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_index_and_weight = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interp = R"doc(Generic linear interpolation with binary search.)doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interp_state_vector = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_P_from_nu = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_T_from_nu = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_V = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_V_from_mach = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_at_index = R"doc(Linear interpolation over vals with a specified weight.)doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_dnu = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_gamma_s_from_mach = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_gamma_s_from_nu = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_h = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_h_from_mach = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_h_from_nu = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_mach = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_mach_from_V = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_nu = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_nu_from_mach = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_state_at_index = R"doc(Linear interpolation over vals with a specified weight.)doc";

static const char *__doc_Goddard_PrandtlMeyerTable_interpolate_state_from_mach = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_is_built = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_is_equilibrium = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_machs = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_nus = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_pressures = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_states = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_temperatures = R"doc()doc";

static const char *__doc_Goddard_PrandtlMeyerTable_velocities = R"doc()doc";

static const char *__doc_Goddard_ReflectedDetonationResult =
R"doc(A detonation reflected as a shock from the closed end of a tube.

The incident detonation leaves its products (state 2) moving toward
the end wall; the reflected shock brings them back to rest (state 5),
at equilibrium.)doc";

static const char *__doc_Goddard_ReflectedDetonationResult_incident =
R"doc(Incident detonation, state 1 to state 2, in the frame of the
detonation. */)doc";

static const char *__doc_Goddard_ReflectedDetonationResult_reflected =
R"doc(Reflected shock, state 2 to state 5, in the reflected-shock frame.
`mach_in` is the reflected-shock Mach number relative to the products,
with their equilibrium sound speed.)doc";

static const char *__doc_Goddard_ReflectedDetonationResult_valid = R"doc(False if the incident detonation or the reflected shock is invalid. */)doc";

static const char *__doc_Goddard_ReflectedShockResult =
R"doc(Incident and reflected normal shocks at the closed end of a shock
tube.

State 1 is the gas at rest ahead of the incident shock, state 2 the
gas behind it, and state 5 the gas brought back to rest behind the
shock reflected from the end wall.)doc";

static const char *__doc_Goddard_ReflectedShockResult_incident = R"doc(Incident shock, state 1 to state 2, in the incident-shock frame. */)doc";

static const char *__doc_Goddard_ReflectedShockResult_reflected =
R"doc(Reflected shock, state 2 to state 5, in the reflected-shock frame.
`mach_in` is the reflected-shock Mach number relative to gas 2.)doc";

static const char *__doc_Goddard_ReflectedShockResult_valid = R"doc(False if either shock is invalid. */)doc";

static const char *__doc_Goddard_RocketCaseParameters = R"doc()doc";

static const char *__doc_Goddard_RocketCaseParameters_combustor_options = R"doc()doc";

static const char *__doc_Goddard_RocketCaseParameters_name = R"doc()doc";

static const char *__doc_Goddard_RocketCaseParameters_nozzle_options = R"doc()doc";

static const char *__doc_Goddard_RocketPerformance = R"doc(Rocket performance of one exit, as CEA reports it. */)doc";

static const char *__doc_Goddard_RocketPerformance_CF = R"doc(Thrust coefficient at optimum expansion [-]. */)doc";

static const char *__doc_Goddard_RocketPerformance_area_ratio = R"doc(Exit area over throat area, Ae/At [-]. */)doc";

static const char *__doc_Goddard_RocketPerformance_cstar = R"doc(Characteristic velocity c* [m/s]. */)doc";

static const char *__doc_Goddard_RocketPerformance_isp =
R"doc(Specific impulse at optimum expansion, as an effective exhaust
velocity [m/s]. */)doc";

static const char *__doc_Goddard_RocketPerformance_isp_s =
R"doc(Specific impulse at optimum expansion in seconds, `isp /
STANDARD_GRAVITY` [s]. */)doc";

static const char *__doc_Goddard_RocketPerformance_ivac = R"doc(Vacuum specific impulse, as an effective exhaust velocity [m/s]. */)doc";

static const char *__doc_Goddard_RocketPerformance_ivac_s = R"doc(Vacuum specific impulse in seconds, `ivac / STANDARD_GRAVITY` [s]. */)doc";

static const char *__doc_Goddard_RocketPerformance_mach_number = R"doc(Exit Mach number [-]. */)doc";

static const char *__doc_Goddard_RocketPerformance_pressure_ratio = R"doc(Chamber (or stagnation) pressure over exit pressure, Pc/Pe [-]. */)doc";

static const char *__doc_Goddard_RocketProblem = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_OF_ratios = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_chemistry = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_combustor_type = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_expansion_ratios = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_expansion_type = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_finite_area_chambers =
R"doc(Finite-area chambers, indexed like `nozzle_states`; empty for
`INFINITE_AREA`. */)doc";

static const char *__doc_Goddard_RocketProblemCaseResult_frozen_NFZ =
R"doc(Freezing station of a FROZEN nozzle; see `NozzleOptions::frozen_NFZ`.
*/)doc";

static const char *__doc_Goddard_RocketProblemCaseResult_inlet_states = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_mixture_type =
R"doc(How `OF_ratios` are interpreted; see `CombustorOptions::mixture_type`.
*/)doc";

static const char *__doc_Goddard_RocketProblemCaseResult_nozzle_states = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_pressures = R"doc()doc";

static const char *__doc_Goddard_RocketProblemCaseResult_process = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_2 = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_chemistry = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_combustor_type = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_expansion_ratios = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_expansion_type = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_frozen_NFZ =
R"doc(Freezing station of a FROZEN nozzle; see `NozzleOptions::frozen_NFZ`.
*/)doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_mass_flux =
R"doc(Mass flux mdot/Ac [kg/(m^2 s)] of each operating point, row-major over
(of_index, pressure_index); empty for `CombustorType::INFINITE_AREA`.
Not derivable from a station's `ThermodynamicState`, so it is captured
here for the report.)doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_mixture_type = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_of_ratios = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_pressures = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_CaseMeta_process = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_RocketProblemResults =
R"doc(Build the flat station list from the results of every case.

Parameter ``case_results``:
    Per-case combustion and nozzle states, consumed by this call.

Parameter ``gas``:
    Product mixture the states were computed with. It carries the
    candidate condensed species, so each station's state vector is
    read back with the condensed amounts it was saved with.)doc";

static const char *__doc_Goddard_RocketProblemResults_calculate_performance = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_case_names = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_chamber =
R"doc(Chamber station of one operating point. For a finite-area combustor,
the injector face.

Parameter ``of_index``:
    Index into the case's mixture ratios [-].

Parameter ``pressure_index``:
    Index into the case's chamber pressures [-].

Parameter ``case_name``:
    Case to read. May be empty when the results hold a single case.

Throws:
    std::runtime_error if the case or the station does not exist.)doc";

static const char *__doc_Goddard_RocketProblemResults_combustion_end =
R"doc(Combustion-end station of a finite-area combustor.

Parameter ``of_index``:
    Index into the case's mixture ratios [-].

Parameter ``pressure_index``:
    Index into the case's chamber pressures [-].

Parameter ``case_name``:
    Case to read. May be empty when the results hold a single case.

Throws:
    std::runtime_error if the case uses an infinite-area combustor, or
    if the case or the station does not exist.)doc";

static const char *__doc_Goddard_RocketProblemResults_exits =
R"doc(Exit stations of one operating point, ordered by expansion index.

Parameter ``of_index``:
    Index into the case's mixture ratios [-].

Parameter ``pressure_index``:
    Index into the case's chamber pressures [-].

Parameter ``case_name``:
    Case to read. May be empty when the results hold a single case.

Throws:
    std::runtime_error if the case does not exist.)doc";

static const char *__doc_Goddard_RocketProblemResults_m_case_meta = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_m_gas = R"doc(Product mixture used to read the stored station states back. */)doc";

static const char *__doc_Goddard_RocketProblemResults_m_stations = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_of_ratios = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_performance =
R"doc(Rocket performance of one operating point at one exit station, from
its stagnation, throat and exit states (see `calculate_performance`).

Parameter ``of_index``:
    Index into the case's mixture ratios [-].

Parameter ``pressure_index``:
    Index into the case's chamber pressures [-].

Parameter ``exit_index``:
    Index into the exit stations of the operating point [-].

Parameter ``case_name``:
    Case to read. May be empty when the results hold a single case.

Throws:
    std::runtime_error if the case or a station does not exist, or if
    `exit_index` is out of range.)doc";

static const char *__doc_Goddard_RocketProblemResults_report = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_resolve_case = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_stagnation =
R"doc(Stagnation state that the nozzle expands from: the "inf" state of a
finite-area combustor, or the chamber state of an infinite-area
combustor.

Parameter ``of_index``:
    Index into the case's mixture ratios [-].

Parameter ``pressure_index``:
    Index into the case's chamber pressures [-].

Parameter ``case_name``:
    Case to read. May be empty when the results hold a single case.

Throws:
    std::runtime_error if the case or the station does not exist.)doc";

static const char *__doc_Goddard_RocketProblemResults_stations = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_stations_of_type = R"doc()doc";

static const char *__doc_Goddard_RocketProblemResults_throat =
R"doc(Throat station of one operating point.

Parameter ``of_index``:
    Index into the case's mixture ratios [-].

Parameter ``pressure_index``:
    Index into the case's chamber pressures [-].

Parameter ``case_name``:
    Case to read. May be empty when the results hold a single case.

Throws:
    std::runtime_error if the case or the station does not exist.)doc";

static const char *__doc_Goddard_RocketProblem_RocketProblem =
R"doc(Parameter ``chem_params``:
    Thermodynamic data, reactants and mixture ratios.

Parameter ``cases``:
    Combustor and nozzle settings of each case.

Parameter ``phase_name``:
    Phase of `chem_params.thermo_file` to build the products from;
    empty selects the first phase in the file.

Parameter ``transport``:
    Not implemented; `solve` raises NotImplementedError if true.

Parameter ``ionized_species``:
    Not implemented; `solve` raises NotImplementedError if true.

Parameter ``trace``:
    Not implemented; `solve` raises NotImplementedError if non-zero.)doc";

static const char *__doc_Goddard_RocketProblem_chemical_params = R"doc()doc";

static const char *__doc_Goddard_RocketProblem_include_ionized_species = R"doc(Not implemented: `solve` raises `NotImplementedError` if true. */)doc";

static const char *__doc_Goddard_RocketProblem_include_transport = R"doc(Not implemented: `solve` raises `NotImplementedError` if true. */)doc";

static const char *__doc_Goddard_RocketProblem_m_condensed_prototype =
R"doc(Product gas carrying the candidate condensed species, built once so
that every `Gas` handed to a solver shares the same set. Empty when no
`condensed_file` was given.)doc";

static const char *__doc_Goddard_RocketProblem_m_sln = R"doc()doc";

static const char *__doc_Goddard_RocketProblem_problem_cases = R"doc()doc";

static const char *__doc_Goddard_RocketProblem_product_gas =
R"doc(Product `Gas` for one solver stage. Every returned `Gas` references
the same `Solution` and the same candidate condensed species set, so
attaching candidates once in the constructor is enough for the
combustor and the nozzle to see them. The chemistry mode is left at
its default: the combustor does not read it, and the nozzle and the
results set their own.)doc";

static const char *__doc_Goddard_RocketProblem_reactant_gas =
R"doc(Reactant stream `Gas` built from `ChemicalParameters::reactant_file`
and set to `state`. */)doc";

static const char *__doc_Goddard_RocketProblem_solution = R"doc()doc";

static const char *__doc_Goddard_RocketProblem_solve = R"doc()doc";

static const char *__doc_Goddard_RocketProblem_thermo = R"doc()doc";

static const char *__doc_Goddard_RocketProblem_trace_cutoff = R"doc(Not implemented: `solve` raises `NotImplementedError` if non-zero. */)doc";

static const char *__doc_Goddard_RocketStation = R"doc()doc";

static const char *__doc_Goddard_RocketStation_area_ratio =
R"doc(A/A_t [-]: 0 for CHAMBER and STAGNATION, A_c/A_t for COMBUSTION_END, 1
for THROAT, Ae/At for EXIT (also for pressure-ratio exits). */)doc";

static const char *__doc_Goddard_RocketStation_case_name = R"doc()doc";

static const char *__doc_Goddard_RocketStation_expansion_index = R"doc()doc";

static const char *__doc_Goddard_RocketStation_of_index = R"doc()doc";

static const char *__doc_Goddard_RocketStation_pressure_index = R"doc()doc";

static const char *__doc_Goddard_RocketStation_thermo =
R"doc(Mixture state of the station. Its `stagnation_enthalpy` [J/kg] is the
enthalpy of the stagnation state of the station's operating point (see
`RocketProblemResults::stagnation`), which the adiabatic expansion
conserves.)doc";

static const char *__doc_Goddard_RocketStation_type = R"doc()doc";

static const char *__doc_Goddard_ShockResult =
R"doc(Jump conditions across one normal shock, in the frame of that shock.

Ratios are downstream over upstream. For an oblique shock this is the
normal component of the flow, except `total_pressure_ratio`, which
uses the full velocity.

An invalid result (subsonic inflow, gamma <= 1) has `valid = false`
and every other field set to -1.)doc";

static const char *__doc_Goddard_ShockResult_density_ratio =
R"doc(Density ratio rho_out/rho_in [-]. It is also the ratio of upstream to
downstream normal velocity. */)doc";

static const char *__doc_Goddard_ShockResult_mach_in =
R"doc(Upstream Mach number relative to the shock [-]. For the incident shock
of `ShockSolver` it uses the frozen sound speed of the pre-shock
state.)doc";

static const char *__doc_Goddard_ShockResult_mach_out =
R"doc(Downstream Mach number relative to the shock [-], with the sound speed
of the solver's chemistry. */)doc";

static const char *__doc_Goddard_ShockResult_static_pressure_ratio = R"doc(Static pressure ratio P_out/P_in [-]. */)doc";

static const char *__doc_Goddard_ShockResult_static_temperature_ratio = R"doc(Static temperature ratio T_out/T_in [-]. */)doc";

static const char *__doc_Goddard_ShockResult_total_pressure_ratio = R"doc(Stagnation pressure ratio P0_out/P0_in [-]. */)doc";

static const char *__doc_Goddard_ShockResult_valid = R"doc(False if the inputs admit no shock. */)doc";

static const char *__doc_Goddard_ShockSolver =
R"doc(Normal, reflected and oblique shocks in a real gas, following Gordon &
McBride, NASA RP-1311 Part I, chapter 7.

The chemistry of the `Gas` selects the model behind the shock:
PERFECT_GAS uses the perfect-gas relations with `Gas::gamma_s()`,
FROZEN holds the composition fixed, and EQUILIBRIUM brings the post-
shock gas to chemical equilibrium. KINETIC chemistry and condensed
species are not supported.

The pre-shock state is the state of the `Gas` at construction. It is
taken as given: it need not be at equilibrium (e.g. an unburned fuel-
oxidizer mixture), and its sound speed is always the frozen one. Mach-
number inputs are converted to velocity with that frozen sound speed.

@note An equilibrium shock in an exothermic mixture only exists above
the Chapman-Jouguet detonation speed; below it the solver throws
`ConvergenceError`. Above it the shock is an overdriven detonation,
which `DetonationSolver` also computes. Oblique shocks in such a
mixture are oblique detonations and are not modelled.)doc";

static const char *__doc_Goddard_ShockSolver_ShockSolver =
R"doc(Parameter ``gas``:
    pre-shock gas and chemistry model. The solver works on its own
    copy of `gas`; the caller's `Gas` is not modified.

Parameter ``options``:
    Newton tolerance on the log pressure and temperature ratios, and
    iteration limit

Throws:
    std::invalid_argument for KINETIC chemistry)doc";

static const char *__doc_Goddard_ShockSolver_apply_perfect_gas_jump =
R"doc(Apply a perfect-gas jump to the gas, which must hold the pre-shock
state: T and P are multiplied by the given ratios [-] and the
composition is kept. An invalid jump leaves the gas unchanged.)doc";

static const char *__doc_Goddard_ShockSolver_m_chemistry = R"doc()doc";

static const char *__doc_Goddard_ShockSolver_m_gas = R"doc()doc";

static const char *__doc_Goddard_ShockSolver_m_options = R"doc()doc";

static const char *__doc_Goddard_ShockSolver_m_post_shock_state = R"doc()doc";

static const char *__doc_Goddard_ShockSolver_m_pre_shock_state = R"doc()doc";

static const char *__doc_Goddard_ShockSolver_max_deflection =
R"doc(Largest deflection an attached oblique shock can produce.

`post_shock_state()` is then the state behind the shock at maximum
deflection.

Parameter ``mach``:
    upstream Mach number, with the frozen pre-shock sound speed [-]

Returns:
    maximum deflection angle [rad]; NaN if mach < 1)doc";

static const char *__doc_Goddard_ShockSolver_normal_shock =
R"doc(Normal shock.

Parameter ``mach``:
    upstream Mach number, with the frozen pre-shock sound speed [-]

Returns:
    jump conditions; invalid if mach < 1

Throws:
    std::invalid_argument if condensed species are present before or
    form behind the shock

Throws:
    ConvergenceError if the jump conditions cannot be solved)doc";

static const char *__doc_Goddard_ShockSolver_normal_shock_from_velocity =
R"doc(Normal shock from the upstream velocity.

Parameter ``velocity``:
    upstream velocity relative to the shock [m/s]

Returns:
    jump conditions; invalid if the velocity is below the frozen pre-
    shock sound speed)doc";

static const char *__doc_Goddard_ShockSolver_oblique_shock_from_deflection =
R"doc(Oblique shock with a given deflection.

For real gases the maximum deflection is found by golden-section
search, and the wave angle by bisection on the weak or strong branch.

Parameter ``mach``:
    upstream Mach number, with the frozen pre-shock sound speed [-]

Parameter ``deflection``:
    flow deflection [rad]

Parameter ``weak``:
    true for the weak solution, false for the strong one

Returns:
    oblique shock; invalid if the deflection exceeds the maximum
    deflection)doc";

static const char *__doc_Goddard_ShockSolver_oblique_shock_from_wave_angle =
R"doc(Oblique shock with a given wave angle.

Parameter ``mach``:
    upstream Mach number, with the frozen pre-shock sound speed [-]

Parameter ``wave_angle``:
    wave angle [rad]

Returns:
    oblique shock; invalid if the normal Mach number is below 1)doc";

static const char *__doc_Goddard_ShockSolver_post_shock_state =
R"doc(Gas at the state behind the most recently solved shock.

For PERFECT_GAS chemistry it is the pre-shock gas with its temperature
and pressure multiplied by the perfect-gas jump ratios, and its
composition unchanged.

The result is an independent copy: later solves do not change it.)doc";

static const char *__doc_Goddard_ShockSolver_pre_shock_sound_speed = R"doc(Frozen sound speed of the pre-shock state [m/s]. */)doc";

static const char *__doc_Goddard_ShockSolver_pre_shock_state = R"doc(Gas at the pre-shock state: an independent copy. */)doc";

static const char *__doc_Goddard_ShockSolver_reflected_shock =
R"doc(Incident and reflected shocks in a shock tube, both with the solver's
chemistry.

Parameter ``mach``:
    incident shock Mach number, with the frozen pre-shock sound speed
    [-]

Returns:
    both jumps; `post_shock_state()` is then state 5)doc";

static const char *__doc_Goddard_ShockSolver_reflected_shock_2 =
R"doc(Incident and reflected shocks in a shock tube, each with its own
chemistry.

Parameter ``mach``:
    incident shock Mach number, with the frozen pre-shock sound speed
    [-]

Parameter ``incident_chemistry``:
    FROZEN or EQUILIBRIUM, for state 2

Parameter ``reflected_chemistry``:
    FROZEN or EQUILIBRIUM, for state 5

Throws:
    std::invalid_argument for any other chemistry, or if the solver's
    gas is PERFECT_GAS)doc";

static const char *__doc_Goddard_ShockSolver_reflected_shock_from_velocity =
R"doc(Incident and reflected shocks in a shock tube, from the incident shock
speed.

Parameter ``velocity``:
    incident shock speed relative to the gas at rest ahead of it [m/s])doc";

static const char *__doc_Goddard_ShockSolver_reflected_shock_from_velocity_2 =
R"doc(Incident and reflected shocks in a shock tube, from the incident shock
speed, each with its own chemistry.

Parameter ``velocity``:
    incident shock speed relative to the gas at rest ahead of it [m/s]

Parameter ``incident_chemistry``:
    FROZEN or EQUILIBRIUM, for state 2

Parameter ``reflected_chemistry``:
    FROZEN or EQUILIBRIUM, for state 5)doc";

static const char *__doc_Goddard_ShockSolver_reset = R"doc(Restore the pre-shock state and the construction chemistry. */)doc";

static const char *__doc_Goddard_ShockSolver_store_post_shock_state =
R"doc(Record the current state as the post-shock state and restore the
construction chemistry. */)doc";

static const char *__doc_Goddard_SolverOptions = R"doc()doc";

static const char *__doc_Goddard_SolverOptions_abstol = R"doc()doc";

static const char *__doc_Goddard_SolverOptions_max_iterations = R"doc()doc";

static const char *__doc_Goddard_StationType = R"doc()doc";

static const char *__doc_Goddard_StationType_CHAMBER = R"doc(< Chamber state. For a finite-area combustor, the injector face.)doc";

static const char *__doc_Goddard_StationType_COMBUSTION_END = R"doc(< Finite-area combustor only: end of the constant-area chamber.)doc";

static const char *__doc_Goddard_StationType_EXIT = R"doc()doc";

static const char *__doc_Goddard_StationType_STAGNATION = R"doc(< Finite-area combustor only: stagnation state "inf" at P_inf.)doc";

static const char *__doc_Goddard_StationType_THROAT = R"doc()doc";

static const char *__doc_Goddard_ThermoArray =
R"doc(Wrapper around `Cantera::SolutionArray` with a higher-level API for
broadcasting thermodynamic operations. Supports up to 3D arrays.

If the third dimension is used, it is assumed to be used to set a
composition.)doc";

static const char *__doc_Goddard_ThermoArray_HP = R"doc(Set the enthalpy and pressure of the array.)doc";

static const char *__doc_Goddard_ThermoArray_HPX =
R"doc(Set the enthalpy, pressure, and mole fractions (HPX) or mass fractions
(HPY) of the array. The composition matrix columns should represent
species and the rows should represent distinct compositions.)doc";

static const char *__doc_Goddard_ThermoArray_HPY = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_SH = R"doc(Set the entropy and the enthalpy of the array.)doc";

static const char *__doc_Goddard_ThermoArray_SP = R"doc(Set the entropy and pressure of the array.)doc";

static const char *__doc_Goddard_ThermoArray_SPX =
R"doc(Set the entropy, pressure, and mole fractions (SPX) or mass fractions
(SPY) of the array. The composition matrix columns should represent
species and the rows should represent distinct compositions.)doc";

static const char *__doc_Goddard_ThermoArray_SPY = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_TD = R"doc(Set the temperature and density of the array.)doc";

static const char *__doc_Goddard_ThermoArray_TP = R"doc(Set the temperature and pressure of the array.)doc";

static const char *__doc_Goddard_ThermoArray_TPX =
R"doc(Set the temperature, pressure, and mole fractions of the array. The
mole fraction matrix columns should represent species and the rows
should represent distinct compositions.)doc";

static const char *__doc_Goddard_ThermoArray_TPY = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_TV = R"doc(Set the temperature [K] and specific volume [m^3/kg] of the array.)doc";

static const char *__doc_Goddard_ThermoArray_ThermoArray =
R"doc(NOTE: `ThermoArray` is a straightforward extension of Cantera'
`SolutionArray`. The most efficient solution, in terms of performance
and repeated code, would be to make it a derived class of
`SolutionArray`. Unfortunately, as of 2025-04-19 all of
`SolutionArray`'s constructors are private, which makes it impossible
to inherit from. The only way to construct a `SolutionArray` class is
to call the `create` method.

The alternative for now is to implement `ThermoArray` as a wrapper
class that contains a pointer to a `SolutionArray`.

Entries are stored with the first dimension varying fastest: the entry
at indices (i, j, k) of an array with shape (n0, n1, n2) is at flat
location `i + j*n0 + k*n0*n1` (see `flat_index`).

The array holds its own copy of the `Solution` passed to the
constructor, so changes to that `Solution` made elsewhere (e.g.
through a `Gas` or `Nozzle` sharing it) do not affect the stored
states. Copies of a `ThermoArray` share the same storage.)doc";

static const char *__doc_Goddard_ThermoArray_ThermoArray_2 =
R"doc(Create an array with the given shape. Every entry is initialized to
the current state of `sol`. An empty `shape` leaves the shape unset;
the first setter call then sets it.)doc";

static const char *__doc_Goddard_ThermoArray_ThermoArray_3 =
R"doc(Create an array from a `Gas`, carrying over its candidate condensed
species.

The array stores its own clone of the gas's condensed species set, so
the candidate list and the phase objects are independent of the `Gas`
it was built from. Every entry additionally carries its own amount of
each candidate, initialized to the amounts of `gas`.)doc";

static const char *__doc_Goddard_ThermoArray_UV = R"doc(Set the internal energy and the volume of the array.)doc";

static const char *__doc_Goddard_ThermoArray_check_dimensionality = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_check_ndim = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_condensed_species_names = R"doc(Names of the candidate condensed species, in candidate order. */)doc";

static const char *__doc_Goddard_ThermoArray_enthalpy_mass = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_enthalpy_mole = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_entropy_mass = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_entropy_mole = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_equilibrate =
R"doc(Equilibrate every entry, holding the two properties named by `XY`
constant.

When the array carries condensed species each location is equilibrated
through a `Gas` built on the array's phase and candidate set, so the
condensed amounts take part and are written back; `rtol` and
`max_steps` are forwarded to the Gibbs solver and the remaining
arguments are ignored.

Throws:
    std::invalid_argument if `solver` is "vcs" and the array carries
    condensed species.)doc";

static const char *__doc_Goddard_ThermoArray_flat_index =
R"doc(Flat storage location of the entry at indices (i, j, k). Indices
beyond the array's number of dimensions must be zero.)doc";

static const char *__doc_Goddard_ThermoArray_get_condensed_moles =
R"doc(Amounts of the candidate condensed species at flat location `loc`
[kmol per kg of mixture], in candidate order. Empty for a gas-only
array.)doc";

static const char *__doc_Goddard_ThermoArray_get_state =
R"doc(State vector of the entry at flat location `loc`: the Cantera state of
the gas phase followed by one entry per candidate condensed species
(kmol per kg of mixture), exactly as `Gas::save_state()` lays it out.)doc";

static const char *__doc_Goddard_ThermoArray_internal_energy_mass = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_internal_energy_mole = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_is_shape_set = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_m_condensed = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_m_condensed_moles = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_m_shape_is_set = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_m_solution = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_m_states = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_mean_molecular_weight = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_ndim = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_num_condensed =
R"doc(Number of candidate condensed species carried by the array. Zero for a
gas-only array. */)doc";

static const char *__doc_Goddard_ThermoArray_pressure = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_reshape = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_retrieve_thermo_data = R"doc(retrieve thermodynamic state values from SolutionArray.)doc";

static const char *__doc_Goddard_ThermoArray_set_condensed_moles =
R"doc(Set the amounts of the candidate condensed species at flat location
`loc` [kmol per kg of mixture].

Throws:
    std::invalid_argument if `moles` does not have one entry per
    candidate.)doc";

static const char *__doc_Goddard_ThermoArray_set_state =
R"doc(Set the entry at flat location `loc` from a state vector of the same
phase.

Accepts either the extended length returned by `get_state()` or the
bare Cantera state length, in which case the condensed amounts of that
entry are set to zero.

Throws:
    std::invalid_argument if the vector has neither length.)doc";

static const char *__doc_Goddard_ThermoArray_shape = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_size = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_temperature =
R"doc(Property getters. The result has shape (n0, 1) for a 1-D array and
(n0, n1) for 2-D and 3-D arrays, where element (i, j) is the entry at
`flat_index(i, j, slice)`. `slice` selects the index along the third
dimension and must be 0 for arrays with fewer than 3 dimensions.

The values are those of the gas phase alone: condensed species are not
included in the energies, entropies or molecular weight.)doc";

static const char *__doc_Goddard_ThermoArray_update_states = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_update_states_2 = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_update_states_3 =
R"doc(Updates thermodynamic state by broadcasting a function `f` with values
in `var1` and `var2`.

This template function is defined in the source file because it is a
private method that is only used within the same source file.)doc";

static const char *__doc_Goddard_ThermoArray_update_states_with_composition = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_update_states_with_composition_2 = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_update_states_with_mass_composition = R"doc()doc";

static const char *__doc_Goddard_ThermoArray_update_states_with_mole_composition = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState = R"doc(Convenience class containing relevant thermodynamic results)doc";

static const char *__doc_Goddard_ThermodynamicState_composition = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_cp = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_density = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_dlV_dlP_T = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_dlV_dlT_P = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_enthalpy = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_entropy = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_gamma_s = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_gas_mass_fraction =
R"doc(Mass fraction of the gas phase in the mixture [-]. 1 with no condensed
phase present. */)doc";

static const char *__doc_Goddard_ThermodynamicState_gibbs = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_internal_energy = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_mixture_molecular_weight =
R"doc(Mixture molecular weight [kg/kmol], CEA's "MW": one kg of mixture
divided by the moles of gas *plus* condensed species it holds. Equal
to `molecular_weight` (CEA's "M" = 1/n, which counts the gas alone)
when no condensed phase is present.)doc";

static const char *__doc_Goddard_ThermodynamicState_molecular_weight = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_pinned_transition =
R"doc(True when the state sits exactly at a condensed phase transition with
both polymorphs present and the chemistry is shifting, so the
equilibrium specific heat is infinite. CEA prints a specific heat of
zero for such a station.)doc";

static const char *__doc_Goddard_ThermodynamicState_pressure = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_speed_of_sound = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_stagnation_enthalpy = R"doc()doc";

static const char *__doc_Goddard_ThermodynamicState_temperature = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_2 = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_H_stagnation = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_P_inlet = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_S_inlet = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_area_ratio = R"doc(Area ratio A/A_t [-]: 1 by definition. */)doc";

static const char *__doc_Goddard_ThroatCondition_dlV_dlP_T = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_dlV_dlT_P = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_gamma_s = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_mach =
R"doc(Mach number at the throat [-]: 1 to within the throat solve's
tolerance. */)doc";

static const char *__doc_Goddard_ThroatCondition_pinned_transition =
R"doc(True when the throat sits exactly at a condensed phase transition, so
its temperature is pinned at the transition temperature [K] and both
polymorphs coexist. The equilibrium specific heat is then infinite and
`gamma_s` is -1 / `dlV_dlP_T`; `dlV_dlT_P` is infinite. Always false
for a mixture without condensed phases.)doc";

static const char *__doc_Goddard_ThroatCondition_speed_of_sound = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_state = R"doc()doc";

static const char *__doc_Goddard_ThroatCondition_thermo = R"doc(Mixture state at the throat; see `NozzleStation::thermo`. */)doc";

static const char *__doc_Goddard_ThroatCondition_velocity = R"doc(Flow velocity at the throat [m/s]. */)doc";

static const char *__doc_Goddard_ThrustCoefficient =
R"doc(Thrust performance obtained by integrating over a solved exit plane.
*/)doc";

static const char *__doc_Goddard_ThrustCoefficient_Cf = R"doc(< Thrust coefficient at the specified ambient pressure.)doc";

static const char *__doc_Goddard_ThrustCoefficient_Cf_vacuum = R"doc(< Vacuum thrust coefficient.)doc";

static const char *__doc_Goddard_ThrustCoefficient_momentum_thrust = R"doc(< Momentum component, normalized by p0 * A_throat.)doc";

static const char *__doc_Goddard_ThrustCoefficient_pressure_thrust = R"doc(< Pressure component, normalized by p0 * A_throat.)doc";

static const char *__doc_Goddard_average_angle = R"doc()doc";

static const char *__doc_Goddard_average_cminus_angle = R"doc()doc";

static const char *__doc_Goddard_average_cplus_angle = R"doc()doc";

static const char *__doc_Goddard_chapman_jouguet_detonation =
R"doc(Chapman-Jouguet detonation in a calorically perfect gas that releases
heat.

The CJ Mach number is M_CJ = sqrt(H) + sqrt(H + 1), with H = (gamma^2
- 1) Q / (2 gamma).

Parameter ``gamma``:
    ratio of specific heats of reactants and products [-]

Parameter ``heat_release``:
    heat release q/(R T1), with R the specific gas constant [-]

Returns:
    jump conditions; invalid if gamma <= 1 or heat_release < 0)doc";

static const char *__doc_Goddard_characteristic_intersection_with_angle =
R"doc(Get the coordinates of a downstream characteristic, based on the
intersection of the characteristics of two upstream parent points and
specified characteristic angles.

Parameter ``p1``:
    Point along C- characteristic

Parameter ``p2``:
    Point along C+ characteristic

Parameter ``angle1``:
    C- characteristic angle

Parameter ``angle2``:
    C+ characteristic angle)doc";

static const char *__doc_Goddard_characteristic_isentropic_PT_from_parent =
R"doc(Set temperature and pressure based on isentropic relations and
thermodynamic state of upstream characteristic node.)doc";

static const char *__doc_Goddard_compute_thrust_coefficient =
R"doc(Compute thrust coefficient by integrating over the MoC exit plane.

Uses the relation rho*V^2 = gamma_s * p * M^2 which holds for all
chemistry types. The integrand at each exit plane point is: f = p *
(gamma_s * M^2 * cos^2(theta) + 1). The integral is normalized by p0
times the throat area built from `result.throat_radius` (2 r_t for
planar flow, pi r_t^2 axisymmetric), so the coefficient does not
depend on the length unit.

Parameter ``result``:
    MoC solution result (must have populated exit_plane with gamma_s)

Parameter ``flow_type``:
    PLANAR or AXISYMMETRIC (determines integration measure)

Parameter ``ambient_pressure_ratio``:
    p_amb / p0 (0 for vacuum)

Throws:
    std::invalid_argument if the exit plane has fewer than 2 points.)doc";

static const char *__doc_Goddard_detonation =
R"doc(Overdriven or under-driven detonation in a calorically perfect gas
that releases heat.

With no heat release the overdriven branch is a normal shock and the
under-driven branch a vanishing wave.

Parameter ``drive_factor``:
    wave speed over the CJ speed [-]

Parameter ``gamma``:
    ratio of specific heats of reactants and products [-]

Parameter ``heat_release``:
    heat release q/(R T1), with R the specific gas constant [-]

Parameter ``branch``:
    root of the Rayleigh line and the Hugoniot

Returns:
    jump conditions; invalid if drive_factor < 1, gamma <= 1 or
    heat_release < 0)doc";

static const char *__doc_Goddard_find_like_characteristic_crossings =
R"doc(Scan a net for same-family characteristic crossings (see
MocCrossings).

Compares every pair of same-family characteristic segments whose
x-ranges overlap, rather than only chains that are neighbours on the
front: the net carries no chain-adjacency relation, and coalescence is
a statement about any two characteristics of a family, not only about
neighbours. Segments sharing an endpoint are excluded -- chains
legitimately meet at reflection points.)doc";

static const char *__doc_Goddard_get_cpR_vector =
R"doc(Standard-state molar heat capacity of each gas species divided by R
[-]. */)doc";

static const char *__doc_Goddard_get_enthalpyRT_vector =
R"doc(Standard-state molar enthalpy of each gas species divided by R*T [-].
*/)doc";

static const char *__doc_Goddard_get_equilibrium_gamma =
R"doc(Isentropic exponent -(d log P / d log V)_s [-] of the gas phase at
equilibrium. */)doc";

static const char *__doc_Goddard_get_equilibrium_gamma_2 =
R"doc(Isentropic exponent -(d log P / d log V)_s [-] of the mixture at
equilibrium. */)doc";

static const char *__doc_Goddard_get_frozen_properties =
R"doc(Frozen-composition expansion properties: gas composition and condensed
amounts are held fixed, so the volume derivatives are +/-1 and the
specific heats are the mixture's frozen ones. Condensed species are
incompressible and so contribute equally to cp and cv.)doc";

static const char *__doc_Goddard_get_mole_vector = R"doc(Amount of each gas species [kmol per kg of gas]. */)doc";

static const char *__doc_Goddard_get_stoichiometric_coeffs =
R"doc(@name Gas-only helpers

These read the gas phase alone. Amounts are per kg of *gas*, so a
mixture carrying condensed products must scale them by
`Gas::gas_mass_fraction()` to reach a per kg of mixture basis. @{

Get matrix of stoichiometric coefficients of the species contained in
the `ThermoPhase` object.

Returns:
    2D Eigen array of stoichiometric coefficients. Rows represent
    species, while columns represent elements. Ordering is the same as
    the data file used to generated the `ThermoPhase` object.)doc";

static const char *__doc_Goddard_get_thermo_equilibrium_derivatives = R"doc()doc";

static const char *__doc_Goddard_get_thermo_equilibrium_derivatives_2 =
R"doc(@name Mixture-aware overloads

These take a `Gas`, so they account for condensed products: amounts
are per kg of mixture and the Gordon & McBride system carries one
extra unknown and one extra row per condensed species that is present.
For a `Gas` without condensed phases they reduce exactly to the
`Cantera::ThermoPhase` overloads above.

At a pinned phase transition (`Gas::at_phase_transition()` with both
polymorphs present) the two polymorphs share one element row, so they
are merged into a single condensed unknown and only the pressure block
is solved. The temperature derivatives are then undefined and are
reported as NaN, `spec_heat_p`, `spec_heat_v` and `dlV_dlT_P` are
infinite, and `gamma_s = -1 / dlV_dlP_T`. @{)doc";

static const char *__doc_Goddard_get_thermo_equilibrium_properties = R"doc()doc";

static const char *__doc_Goddard_get_thermo_equilibrium_properties_2 = R"doc()doc";

static const char *__doc_Goddard_get_thermo_equilibrium_properties_3 = R"doc()doc";

static const char *__doc_Goddard_get_thermo_equilibrium_properties_4 = R"doc()doc";

static const char *__doc_Goddard_ideal_gas_D_to_P =
R"doc(Get the pressure of an ideal gas from its density, temperature, and
molar mass.

Parameter ``D``:
    density in kg/m3

Parameter ``T``:
    temperature in K

Parameter ``molar_mass``:
    molar mass in kg/kmol

Returns:
    Pressure in Pa)doc";

static const char *__doc_Goddard_ideal_gas_P_to_D =
R"doc(Get the density of an ideal gas from its pressure, temperature, and
molar mass.

Parameter ``P``:
    pressure in Pa

Parameter ``T``:
    temperature in K

Parameter ``molar_mass``:
    molar mass in kg/kmol

Returns:
    Density in kg/m3)doc";

static const char *__doc_Goddard_mach_from_prandtl_meyer = R"doc()doc";

static const char *__doc_Goddard_normal_shock =
R"doc(Normal shock in a calorically perfect gas.

Parameter ``mach``:
    upstream Mach number [-]

Parameter ``gamma``:
    ratio of specific heats [-]

Returns:
    jump conditions; invalid if mach < 1 or gamma <= 1)doc";

static const char *__doc_Goddard_oblique_shock_deflection_angle =
R"doc(Flow deflection of an oblique shock with a given wave angle, in a
calorically perfect gas.

Parameter ``mach``:
    upstream Mach number [-]

Parameter ``wave_angle``:
    wave angle [rad]

Parameter ``gamma``:
    ratio of specific heats [-]

Returns:
    deflection angle [rad]; negative if the wave angle is below the
    Mach angle)doc";

static const char *__doc_Goddard_oblique_shock_from_deflection =
R"doc(Oblique shock with a given deflection, in a calorically perfect gas.

Parameter ``mach``:
    upstream Mach number [-]

Parameter ``deflection_angle``:
    flow deflection [rad]

Parameter ``gamma``:
    ratio of specific heats [-]

Parameter ``weak``:
    true for the weak solution, false for the strong one

Returns:
    oblique shock; invalid if the deflection exceeds the maximum
    deflection)doc";

static const char *__doc_Goddard_oblique_shock_from_wave_angle =
R"doc(Oblique shock with a given wave angle, in a calorically perfect gas.

Parameter ``mach``:
    upstream Mach number [-]

Parameter ``wave_angle``:
    wave angle [rad]

Parameter ``gamma``:
    ratio of specific heats [-]

Returns:
    oblique shock; invalid if the wave angle is below the Mach angle)doc";

static const char *__doc_Goddard_oblique_shock_max_deflection =
R"doc(Largest deflection an attached oblique shock can produce, in a
calorically perfect gas.

Parameter ``mach``:
    upstream Mach number [-]

Parameter ``gamma``:
    ratio of specific heats [-]

Returns:
    maximum deflection angle [rad]; NaN if mach < 1 or gamma <= 1)doc";

static const char *__doc_Goddard_oblique_shock_max_deflection_wave_angle =
R"doc(Wave angle at which the deflection of an oblique shock is largest, in
a calorically perfect gas (NACA 1135, eq. 168).

Parameter ``mach``:
    upstream Mach number [-]

Parameter ``gamma``:
    ratio of specific heats [-]

Returns:
    wave angle [rad]; NaN if mach < 1 or gamma <= 1)doc";

static const char *__doc_Goddard_oblique_shock_wave_angle =
R"doc(Weak and strong wave angles of an oblique shock with a given
deflection, in a calorically perfect gas.

Parameter ``mach``:
    upstream Mach number [-]

Parameter ``deflection_angle``:
    flow deflection [rad]

Parameter ``gamma``:
    ratio of specific heats [-]

Returns:
    (weak, strong) wave angles [rad]. Both are NaN if the deflection
    is negative or exceeds the maximum deflection (a detached shock),
    or if mach < 1 or gamma <= 1. At zero deflection they are the Mach
    angle and pi/2.)doc";

static const char *__doc_Goddard_prandtl_meyer = R"doc()doc";

static const char *__doc_Goddard_prandtl_meyer_derivative = R"doc()doc";

static const char *__doc_Goddard_reflected_detonation =
R"doc(Detonation and its reflection from the closed end of a tube, in a
calorically perfect gas that releases heat.

Parameter ``drive_factor``:
    wave speed over the CJ speed [-]

Parameter ``gamma``:
    ratio of specific heats of reactants and products [-]

Parameter ``heat_release``:
    heat release q/(R T1), with R the specific gas constant [-]

Parameter ``branch``:
    root of the Rayleigh line and the Hugoniot

Returns:
    incident detonation and reflected shock; invalid if the detonation
    is invalid)doc";

static const char *__doc_Goddard_reflected_shock =
R"doc(Incident shock and its reflection from the closed end of a shock tube,
in a calorically perfect gas.

Parameter ``mach``:
    incident shock Mach number relative to the gas at rest ahead of it
    [-]

Parameter ``gamma``:
    ratio of specific heats [-]

Returns:
    both jumps; invalid if mach < 1 or gamma <= 1)doc";

static const char *__doc_Goddard_summarize_front_shear =
R"doc(Reduce a pass-diagnostics history to the front-shear summary above.

Parameter ``pass_diagnostics``:
    Per-pass records, in order, from MocResult.

Parameter ``window``:
    Number of passes averaged at each end of the march.

Returns:
    The summary; `valid` is false when fewer than 2*window passes were
    recorded.)doc";

static const char *__doc_Goddard_to_string =
R"doc(Human-readable name for a MocStepLimiter, for log/diagnostic messages.
*/)doc";

static const char *__doc_Goddard_to_string_2 = R"doc(Human-readable name for a MocErrorCode, for log/diagnostic messages.)doc";

static const char *__doc_Goddard_validate_moc_options =
R"doc(Check a MocOptions for self-consistency, throwing
std::invalid_argument on any value outside its usable range.

Called from the MocNozzle constructors and again at the top of
solve(), since MocNozzle::options is public and can be changed after
construction. This is programmer error, not a numerical failure, so it
throws rather than returning a MocFailure -- solve()'s "never throws
for numerical failures" contract is unaffected.

theta_max is mode-dependent and deliberately not checked: it is
legitimately left unset in ANALYSIS mode.

Throws:
    NotImplementedError for MocMode::DESIGN_CENTERLINE, a non-zero
    exit_mach, or a non-default
    geometry.upstream_wall_curvature_radius.)doc";

#if defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

