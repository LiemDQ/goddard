import numpy as np
import cantera as ct
import pint

def speed_of_sound(gas: ct.Mixture, gamma: float) -> float: 
    return np.sqrt(ct.gas_constant * gas.T * gamma / gas.mean_molecular_weight)

def get_thermo_derivatives(gas):
    '''
    Gets thermo derivatives based on shifting equilibrium. 

    Taken from https://kyleniemeyer.github.io/rocket-propulsion/thermochemistry/cea_cantera.html. Cantera calculates
    thermodynamic properties, such as the specific heat, based on an assumption of frozen (i.e. invariant) composition, whereas
    CEA assumes the equilibrium changes instantaneously. To reconcile the two, the derivatives of     
    '''
    # unknowns for system with no condensed species:
    # dpi_i_dlogT_P (# elements)
    # dlogn_dlogT_P
    # dpi_i_dlogP_T (# elements)
    # dlogn_dlogP_T
    # total unknowns: 2*n_elements + 2

    num_var = 2 * gas.n_elements + 2
    

    coeff_matrix = np.zeros((num_var, num_var))
    right_hand_side = np.zeros(num_var)

    tot_moles = 1.0 / gas.mean_molecular_weight
    # print(f"Mean molecular weight: {gas.mean_molecular_weight}")
    moles = gas.X * tot_moles

    condensed = False

    # indices
    idx_dpi_dlogT_P = 0
    idx_dlogn_dlogT_P = idx_dpi_dlogT_P + gas.n_elements
    idx_dpi_dlogP_T = idx_dlogn_dlogT_P + 1
    idx_dlogn_dlogP_T = idx_dpi_dlogP_T + gas.n_elements

    # construct matrix of elemental stoichiometric coefficients
    stoich_coeffs = np.zeros((gas.n_elements, gas.n_species))
    for i, elem in enumerate(gas.element_names):
        for j, sp in enumerate(gas.species_names):
            stoich_coeffs[i,j] = gas.n_atoms(sp, elem)

    # equations for derivatives with respect to temperature
    # first n_elements equations
    for k in range(gas.n_elements):
        for i in range(gas.n_elements):
            coeff_matrix[k,i] = np.sum(stoich_coeffs[k,:] * stoich_coeffs[i,:] * moles)
        coeff_matrix[k, gas.n_elements] = np.sum(stoich_coeffs[k,:] * moles)
        right_hand_side[k] = -np.sum(stoich_coeffs[k,:] * moles * gas.standard_enthalpies_RT)

    # skip equation relevant to condensed species

    for i in range(gas.n_elements):
        coeff_matrix[gas.n_elements, i] = np.sum(stoich_coeffs[i, :] * moles)
    right_hand_side[gas.n_elements] = -np.sum(moles * gas.standard_enthalpies_RT)

    # equations for derivatives with respect to pressure

    for k in range(gas.n_elements):
        for i in range(gas.n_elements):
            coeff_matrix[gas.n_elements+1+k,gas.n_elements+1+i] = np.sum(stoich_coeffs[k,:] * stoich_coeffs[i,:] * moles)
        coeff_matrix[gas.n_elements+1+k, 2*gas.n_elements+1] = np.sum(stoich_coeffs[k,:] * moles)
        right_hand_side[gas.n_elements+1+k] = np.sum(stoich_coeffs[k,:] * moles)

    for i in range(gas.n_elements):
        coeff_matrix[2*gas.n_elements+1, gas.n_elements+1+i] = np.sum(stoich_coeffs[i, :] * moles)
    right_hand_side[2*gas.n_elements+1] = np.sum(moles)

    # print("DEBUG:")
    # print(f"{coeff_matrix}")
    # print("Moles")
    # print(f"{moles}")
    # print("Enthalpies")
    # print(f"{gas.standard_enthalpies_RT}")
    # print("------------------")
    
    derivs = np.linalg.solve(coeff_matrix, right_hand_side)

    dpi_dlogT_P = derivs[idx_dpi_dlogT_P : idx_dpi_dlogT_P + gas.n_elements]
    dlogn_dlogT_P = derivs[idx_dlogn_dlogT_P]
    dpi_dlogP_T = derivs[idx_dpi_dlogP_T]
    dlogn_dlogP_T = derivs[idx_dlogn_dlogP_T]

    # dpi_dlogP_T is not used
        
    return dpi_dlogT_P, dlogn_dlogT_P, dlogn_dlogP_T

def _get_thermo_properties(gas, dpi_dlogT_P, dlogn_dlogT_P, dlogn_dlogP_T):
    '''Calculates specific heats, volume derivatives, and specific heat ratio.
    
    Based on shifting equilibrium for mixtures.
    '''
    
    tot_moles = 1.0 / gas.mean_molecular_weight
    moles = gas.X * tot_moles
    
    # construct matrix of elemental stoichiometric coefficients
    stoich_coeffs = np.zeros((gas.n_elements, gas.n_species))
    for i, elem in enumerate(gas.element_names):
        for j, sp in enumerate(gas.species_names):
            stoich_coeffs[i,j] = gas.n_atoms(sp, elem)
    
    # print("DEBUG PROPERTIES:")
    # print("Moles")
    # print(f"{moles}")
    # print("Enthalpies")
    # print(f"{gas.standard_enthalpies_RT}")
    # print("Heat capacities")
    # print(f"{gas.standard_cp_R}")
    # print(f"Volume: {gas.v}")
    # print(f"R: {ct.gas_constant}")
    # print("------------------")

    spec_heat_p = ct.gas_constant * (
        np.sum([dpi_dlogT_P[i] * 
                np.sum(stoich_coeffs[i,:] * moles * gas.standard_enthalpies_RT) 
                for i in range(gas.n_elements)
                ]) +
        np.sum(moles * gas.standard_enthalpies_RT) * dlogn_dlogT_P +
        np.sum(moles * gas.standard_cp_R) +
        np.sum(moles * gas.standard_enthalpies_RT**2)
        )
    
    dlogV_dlogT_P = 1 + dlogn_dlogT_P
    dlogV_dlogP_T = -1 + dlogn_dlogP_T
    
    spec_heat_v = (
        spec_heat_p + gas.P * gas.v / gas.T * dlogV_dlogT_P**2 / dlogV_dlogP_T
        )

    gamma = spec_heat_p / spec_heat_v
    gamma_s = -gamma/dlogV_dlogP_T
    print(f"{gamma_s}, {spec_heat_p}")
    return dlogV_dlogT_P, dlogV_dlogP_T, spec_heat_p, gamma_s

def get_thermo_properties(gas):
    derivs = get_thermo_derivatives(gas)
    print(derivs)
    return _get_thermo_properties(gas, *derivs)
