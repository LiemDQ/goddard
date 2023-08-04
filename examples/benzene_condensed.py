import cantera as ct

h2o2 = ct.Species.list_from_file("h2o2.yaml")
O2 = h2o2[3]

condensed = ct.Species.list_from_file("nasa_condensed.yaml")
condensed_dict = {species.name:species for species in condensed}

benzene = condensed_dict["C6H6(L)"]

fuel = ct.Solution(thermo='IdealGas', species = [benzene])
oxidizer = ct.Solution(thermo='IdealGas', species = [O2])

# set temperature and pressure
pressure = 2e7 #200 bar
temperature = 293.15 #room temperature inlet
fuel.TP = (temperature, pressure)
oxidizer.TP = (temperature, pressure)

# generate list of all possible products by comparing species elements to reactant elements
nasa_gas = ct.Species.list_from_file("nasa_gas.yaml")
reaction_elements = ['C', 'H', 'O']

product_species = ['CO','CO2','COOH','H', 'H2', 'HO2', 'HCO', 'HCOOH', 'H2O', 'H2O2', 'O', 'OH', 'O2']

# product_list = [species for species in nasa_gas if species.name in product_species]
product_list = [species for species in nasa_gas if all(x in reaction_elements for x in species.composition)]
products = ct.Solution(thermo='IdealGas', species=product_list)

# set up combustion mixture
o_f_ratio = 2.6 #ratio by weight

molar_ratio = o_f_ratio/(oxidizer.mean_molecular_weight / fuel.mean_molecular_weight)
moles_O2 = molar_ratio/(1+molar_ratio)
moles_fuel = 1-moles_O2

mixture = ct.Mixture([(fuel, moles_fuel), (oxidizer, moles_O2), (products, 0)])

# simulate complete combustion
mixture.equilibrate("HP")

#print report
mixture()