# Heat Capacity Ratio


One of the most important values in gas dynamics is the heat capacity ratio $\gamma$, which is defined as the ratio of the isobaric heat capacity to the isochoric heat capacity:
$$
\gamma = \frac{C_{p}}{C_{v}}
$$
Unfortunately, the way $\gamma$ is calculated in Cantera (and most other equilibrium solvers) is incorrect for gas mixtures at chemical equilibrium, and using the standard values leads to significant error when computing thermodynamic expansion or compression processes. 

## A physical intuition for the heat capacity ratio
Due to the heat capacity ratio's fundamental importance in gas dynamics, it's ubiquitous in compressible flow equations. As such, it's worth briefly discussing what the heat capacity ratio actually physically represents. Otherwise, the flow equations will seem like tedious algebraic manipulation with little physical insight. 

The heat capacity ratio is above all linked to isentropic processes, with the well-known relations
$$
\frac{T_{2}}{T_{1}} 
= \left(\frac{P_{2}}{P_{1}}\right)^{\frac{\gamma-1}{\gamma}}
= \left(\frac{V_{1}}{V_{2}}\right)^{\gamma-1}
$$

This is a strange and unintuitive result. An isentropic process, by definition, must have *no heat transfer*. So why do heat capacities show up here?

---

At constant volume, when heat is added to a gas, all of the energy goes into the kinetic motions of the gas: the particle velocity, and the vibrations and gyrations of its various bonds and atoms with respect to one another. The isochoric heat capacity represents the amount of heat required to change the average kinetic velocity of the gas particles by a certain amount; we define this as temperature.
$$
C_{v}= \left(\frac{\partial U}{\partial T}\right)_{V}
$$
At constant pressure, the gas is free to expand and contract. So when heat is added to it, in addition to increasing its kinetic velocity and its bond vibrations, it also expands, performing work against its surroundings ($W = P\Delta V)$. 

$$
C_{p}= \left(\frac{\partial H}{\partial T}\right)_{P}
$$
If the gas is performing work on its surroundings then it will require more energy to increase the temperature than if it were at constant volume. 

The heat capacity ratio is therefore **a measure of how much thermal energy goes towards changing the volume of the gas rather than changing the temperature (internal energy) as heat is added or removed at constant pressure**. 

When viewed through this lens, certain opaque mathematical expressions become more intuitive. For example. the common expression $\frac{\gamma-1}{\gamma}$ that appears in e.g. the [turbine work equation](https://www.grc.nasa.gov/www/k-12/airplane/powtrbth.html) is simply the fraction of energy that goes into volume change. 

So while the name "heat capacity ratio" is correct and descriptive in a literal sense, it's not very descriptive in a physical sense. It would be like calling electrical resistance the "voltage-current ratio" based on Kirchoff's voltage law.

Since compressible flows, almost by definition, involve in some manner the expansion, contraction and compression of a fluid body against its surroundings, it's easy to see why the heat capacity ratio is of relevance here. 

## $C_p$ of reactive flows

The thermodynamics of reacting flows are simply an extension of the reasoning applied in the previous section. When the thermodynamic state (temperature and pressure) of a mixture of reacting chemicals shifts, this also results in a shift in the equilibrium composition of the mixture, as per Le Chatelier's principle. Chemical reactions can be exothermic or endothermic, releasing or absorbing more energy, and they can increase or decrease the amount of particles in the gas (e.g. by combining or splitting molecules). This means that a reacting gas that is heated may expand significantly more or less than a gas with a fixed composition, depending on its specific chemical composition. 

The full heat capacity equation, including the contribution of chemical reactions to the specific heat is given by

\begin{align}
C_{p,e} &= C_{p,f} + C_{p,r} \\
&= \sum_{j=1}^{N_s} n_j C_{p,j}^{\circ} + \sum_{j=1}^{N_g} n_j \frac{H_j^{\circ}}{T} \left( \frac{\partial \log n_j}{\partial \log T}\right)_P + \sum_{j=N_g+1}^{N_s} \frac{H_j^{\circ}}{T} \left( \frac{\partial n_j}{\partial \log T}\right)_P
\end{align}

This, however, leaves open the question of how to actually calculate the thermodynamic derivatives for a given thermodynamic state. 

This involves minimizing the chemical potential of the system, while ensuring that the mass and energy of the system are conserved. This results in a linear system of equations:


\begin{align}
\sum_{i=1}^{N_e} \sum_{j=1}^{N_g} a_{kj} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \log T}\right)_P + \sum_{j=N_g+1}^{N_s} a_{ij} \left( \frac{\partial n_j}{\partial \log T}\right)_P + \sum_{j=1}^{N_g} a_{kj} n_j \left( \frac{\partial \log n}{\partial \log T} \right)_P &= -\sum_{j=1}^{N_g} \frac{a_{kj} n_j H_j^{\circ}}{RT} \;, \quad k=1, \ldots, {N_e} \\
\sum_{i=1}^{N_e} a_{ij} \left( \frac{\partial \pi_i}{\partial \log T}\right)_P &= - \frac{H_j^{\circ}}{RT} \;, \quad j = N_g + 1, \ldots, N_s \\
\sum_{i=1}^{N_e} \sum_{j=1}^{N_g} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \log T}\right)_P &= -\sum_{j=1}^{N_g} \frac{n_j H_j^{\circ}}{RT} \\
\sum_{i=1}^{N_e} \sum_{j=1}^{N_g} a_{kj} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \log P}\right)_T + \sum_{j=N_g + 1}^{N_s} a_{kj} \left( \frac{\partial n_j}{\partial \log P}\right)_T + \sum_{j=1}^{N_g} a_{ij} n_j \left( \frac{\partial \log n}{\partial \log P}\right)_T &= \sum_{j=1}^{N_g} a_{kj} n_j \;, \quad k=1, \ldots, {N_e} \\
\sum_{i=1}^{N_e} a_{ij} \left( \frac{\partial \pi_i}{\partial \log P}\right)_T &= 0 \;, \quad j = N_g + 1, \ldots, N_s \\
\sum_{i=1}^{N_e} \sum_{j=1}^{N_g} a_{ij} n_j \left( \frac{\partial \pi_i}{\partial \log P}\right)_T &= \sum_{j=1}^{N_g} n_j \;,
\end{align}


where the unknowns are the thermodynamic derivatives. 

## $C_v$ of reacting flows
Does the isochoric heat capacity $C_v$ also require a similar treatment? Even if the gas can't expand, surely the heat absorbed or released from chemical reactions could affect its heat capacity? Isn't comparing the reactive heat capacity to the frozen heat capacity an apples-to-oranges comparison?

The answer is no, and the reasoning lies, once again, in our interpretation of $\gamma$. It is a ratio, but we're not interested in the heat capacities *per se*, but rather **how much energy goes towards gas expansion in the isobaric case compared to a reference case with no expansion.** Therefore $C_{v,f}$ is merely a reference value, and to a certain extent it does not matter what it is as long as this reference is kept consistent. Pragmatically, however, $C_{v,f}$ is very easy to calculate, and its use simplifies the equations considerably. 