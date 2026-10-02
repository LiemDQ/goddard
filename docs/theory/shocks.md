# Shocks

## Normal shocks

The relationship between the thermodynamic properties before and after a shock wave are governed by the three conversation laws:

\begin{align}
\rho_{1}u_{1} &= \rho_{2}u_{2}\equiv m \quad &\text{(mass)}\\
\rho_{1}u_{1}^{2} + P_{1} &= \rho_{2}u_{2}^{2} + P_{2}\quad &\text{(momentum)} \\
h_{1} + \frac{1}{2}u_{1}^{2} &= h_{2} + \frac{1}{2}u_{2}^{2} \quad &\text{(energy)}\\
\end{align}

### Hugoniot equation

The Hugoniot equation is:
$$
h_2-h_1 = (P_2 - P_1)\left(\frac{1}{2}\right)\left(\frac{1}{\rho_2}+\frac{1}{\rho_1}\right)
$$

Some key findings:

- Enthalpy change equals pressure difference times mean volume.
- Enthalpy change is independent of wave speed and velocity.
- Enthalpy change is independent of equation of state.

For a calorically perfect gas, we can write $h_2-h_1 = -q + c_p(T_2-T_1)$ which allows us to write the Hugoniot equation only in terms of the pressures and densities:
$$
\left(\frac{\gamma}{\gamma-1}\right)\left(\frac{P_{2}}{\rho_{2}} - \frac{P_{1}}{\rho_{1}}\right) - \frac{1}{2}\left(\frac{1}{\rho_{2}} + \frac{1}{\rho_{1}}\right)(P_{2}- P_{1}) = q
$$
where $P$ is the *static* pressure, $q$ is the heat transferred across the jump and $\gamma$ is the heat capacity ratio. 

### Rankine-Hugoniot relation

Combining the Hugoniot equation with the Rayleigh line allows us to obtain a relation between the pressure and density:

$$
\frac{P_2}{P_1} = \frac{\frac{\gamma+1}{\gamma-1} \frac{\rho_2}{\rho_1} - 1}{\frac{\gamma+1}{\gamma-1} - \frac{\rho_2}{\rho_1}}
$$

### Normal shock relations for real gases

In the general case, we cannot assume that $\gamma$ is constant nor that the gas is calorically perfect. Under these circumstances, a numerical solution is required using the conservation equations (mass, momentum, energy), as explained in Gordon & McBride[^1].

Combining these equations gives us the following general expressions:

$$
\frac{P_{2}}{P_{1}}= 1 - \frac{\rho_{1}u_{1}^{2}}{P_{1}}\left(\frac{\rho_{1}}{\rho_{2}}-1\right) = P^{*}
$$

$$
h_{2}= h_{1} + \frac{u_{1}^{2}}{2}\left[1-\left(\frac{\rho_{1}}{\rho_{2}}\right)^{2}\right] = h^{*}
$$

We can solve for the conditions across the shock in terms of the logarithm of the pressure and temperature ratios $(\frac{P_{2}}{P_{1}}, \frac{T_{2}}{T_{1}})$, using the Newton-Raphson method:

$$
\frac{\partial \left(P^{*}- \frac{P_{2}}{P_{1}}\right)}{\partial \ln \frac{P_{2}}{P_{1}}}\Delta \ln \frac{P_{2}}{P_{1}}
+ \frac{\partial \left(P^{*}- \frac{P_{2}}{P_{1}}\right)}{\partial \ln \frac{T_{2}}{T_{1}}}\Delta \ln \frac{T_{2}}{T_{1}}
= \frac{P_{2}}{P_{1}}- P^{*}
$$

$$
\frac{\partial \left(\frac{h^{*}-h_{2}}{R}\right)}{\partial \ln \frac{P_{2}}{P_{1}}}\Delta \ln \frac{P_{2}}{P_{1}}
+ \frac{\partial \left(\frac{h^{*}-h_{2}}{R}\right)}{\partial \ln \frac{T_{2}}{T_{1}}}\Delta \ln \frac{T_{2}}{T_{1}}
= \frac{h^{*}-h_{2}}{R}
$$

The partial derivatives all have exact expressions. With $\rho_1/\rho_2 = \rho_1 V_2$, they need only the volume derivatives $(\partial \ln V/\partial \ln T)_P$ and $(\partial \ln V/\partial \ln P)_T$, the specific heat $c_p$, and

$$
\left(\frac{\partial h}{\partial \ln P}\right)_T = \frac{R T}{\mathcal{M}}\left[1 - \left(\frac{\partial \ln V}{\partial \ln T}\right)_P\right],
$$

which follows from $dh = T\,ds + v\,dP$ and a Maxwell relation. $\mathcal{M}$ is the mixture molar mass.

The same equations hold for frozen and equilibrium chemistry. They differ only in how each Newton iterate is evaluated:

| | frozen | equilibrium |
|---|---|---|
| trial state at $(T_2, P_2)$ | composition held fixed | equilibrium at $(T_2, P_2)$ |
| $(\partial \ln V/\partial \ln T)_P$ | $1$ | $1 + (\partial \ln n/\partial \ln T)_P$ |
| $(\partial \ln V/\partial \ln P)_T$ | $-1$ | $-1 + (\partial \ln n/\partial \ln P)_T$ |
| $c_p$ | frozen | equilibrium |

The perfect-gas relations with the frozen $\gamma_1$ give the initial guess. An equilibrium shock starts from the frozen solution. Steps are limited as in Gordon & McBride.

The pre-shock state 1 is taken as given, and it need not be at equilibrium: in a shock tube it is often an unburned mixture that burns behind the shock. Its sound speed, and so the Mach number $M_1 = u_1/a_1$, is therefore the frozen one, as in CEA. Downstream Mach numbers use the sound speed of the post-shock chemistry, and stagnation pressures come from isentropic deceleration with that chemistry.

An equilibrium shock in an exothermic mixture only exists above the Chapman-Jouguet detonation speed. Below it the Newton iteration fails to converge.

## Reflected shocks

In a shock tube, the incident shock sets the gas behind it (state 2) moving toward the closed end at

$$
u_p = u_1\left(1 - \frac{\rho_1}{\rho_2}\right).
$$

The shock reflected from the end wall brings it back to rest (state 5). Let $W_R$ be the speed of the reflected shock and $\rho_{52} = \rho_5/\rho_2$. In the reflected-shock frame gas 2 enters at $u_p + W_R$ and gas 5 leaves at $W_R$, so mass conservation gives $W_R = u_p/(\rho_{52} - 1)$. Momentum and energy conservation then become

$$
\frac{P_5}{P_2} = 1 + \frac{\mathcal{M}_2 u_p^2}{R T_2}\,\frac{\rho_{52}}{\rho_{52} - 1},
\qquad
h_5 - h_2 = \frac{u_p^2}{2}\,\frac{\rho_{52} + 1}{\rho_{52} - 1},
$$

which are solved with the same Newton method in $\ln(P_5/P_2)$ and $\ln(T_5/T_2)$. The reflected-shock Mach number relative to gas 2 is $M_R = (u_p + W_R)/a_2$. The incident and reflected shocks can use different chemistry, e.g. a frozen incident shock followed by an equilibrium reflected shock.

For a calorically perfect gas, $M_R$ follows from the incident shock Mach number $M_s$:

$$
\frac{M_R}{M_R^2 - 1} = \frac{M_s}{M_s^2 - 1}\sqrt{1 + \frac{2(\gamma - 1)}{(\gamma + 1)^2}\left(M_s^2 - 1\right)\left(\gamma + \frac{1}{M_s^2}\right)}
$$

## Oblique shocks

When an object moves supersonically in a gas, pressure disturbances originating from it form a front traveling away from it. This line of disturbances is called a Mach wave, and the direction of the wave relative to the direction of travel of the object is called the Mach angle $\mu$. 

If the speed of sound is known, the Mach angle is simply determined by the local Mach number $M$ as
$$
\mu = \sin^{-1} \frac{1}{M}
$$
A Mach wave is a limiting case of a more general phenomenon of oblique shocks, which are stronger pressure disturbances that occur when supersonic flow is forced to "turn into itself". Moving through an oblique shock results in a change in the direction of flow. 


### $\theta$-$\beta$-$M$ relation 
From the trigonometric relations, and the normal shock relations for a [[Calorically perfect gas]], we have
$$
\tan \theta = 2 \cot \beta \left[\frac{M_{1}^{2}\sin^{2}\beta - 1}{M_{1}^{2}(\gamma + \cos 2\beta) + 2}\right]
$$
Which is called the $\theta\text{-}\beta\text{-}M$ relation and specifies the turn angle $\theta$ as a unique function of the Mach number and shock angle $\beta$. 

The opposite relation, $\beta$ as a function of $\theta$, also exits, but it is more elaborate and has two solutions:
$$
\tan \beta = \frac{M^{2}- 1 + 2\lambda \cos[(4\pi\delta + \arccos\chi)/3]}{3(1 + \frac{\gamma-1}{2}M^{2})\tan\theta}
$$
Where $\delta = 0$  for the strong shock solution and $\delta = 1$ for the weak shock, and where:
$$
\lambda = \left[(M^{2}-1)^{2}- 3\left(1 + \frac{\gamma-1}{2}M^{2}\right)\left(1 + \frac{\gamma+1}{2}M^{2}\right)\tan^{2}\theta\right]^{\frac{1}{2}}
$$
$$
\chi = \frac{1}{\lambda^{3}} \left[(M^{2}-1)^{3}- 9\left(1 + \frac{\gamma-1}{2}M^{2}\right)\left(1 + \frac{\gamma-1}{2}M^{2} + \frac{\gamma+1}{4}M^{4}\right)\tan^{2}\theta\right]
$$

### Real gases

If the gas is not a perfect gas, then the previous equations are inapplicable. In this case the general equation is:
$$
\frac{\tan(\beta - \theta)}{\tan\beta} = \frac{u_{2}}{u_{1}}
$$

$\frac{u_{2}}{u_{1}}$ can be obtained from solving the normal shock relations. Solving for $\theta$:
$$
\theta = \beta - \arctan\left(\frac{u_{2}}{u_{1}}\tan\beta\right )
$$
For any chemistry, the tangential velocity is conserved and mass conservation gives $u_{n2}/u_{n1} = \rho_1/\rho_2$, so $\frac{u_2}{u_1}$ above is the density ratio of the normal shock at $u_{n1} = u_1 \sin\beta$.

The opposite situation where $\beta$ must be obtained from $\theta$ is more difficult, as we need the shock angle to determine the normal Mach number $M_{n1}$. In this situation, numerical iteration is required. $\theta(\beta)$ is zero at the Mach angle $\beta = \mu_1$ and at $\beta = \pi/2$, and positive with a single maximum $\theta_{max}$ in between. The solver:

1. finds $\beta_{max}$ and $\theta_{max}$ by golden-section search on $(\mu_1, \pi/2)$. A deflection above $\theta_{max}$ has no attached shock (the shock detaches), and the result is invalid;
2. bisects on $[\mu_1, \beta_{max}]$ for the weak shock, or on $[\beta_{max}, \pi/2]$ for the strong shock. The deflection at the bracket ends is known, so no normal shock is solved at $M_{n1} = 1$.

For a calorically perfect gas the maximum deflection has a closed form (NACA 1135, eq. 168):

$$
\sin^2\beta_{max} = \frac{(\gamma+1)M^2 - 4 + \sqrt{(\gamma+1)\left[(\gamma+1)M^4 + 8(\gamma-1)M^2 + 16\right]}}{4\gamma M^2}
$$

Oblique shocks in an exothermic mixture are oblique detonations and are not modelled.

## Limitations

Shocks are gas-phase only: condensed species are rejected, before or behind the shock. Finite-rate (KINETIC) chemistry is not supported.


[^1]: Gordon, Sanford, and Bonnie J. McBride. 1994. _Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. Part 1: Analysis_. [https://ntrs.nasa.gov/citations/19950013764](https://ntrs.nasa.gov/citations/19950013764).