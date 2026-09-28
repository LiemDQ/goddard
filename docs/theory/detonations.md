# Detonations

A detonation is a shock wave followed by a reaction zone that releases heat. Seen from the wave, the unburned gas (state 1) enters at the wave speed $u_1$ and the burned gas (state 2) leaves at $u_2$. Goddard treats the wave as a discontinuity with the products at chemical equilibrium. The jump conditions are those of a [normal shock](shocks.md):

\begin{align}
\rho_{1}u_{1} &= \rho_{2}u_{2} \quad &\text{(mass)}\\
P_{1} + \rho_{1}u_{1}^{2} &= P_{2} + \rho_{2}u_{2}^{2} \quad &\text{(momentum)} \\
h_{1} + \tfrac{1}{2}u_{1}^{2} &= h_{2} + \tfrac{1}{2}u_{2}^{2} \quad &\text{(energy)}
\end{align}

The difference is that $h_1$ includes the chemical energy of the unburned gas. The unburned gas is taken as given, so its sound speed $a_1$ and the wave Mach number $M_1 = u_1/a_1$ are frozen.

## Rayleigh line and Hugoniot

With specific volume $v = 1/\rho$ and mass flux $m = \rho_1u_1$, mass and momentum give the Rayleigh line, and adding energy gives the Hugoniot:

$$
P_2 - P_1 = m^2(v_1 - v_2), \qquad h_2 - h_1 = \tfrac{1}{2}(P_2 - P_1)(v_1 + v_2).
$$

The Hugoniot of the equilibrium products does not pass through state 1: it lies above it by the heat of reaction. How the Rayleigh line meets it depends on the wave speed:

- Above the Chapman-Jouguet (CJ) speed $u_{CJ}$ the line cuts the Hugoniot twice. The upper root is the **overdriven** (strong) detonation, with products subsonic relative to the wave. The lower root is the **under-driven** (weak) detonation, with supersonic products.
- At $u_{CJ}$ the line is tangent to the Hugoniot. Tangency is equivalent to the Jouguet condition $u_2 = a_2$, with $a_2$ the equilibrium sound speed of the products.
- Below $u_{CJ}$ there is no steady solution.

The drive factor $f = u_1/u_{CJ}$ is therefore at least 1 for both branches, and the two branches merge at $f = 1$. A free, self-sustained detonation travels at the CJ speed. An overdriven detonation needs support from behind, e.g. a piston. The under-driven root is a solution of the jump conditions, but a steady ZND wave cannot reach it through a leading shock.

The CJ point is also the minimum of entropy along the Hugoniot.

## Chapman-Jouguet detonations

Goddard follows Gordon & McBride[^1], chapter 8. With $u_2 = a_2$, $a_2^2 = \gamma_s P_2/\rho_2$ and $r = \rho_2/\rho_1$, momentum and energy become

$$
\frac{P_1}{P_2} = 1 + \gamma_s(1 - r), \qquad
h_2 = h_1 + \frac{\gamma_s R T_2}{2\mathcal{M}_2}\left(r^2 - 1\right),
$$

where $\gamma_s$ is the equilibrium isentropic exponent of the products, $\mathcal{M}$ the molar mass and $R$ the universal gas constant. The velocity has been eliminated: the CJ speed follows from the solution as $u_{CJ} = r\,a_2$.

The unknowns are $x_P = \ln(P_2/P_1)$ and $x_T = \ln(T_2/T_1)$, and the products are equilibrated at each trial $(T_2, P_2)$. The residuals are

$$
f_P = \frac{P_1}{P_2} - 1 + \gamma_s(r - 1), \qquad
f_h = \frac{h_2 - h_1}{R} - \frac{\gamma_s T_2}{2\mathcal{M}_2}\left(r^2 - 1\right).
$$

Gordon & McBride use the reciprocal pressure ratio because the Newton iteration has fewer problems in that form. In the Jacobian, $\gamma_s$ is held constant, as in Gordon & McBride. With $\delta_T = (\partial \ln V/\partial \ln T)_P$, $\delta_P = (\partial \ln V/\partial \ln P)_T$ and $c_p$ all equilibrium values, and using $\partial \ln r = -\partial \ln V$ and $R T/\mathcal{M} = PV$:

$$
\frac{\partial f_P}{\partial x_P} = -\frac{P_1}{P_2} - \gamma_s r\,\delta_P, \qquad
\frac{\partial f_P}{\partial x_T} = -\gamma_s r\,\delta_T,
$$

$$
\frac{\partial f_h}{\partial x_T} = \frac{T_2 c_p}{R} + \frac{\gamma_s T_2}{2\mathcal{M}_2}\left(r^2 + 1\right)\delta_T, \qquad
\frac{\partial f_h}{\partial x_P} = \frac{T_2}{\mathcal{M}_2}\left(1 - \delta_T\right) - \frac{\gamma_s T_2}{2\mathcal{M}_2}\left[\left(r^2 - 1\right) - \left(r^2 + 1\right)\delta_P\right].
$$

The step limiter and the convergence test are those of the shock solver.

### Initial estimate

The initial estimate follows Gordon & McBride, section 8.3:

1. Assume $(P_2/P_1)_0 = 15$. Find the flame temperature $\tilde T$ of an equilibrium at that pressure with $h_2 = h_1 + \tfrac34 (R T_1/\mathcal{M}_1)(P_2/P_1)_0$.
2. Take $\mathcal{M}_2$, $\gamma_s$ and $c_p$ at that state and apply three times:

$$
\alpha = \frac{T_1}{T_2}\frac{\mathcal{M}_2}{\mathcal{M}_1}, \qquad
\frac{P_2}{P_1} = \frac{1 + \gamma_s}{2\gamma_s\alpha}\left[1 + \sqrt{1 - \frac{4\gamma_s\alpha}{(1+\gamma_s)^2}}\right], \qquad r = \alpha\frac{P_2}{P_1},
$$

$$
\frac{T_2}{T_1} = \frac{\tilde T}{T_1} - \frac34\frac{R}{\mathcal{M}_1 c_p}\left(\frac{P_2}{P_1}\right)_0 + \frac{R\gamma_s}{2\mathcal{M}_1 c_p}\frac{r^2 - 1}{r}\frac{P_2}{P_1}.
$$

The first is the pressure root of the momentum equation. The second is the energy equation, linearized about $\tilde T$.

A mixture that releases no heat has no detonation, so Goddard first checks that an equilibrium at $(h_1, P_1)$ is hotter than $T_1$. If it is not, the result is invalid.

## Overdriven and under-driven detonations

At a given wave speed the jump conditions are exactly those of an incident shock, with the products at equilibrium. Goddard solves them with the shock solver's Newton iteration. The initial guess chooses the root: a one-gamma model fitted to the CJ solution.

For a calorically perfect gas (see below), the compression $s = 1 - \rho_1/\rho_2$ at Mach number $M$ is $s = \big[(M^2 - 1) \pm \sqrt{(M^2-1)^2 - bM^2}\big]/(aM^2)$. The fit makes $a$ and $b$ reproduce the real CJ point, $M_{CJ} = u_{CJ}/a_1$ with compression $s_{CJ}$:

$$
a = \frac{M_{CJ}^2 - 1}{s_{CJ} M_{CJ}^2}, \qquad b = \frac{(M_{CJ}^2 - 1)^2}{M_{CJ}^2}.
$$

The guess at $M = f M_{CJ}$ takes the + root for the overdriven branch and the − root for the under-driven one. Pressure comes from the Rayleigh line, $P_2/P_1 = 1 + (\rho_1u_1^2/P_1)\,s$, and temperature from the ideal-gas law with the molar mass of the CJ products. The guess is exact at $f = 1$ and separates the branches as $\sqrt{f - 1}$, like the true roots.

After the solve, the Mach number of the products is checked against 1 to confirm the requested branch. Very close to $f = 1$ the two roots are closer than the solver can separate, and the solver throws `ConvergenceError` rather than return the other root. For stoichiometric H2/O2 this happens below about $f = 1 + 10^{-5}$.

Above the CJ speed, an equilibrium normal shock computed by `ShockSolver` is the overdriven detonation. The shock solver reaches it from the frozen shock.

## ZND structure: the von Neumann state

In the ZND model, the detonation is a non-reacting shock followed by the reaction zone. Across the shock the composition is frozen, which gives the von Neumann spike: pressure and density well above the CJ values. Every `DetonationResult` reports this frozen shock at the same wave speed in `von_neumann`.

## Reflected detonations

When a detonation reaches the closed end of a tube, its products move toward the wall at

$$
u_p = u_1\left(1 - \frac{\rho_1}{\rho_2}\right).
$$

A reflected shock brings them back to rest. It is solved as a [reflected shock](shocks.md#reflected-shocks) in the products, with equilibrium on both sides. The initial guess is the perfect-gas reflected shock with the $\gamma_s$ of the products. Its Mach number relative to the products is $M_R = c + \sqrt{c^2 + 1}$, with $c = (\gamma_s + 1)u_p/(4a_2)$.

## Calorically perfect gas

In the one-gamma model, reactants and products share $\gamma$ and molar mass, and the reaction releases heat $q$ per unit mass. With $Q = q/(R_s T_1)$ ($R_s$ the specific gas constant) and $\mathcal{H} = (\gamma^2 - 1)Q/(2\gamma)$, the Rayleigh line and the Hugoniot give a quadratic in the compression:

$$
(\gamma+1)M_1^2 s^2 - 2\left(M_1^2 - 1\right)s + \frac{4\mathcal{H}}{\gamma + 1} = 0,
\qquad
s_\pm = \frac{\left(M_1^2 - 1\right) \pm \sqrt{\left(M_1^2 - 1\right)^2 - 4\mathcal{H}M_1^2}}{(\gamma + 1)M_1^2}.
$$

The CJ Mach number makes the discriminant vanish:

$$
M_{CJ} = \sqrt{\mathcal{H}} + \sqrt{\mathcal{H} + 1}.
$$

The remaining ratios follow from $s$:

$$
\frac{P_2}{P_1} = 1 + \gamma M_1^2 s, \qquad \frac{\rho_2}{\rho_1} = \frac{1}{1 - s}, \qquad \frac{T_2}{T_1} = (1 - s)\frac{P_2}{P_1}, \qquad M_2^2 = M_1^2\,(1 - s)\frac{P_1}{P_2}.
$$

With $Q = 0$ the + root is the normal shock and the − root the vanishing wave.

## Limitations

- Detonations are gas-phase only. Condensed species in the unburned gas or the products are rejected.
- The unburned gas must have FROZEN or EQUILIBRIUM chemistry. Use the free functions for a perfect gas.
- Oblique detonations and detonation polars are not modelled.

[^1]: Gordon, Sanford, and Bonnie J. McBride. 1994. _Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. Part 1: Analysis_. NASA RP-1311. [https://ntrs.nasa.gov/citations/19950013764](https://ntrs.nasa.gov/citations/19950013764).
