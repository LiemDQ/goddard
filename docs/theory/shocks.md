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

The partial derivatives all have exact expressions. 

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
The opposite situation where $\beta$ must be obtained from $\theta$ is more difficult, as we need the shock angle to determine the normal Mach number $M_{n1}$. In this situation, numerical iteration is required. Using product-to-sum identities, we have
$$
\frac{\tan(\beta - \theta)}{\tan\beta} = \frac{\sin(2\beta-\theta) - \sin\theta}{\sin(2\beta-\theta) + \sin\theta}
$$
This is monotonically increasing in $\sin(2\beta-\theta)$, so the maximum value occurs when $\sin(2\beta-\theta) = 1$. Which geometrically implies that the two solutions for $\beta$ must be bracketed by
$$
\beta_{peak} = \frac{\pi}{4}+ \frac{\theta}{2}
$$
As a result, the weak shock is bracketed by $[0, \beta_{peak}]$ while the strong shock is bracketed by $[\beta_{peak}, \pi/2]$. With these brackets, the bisection method can be straightforwardly done.


[^1]: Gordon, Sanford, and Bonnie J. McBride. 1994. _Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications. Part 1: Analysis_. [https://ntrs.nasa.gov/citations/19950013764](https://ntrs.nasa.gov/citations/19950013764).