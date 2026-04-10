# Numerical Algorithms

## Stagnation pressure of a real gas

For a generic gas with static pressure $P_1$, velocity $u_1$, entropy $S_1$, and enthalpy function $H(S,P)$, the stagnation pressure can be determined from the following:

1. The stagnation pressure is related to the static pressure through an isentropic deceleration of the flow, so $S_{0} = S_1$.
2. Through energy conservation, the stagnation enthalpy must be $H_0 = H_1 + \frac{u_{1}^{2}}{2}$.
3. Define $f(P) \equiv H(S_1, P) - H_0$; the root of $f$ is found at the stagnation pressure $P_0$. 
4. We can use Newton's method to find the root. The iteration is: $P_{k+1} = P_{k} - \frac{f(P_{k})}{(df/dP)_k}$.
5. Expanding the derivative, we obtain $\frac{df}{dP} = \left(\frac{\partial h}{\partial P}\right)_{S} = \frac{1}{\rho}$ where the second equality is obtained from the thermodynamic definition of enthalpy.
6. Therefore, the iteration is: $P_{k+1} = P_{k} - \rho(S,P_{k}) \cdot f(P_{k})$
