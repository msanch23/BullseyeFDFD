# Eigenfrequency & Driven Simulations

We can now simplify the notation to 
$$\begin{equation}\textbf{C}_e\mathbf{e}=k_0[\mu]\mathbf{h}\end{equation}$$
$$\begin{equation}\textbf{C}_h\mathbf{h}=\mathbf{\tilde{j}}+k_0[\epsilon]\mathbf{e}\end{equation}$$

Solving Equation 10 for $\mathbf{h}$ yields $\mathbf{h}=k^{-1}_0[\mu]^{-1}\mathbf{C}_e\mathbf{e}$

Plugging in $\mathbf{h}$ into Equation 11 gives:
$$\begin{equation}\begin{aligned}
\textbf{C}_h(k^{-1}_0[\mu]^{-1}\mathbf{C}_e\mathbf{e})&=\mathbf{\tilde{j}}+k_0[\epsilon]\mathbf{e}\\
([\epsilon]^{-1}\mathbf{C}_h[\mu]^{-1}\mathbf{C}_e-k^2_0)\mathbf{e}&=k_0[\epsilon]^{-1}\mathbf{\tilde{j}}\\
(\mathbf{A}-k^2_0\mathbf{I})\mathbf{e}&=\mathbf{b}
\end{aligned}\end{equation}$$

## Driven Simulation
For driven simulations, $\mathbf{b}\ne\mathbf{0}$. Then we need to solve:
$$\begin{equation}
\mathbf{e}=(\mathbf{A}-k_0^2\mathbf{I})\mathbin{\backslash}\mathbf{b}
\end{equation}$$
Note that we need three ingredients:

1. $\mathbf{A}$; a matrix constructed above - enforces Maxwell's curl equations.
1. $k_0^2$; the squared free-space wavenumber of the driving source.
1. $\mathbf{b}$; a vector of the driving source - used to implement the source's shape and polarization.

to calculate a single vector $\mathbf{e}$ - the steady-state electric field of the configured driven simulation.

## Eigenfrequency Simulation
For eigenfrequency simulations, $\mathbf{b}=\mathbf{0}$. Then we need to solve:
$$\begin{equation}
\mathbf{A}\mathbf{e}=k_0^2\ \mathbf{e}
\end{equation}$$
Note that we only need one ingredient, $\mathbf{A}$, which already contains the permittivity distribution, $[\epsilon]$, that describes the structure/device. Then we use an eigenvalue/eigenvector solver to find the $k_0^2$ and $\mathbf{e}$ that satisfy the eigenvalue formulation.

The resulting eigenvalue, $k_0^2$, is a complex number. It can be used to extract the wavelength, $\lambda = \frac{2\pi}{\mathrm{Re}(k_0)}$, and the optical quality factor, $Q = \frac{\mathrm{Re}(k_0)}{2|\mathrm{Im}(k_0)|}$.

The resulting eigenvector, $\mathbf{e}$, is the electric field of the device. It can be used to derive quantities like modal volume, $V$, and far-field transformations for calculating the collection efficiency, $\eta$.
