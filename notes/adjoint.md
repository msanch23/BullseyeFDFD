# Eigenvalue Adjoint Sensitivity
We want to know how the resonance, $\lambda$, changes with respect to the permittivity distribution, $\epsilon$. i.e. we want $\frac{\partial\lambda}{\partial\epsilon}$.

Recall, that $\lambda = \frac{2\pi}{\mathrm{Re}(k_0)}$. Then,
$$\begin{equation}\begin{aligned}
\frac{\partial\lambda}{\partial\epsilon}=2\pi\frac{\partial}{\partial\epsilon}\mathrm{Re}(k_0)^{-1}&=-\frac{2\pi}{\mathrm{Re}(k_0)^2}\frac{\partial}{\partial\epsilon}\mathrm{Re}(k_0)\\
&=-\frac{2\pi}{\mathrm{Re}(k_0)^2}\mathrm{Re}\left(\frac{\partial k_0}{\partial\epsilon}\right)
\end{aligned}\end{equation}$$

We also want to know how the optical quality factor, $Q$, changes with respect to the permittivity distribution, $\epsilon$. i.e. we want $\frac{\partial Q}{\partial\epsilon}$.

Recall, that $Q = \frac{\mathrm{Re}(k_0)}{2|\mathrm{Im}(k_0)|}$. Then,
$$\begin{equation}\begin{aligned}
\frac{\partial Q}{\partial\epsilon}&=\frac{1}{2}\left[\frac{1}{|\mathrm{Im}(k_0)|}\frac{\partial}{\partial\epsilon}\mathrm{Re}(k_0)-\frac{\mathrm{Re}(k_0)}{\mathrm{Im}(k_0)^2}\left(\frac{\partial}{\partial\epsilon}\mathrm{Im}(k_0)\right)\frac{\partial}{\partial\mathrm{Im}(k_0)}|\mathrm{Im}(k_0)|\right]\\
&=\frac{1}{2}\left[\frac{1}{|\mathrm{Im}(k_0)|}\frac{\partial}{\partial\epsilon}\mathrm{Re}(k_0)-\frac{\mathrm{Re}(k_0)}{\mathrm{Im}(k_0)^2}\left(\frac{\partial}{\partial\epsilon}\mathrm{Im}(k_0)\right)\mathrm{sgn(\mathrm{Im}(k_0))}\right]\\
&=\frac{1}{2}\left[\frac{1}{|\mathrm{Im}(k_0)|}\mathrm{Re}\left(\frac{\partial k_0}{\partial\epsilon}\right)-\frac{\mathrm{Re}(k_0)}{\mathrm{Im}(k_0)^2}\mathrm{Im}\left(\frac{\partial k_0}{\partial\epsilon}\right)\mathrm{sgn(\mathrm{Im}(k_0))}\right]\\
\end{aligned}\end{equation}$$

Notice that both of them depend on $\frac{\partial k_0}{\partial\epsilon}$. 

We can derive that by recalling that $k_0 = \sqrt{k_0^2}$. Then, by letting $k_0=u^{1/2}$, we have:
$$\begin{equation}\begin{aligned}
\frac{\partial k_0}{\partial\epsilon} = \frac{\partial u^{1/2}}{\partial\epsilon} = \frac{u^{-1/2}}{2}\frac{\partial u}{\partial\epsilon}=\frac{1}{2\sqrt{u}}\frac{\partial u}{\partial\epsilon} = \frac{1}{2k_0}\frac{\partial k_0^2}{\partial\epsilon}
\end{aligned}\end{equation}$$

Thus, $\frac{\partial\lambda}{\partial\epsilon}$ and $\frac{\partial Q}{\partial\epsilon}$ depend on $\frac{\partial k_0}{\partial\epsilon}$, which depends on $\frac{\partial k_0^2}{\partial\epsilon}$.

To find $\frac{\partial k_0^2}{\partial\epsilon}$, recall the eigenfrequency formulation $\mathbf{A}\mathbf{e}=k_0^2\mathbf{e}$. Then,
$$\begin{equation}\begin{aligned}
\frac{\partial}{\partial\epsilon}(\mathbf{A}\mathbf{e}) &= \frac{\partial}{\partial\epsilon}(k_0^2\mathbf{e})\\
\left(\frac{\partial\mathbf{A}}{\partial\epsilon}\right)\mathbf{e}+\mathbf{A}\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right) &= \left(\frac{\partial k_0^2}{\partial\epsilon}\right)\mathbf{e}+k_0^2\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right)
\end{aligned}\end{equation}$$

Here we introduce the left eigenvector, $\psi^\dagger$, which satisfies $\psi^\dagger\mathbf{A} = \psi^\dagger k_0^2$ and $\mathbf{A}^\dagger\psi=\bar{k}_0^2\psi$.

Multiplying from the left gives:
$$\begin{equation}\begin{aligned}
\psi^\dagger\left(\frac{\partial\mathbf{A}}{\partial\epsilon}\right)\mathbf{e}+\psi^\dagger\mathbf{A}\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right) &= \psi^\dagger\left(\frac{\partial k_0^2}{\partial\epsilon}\right)\mathbf{e}+\psi^\dagger k_0^2\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right)
\end{aligned}\end{equation}$$

Notice that the second and fourth term cancel out, so we're left with 
$$\begin{equation}\begin{aligned}
\psi^\dagger\left(\frac{\partial\mathbf{A}}{\partial\epsilon}\right)\mathbf{e}&= \psi^\dagger\left(\frac{\partial k_0^2}{\partial\epsilon}\right)\mathbf{e}
\end{aligned}\end{equation}$$

Solving for $\frac{\partial k_0^2}{\partial\epsilon}$:
$$\begin{equation}\begin{aligned}
\frac{\partial k_0^2}{\partial\epsilon}=\frac{\psi^\dagger\left(\frac{\partial\mathbf{A}}{\partial\epsilon}\right)\mathbf{e}}{\psi^\dagger\mathbf{e}}
\end{aligned}\end{equation}$$

Note we have everything needed for $\frac{\partial k_0^2}{\partial\epsilon}$:

1. $\mathbf{e}$ is obtained from the forward eigensolve $\mathbf{A}\mathbf{e}=k_0^2\mathbf{e}$
1. $\psi^\dagger$ is obtained from a second eigensolve $\mathbf{A}^\dagger\psi=\bar{k}_0^2\psi$, then conjugate transposing $\psi$
1. $\frac{\partial \mathbf{A}}{\partial\epsilon}$ is where the geometry enters:

$$\begin{equation}\begin{aligned}
\frac{\partial \mathbf{A}}{\partial\epsilon}&=\frac{\partial}{\partial\epsilon}\left([\epsilon]^{-1}\mathbf{C}_h[\mu]^{-1}\mathbf{C}_e\right)\\
&=\left(\frac{\partial [\epsilon]^{-1}}{\partial\epsilon}\right)\mathbf{C}_h[\mu]^{-1}\mathbf{C}_e\\
&=-[\epsilon]^{-1}\frac{\partial [\epsilon]}{\partial \epsilon}[\epsilon]^{-1} \mathbf{C}_h[\mu]^{-1}\mathbf{C}_e\\
&=-[\epsilon]^{-1}\frac{\partial [\epsilon]}{\partial \epsilon}\mathbf{A}
\end{aligned}\end{equation}$$

Plugging back into our expression for $\frac{\partial k_0^2}{\partial\epsilon}$, we get:
$$\begin{equation}\begin{aligned}
\frac{\partial k_0^2}{\partial\epsilon}&=\frac{\psi^\dagger\left(-[\epsilon]^{-1}\frac{\partial [\epsilon]}{\partial \epsilon}\mathbf{A}\right)\mathbf{e}}{\psi^\dagger\mathbf{e}}\\
&=\frac{-\psi^\dagger[\epsilon]^{-1}\frac{\partial [\epsilon]}{\partial \epsilon}\mathbf{A}\mathbf{e}}{\psi^\dagger\mathbf{e}}
\end{aligned}\end{equation}$$

Note, that since $[\epsilon]$ is a diagonal matrix, $\frac{\partial [\epsilon]}{\partial \epsilon}$ selects each diagonal element $i$ - in our case this is each cell/pixel of the permittivity. Thus, we have:
$$\begin{equation}\begin{aligned}
\frac{\partial k_0^2}{\partial\epsilon_i}=\frac{-\left(\psi^\dagger[\epsilon]^{-1}\right)_i\left(\mathbf{A}\mathbf{e}\right)_i}{\psi^\dagger\mathbf{e}}
\end{aligned}\end{equation}$$

The contribution of each pixel $i$ to the gradient of $k_0^2$ with respect to $\epsilon$.

We can make use of $\mathbf{A}\mathbf{e}=k_0^2\mathbf{e}$ to show explicitly that $\frac{\partial k_0^2}{\partial\epsilon_i}$ is the (weighted) adjoint eigenvector (adjoint electric field) times the forward eigenvector (electric field).
$$\begin{equation}\begin{aligned}
\frac{\partial k_0^2}{\partial\epsilon_i}&=\frac{-\left(\psi^\dagger[\epsilon]^{-1}\right)_i\left(k_0^2\mathbf{e}\right)_i}{\psi^\dagger\mathbf{e}}\\
&=-k_0^2\frac{\left(\psi^\dagger[\epsilon]^{-1}\right)_i\left(\mathbf{e}\right)_i}{\psi^\dagger\mathbf{e}}
\end{aligned}\end{equation}$$

# Eigenvector Adjoint Sensitivity
From the eigenvector (the electric field), $\mathbf{e}$, we want to calculate quantities like modal volume, $V$, flux, far-fields, and collection efficiency, $\eta$. Note that the eigenvector depends on the permittivity distribution, $\epsilon$, so we can write it as $\mathbf{e}(\epsilon)$. Then, a quantity derived from $\mathbf{e}(\epsilon)$ could be written as $g(\mathbf{e}(\epsilon),\epsilon)$ since it depends on both the permittivity distribution, $\epsilon$, and the eigenvector $\mathbf{e}(\epsilon)$.

For the adjoint sensitivity, we want to know how a quantity, $g(\mathbf{e}(\epsilon),\epsilon)$, changes with respect to the permittivity distribution, $\epsilon$. i.e. we want $\frac{d g}{d \epsilon}$. (Function $g$ arguments are left out going forward). Then, we have:
$$\begin{equation}\begin{aligned}
\frac{dg}{d\epsilon}&=\frac{\partial g}{\partial\epsilon} + \frac{\partial g}{\partial \mathbf{e}}\frac{\partial \mathbf{e}}{\partial \epsilon} + \frac{\partial g}{\partial \bar{\mathbf{e}}}\frac{\partial \bar{\mathbf{e}}}{\partial \epsilon} 
\end{aligned}\end{equation}$$

Note the three terms because $\mathbf{e}$ is complex. But because $g$ and $\epsilon$ are real, the expression collapses to:
$$\begin{equation}\begin{aligned}
\frac{dg}{d\epsilon}&=\frac{\partial g}{\partial\epsilon} + 2\mathrm{Re}\left(\frac{\partial g}{\partial \mathbf{e}}\frac{\partial \mathbf{e}}{\partial \epsilon}\right)
\end{aligned}\end{equation}$$

As in the eigenvalue adjoint sensitivity formulation, to find $\frac{\partial g}{\partial \mathbf{e}}\frac{\partial \mathbf{e}}{\partial \epsilon} $, recall the eigenfrequency formulation $\mathbf{A}\mathbf{e}=k_0^2\mathbf{e}$. Then,
$$\begin{equation}\begin{aligned}
\frac{\partial}{\partial\epsilon}(\mathbf{A}\mathbf{e}) &= \frac{\partial}{\partial\epsilon}(k_0^2\mathbf{e})\\
\left(\frac{\partial\mathbf{A}}{\partial\epsilon}\right)\mathbf{e}+\mathbf{A}\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right) &= \left(\frac{\partial k_0^2}{\partial\epsilon}\right)\mathbf{e}+k_0^2\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right)
\end{aligned}\end{equation}$$
Let's group the terms:
$$\begin{equation}\begin{aligned}
\mathbf{A}\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right)-k_0^2\left(\frac{\partial\mathbf{e}}{\partial\epsilon}\right) &= \left(\frac{\partial k_0^2}{\partial\epsilon}\right)\mathbf{e} - \left(\frac{\partial\mathbf{A}}{\partial\epsilon}\right)\mathbf{e}\\
\left[\mathbf{A}-k_0^2\mathbf{I}\right]\frac{\partial\mathbf{e}}{\partial\epsilon}&=\left[\frac{\partial k_0^2}{\partial\epsilon}\mathbf{I}-\frac{\partial\mathbf{A}}{\partial\epsilon}\right]\mathbf{e}
\end{aligned}\end{equation}$$
Notice that we have all the terms on the right side from the eigenvalue adjoint sensitivity!
However, $\left[\mathbf{A}-k_0^2\mathbf{I}\right]$ is singular. This is not an issue because our FOM $g$ is scale-invariant ($\frac{\partial g}{\partial \mathbf{e}}\mathbf{e}=0$).

We can play a nice trick by multiplying by some vector $v^\dagger$, and choosing it such that $v^\dagger\left[\mathbf{A}-k_0^2\mathbf{I}\right]= \frac{\partial g}{\partial \mathbf{e}}$.

Then we have:
$$\begin{equation}\begin{aligned}
v^\dagger\left[\mathbf{A}-k_0^2\mathbf{I}\right]\frac{\partial\mathbf{e}}{\partial\epsilon}&=v^\dagger\left[\frac{\partial k_0^2}{\partial\epsilon}\mathbf{I}-\frac{\partial\mathbf{A}}{\partial\epsilon}\right]\mathbf{e}\\
\frac{\partial g}{\partial \mathbf{e}}\frac{\partial\mathbf{e}}{\partial\epsilon}&=v^\dagger\left[\frac{\partial k_0^2}{\partial\epsilon}\mathbf{I}-\frac{\partial\mathbf{A}}{\partial\epsilon}\right]\mathbf{e}
\end{aligned}\end{equation}$$

So then $\frac{d g}{d \epsilon}$ becomes:
$$\begin{equation}\begin{aligned}
\frac{dg}{d\epsilon}&=\frac{\partial g}{\partial\epsilon} + 2\mathrm{Re}\left(v^\dagger\left[\frac{\partial k_0^2}{\partial\epsilon}\mathbf{I}-\frac{\partial\mathbf{A}}{\partial\epsilon}\right]\mathbf{e}\right)
\end{aligned}\end{equation}$$
Now instead of having to find $\frac{\partial e}{\partial\epsilon}$ (expensive) we just need $v$!

From $v^\dagger\left[\mathbf{A}-k_0^2\mathbf{I}\right]= \frac{\partial g}{\partial \mathbf{e}}$ we can take the conjugate transpose to get $\left[\mathbf{A}-k_0^2\mathbf{I}\right]^\dagger v= \left(\frac{\partial g}{\partial \mathbf{e}}\right)^\dagger$.

We've taken this as far as we can go. To find $v$ we will need an adjoint source $\left(\frac{\partial g}{\partial \mathbf{e}}\right)^\dagger$ for **each** FOM.
## Modal Volume, $V$
We want to find how the modal volume, $V$ changes with respect to the permittivity distribution, $\epsilon$. i.e. we want $\frac{d V}{d\epsilon}$.

From above,
$$\begin{equation}\begin{aligned}
\frac{d V}{d\epsilon}=\frac{\partial V}{\partial\epsilon} + 2\mathrm{Re}\left(v^\dagger\left[\frac{\partial k_0^2}{\partial\epsilon}\mathbf{I}-\frac{\partial\mathbf{A}}{\partial\epsilon}\right]\mathbf{e}\right)
\end{aligned}\end{equation}$$
and
$$\begin{equation}\begin{aligned}
\left[\mathbf{A}-k_0^2\mathbf{I}\right]^\dagger v= \left(\frac{\partial V}{\partial e}\right)^\dagger
\end{aligned}\end{equation}$$

Notice that we only need $\frac{\partial V}{\partial\epsilon}$, and $\left(\frac{\partial V}{\partial e}\right)^\dagger$, both of which depend on our definition of $V$. 

We calculate $V$ as
$$\begin{equation}\begin{aligned}
V=\frac{\int\epsilon|\mathbf{E}|^2\mathrm{dV}}{\mathrm{max}(\epsilon|\mathbf{E}|^2)}
\end{aligned}\end{equation}$$

when discretized we have

$$\begin{equation}\begin{aligned}
V=\frac{\sum_i\left[\sum_{c\in\rho,\phi,z}\epsilon_{c,i}|\mathbf{E}_{c,i}|^2\right]\Delta V_i}{\mathrm{max_i}(\sum_{c\in\rho,\phi,z}\epsilon_{c,i}|\mathbf{E}_{c,i}|^2)}
\end{aligned}\end{equation}$$

If we let

$$\begin{equation}\begin{aligned}
W_i=\sum_{c\in\rho,\phi,z}\epsilon_{c,i}|\mathbf{E}_{c,i}|^2
\end{aligned}\end{equation}$$
$$\begin{equation}\begin{aligned}
U=\sum_i W_i
\Delta V_i\end{aligned}\end{equation}$$
$$\begin{equation}\begin{aligned}
W_{\mathrm{max}}=\mathrm{max}_i\left(\sum_{c\in\rho,\phi,z}\epsilon_{c,i}|\mathbf{E}_{c,i}|^2\right)
\end{aligned}\end{equation}$$

Then

$$\begin{equation}\begin{aligned}
V=\frac{U}{W_{\mathrm{max}}}
\end{aligned}\end{equation}$$

Note that $\mathbf{e} = \mathbf{E}$. So for $\frac{\partial V}{\partial \mathbf{E}}$, we use the quotient rule:

$$\begin{equation}\begin{aligned}
\frac{\partial V}{\partial E_{c,j}}&=\frac{W_\mathrm{max}\frac{\partial U}{\partial \mathbf{E}_{c,j}}-U\frac{\partial W_{\mathrm{max}}}{\partial \mathbf{E}_{c,j}}}{(W_{\mathrm{max}})^2}\\
&=\frac{1}{W_{\mathrm{max}}}\frac{\partial U}{\partial \mathbf{E}_{c,j}}-\frac{U}{(W_{\mathrm{max}})^2}\frac{\partial W_{\mathrm{max}}}{\partial \mathbf{E}_{c,j}}
\end{aligned}\end{equation}$$

For $\frac{\partial U}{\partial \mathbf{E}_{c,j}}$
$$\begin{equation}\begin{aligned}
\frac{\partial U}{\partial \mathbf{E}_{c,j}} &= \frac{\partial}{\partial \mathbf{E}_{c,j}}\left(\sum_i W_i \Delta V_i\right)\\
&=\frac{\partial}{\partial \mathbf{E}_{c,j}}\left(\sum_i \left[\sum_{c\in\rho,\phi,z}\epsilon_{c,i}|\mathbf{E}_{c,i}|^2\right] \Delta V_i\right)\\
&=\sum_i\epsilon_{c,i}\Delta V_i \frac{\partial}{\partial \mathbf{E}_{c,j}}|E_{c,i}|^2\\
&=\sum_i\epsilon_{c,i}\Delta V_i \frac{\partial}{\partial \mathbf{E}_{c,j}}E_{c,i}\bar{E}_{c,i}\\
&=\sum_i\epsilon_{c,i}\Delta V_i \bar{E}_{c,i}\frac{\partial}{\partial \mathbf{E}_{c,j}}E_{c,i}\\
&=\sum_i\epsilon_{c,i}\Delta V_i \bar{E}_{c,i}\delta_{ij}\\
&=\epsilon_{c,j}\Delta V_j \bar{E}_{c,j}\\
\end{aligned}\end{equation}$$

For $\frac{\partial W_{\mathrm{max}}}{\partial \mathbf{E}_{c,j}}$, note that the max is attained at cell $i^*$, then
$$\begin{equation}\begin{aligned}
\frac{\partial W_{\mathrm{max}}}{\partial \mathbf{E}_{c,j}} = \frac{\partial}{\partial \mathbf{E}_{c,j}}\sum_i W_{i}\delta_{i,i^*}
&= \sum_i\left(\frac{\partial}{\partial \mathbf{E}_{c,j}} W_{i}\right) \delta_{i,i^*}\\
&= \sum_i\left(\frac{\partial}{\partial \mathbf{E}_{c,j}} \sum_{c\in\rho,\phi,z}\epsilon_{c,i}|\mathbf{E}_{c,i}|^2\right) \delta_{i,i^*}\\
&= \sum_i\sum_{c\in\rho,\phi,z}\epsilon_{c,i}\left(\frac{\partial}{\partial \mathbf{E}_{c,j}}|\mathbf{E}_{c,i}|^2\right) \delta_{i,i^*}\\
&= \sum_i\sum_{c\in\rho,\phi,z}\epsilon_{c,i} \bar{\mathbf{E}}_{c,i}\delta_{ij} \delta_{i,i^*}\\
&= \epsilon_{c,j} \bar{\mathbf{E}}_{c,j} \delta_{j,i^*}\\
\end{aligned}\end{equation}$$

Thus, our expression for $\frac{\partial V}{\partial E_{c,j}}$ is now
$$\begin{equation}\begin{aligned}
\frac{\partial V}{\partial E_{c,j}} = \frac{\epsilon_{c,j}\Delta V_j \bar{E}_{c,j}}{W_{\mathrm{max}}} - \frac{U\epsilon_{c,j} \bar{\mathbf{E}}_{c,j} \delta_{j,i^*}}{(W_{\mathrm{max}})^2}
\end{aligned}\end{equation}$$

Using $U = V\ W_{\mathrm{max}}$, we have
$$\begin{equation}\begin{aligned}
\frac{\partial V}{\partial E_{c,j}} &= \frac{\epsilon_{c,j} \bar{E}_{c,j}\Delta V_j}{W_{\mathrm{max}}} - \frac{V\epsilon_{c,j}\bar{\mathbf{E}}_{c,j}\delta_{j,i^*}}{W_{\mathrm{max}}}\\
&= \frac{\epsilon_{c,j} \bar{E}_{c,j}}{W_{\mathrm{max}}}\left(\Delta V_j-V\delta_{j,i^*}\right)
\end{aligned}\end{equation}$$

So the adjoint source for $V$ is the conjugated field $\bar{E}_{c,j}$ weighted by $\frac{\epsilon_{c,j} \Delta V_j}{W_{\mathrm{max}}}$ minus a spike at a single cell weighted by $V/W_{\mathrm{max}}$.

Since our definition of $V$ includes $\epsilon$, $\frac{\partial V}{\partial \epsilon}\ne 0$.

We still need $\frac{\partial V}{\partial \epsilon}$. It follows a very similar derivation, and ultimately looks like
$$\begin{equation}\begin{aligned}
\frac{\partial V}{\partial\epsilon_{c,j}}=\frac{|E_{c,j}|^2}{W_{\max}}\big(\Delta V_j-V\delta_{j,i^*}\big)
\end{aligned}\end{equation}$$