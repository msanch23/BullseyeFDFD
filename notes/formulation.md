# Cylindrical FDFD Formulation
The most general form of Maxwell's equations are the time-domain integral form:
$$\begin{equation}
\begin{array}{cc}
\displaystyle \oiint_\mathrm{S}\vec{D}(t)\cdot d\vec{s} = \iiint_\mathrm{V}\rho_\mathrm{v}(t)\,dv
& \displaystyle \oiint_\mathrm{S}\vec{B}(t)\cdot d\vec{s} = 0 \\[1em]
\displaystyle \oint_\mathrm{L}\vec{E}(t)\cdot d\vec{l} = -\iint_\mathrm{S}\frac{\partial\vec{B}(t)}{\partial t}\cdot d\vec{s}
& \displaystyle \oint_\mathrm{L}\vec{H}(t)\cdot d\vec{l} = \iint_\mathrm{S}\left[\vec{J}(t)+\frac{\partial\vec{D}(t)}{\partial t}\right]\cdot d\vec{s}
\end{array}
\end{equation}$$

For Finite-Difference Frequency-Domain (FDFD) simulations, Maxwell's Equations will be handled exclusively in the frequency-domain differential form. These are derived by Fourier transforming Maxwell's equations and applying Stokes' theorem and the divergence theorem.

$$\begin{equation}
\begin{array}{cc} 
\displaystyle \nabla\cdot\vec{D} = \rho_\mathrm{v} 
&
\displaystyle \nabla\cdot\vec{B} = 0
\\[1em]
\displaystyle \nabla\times\vec{E}=-j\omega\vec{B}
&
\displaystyle \nabla\times\vec{H}=\vec{J} + j\omega\vec{D}
\end{array}
\end{equation}$$

Using the constitutive relations ($\vec{D}=\epsilon \vec{E}$, $\vec{B}=\mu \vec{H}$), we can rewrite the equations as:

$$\begin{equation}
\begin{array}{cc} 
\displaystyle \nabla\cdot\left(\epsilon\vec{E}\right) = \rho_\mathrm{v}
&
\displaystyle \nabla\cdot\left(\mu\vec{H}\right) = 0
\\[1em]
\displaystyle \nabla\times\vec{E}=-j\omega\mu\vec{H}
&
\displaystyle \nabla\times\vec{H}=\vec{J} + j\omega\epsilon\vec{E}
\end{array}
\end{equation}$$

Note that $\epsilon=\epsilon_0\epsilon_r$ and $\mu=\mu_0\mu_r$ are complex tensors.

Finally, we normalize all of the functions and parameters so that they are all the same order of magnitude by letting $\vec{\tilde{H}}=-j\eta_0\vec{H}$, where $\eta_0=\sqrt{\mu_0/\epsilon_0}$ is the free-space impedance. We also divide divergence equations by $\epsilon_0$ and $\mu_0$ respectively:

$$\begin{equation}
\begin{array}{cc} 
\displaystyle \nabla\cdot\left(\epsilon_r\vec{E}\right) = \rho_\mathrm{v}/\epsilon_0
&
\displaystyle \nabla\cdot\left(\mu_r\vec{\tilde{H}}\right) = 0
\\[1em]
\displaystyle \nabla\times\vec{E}=k_0\mu_r\vec{\tilde{H}}
&
\displaystyle \nabla\times\vec{\tilde{H}}=-j\eta_0\vec{J} + k_0\epsilon_r\vec{E}\\
& =\vec{\tilde{J}}+k_0\epsilon_r\vec{E}
\end{array}
\end{equation}$$

Note that we made use of the free-space wavenumber, defined as $k_0=\omega\sqrt{\mu_0\epsilon_0}=\omega/c$.

## Cylindrical Coordinate Expansion 
Going forward, we focus primarily on the curl equations.

Recall the curl equation in cylindrical coordinates, $\nabla\times\vec{A}=
\left(\frac{1}{\rho}\frac{\partial A_z}{\partial \phi} - \frac{\partial A_\phi}{\partial z}\right)\hat{\rho} + 
\left(\frac{\partial A_\rho}{\partial z} - \frac{\partial A_z}{\partial \rho}\right)\hat{\phi}+
\frac{1}{\rho}\left(\frac{\partial(\rho A_\phi)}{\partial \rho} - \frac{\partial A_\rho}{\partial \phi}\right)\hat{z}$, then the two curl equations can be expanded like so:

$$\begin{equation}
\begin{array}{cc} 
\displaystyle \frac{1}{\rho}\frac{\partial E_z}{\partial\phi} - \frac{\partial E_\phi}{\partial z} = k_0\mu_{\rho\rho}\tilde{H}_\rho + k_0\mu_{\rho\phi}\tilde{H}_\phi + k_0\mu_{\rho z}\tilde{H}_z
\quad & \quad
\displaystyle \frac{1}{\rho}\frac{\partial \tilde{H}_z}{\partial\phi} - \frac{\partial \tilde{H}_\phi}{\partial z} = \tilde{J}_\rho + k_0\epsilon_{\rho\rho}E_\rho + k_0\epsilon_{\rho\phi}E_\phi + k_0\epsilon_{\rho z}E_z
\\[2em]
\displaystyle \frac{\partial E_\rho}{\partial z} - \frac{\partial E_z}{\partial\rho} = k_0\mu_{\phi\rho}\tilde{H}_\rho + k_0\mu_{\phi\phi}\tilde{H}_\phi + k_0\mu_{\phi z}\tilde{H}_z
\quad & \quad
\displaystyle \frac{\partial \tilde{H}_\rho}{\partial z} - \frac{\partial \tilde{H}_z}{\partial\rho} = k_0\epsilon_{\phi\rho}E_\rho + \tilde{J}_\phi + k_0\epsilon_{\phi\phi}E_\phi + k_0\epsilon_{\phi z}E_z
\\[2em]
\displaystyle \frac{1}{\rho}\left(\frac{\partial \left(\rho E_\phi\right)}{\partial\rho} - \frac{\partial E_\rho}{\partial \phi}\right) = k_0\mu_{z\rho}\tilde{H}_\rho + k_0\mu_{z\phi}\tilde{H}_\phi + k_0\mu_{zz}\tilde{H}_z
\quad & \quad
\displaystyle \frac{1}{\rho}\left(\frac{\partial \left(\rho \tilde{H}_\phi\right)}{\partial\rho} - \frac{\partial \tilde{H}_\rho}{\partial \phi}\right) = k_0\epsilon_{z\rho}E_\rho + k_0\epsilon_{z\phi}E_\phi + \tilde{J}_z + k_0\epsilon_{zz}E_z
\end{array}
\end{equation}$$

This is the final form of Maxwell's equations that we will be using and solving.

## Matrix Notation
We can cast the Maxwell's equations into matrix form to facilitate the application of finite-differences. We do this by sorting the left side by the field component.

For the electric field we have:
$$\begin{equation}\begin{bmatrix}
0 & -\partial^{E}_z & \rho^{-1}\partial^{E}_\phi
\\[1em]
\partial^{E}_z & 0 & -\partial^{E}_\rho
\\[1em]
-\rho^{-1}\partial^{E}_\phi & \rho^{-1}\partial^{E}_\rho\rho & 0
\end{bmatrix}
\begin{bmatrix}
E_\rho\\[1em]E_\phi\\[1em]E_z
\end{bmatrix}
=
k_0\begin{bmatrix}
\mu_{\rho\rho} & \mu_{\rho\phi} & \mu_{\rho z}\\[1em]
\mu_{\phi\rho} & \mu_{\phi\phi} & \mu_{\phi z}\\[1em]
\mu_{z\rho} & \mu_{z\phi} & \mu_{z z}
\end{bmatrix}
\begin{bmatrix}
\tilde{H}_\rho\\[1em]\tilde{H}_\phi\\[1em]\tilde{H}_z
\end{bmatrix}\end{equation}$$

For the magnetic field we have:

$$\begin{equation}\begin{bmatrix}
0 & -\partial^{\tilde{H}}_z & \rho^{-1}\partial^{\tilde{H}}_\phi
\\[1em]
\partial^{\tilde{H}}_z & 0 & -\partial^{\tilde{H}}_\rho
\\[1em]
-\rho^{-1}\partial^{\tilde{H}}_\phi & \rho^{-1}\partial^{\tilde{H}}_\rho\rho & 0
\end{bmatrix}
\begin{bmatrix}
\tilde{H}_\rho\\[1.2em]\tilde{H}_\phi\\[1.2em]\tilde{H}_z
\end{bmatrix}
=
\begin{bmatrix}
\tilde{J}_\rho\\[1em]\tilde{J}_\phi\\[1em]\tilde{J}_z
\end{bmatrix} 
+ 
k_0\begin{bmatrix}
\epsilon_{\rho\rho} & \epsilon_{\rho\phi} & \epsilon_{\rho z}\\[1em]
\epsilon_{\phi\rho} & \epsilon_{\phi\phi} & \epsilon_{\phi z}\\[1em]
\epsilon_{z\rho} & \epsilon_{z\phi} & \epsilon_{z z}
\end{bmatrix}
\begin{bmatrix}
E_\rho\\[1em]E_\phi\\[1em]E_z
\end{bmatrix}\end{equation}$$

So far we have kept the derivation general and explicit. Moving forward we will assume isotropic materials; all non-diagonal elements of $\epsilon$ and $\mu$ are set to 0.

## Azimuthal Mode Decomposition
If we assume the geometry is axisymmetric, the structure couples only fields of the same azimuthal order.
Therefore we can expand each field as a single azimuthal harmonic: $E(\rho,\phi,z)=E(\rho,z)e^{jm\phi}$, which means $\partial_\phi = jm$, and collapses the simulation domain from 3D to 2D in the $(\rho, z)$ plane.

Here $m$ is the integer azimuthal mode number (the number of oscillations the field undergoes around the axis of symmetry). Its value sets which components survive on the axis $\rho=0$, where regularity requires the $z$-component to scale as $\rho^{|m|}$ and the transverse components as $\rho^{|m\pm1|}$:

$$\begin{array}{ll}
m=0: & E_z\ne0,\quad E_\phi=H_\rho=0 \\[0.4em]
m=\pm1: & E_\phi=H_\rho\ne0,\quad E_z=0 \\[0.4em]
|m|>1: & E_\phi=E_z=H_\rho=0
\end{array}$$

Applying the isotropic assumption and azimuthal mode decomposition, our equations currently look like:
$$\begin{equation}
\begin{bmatrix}
0 & -\partial^{E}_z & \frac{jm}{\rho}
\\[1em]
\partial^{E}_z & 0 & -\partial^{E}_\rho
\\[1em]
-\frac{jm}{\rho} & \rho^{-1}\partial^{E}_\rho\rho & 0
\end{bmatrix}
\begin{bmatrix}
E_\rho\\[1em]E_\phi\\[1em]E_z
\end{bmatrix}
=
k_0\begin{bmatrix}
\mu_{\rho\rho} & 0 & 0\\[1em]
0 & \mu_{\phi\phi}& 0\\[1em]
0 & 0 & \mu_{zz}
\end{bmatrix}
\begin{bmatrix}
\tilde{H}_\rho\\[1em]\tilde{H}_\phi\\[1em]\tilde{H}_z
\end{bmatrix}\end{equation}$$

$$\begin{equation}
\begin{bmatrix}
0 & -\partial^{\tilde{H}}_z & \frac{jm}{\rho}
\\[1em]
\partial^{\tilde{H}}_z & 0 & -\partial^{\tilde{H}}_\rho
\\[1em]
-\frac{jm}{\rho} & \rho^{-1}\partial^{\tilde{H}}_\rho\rho & 0
\end{bmatrix}
\begin{bmatrix}
\tilde{H}_\rho\\[1.2em]\tilde{H}_\phi\\[1.2em]\tilde{H}_z
\end{bmatrix}
=
\begin{bmatrix}
\tilde{J}_\rho\\[1em]\tilde{J}_\phi\\[1em]\tilde{J}_z
\end{bmatrix}
+
k_0\begin{bmatrix}
\epsilon_{\rho\rho} & 0 & 0\\[1em]
0 & \epsilon_{\phi\phi} & 0\\[1em]
0 & 0 & \epsilon_{zz}
\end{bmatrix}
\begin{bmatrix}
E_\rho\\[1em]E_\phi\\[1em]E_z
\end{bmatrix}\end{equation}$$

We show how to solve these decomposed equations in [`system-assembly.md`](system-assembly.md) by casting them into a discrete linear system, which can be solved for either the eigenfrequency or the driven case.