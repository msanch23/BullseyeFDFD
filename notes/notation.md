# Notation

Canonical symbol table for the theory notes. Every symbol is defined once here;
other notes link to this file instead of redefining. Meanings are collated from
the linked notes. Rows marked *(standard)* are conventional EM notation not
spelled out in-text — the gloss is the usual one.

## Coordinates & operators

| Symbol | Meaning | Defined in |
|---|---|---|
| $\rho,\ \phi,\ z$ | Cylindrical coordinates (radial, azimuthal, axial) | [formulation.md](formulation.md) |
| $\hat\rho,\ \hat\phi,\ \hat z$ | Cylindrical unit vectors *(standard)* | [formulation.md](formulation.md) |
| $\nabla\times$ | Curl operator *(standard)* | [formulation.md](formulation.md) |
| $\partial_\phi = jm$ | Azimuthal derivative under mode decomposition | [formulation.md](formulation.md) |
| $\partial^{E}_\rho,\ \partial^{\tilde H}_\rho,\ \dots$ | Discrete derivative operators on the $E$-grid / $\tilde H$-grid | [formulation.md](formulation.md) |

## Fields & sources

| Symbol | Meaning | Defined in |
|---|---|---|
| $\vec E,\ \vec H$ | Electric, magnetic field *(standard)* | [formulation.md](formulation.md) |
| $\vec D,\ \vec B$ | Electric, magnetic flux density *(standard)* | [formulation.md](formulation.md) |
| $\vec J$ | Electric current density *(standard)* | [formulation.md](formulation.md) |
| $\rho_\mathrm{v}$ | Volume charge density *(standard)* | [formulation.md](formulation.md) |
| $\vec{\tilde H} = -j\eta_0\vec H$ | Normalized magnetic field | [formulation.md](formulation.md) |
| $\vec{\tilde J} = -j\eta_0\vec J$ | Normalized current density | [formulation.md](formulation.md) |
| $E_\rho, E_\phi, E_z$ (etc.) | Field components | [formulation.md](formulation.md) |
| $m$ | Integer azimuthal mode number | [formulation.md](formulation.md) |

## Material parameters & constants

| Symbol | Meaning | Defined in |
|---|---|---|
| $\epsilon = \epsilon_0\epsilon_r$ | Permittivity (complex tensor) | [formulation.md](formulation.md) |
| $\mu = \mu_0\mu_r$ | Permeability (complex tensor) | [formulation.md](formulation.md) |
| $\epsilon_r,\ \mu_r$ | Relative permittivity, permeability | [formulation.md](formulation.md) |
| $\epsilon_0,\ \mu_0$ | Free-space permittivity, permeability *(standard)* | [formulation.md](formulation.md) |
| $\epsilon_{\rho\rho},\ \mu_{\rho\rho},\ \dots$ | Material tensor components | [formulation.md](formulation.md) |
| $\eta_0 = \sqrt{\mu_0/\epsilon_0}$ | Free-space impedance | [formulation.md](formulation.md) |
| $k_0 = \omega\sqrt{\mu_0\epsilon_0} = \omega/c$ | Free-space wavenumber | [formulation.md](formulation.md) |
| $\omega$ | Angular frequency *(standard)* | [formulation.md](formulation.md) |
| $c$ | Speed of light in vacuum *(standard)* | [formulation.md](formulation.md) |
| $j$ | Imaginary unit *(standard)* | [formulation.md](formulation.md) |

## Discrete linear system

| Symbol | Meaning | Defined in |
|---|---|---|
| $[\epsilon],\ [\mu]$ | Material tensors as matrices | [system-assembly.md](system-assembly.md) |
| $\mathbf C_e,\ \mathbf C_h$ | Discrete curl operators | [system-assembly.md](system-assembly.md) |
| $\mathbf e,\ \mathbf h$ | Vectorized $E$, $\tilde H$ fields | [system-assembly.md](system-assembly.md) |
| $\mathbf{\tilde j}$ | Vectorized source | [system-assembly.md](system-assembly.md) |
| $\mathbf A = [\epsilon]^{-1}\mathbf C_h[\mu]^{-1}\mathbf C_e$ | System matrix | [system-assembly.md](system-assembly.md) |
| $\mathbf b = k_0[\epsilon]^{-1}\mathbf{\tilde j}$ | Source vector | [system-assembly.md](system-assembly.md) |
| $\mathbf I$ | Identity matrix *(standard)* | [system-assembly.md](system-assembly.md) |

## Derived quantities

| Symbol | Meaning | Defined in |
|---|---|---|
| $\lambda = 2\pi/\mathrm{Re}(k_0)$ | Wavelength | [system-assembly.md](system-assembly.md) |
| $Q = \mathrm{Re}(k_0)/2\lvert\mathrm{Im}(k_0)\rvert$ | Quality factor (modal loss) | [system-assembly.md](system-assembly.md) |
| $V$ | Modal volume | [system-assembly.md](system-assembly.md) |
| $\eta$ | Collection efficiency | [system-assembly.md](system-assembly.md) |

> **Symbol clash:** $\eta$ (collection efficiency) vs $\eta_0$ (free-space impedance) are unrelated despite the shared letter.
