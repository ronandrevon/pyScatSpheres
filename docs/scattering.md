{% set figs= '/figures/' %}

# Scattering

## The scalar wave equation

The [wave equation](https://en.wikipedia.org/wiki/Wave_equation) is commonly found in various areas of physics such as [electromagnetics](https://en.wikipedia.org/wiki/Electromagnetic_wave_equation), electric transmission lines, [acoustics](https://en.wikipedia.org/wiki/Acoustic_wave_equation#Derivation), solid mechanics.

The scalar wave equation is established in linear isotropic medium as :

\begin{equation}
  \Big( \grad^2 -\frac{1}{c(\bb r,t)^2}\dP_t^2 \Big) u = 0
\end{equation}

where $u$ is a scalar physical quantity that propagates such as the pressure for acoustic waves, the electric potential in electromagnetics, the strain in solid mechanics. The quantity $c$ is the speed of the wave which depends on the medium.

This equation can be written in harmonic form using $u(\bb r,t)=u(\bb r)e^{j\omega t}$ (further assuming a homogeneous medium) and becomes known as the Helmholtz equation  :  

\begin{equation}
  \Big( \grad^2 + k^2 \Big) u = 0
\end{equation}
where $k=\omega/c(\omega)$ is known as the wave number. Using $\omega=2\pi/T$ and writing $k=2\pi/\lambda$ this relationship becomes $\lambda=cT$. This is commonly known as the dispersion relation of the wave which relates its period $T$ to its wavelength $\lambda$ through the speed of the wave.

In 1D the family of solutions to this equation are  $e^{\pm jk z}$. As a result the solutions to the time dependent equation are $e^{j(\omega t\mp kz)}$ which indeed represent forward and backward propagating waves if $k$ is real.


## Schroedinger's equation

The time dependent [Schroedinger's equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation) reads :
\begin{equation}
  j\hbar\dP_t\Psi(t,\bb r) = \Big(-\frac{\hbar^2}{2m}\grad^2 + V(t,\bb r)\Big)\Psi(t,\bb r)
\end{equation}

where $V(t,\bb r)$ is the potential energy.

The time independent equation is obtained by assuming the time dependent part of the wave function $\Psi(t,\bb r)=\Psi(\bb r)e^{jE/\hbar t}$ where $E$ is the energy of the system.

Injecting above and further assuming time independent potential yields the following eigen value problem :

\begin{equation}
  \Big(-\frac{\hbar^2}{2m}\grad^2 + V(\bb r)\Big)\Psi(\bb r) = E\Psi(\bb r)
\end{equation}

Assuming a constant potential energy $V$ and setting $k^2=2m/\hbar^2\left(E-V\right)$, the above equation becomes the Helmholtz equation.

In the case of electron-atom interaction, the charge of the nucleus is positive concentrated at the centre of the atom. The electron cloud is a distributed negative charge. The electrostatic potential $\varphi$ of the atom is therefore locally positive.
The potential energy $V=-e\varphi$ of an incident electron due to the presence of an atom is negative since the Coulomb force is attractive.

<!-- We use the convention $V_0=e\varphi$ so that :

\begin{equation}
  k^2=2m/\hbar^2\left(E+V_0\right)
\end{equation}

Note that the electron would be in a bound state for energies $E\le-V_0$ and in a bound -->




## Links
- [Schroedinger's equation](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation)
- [wave equation](https://en.wikipedia.org/wiki/Wave_equation)
