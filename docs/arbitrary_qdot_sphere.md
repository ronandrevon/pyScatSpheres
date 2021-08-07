# Linear array N constant potential spheres
{% set figs='figures/qdotSphereArray' %}

## Formulation
### Scattered field
The scattered wave function inside and outside of the $p$th sphere, $p=1..N$, is expressed as :
\begin{eqnarray}
  f_p^{(in )}(\bb r_p) &=& \sum_{l=0}^{\infty}j_l(k_pr_p)\sum_{m=-l}^{m=l}a_{p;lm}Y_l^m(\theta_p,\phi_p)\\
  f_p^{(out)}(\bb r_p) &=& \sum_{l=0}^{\infty}\hl(k_0r_p)\sum_{m=-l}^{m=l}b_{p;lm}Y_l^m(\theta_p,\phi_p)\\
  k_p &=& k_0\sqrt{1+\frac{V_p}{E}}\\
  k_0 &=& \sqrt{\frac{2m_eE}{\hbar^2}}\\
\end{eqnarray}

where $E_0$, $k_0$ are the energy and wave number of the incident wave, $k_p$, $V_p$ the wave number and constant potential inside the sphere,
$j_l$ and $\hl$ are the spherical Bessel and Hankel functions of the first kind, $Y_l^m$ are the spherical harmonics with azimuthal order $m$.


### Continuity relations

The coefficients $a_{p;lm}$, $b_{p;lm}$ are found imposing the continuity of the wave function and its gradient at the surface of the $p$th sphere.
Since spherical harmonics are used, the continuity of the derivative of the radial part is sufficient to fulfil this condition :
\begin{eqnarray}
      \Big(\sum_{q=1}^N f_q^{(out)}+f^{(i)}\Big)\big|_{r_p=a_p}
  &=& \Big(f_p^{(in)}\Big)\big|_{r_p=a_p} \\
      \dP_{r_p}\Big(\sum_{q=1}^N f_q^{(out)}+f^{(i)}\Big)\big|_{r_p=a_p}
  &=& \dP_{r_p}\Big(f_p^{(in)}\Big)\big|_{r_p=a_p} \\
\end{eqnarray}

Using the translational additions theorem and the orthogonality of the spherical harmonics $\int_{\Omega} Y_{l}^{m}Y_{l'}^{m'*}d\Omega=\delta_{l,l'}\delta_{m,m'}$,
the following linear system yields the unknown coefficients :
\begin{eqnarray}
  a_{p;lm}j_l(k_pa_p)
    &-& b_{p;lm}\hl(k_0a_p)
    - j_l(k_0a_p)\sum_{q\neq p}^{N} \\
      &&\sum_{\nu=0}^{\infty}\sum_{\mu=-\nu}^{\mu=\nu}
      a_{l,m,\nu,\mu}^{(out-in)}(k_0d_{pq},\theta_{pq},\phi_{pq})b_{q;\nu\mu}
      = e^{jk_0d_p\zeta_p}c_{lm}j_l(k_0a_p) \\
  k_p a_{p;lm} \dP_{\rho}j_l(k_pa_p)
    &-& k_0b_{p;lm}\dP_{\rho}\hl(k_0a_p) - k_0 \dP_{\rho} j_l(k_0a_p)\sum_{q\neq p}^{N} \\
    &&\sum_{\nu=0}^{\infty}\sum_{\mu=-\nu}^{\mu=\nu}
      a_{l,m,\nu,\mu}^{(out-in)}(k_0d_{pq},\theta_{pq},\phi_{pq})b_{q;\nu\mu}
      = k_0e^{jk_0d_p\zeta_p}c_{lm}\dP_{\rho}j_l(k_0a_p) \\
\end{eqnarray}

where :
\begin{eqnarray}
  c_{lm}  &=& 4\pi j^l Y_l^{m*}(\alpha,\pi/2)\\
  \zeta_p &=& \sin(\Theta_p)\sin(\Phi_p)\sin(\alpha)+
                \cos(\Theta_p)\cos(\alpha)\\
\end{eqnarray}

where $0\le\alpha\le\pi$ is the angle of the incident wave with respect to $z$ in the $(y,z)$ plane and $d_p,\Theta_p,\Phi_p$ the centre of sphere $p$ in the global coordinate system.

The coupling coefficients involve the translational addition theorem coefficients :
\begin{eqnarray}
  a_{lm;\nu\mu}^{(out-in)}(kd_{pq},\theta_{pq},\phi_{pq})
    = 4\pi\sum_{q=|l-\nu|}^{l+\nu}&&(-j)^{l-\nu-q} \hl(kd_{pq})\\
      && Y_q^{m-\mu}(\theta_{pq},\phi_{pq})\cc G(l,\nu,q, m,-\mu,-m+\mu)
\end{eqnarray}

<!-- ### Linear system
The linear system can also be written :
\begin{equation}
  \big(\bb P - \bb T \big)\bb A = \bb L
\end{equation}

where $\bb A=(\bb a_{pl}, \bb b_{pl})$ is the unknown vector,
$\bb P$ the matrix of each individual uncoupled sphere :
\begin{equation}
\bb P = \left[
  \begin{array}{ccccccc}
  j_0(k_1a_1)       &\bb 0&  0               &-h_0(k_0a_1)    &\bb 0&  0              \\
                    & ... &                  &                & ... &                 \\
    0               &\bb 0&j_M(k_Na_N)       &  0             &\bb 0&-h_M(k_0a_N)     \\
  n_1j_0^{'}(k_1a_1)&\bb 0&  0               &-h_0^{'}(k_0a_1)&\bb 0& 0               \\
                    & ... &                  &                & ... &                 \\
    0               &\bb 0&n_Nj^{'}_M(k_Na_N)&  0             &\bb 0&-h_M^{'}(k_0a_N) \\
  \end{array}\right]
\end{equation}

, $\bb T$ is the cross-coupling matrix and :
\begin{eqnarray}
\bb T &=& \left[
  \begin{array}{cc}
    \bb 0 & \bb T_p     \\
    \bb 0 & \bb T_p^{'}
  \end{array}\right] \\
\bb T_p     &=& \bb j_l(k_0a_p)     \sum_{q\neq p}^{N}\bb A_{pq}^{(out-in)}\\
\bb T_p^{'} &=& \bb j_l^{'}(k_0a_p) \sum_{q\neq p}^{N}\bb A_{pq}^{(out-in)}
\end{eqnarray}
and :
\begin{equation}
\bb A_{pq}^{(out-in)} = \left[
  \begin{array}{cccc}
  0 & .. & a_{M0;00}^{(out-in)}(k_0d_{pq},\theta_{pq}) \\
  .. & 0 & ..\\
  a_{00;M0}^{(out-in)}(k_0d_{pq},\theta_{pq}) & .. & 0 \\
  \end{array}\right]
\end{equation}

and $\bb L$ the incident wave :
\begin{equation}
\bb L =
  \left[\begin{array}{c}
    c_0e^{jk_0d_1}j_0(k_0a_1)\\...\\c_Me^{jk_0d_N}j_M(k_0a_N)\\c_0e^{jk_0d_1}j_0^{'}(k_0a_1)\\...\\c_Me^{jk_0d_N}j_M^{'}(k_0a_N)
  \end{array}\right]
\end{equation}

where $z_l^{'}=\dP_{\rho}z_l(\rho)$.

### Far field scattering
In the far field,
$\hl(k_0r_p)\approx (-j)^{l+1}\frac{e^{jk_0rp}}{k_0r_p}$, $r_p=r-d_p\cos\theta$, $\theta_p=\theta$
so the scattering amplitude from the $pth$ sphere $f_p(\theta)$  can be written :  

\begin{equation}
  f_p(\theta) = \sum_{l=0}^{\infty} (-j)^{l+1}b_{pl}Y_l^0(\theta)
\end{equation}
where we have used the notation $f_p(r_p,\theta)=\frac{e^{jk_0r_p}}{k_0r_p}f_p(\theta)$.

The total scattering amplitude is the sum of the contribution from all individual spheres.
\begin{eqnarray}
  f(\theta) &=& \sum_{p=1}^{N} f_p(\theta)e^{-jk_0d_p\cos\theta} \\
            &=& \sum_{l=1}^{\infty}(-j)^{l+1}Y_l^0(\theta)
              \sum_{p=1}^{N} b_{pl}e^{-jk_0d_p\cos\theta} \\
\end{eqnarray}
where we have used the notation $f(r,\theta) = \frac{e^{jk_0r}}{k_0r}f(\theta)$.

The normalized differential scattering cross section is therefore :
\begin{equation}
  \frac{\sigma(\theta)}{\pi a_p^2} = \frac{4|f(\theta)|^2}{\left(k_0a_p\right)^2}
\end{equation}

where we have used the definition $\sigma=4\pi r^2 \Bigg |\frac{f(r,\theta)|}{f^{(i)}(r,\theta)}\Bigg|^2$. -->






## Single sphere 
### Analytical solution

In the case of a single sphere, the analytical solution would be found as :

\begin{eqnarray}
  j_l(k_pa_p) a_{p;lm} - \hl(k_0 a_p) b_{p;lm}
    &=& j_l(k_0a_p) c_{lm}e^{jk_0 a_d\zeta_p} \\
  n_p\dP_{\rho} j_l(k_pa_p) a_{p;lm} - \dP_{\rho} \hl(k_0 a_p) b_{p;lm}
    &=& \dP_{\rho} j_l(k_0a_p) c_{lm}e^{jk_0d_p\zeta_p} \\
\end{eqnarray}

which would be solved as :
\begin{eqnarray}
  a_{p;lm} &=& c_{lm}e^{jk_0d_p\zeta_p}\frac{
    -h_l^{'}(k_0a_p)j_l(k_0a_p) + h_l(k_0a_p)j_l(k_0a_p)^{'}}{n_pj_l^{'}(k_pad_p)h_l(k_0a_p) - j_l(k_pa_p)h_l^{'}(k_0a_p)} \\
  b_{p;lm} &=& c_{lm}e^{jk_0d_p\zeta_p}\frac{
    -n_pj_l^{'}(k_pa_p)j_l(k_0a_p) + j_l(k_pa_p)j_l^{'}(k_0a_p)}{n_pj_l^{'}(k_pa_p)h_l(k_0a_p) - j_l(k_pa_p)h_l^{'}(k_0a_p)} \\
  c_{lm}  &=& 4\pi j^l Y_l^{m*}(\alpha,\pi/2)\\
  \zeta_p &=& \sin(\Theta_p)\sin(\Phi_p)\sin(\alpha)+
                \cos(\Theta_p)\cos(\alpha)\\    
\end{eqnarray}

where $z_l^{'}=\dP_{\rho}z_l(\rho)$ for $z_l=j_l,h_l^{(1)}$.

Therefore $a_{p;lm}=a_{p;l}Y_l^{m*}(\alpha)/Y_l^{0}(0)$ and
$b_{p;lm}=b_{p;l}Y_l^{m*}(\alpha)/Y_l^{0}(0)$.
<!--
Using $d_p=0$ the continuity at the interface should be $f_p^{(in )}(\bb a_p)=f_p^{(out )}(\bb a_p)+e^{jk(z\cos\alpha+y\sin\alpha)}$

\begin{eqnarray}
  f_p^{(in )}(\bb a_p)
    &=& \sum_{l=0}^{\infty}j_l(k_pa_p)\sum_{m=-l}^{m=l}
      a_{p;lm}Y_l^m(\theta_p,\phi_p)\\
    &=&\sum_{l=0}^{\infty}j_l(k_pa_p) a_{p;l0}\sum_{m=-l}^{m=l}
      \frac{Y_l^{m}(\alpha)}{Y_l^{0}(0)} Y_l^m(\theta_p,\phi_p)\\
\end{eqnarray} -->
