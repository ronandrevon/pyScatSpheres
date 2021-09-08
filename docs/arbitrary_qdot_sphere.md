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
      a_{\nu,\mu;l,m}^{(out-in)}(d_{pq},\theta_{pq},\phi_{pq})b_{q;\nu\mu}
      = e^{jk_0d_p\zeta_p}c_{lm}j_l(k_0a_p) \\
  k_p a_{p;lm} \dP_{\rho}j_l(k_pa_p)
    &-& k_0b_{p;lm}\dP_{\rho}\hl(k_0a_p) - k_0 \dP_{\rho} j_l(k_0a_p)\sum_{q\neq p}^{N} \\
    &&\sum_{\nu=0}^{\infty}\sum_{\mu=-\nu}^{\mu=\nu}
      a_{\nu,\mu;l,m}^{(out-in)}(d_{pq},\theta_{pq},\phi_{pq})b_{q;\nu\mu}
      = k_0e^{jk_0d_p\zeta_p}c_{lm}\dP_{\rho}j_l(k_0a_p) \\
\end{eqnarray}

where :
\begin{eqnarray}
  c_{lm}  &=& 4\pi j^l Y_l^{m*}(\alpha,\pi/2)\\
  d_{pq}  &=& |\bb d_q-\bb d_p| \\
  \theta_{pq} &=& \arccos\Big(\frac{z_q-z_p}{d_{pq} }\Big)\\
  \phi_{pq}   &=& \arctan\Big(\frac{y_q-y_p}{x_q-x_p}\Big)\\
  \zeta_p     &=& \sin(\Theta_p)\sin(\Phi_p)\sin(\alpha)+
                \cos(\Theta_p)\cos(\alpha)\\
\end{eqnarray}

where $0\le\alpha\le\pi$ is the angle of the incident wave with respect to $z$ in the $(y,z)$ plane and $d_p,\Theta_p,\Phi_p$ the centre of sphere $p$ in the global coordinate system.

The coupling coefficients involve the translational addition theorem coefficients :
\begin{eqnarray}
  a_{\nu,\mu;l,m}^{(out-in)}(d_{pq},\theta_{pq},\phi_{pq})
    = 4\pi\sum_{q=|l-\nu|}^{l+\nu}  
      &&(-j)^{l-\nu-q}h_q^{(1)}(kd_{pq}) Y_q^{m-\mu}(\theta_{pq},\phi_{pq})\\
       &&(-1)^m \cc G(l,\nu,q, m,-\mu,\mu-m))
\end{eqnarray}

where
$\cc G(l_1,l_2,l_3,m_1,m_2,m_3) =
\int_{\Omega}Y_{l_1}^{m_1}Y_{l_2}^{m_2}Y_{l_3}^{m_3}d\Omega$
are the [Gaunt coefficients](https://doc.sagemath.org/html/en/reference/functions/sage/functions/wigner.html) where we have used :
\begin{equation}
  \int_{\Omega}Y_{l}^{m}Y_{\nu}^{\mu*}Y_{q}^{m-\mu~*}d\Omega =
(-1)^{\mu-m-\mu}\int_{\Omega}Y_{l}^{m}Y_{\nu}^{-\mu}Y_{q}^{-m+\mu}d\Omega
\end{equation}

### Alternative expression
We can rewrite the equations by doing
$   h_l^{'}(k_0a_p)eqa_{p;lm} - h_l(k_0a_p)eqb_{p;lm}$ and
$n_pj_l^{'}(k_pa_p)eqa_{p;lm} - j_l(k_pa_p)eqb_{p;lm}$
which gives :
\begin{eqnarray}
  a_{p;lm} &=& u_{p;l}e^{jk_0d_p\zeta_p}c_{lm} + u_{p;l}\sum_{q\neq p}^{N}
      \sum_{\nu=0}^{\infty}\sum_{\mu=-\nu}^{\mu=\nu}
      a_{\nu,\mu;l,m}^{(out-in)}(d_{pq},\theta_{pq},\phi_{pq})b_{q;\nu\mu}\\
  b_{p;lm} &=& v_{p;l}e^{jk_0d_p\zeta_p}c_{lm} + v_{p;l} \sum_{q\neq p}^{N}
      \sum_{\nu=0}^{\infty}\sum_{\mu=-\nu}^{\mu=\nu}
      a_{\nu,\mu;l,m}^{(out-in)}(d_{pq},\theta_{pq},\phi_{pq})b_{q;\nu\mu}\\
\end{eqnarray}

where :  
\begin{eqnarray}
  u_{p;l} &=& \frac{h_l^{'}(k_0a_p)j_l(k_0a_p) - h_l(k_0a_p)j_l(k_0a_p)^{'}}
    {j_l(k_pa_p)h_l^{'}(k_0a_p)-n_pj_l^{'}(k_pad_p)h_l(k_0a_p)} \\
  v_{p;l} &=& \frac{n_pj_l^{'}(k_pa_p)j_l(k_0a_p) - j_l(k_pa_p)j_l^{'}(k_0a_p)}
    {j_l(k_pa_p)h_l^{'}(k_0a_p)-n_pj_l^{'}(k_pa_p)h_l(k_0a_p)} \\
\end{eqnarray}


### Linear system
The linear system can also be written :
\begin{equation}
  \big(\bb I - \bb T \big)\bb A = \bb L
\end{equation}

where $\bb A=(\bb a_{pl}, \bb b_{pl})$ is the unknown vector,

$\bb T$ is the cross-coupling matrix and :
\begin{eqnarray}
\bb T &=& \left[
  \begin{array}{cc}
    \bb 0 & \bb T_u  \\
    \bb 0 & \bb T_v
  \end{array}\right] \\
\bb T_u &=& \left[
  \begin{array}{cccc}
    \bb 0        & ..~\bb T_{u;1p}~.. & ..~\bb T_{u;1q}~.. & \bb T_{u;1N}\\
    \bb T_{u;p1} & ..~\bb 0       ~.. & ..~\bb T_{u;pq}~.. & \bb T_{u;pN}\\
    \bb T_{u;q1} & ..~\bb T_{u;qp}~.. & ..~\bb 0       ~.. & \bb T_{u;qN}\\
    \bb T_{u;N1} & ..~\bb T_{u;Np}~.. & ..~\bb T_{u;Nq}~.. & \bb 0       \\
  \end{array}\right] \\
  \bb T_{u;pq} &=& \bb u_{p;l}\otimes\bb A_{pq}^{(out-in)}\\
  \bb T_{v;pq} &=& \bb v_{p;l}\otimes\bb A_{pq}^{(out-in)}\\
  \bb u_{p;l} &=& \left[ \bb u_{p;0},\bb u_{p;1},\bb u_{p;1},\bb u_{p;1},\bb u_{p;2}~... \right]^{T}\\
  \bb v_{p;l} &=& \left[ \bb v_{p;0},\bb v_{p;1},\bb v_{p;1},\bb v_{p;1},\bb u_{p;2}~... \right]^{T}\\
  \bb A_{pq}^{(out-in)} &=& \big(a^{(out-in)}(k_0d_{pq},\theta_{pq},\phi_{pq})\big)_{l,m}
\end{eqnarray}

and $\bb L$ the incident wave coefficients:
\begin{eqnarray}
\bb L   &=& \left[\bb u_{p;l}\bb L_0,\bb v_{p;l}\bb L_0\right]\\
\bb L_0 &=& \left[\bb c_{lm}e^{jk_0d_1\zeta_1},..~
  \bb c_{lm}e^{jk_0d_p\zeta_p}..~,\bb c_{lm}e^{jk_0d_N\zeta_N}\right]^T\\
\end{eqnarray}

where $\square^T$ denotes transpose and $z_l^{'}=\dP_{\rho}z_l(\rho)$.

### Far field scattering
In the far field,
$\hl(k_0r_p)\approx (-j)^{l+1}\frac{e^{jk_0rp}}{k_0r_p}$, $r_p=r-d_p\cos(\theta)$, $\theta_p=\theta$,$\phi_p=\phi$
so the scattering amplitude from the $pth$ sphere $f_p(\theta-\Theta_p,\phi)$  can
be written :  

<!-- \begin{eqnarray}
  r^2 &=& r_p^2 + d_p^2 + 2r_pd_p\cos\left(\theta-\Theta_p\right)\\
      &\approx& r_p^2\big(1+2d_p/r_p\cos\left(\theta-\Theta_p\right)\big)\\
  r   &\approx& r_p + d_p\cos\left(\theta-\Theta_p\right)\\
\end{eqnarray} -->

\begin{equation}
  f_p(\theta,\phi) = \sum_{l=0}^{\infty}\sum_{m=-l}^{l} (-j)^{l+1}b_{p;lm}Y_l^m(\theta,\phi)
\end{equation}
where we have used the notation $f_p(r_p,\theta,\phi)=\frac{e^{jk_0r_p}}{k_0r_p}f_p(\theta,\phi)$.

The total scattering amplitude is the sum of the contribution from all individual spheres.
\begin{eqnarray}
  f(\theta,\phi) &=& \sum_{p=1}^{N} f_p(\theta,\phi)e^{-jk_0d_p\cos(\theta-\Theta_p)} \\
                 &=& \sum_{l=1}^{\infty}(-j)^{l+1}\sum_{m=-l}^{l}Y_l^m(\theta,\phi)
                      \sum_{p=1}^{N} b_{p;lm}e^{-jk_0d_p\cos(\theta-\Theta_p)} \\
\end{eqnarray}
where we have used the notation $f(r,\theta,\phi) = \frac{e^{jk_0r}}{k_0r}f(\theta,\phi)$.

The normalized differential scattering cross section is therefore :
\begin{equation}
  \frac{\sigma(\theta,\phi)}{\pi a_p^2} = \frac{4|f(\theta,\phi)|^2}{\left(k_0a_p\right)^2}
\end{equation}

where we have used the definition $\sigma=4\pi r^2 \Bigg |\frac{f(r,\theta,\phi)|}{f^{(i)}(r,\theta,\phi)}\Bigg|^2$.






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
