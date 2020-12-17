# Linear array N constant potential spheres
{% set figs='/figures/qdotSphereArray' %}
{% set figs1='/figures/qdotSphereSingle_' %}
{% set figs2='/figures/qdotSphereArray2_' %}

## Formulation
### Scattered field
The scattered wave function inside and outside of the $p$th sphere, $p=1..N$, is expressed as :
\begin{eqnarray}
  f_p^{(in )}(\bb r_p) &=& \sum_{l=0}^{\infty}a_{pl}j_l(k_pr_p)Y_l^0(\theta_p)\\
  f_p^{(out)}(\bb r_p) &=& \sum_{l=0}^{\infty}b_{pl}\hl(k_0r_p)Y_l^0(\theta_p)\\
  k_p &=& k_0\sqrt{1+\frac{V_p}{E}}\\
  k_0 &=& \sqrt{\frac{2m_eE}{\hbar^2}}\\
\end{eqnarray}

where $E_0$, $k_0$ are the energy and wave number of the incident wave, $k_p$, $V_p$ the wave number and constant potential inside the sphere,
$j_l$ and $\hl$ are the spherical Bessel and Hankel functions of the first kind, $Y_l^0$ are the spherical harmonics with azimuthal order 0.


### Continuity relations

The coefficients $a_{p;l}$, $b_{p;l}$ are found imposing the continuity of the wave function and its gradient at the surface of the $p$th sphere.
Since spherical harmonics are used, the continuity of the derivative of the radial part is sufficient to fulfil this condition :
\begin{eqnarray}
      \Big(\sum_{q=1}^N f_q^{(out)}+f^{(i)}\Big)\big|_{r_p=a_p}
  &=& \Big(f_p^{(in)}\Big)\big|_{r_p=a_p} \\
      \dP_r\Big(\sum_{q=1}^N f_q^{(out)}+f^{(i)}\Big)\big|_{r_p=a_p}
  &=& \dP_r\Big(f_p^{(in)}\Big)\big|_{r_p=a_p} \\
\end{eqnarray}

Using the orthogonality of the spherical harmonics.
solving the following linear system yields the unknown coefficients :
\begin{eqnarray}
  \sum_{l=0}^{\infty} a_{pl}j_l(k_pa_p)
    &-& b_{pl}\hl(k_0a_p)\\
    &-& \sum_{q\neq p}^{N} \sum_{n=0}^{\infty} a_{l0,n0}(k_0d_{pq},\theta_{pq})b_{qn}j_n(k_0a_p)
      = e^{jk_0d_p}c_lj_l(k_0a_p) \\
  k_p\sum_l a_{pl}\dP_{\rho}j_l(k_pa_p)
    &-& k_0\dP_{\rho}b_{pl}\hl(k_0a_p) \\
    &-& k_0 \sum_{q\neq p}^{N} \sum_{n} a_{l0,n0}(k_0d_{pq},\theta_{pq})b_{qn}\dP_{\rho}j_n(k_0a_p)
      = k_0e^{jk_0d_p}c_l\dP_{\rho}j_l(k_0a_p)
\end{eqnarray}

where :
\begin{eqnarray}
  c_l  &=& j^l(2l+1)\sqrt{\frac{4\pi}{2l+1}}\\
  a_{l0;n 0}^{out-in}(kd_{pq},\theta_{pq})
    &=& 4\pi\sum_{q=|l-n|}^{l+n}(-j)^{l-n-q} \hl(kd_{pq})Y_l^0(\theta_{pq})\cc G(l,n,q,0,0,0)
\end{eqnarray}

### Linear system
The linear system can also be written :
\begin{equation}
  \big(\bb P - \bb T \big)\bb A = \bb L
\end{equation}

where $\bb A=(\bb a_{pl}, \bb b_{pl})$ is the unknown vector, $\bb T$ the coupling matrix and :
\begin{eqnarray}
\bb P &=& \left[
  \begin{array}{ccccccc}
  j_0^1      &\bb 0&  0        &-h_0^1   &\bb 0&  0      \\
             & ... &           &         & ... &         \\
    0        &\bb 0&j_M^N      &  0      &\bb 0&-h_M^N   \\
  n_1j_0^{1'}&\bb 0&  0        &-h_0^{1'}&\bb 0& 0       \\
             & ... &           &         & ... &         \\
    0        &\bb 0&n_Nj_M^{N'}&  0      &\bb 0&-h_M^{N'}\\
  \end{array}\right]
~\mbox{and}~~~
\bb L &=&
  \left[\begin{array}{c}
    e^{jk_0d_1}j_0^1\\...\\e^{jk_0d_N}j_M^N\\e^{jk_0d_1}j_0^{N'}\\...\\e^{jk_0d_N}j_M^{N'}
  \end{array}\right]
\end{eqnarray}

where $j_l^p = j_l(k_pa_p)$, $h_l^p = h_l^{(1)}(k_0a_p)$ and $z_l^{p'}=\dP_{\rho}z_l^p(\rho)$.

### Far field scattering
In the far field,
$\hl(kr_p)\approx (-j)^{l+1}\frac{e^{jkrp}}{kr_p}$, $r_p=r-d_p\cos\theta$, $\theta_p=\theta$
so the scattering amplitude from the $pth$ sphere can be written $\frac{e^{jkr}}{kr} f_p^{(out)}(\theta)$ where :

\begin{equation}
  f_p(\theta) = \sum_{l=0}^{\infty} (-j)^{l+1}b_{pl}Y_l^0(\theta)
\end{equation}

The total scattering amplitude is the sum of the contribution from individual spheres is :
\begin{eqnarray}
  f(\theta) &=& \sum_{p=1}^{N} f_p(\theta)e^{-jkd_p\cos\theta} \\
            &=& \sum_{l=1}^{\infty}(-j)^{l+1}Y_l^0(\theta)
              \sum_{p=1}^{N} b_{pl}e^{-jkd_p\cos\theta} \\
\end{eqnarray}

The normalized differential scattering cross section is therefore :
\begin{equation}
  \frac{\sigma(\theta)}{a_p^2} = \frac{|f(\theta)|^2}{\left(ka_p\right)^2}
\end{equation}




<!--
#######################################################################
                  Single sphere
#######################################################################
 -->
## Single qdot sphere scattering

Using the orthogonality of the spherical harmonics the unknown coefficients are :

\begin{eqnarray}
  j_l(k_pa_p) a_{pl} - \hl(k_0 a_p) b_{pl}
    &=& j_l(k_0a_p) c_{l}e^{jkd_p} \\
  n_p\dP_{\rho} j_l(k_pa_p) a_{pl} - \dP_{\rho} \hl(k_0 a_p) b_{pl}
    &=& \dP_{\rho} j_l(k_0a_p) c_{l}e^{jkd_p} \\
\end{eqnarray}

\begin{eqnarray}
  a_{pl} &=& c_{l}\frac{-h_0^{'}j_0 + h_0j_0^{'}}{n_pj_1^{'}h_0 - j_1h_0^{'}} \\
  b_{pl} &=& c_{l}\frac{-n_pj_1^{'}j_0 + j_1j_0^{'}}{n_pj_1^{'}h_0 - j_1h_0^{'}} \\
  c_l &=& j^l(2l+1)\sqrt{\frac{4\pi}{2l+1}}\\
\end{eqnarray}

where :
$j_1=j_l(k_pa_p)$, $j_0=j_l(k_0a_p)$, $h_0=\hl(k_0a_p)$ and
$j_1^{'}=\dP_{k_pr}j_l(k_pa_p)$, $j_0^{'}=\dP_{k_0r}j_l(k_0a_p)$, $h_0^{'}=\dP_{k_0r}\hl(k_0a_p)$.


 . | a) $i$ | b) $s$ | c) $t$  
-- | ------ | ------ | ------
$\psi(\bb r)$  | [![]({{figs}}1_fi.png)]({{figs}}1_fi.png) | [![]({{figs}}1_fs.png)]({{figs}}1_fs.png) | [![]({{figs}}1_ft.png)]({{figs}}1_ft.png)

 . | d) $i$ | e) $s$ | f) $t$
 -- | ------ | ------ | ------
$\dP_r\psi$  | [![]({{figs}}1_dfi.png)]({{figs}}1_dfi.png) | [![]({{figs}}1_dfs.png)]({{figs}}1_dfs.png) | [![]({{figs}}1_dft.png)]({{figs}}1_dft.png)

a,d) Incident, b,e) scattered and c,f) total wave function $\psi(\bb r)$ and radial derivative $\dP_r\psi(\bb r)$ for a single sphere with $ka=3$, $n_{ref}=1.5$.


### Scattering amplitude and cross section
$k_p=1.0001$ | $k_p=1.1000$ | $k_p=1.5000$
------------ | ------------ | ------------
[![]({{figs}}1_n0_fka.svg)]({{figs}}1_n0_fka.svg) | [![]({{figs}}1_n1_fka.svg)]({{figs}}1_n1_fka.svg) | [![]({{figs}}1_n2_fka.svg)]({{figs}}1_n2_fka.svg)

Scattering amplitude for a few normalized radius $ka$ and potential strength $k_p$.

a) $log(\sigma)$ | b) $lin(\sigma)$
---------------- | ----------------
[![]({{figs}}1_log_ska.svg)]({{figs}}1_log_ska.svg) | [![]({{figs}}1_lin_ska.svg)]({{figs}}1_lin_ska.svg)

Normalized scattering cross section in a) log b) lin scale for a few values of $k_p$ as a function of $ka$.


### Small potential limit : Born approximation

If the potential is small $\epsilon=V_0/E\ll 1$, the first Born approximation is valid so that the scattering amplitude (expressed in $A$) are found from the Fourier transform of the potential :
\begin{eqnarray}
  f(\bb q) &=& \frac{2m_e}{4\pi\hbar^2}
    \int V(\bb r)e^{-j\bb q\cdot\bb r}d\bb r \\
           &=& \frac{k_0^2}{2}\epsilon\int_0^{a}
    \int_0^\pi e^{-jqr\cos\theta}r^2dr\sin\theta d\theta\\
  f(q)  &=& k_0^2\frac{\epsilon}{q}\int_0^{a}\sin(qr)rdr\\
        &=& \epsilon \frac{a^3k_0^2}{q^2a^2}
    \Big(-\cos qa + \frac{\sin qa}{qa}\Big)\\
\end{eqnarray}
where the $q=2k_0\sin\frac{\theta}{2}$ is the momentum transfer wave vector and $\theta$ the scattering angle.

a) $abs(f(\bb qa))$ | b) $f(\theta),\epsilon\ll 1$ | c) $f(\theta),\epsilon\approx 1$
------------- | --------------- | --------------
[![]({{figs1}}fqa.png )]({{figs1}}fqa.png) | [![]({{figs1}}fka1.svg)]({{figs1}}fka1.svg) | [![]({{figs1}}fka2.svg)]({{figs1}}fka2.svg)

a) Fourier transform of the potential in normalized qa momentum space. The circles show the transfer momentum vector for different values of the normalized radius $ka$. b,c) Scattering amplitudes for different values of $ka$ for b) negligible potential $\epsilon=V_0/E\ll 1$,  c) moderate potential $\epsilon\approx 1$. The dots are exact calculation and dashed lines use the Born approximation.

### Multislice approach : weak phase approximation

Using multislice in the weak phase approximation, the wave function in the far field $\psi(\bb q)$ is the Fourier transform of the transmission function $T$ :
\begin{eqnarray}
  \psi_{ms}(\bb q) &=& \hat T(q_x,q_y) =
    \cc F_{2D}\Big[e^{j\sigma V_z(x,y)} \Big] \\
  V_z(x,y) &=& \int_{\cc B} V(x,y,z) dz  =
    \cc F^{-1}\big[\hat V(q_x,q_y,0)\big]\\
  V_z(\rho,\phi)&=& 4\pi V_0\sqrt{a^2-\rho^2} ~\mbox{,}~~ \rho\le a
\end{eqnarray}

where $2\sigma V_0=k_0\epsilon$, $\hat T(q_x,q_y,0)$ is the 2D Fourier transform of the transmission function and $V_z$ the projected potential.

The scattering amplitude can be calculated by removing the forward propagation and using the Fourier transform in polar coordinates, the trans
\begin{eqnarray}
  f(\bb q) &=& \cc F_{2D} \Big[1-e^{j\sigma V_z(x,y)}\Big]\\
           &=& 2\pi \int_0^a \Big(1-e^{j\sigma V_z(\rho)}\Big) J_0(q\rho)\rho d\rho
\end{eqnarray}

 where $J_0$ is the zeroth order Bessel function.

 $ka=10$ | $ka=30$ | $ka=60$
 --------| ------- | --------
 [![]({{figs1}}feka0.svg)]({{figs1}}feka0.svg) | [![]({{figs1}}feka1.svg)]({{figs1}}feka1.svg) | [![]({{figs1}}feka2.svg)]({{figs1}}feka2.svg)

### Multi-shell single sphere scattering
For a multi-shell sphere with constant potential, the radial Fourier transform of the potential can be computed thanks to the linearity of the FT using the expression of the single shell sphere.

For a $N$-shell with radius $r_i$, potential strength $\epsilon_i$:
\begin{eqnarray}
  f(q) &=& \frac{k_0^2}{q}\sum_{i=1}^{N}
    \epsilon_i\int_{r_{i-1}}^{r_i} \sin(qr)rdr\\
       &=& \frac{k_0^2}{q}\sum_{i=1}^{N}\epsilon_i
    \int_{0}^{r_i} \sin(qr)rdr - \epsilon_i\int_{0}^{r_{i-1}} \sin(qr)rdr\\
       &=& \sum_{i=1}^{N}
    \Big(f_i(q)-\frac{\epsilon_{i}}{\epsilon_{i-1}}f_{i-1}(q)\Big)\\
       &=& f_N + \sum_{i=1}^{N}
    \frac{\epsilon_i-\epsilon_{i+1}}{\epsilon_i}f_{i}(q)
\end{eqnarray}

$V_0(r)$ | $f_0(\theta)$ | $V_1(r)$ | $f_1(\theta)$
-------- | ------------- | -------- | -------------
[![]({{figs1}}shells_pot0.png)]({{figs1}}shells_pot0.png) |  [![]({{figs1}}shells_fka0.svg)]({{figs1}}shells_fka0.svg) |  [![]({{figs1}}shells_pot1.png)]({{figs1}}shells_pot1.png) |  [![]({{figs1}}shells_fka1.svg)]({{figs1}}shells_fka1.svg)



<!--
#######################################################################
                  Two sphere
#######################################################################
-->

## 2 qdot spheres scattering

$ka=1.00$ | $ka=5.00$ | $ka=15.0$
--------- | --------- | --------
[![]({{figs2}}0_ft.png)]({{figs2}}0_ft.png) | [![]({{figs2}}1_ft.png)]({{figs2}}1_ft.png) | [![]({{figs2}}2_ft.png)]({{figs2}}2_ft.png)
[![]({{figs2}}0_fka.svg)]({{figs2}}0_fka.svg) | [![]({{figs2}}1_fka.svg)]({{figs2}}1_fka.svg) | [![]({{figs2}}2_fka.svg)]({{figs2}}2_fka.svg)




### Approximate solution
{% set figs2a='/figures/qdotSphereArray2approx_' %}

Evolution as a function of distance  $kd$ for $N=5$, $ka=2$, $n_{ref}=1.2$.

$kd=2.50ka$ | $kd=5.00ka$ | $kd=10.0ka$
----------- | ----------- | ----------
[![]({{figs2a}}kd0_fka.svg)]({{figs2a}}kd0_fka.svg) | [![]({{figs2a}}kd1_fka.svg)]({{figs2a}}kd1_fka.svg) |  [![]({{figs2a}}kd2_fka.svg)]({{figs2a}}kd2_fka.svg)


Evolution as a function of number of spheres $N$ for $ka=2$,$kd=4ka$, $n_{ref}=1.01$.

$N=2$ | $N=5$ | $N=10$
----- | ----- | -----
[![]({{figs2a}}N0_fka.svg)]({{figs2a}}N0_fka.svg) | [![]({{figs2a}}N1_fka.svg)]({{figs2a}}N1_fka.svg) |  [![]({{figs2a}}N2_fka.svg)]({{figs2a}}N2_fka.svg)
