# Linear array N hard spheres
{% set figs='/figures/hardSphereArray' %}

## Scattered field
The scattering amplitude from the $p$th sphere is expressed as :
\begin{equation}
  f_p^{(s)}(\bb r_p) = \sum_{l=0}^{\infty}c_{p;l}\hl(r_p)Y_l^0(\theta_p)
\end{equation}

where $\hl$ is the spherical Hankel function of the first kind, $Y_l^0$ are the spherical harmonics with azimuthal order 0.

The coefficients $c_{p;l}$ are found imposing $\big(\sum_{p} f_p^{(s)}+f^{(i)}\big)\big|_{r_p=a_p}=0$ and using the orthogonality of the spherical harmonics.
This results in solving the system :
\begin{eqnarray}
  c_{p;l} = &&u_l(ka_p)\Bigg[
      e^{jkd_p}b_{l}\sqrt{\frac{4\pi}{2l+1}}\\
      && + \sum_{q=p+1}^{N}\sum_{n=0}^{\infty}a_{l0;n 0}^{out-in}(kd_{pq},\pi) c_{q;n}
         + \sum_{q=1}^{p-1}\sum_{n=0}^{\infty}a_{l0;n 0}^{out-in}(kd_{pq},0  ) c_{q;n}
      \Bigg]
\end{eqnarray}

where $u_l(ka_p)$ are the single sphere continuity coefficients,
$b_l$ are the incident plane wave spherical expansion coefficients and
$a_{l0;n0}^{out-in}(kd_{pq})$ are the out-to-in translational coefficients :
\begin{eqnarray}
  u_l(ka_p) &=& -\frac{j_l(ka_p)}{\hl(ka_p)}\\
  b_l &=& j^l(2l+1)\\
  a_{l0;n 0}^{out-in}(kd_{pq},\theta_{pq}) &=& 4\pi\sum_{q=|l-n|}^{l+n}(-j)^{l-n-q} \hl(kd_{pq})Y_l^0(\theta_{pq})\cc G(l,n,q,0,0,0)
\end{eqnarray}
where $j_l$ is the spherical Bessel function of the first kind and $\cc G$ are the [Gaunt](https://en.wikipedia.org/wiki/3-j_symbol#Relation_to_spherical_harmonics) coefficients.


## Far field scattering
In the far field,
$\hl(kr_p)\approx (-j)^{l+1}\frac{e^{jkrp}}{kr_p}$, $r_p=r-d_p\cos\theta$, $\theta_p=\theta$
so the scattering amplitude from the $pth$ sphere can be written $\frac{e^{jkr}}{kr} f_p(\theta)$ where :

\begin{equation}
  f_p(\theta) = \sum_{l=0}^{\infty} (-j)^{l+1}c_ {p;l}Y_l^0(\theta)
\end{equation}

The total scattering amplitude is the sum of the contribution from individual spheres is :   
\begin{eqnarray}
  f(\theta) &=& \sum_{p=1}^{N} f_p(\theta)e^{-jkd_p\cos\theta} \\
            &=& \sum_{l=1}^{\infty}(-j)^{l+1}Y_l^0(\theta)
              \sum_{p=1}^{N} c_{pl}e^{-jkd_p\cos\theta} \\
\end{eqnarray}


The normalized differential scattering cross section is therefore :
\begin{equation}
  \frac{\sigma(\theta)}{a_p^2} = \frac{|f(\theta)|^2}{\left(ka_p\right)^2}
\end{equation}

The integrated scattering cross section over the unit sphere can be written using the orthonormality of $Y_l^0$ :
\begin{eqnarray}
  \sigma &=&\int_{\cc S} \Big(\sum_p f_p(\theta)e^{-jkd_p\cos\theta}\Big)\Big(\sum_q f_q^{*}(\theta)e^{jkd_q\cos\theta}\Big) d\Omega\\
        &=&\sum_{p=1}^{N}\sum_{l=0}^{\infty}c_{pl}
            \sum_{q=1}^{N}\sum_{n=0}^{\infty}c_{qn}^{*} \int_{\cc S}
            e^{-jkd_{pq}\cos\theta}Y_l^0(\theta) Y_n^{0*}(\theta)d\Omega\\
        &=&\sum_{p=1}^{N}\sigma_p^{(coupled)} +
            \sum_{l=0}^{\infty}\sum_{n=0}^{\infty}
            \sum_{p=1}^{N}\sum_{q>p}^{N}
              2\Re(c_{pl}c_{qn}^{*})\cc I_{lnpq}
\end{eqnarray}
where $\sigma_p^{(coupled)}$ is the scattering expression for the scattering of individual sphere. Note however that since $C_{pl}$ are different for each sphere than of the unperturbed case $\sigma_p^{(coupled)}\neq \sigma_p^{(0)}$.

<!-- Demo for this  -->
<!-- \begin{eqnarray}
  \sigma &=&\int \Big(\sum_p f_p(\theta)\Big)\Big(\sum_q f_q^{*}(\theta) \Big) d\Omega\\
        &=& \sum_{p=1}^{N}\sum_{q=1}^{N}c_{pl}c_{pn}^{*}
            \sum_l\sum_n\int Y_l(\theta)Y_n^{*}(\theta)d\Omega \\
        &=&\sum_{l=0}^{\infty}
            \sum_{p=1}^{N}\sum_{q=1}^{N}c_{pl}c_{ql}^{*} \\
        &=&\sum_{p=1}^{N}\sigma_p +
            \sum_{l=0}^{\infty}\sum_{p=1}^{N}\sum_{q\neq p}c_{pl}c_{ql}^{*} \\
        &=&\sum_{p=1}^{N}\sigma_p +
            2\sum_{l=0}^{\infty}\sum_{p=1}^{N}\sum_{q>p}
              \Re(c_{pl}c_{ql}^{*})
\end{eqnarray} -->



## Single hard sphere scattering

a) $\sigma(ka)$ | b) $f(\theta;ka)$ | c) $abs(h_l(ka))$
--------------- | ----------------- | ---------------
[<img src="{{figs}}1_ska.svg" width="250" />]({{figs}}1_ska.svg) |  [<img src="{{figs}}1_fka.svg" width="250" />]({{figs}}1_fka.svg) |  [<img src="{{figs}}1_hl.svg" width="250" />]({{figs}}1_hl.svg)

d) $c_{pl}(ka=0.5-5)$ | (e) $c_{pl}(ka=5-20)$ | (f) $c_{pl}(ka=30-40)$
--------------------- | --------------------- | ---------------------
[<img src="{{figs}}1_1cpl.svg" width="250" />]({{figs}}1_1cpl.svg) |  [<img src="{{figs}}1_2cpl.svg" width="250" />]({{figs}}1_2cpl.svg) | [<img src="{{figs}}1_cpl.svg" width="250" />]({{figs}}1_cpl.svg)

a) Integrated and b) differential cross sections for a single sphere with normalized radius $ka$.
c) The magnitude of the Hankel functions up to order $l=50$ and $ka=50$.
The modal coefficients $c_{pl}(ka)$ for d) $ka=0.5-5$, e) $ka=5-20$, f) $ka=30-40$.



## 2 hard spheres scattering

Case for 2 identical spheres of normalized radius $ka$ and inter-distance $kd$.

 --          | ka=0.5 | ka=1.0 | ka=2.0 | ka=4.0
------------ | ------ | ------ | ------ | ------
$C_{0l}(kd)$ |[![]({{figs}}_ka0cpl0.svg)]({{figs}}_ka0cpl0.svg) | [![]({{figs}}_ka1cpl0.svg)]({{figs}}_ka1cpl0.svg) | [![]({{figs}}_ka2cpl0.svg)]({{figs}}_ka2cpl0.svg) | [![]({{figs}}_ka3cpl0.svg)]({{figs}}_ka3cpl0.svg)
$C_{1l}(kd)$ |[![]({{figs}}_ka0cpl1.svg)]({{figs}}_ka0cpl1.svg) | [![]({{figs}}_ka1cpl1.svg)]({{figs}}_ka1cpl1.svg) | [![]({{figs}}_ka2cpl1.svg)]({{figs}}_ka2cpl1.svg) | [![]({{figs}}_ka3cpl1.svg)]({{figs}}_ka3cpl1.svg)

Scattering coefficients of sphere 0 and sphere 1 for the coupled(solid lines) and uncoupled(dashed lines) systems with varying inter-sphere distance $kd$.
Different graphs for $ka=0.5,1,2,4$. The different colours on each graph correspond the coefficient of a given mode order.


$\sigma(kd;ka)$ | $f(\theta;ka=2.0)$ | $f(\theta;ka=5.0)$ | $f(\theta;ka=10)$
--------------- | -----------------  | -----------------  | ------------------
[<img src="{{figs}}2_ska.svg"  width="180" />]({{figs}}2_ska.svg) |  [<img src="{{figs}}2_fka0.svg" width="180" />]({{figs}}2_fka0.svg) |  [<img src="{{figs}}2_fka1.svg" width="180" />]({{figs}}2_fka1.svg) |  [<img src="{{figs}}2_fka2.svg" width="180" />]({{figs}}2_fka2.svg)

Scattering amplitudes and total scattering cross section with inter-sphere distance for various normalized radius.


--   | $i$ | $s$ | $t$  
----- | --- | --- | ---
p0    | [![]({{figs}}2_fi_sphere0.png)]({{figs}}2_fi_sphere0.png) | [![]({{figs}}2_fs_sphere0.png)]({{figs}}2_fs_sphere0.png) | [![]({{figs}}2_ft_sphere0.png)]({{figs}}2_ft_sphere0.png)
p1    | [![]({{figs}}2_fi_sphere1.png)]({{figs}}2_fi_sphere1.png) | [![]({{figs}}2_fs_sphere1.png)]({{figs}}2_fs_sphere1.png) |  [![]({{figs}}2_ft_sphere1.png)]({{figs}}2_ft_sphere1.png)
p0+p1 | [![]({{figs}}2_fi.png)]({{figs}}2_fi.png) | [![]({{figs}}2_fs.png)]({{figs}}2_fs.png) | [![]({{figs}}2_ft.png)]({{figs}}2_ft.png)

Incident, scattered and total wave function for the uncoupled $0th$ sphere, uncoupled $1st$ sphere and coupled hard sphere system. Normalized radius $ka=1$ and inter-distance $kd=\pi$.




<!--
#############################################################
          Approximate solution
#############################################################
 -->
## Approximate solution

### Scattering amplitude
Assuming that the scattering amplitude at the $pth$ sphere from the $qth $ sphere is a plane wave with unknown scattering amplitude $c_q$, the scattering function of the $pth$ sphere is :
\begin{eqnarray}
  g_p^{(s)}(\bb r_p)  
    &=&g_{p}^{(s,0)}(r_p,\theta_p)\Big(1 + \sum_{q=1}^{p-1} c_q e^{-jkd_q} \Big) \\
    && g_{p}^{(s,0)}(r_p,\pi-\theta_p)\sum_{q=p+1}^{N}c_q e^{+jkd_q}
\end{eqnarray}

where $g_{p}^{(s,0)}$ is the far field scattering function due to an incident plane wave of the $pth$ sphere located at $d_p$ :
\begin{eqnarray}
  g_{p}^{(s,0)}(r_p,\theta_p) &=& \frac{e^{jkr_p}}{kr_p}f_p(\theta_p) \\
  f_p(\theta_p) &=& \sum_{l=0}^{\infty} (-j)^{l+1} c_{pl}^{(0)}Y_l^0(\theta_p)\\
  c_{pl}^{(0)} &=& \sqrt{\frac{4\pi}{2l+1}}b_l u_l(ka_p)e^{jkd_p}
\end{eqnarray}

### System
The coefficients are found solving the system :
\begin{eqnarray}
  \left(N-1\right)c_p = \sum_{q\neq p}e^{-jkd_{pq}}
     &&\Bigg[g_p(d_{pq},\theta_{pq})\\
    +&& g_p(d_{pq},    \theta_{pq})\sum_{r=1}^{p-1}c_re^{-jkd_r} \\
    +&& g_p(d_{pq},\pi-\theta_{pq})\sum_{r=p+1}^{N}c_re^{jkd_r}
    \Bigg]
\end{eqnarray}
where $\theta_{pq}=0$ for $q>p$ and $\theta_{pq}=\pi$ for $q<p$.

For the case of 2 identical spheres with $d_1=0$ and $d_2=d$,
using $f_1(\theta)=f(\theta)$, $f_2(\theta)=f(\theta)e^{jkd}$,
the system gives :
\begin{eqnarray}
  c_1 &=& \frac{f(0         )}{kd} + \frac{e^{jkd}f(\pi)}{kd}c_2\\
  c_2 &=& \frac{e^{jkd}f(\pi)}{kd} + \frac{e^{jkd}f(\pi)}{kd}c_1\\
\end{eqnarray}


### Far field scattering
{% set figsa='/figures/hardSphereArray2approx_' %}

\begin{eqnarray}
  f(\theta)   &=& \sum_{p=1}^{N} f_p(\theta)e^{-jkd_p\cos\theta} \\
  f_p(\theta) &=& \sum_{l=0}^{\infty} (-j)^{l+1}c_{pl}^{(0)}\Big(
      Y_l^0(\theta)
      + \sum_{q=1}^{p-1} c_qe^{-jkd_q}Y_l^0(\theta)
      + \sum_{q=p+1}^{N} c_qe^{ jkd_q}Y_l^0(\pi-\theta)
    \Big)\\
            &=&\sum_{l=0}^{\infty} (-j)^{l+1}c_{pl}Y_l^0(\theta)\\
  c_{pl}    &=&c_{p;l}^{(0)}\Big(
    1+\sum_{q=1}^{p-1}c_qe^{-jkd_q} + \sum_{q=p+1}^{N}(-1)^l c_qe^{jkd_q} \Big)            
\end{eqnarray}

where the parity of the spherical harmonics
$Y_l^0(\pi-\theta)=(-1)^lY_l^0(\theta)$ has been used.


$f(ka=1.0;kd=20)$ | $f(ka=5.0;kd=40)$ | $c0l(ka=2;kd)$ | $c1l(ka=2;kd)$
----------------- | ----------------- | -------------- | --------------
[<img src="{{figsa}}fka0.svg" width="150" />]({{figsa}}fka1.svg) |  [<img src="{{figsa}}fka1.svg" width="150" />]({{figsa}}fka1.svg) | [<img src="{{figsa}}ka0cpl0.svg" width="150" />]({{figsa}}ka0cpl0.svg) | [<img src="{{figsa}}ka0cpl1.svg" width="150" />]({{figsa}}ka0cpl1.svg)


<!--
#############################################################
          Spherical coordinates appendix
#############################################################
 -->
## apx:Spherical coordinates
### Plane wave expansion at normal incidence
The scalar plane wave travelling along $z$ is expanded upon spherical waves using :
\begin{equation}
  e^{jkz} = e^{jkr\cos\theta} =\sum_{l=0}^{\infty} j^{l}(2l+1) j_l(kr)P_l(\cos\theta)
\end{equation}
where $P_l$ are the Legendre polynomials.

<!-- A plane wave expansion upon spherical waves, with $\lambda=1$ over y,z=[-5,5]. -->

<!-- N=10 | N=20 | N=50 | N=99
---- | ---- | ---- | ----
[![]({{figs}}Exi_sphere1.png)]({{figs}}Exi_sphere1.png) |  [![]({{figs}}Exi_sphere2.png)]({{figs}}Exi_sphere2.png) | [![]({{figs}}Exi_sphere3.png)]({{figs}}Exi_sphere3.png) | [![]({{figs}}Exi_sphere4.png)]({{figs}}Exi_sphere4.png) -->

### Translational addition theorem coefficients
The spherical wave function $\psi_{l,0}^{out}=\hl Y_l^0(\theta_p)$ can be translated from $(r_p,\theta_p,\phi_p)$ to $(r,\theta,\phi)$ where $\bb r=\bb r_p+\bb t$ using :
\begin{eqnarray}
  \psi_{l,0}^{out}(\bb r_p) &=& \sum_{n=0}^{\infty}
    a_{l,0;n,0}^{out-in}(\bb t)\psi_{n,0}^{in}(\bb r)\\
  a_{l,0;n,0}^{out-in}(\bb t) &=& 4\pi\sum_{q=|l-n|}^{l+n}
    \left(-j\right)^{l-n-q}\psi^{out}_{q,0}(\bb t)\cc G(l,n,q,0,0,0)
\end{eqnarray}

Spherical Hankel functions of the first kind in $(r_p,\theta_p,\phi_p)$ and in reference $(r,\theta,\phi)$ using the translational addition theorem with spherical Bessel functions expansion.

$\psi_{0}$ | $\psi_{1}$ | $\psi_{2}$ | $\psi_{3}$  
----------- | ----------- | ----------- | -----------  
[![]({{figs}}_psi0.png)]({{figs}}_psi0.png) | [![]({{figs}}_psi1.png)]({{figs}}_psi1.png) | [![]({{figs}}_psi2.png)]({{figs}}_psi2.png) | [![]({{figs}}_psi3.png)]({{figs}}_psi3.png)
[![]({{figs}}_psi0_t.png)]({{figs}}_psi0_t.png) | [![]({{figs}}_psi1_t.png)]({{figs}}_psi1_t.png) | [![]({{figs}}_psi2_t.png)]({{figs}}_psi2_t.png) | [![]({{figs}}_psi3_t.png)]({{figs}}_psi3_t.png)
