# Spherical coordinates
## links
- [Legendre poly](https://en.wikipedia.org/wiki/Associated_Legendre_polynomials#Reparameterization_in_terms_of_angles)
- [scipy spherical harmonics](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.sph_harm.html#scipy.special.sph_harm)
- [spherical coordinates](https://en.wikipedia.org/wiki/Spherical_coordinate_system#Integration_and_differentiation_in_spherical_coordinates)
- [$\grad\times$ in spherical](https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates)
- [Spherical Bessel and Hankel functions](https://en.wikipedia.org/wiki/Bessel_function#Spherical_Bessel_functions)



## Plane wave with spherical harmonics expansion

### Scalar normal incidence
A scalar plane wave along $z$ is expanded upon spherical waves using :
\begin{equation}
  E_x = e^{jkz} =\sum_{n=0}^{\infty} a_nJ_n(kr)P_n(\cos\theta)
\end{equation}
where $J_n$ are the spherical Bessel functions of the first kind, $P_n$ the associated Legendre polynomials of the first kind and $a_n=j^{-n}(2n+1)$.

A plane wave expansion upon spherical waves, with $\lambda=1$ over y,z=[-5,5].

N=10 | N=20 | N=50 | N=99
---- | ---- | ---- | ----
[![](/figures/Exi_sphere1.png)](/figures/Exi_sphere1.png) |  [![](/figures/Exi_sphere2.png)](/figures/Exi_sphere2.png) | [![](/figures/Exi_sphere3.png)](/figures/Exi_sphere3.png) | [![](/figures/Exi_sphere4.png)](/figures/Exi_sphere4.png)

### Scalar at oblique incidence
A scalar plane wave at an angle $\alpha=20^{\circ}$ is expanded upon spherical harmonics using [wiki](https://en.wikipedia.org/wiki/Plane_wave_expansion):
\begin{equation}
  E_x(r,\theta,\phi) = 4\pi
    \sum_{n=0}^{\infty}\sum_{m=-n}^{n} j^n J_n(kr)Y_n^{m*}(\theta,\phi)Y_n^m(\alpha,0)
\end{equation}

N=10 | N=30 | N=40 | N=50
---- | ---- | ---- | ----
[![](/figures/Exi_alpha_sphere1.png)](/figures/Exi_alpha_sphere1.png) | [![](/figures/Exi_alpha_sphere2.png)](/figures/Exi_alpha_sphere2.png) | [![](/figures/Exi_alpha_sphere3.png)](/figures/Exi_alpha_sphere3.png) | [![](/figures/Exi_alpha_sphere4.png)](/figures/Exi_alpha_sphere4.png)





## Translational addition theorem
- [python sage wigner3j symbols](https://doc.sagemath.org/html/en/reference/functions/sage/functions/wigner.html)
- [translational addition theorem](/articles/addVSH.pdf)

### Scalar
Translates the spherical wave function $\psi_{l,m}$ from reference $O(\bb r)$ to $O^{'}$ where $\bb t=\bb r^{'} - \bb r$.
\begin{eqnarray}
  \psi_{l,m}^{(out)}(\bb r^{'}) &=&
    \sum_{n=0}^{\infty}\sum_{p=-n}^{n}
    a_{l,m;n,p}(\bb t)\psi_{n,p}^{(out)}(\bb r)\\
  a_{l,m;n,p}^{out-out}(\bb t) &=&
    4\pi\sum_{q=|l-n|}^{l+n}
    \left(-i\right)^{l-n-q}\psi^{(in)}_{q,m-p}(\bb t)
    \left(-1\right)^m\cc G(l,n,q,m,-p,-m+p)
\end{eqnarray}
where
$\psi_{l,m} = z_l(r)Y_l^m(\theta,\phi)$ and
$z_l=j_l$ for $\psi^{(in)}$ and
$z_l=h_l^{(1)}$ for $\psi^{(out)}$ and  
$\cc G$ are the Gaunt coefficients.


Application for a few $\psi_{l,m}$ the case $\bb t=d\bb{e_r}= d\left(\cos\theta\bb{e_z} + \sin\theta\bb{e_y}\right)$

Normal Hankel functions of the first kind.

$\psi_{00}$ | $\psi_{10}$ | $\psi_{20}$ | $\psi_{30}$  
----------- | ----------- | ----------- | -----------  
[![](/figures/psi00_out.png)](/figures/psi00_out.png) | [![](/figures/psi10_out.png)](/figures/psi10_out.png) | [![](/figures/psi20_out.png)](/figures/psi20_out.png) | [![](/figures/psi30_out.png)](/figures/psi30_out.png)


Translated Bessel functions of the first kind.

$\psi_{00}$ | $\psi_{10}$ | $\psi_{11}$ | $\psi_{20}$  
----------- | ----------- | ----------- | -----------  
[![](/figures/gaunt00.png)](/figures/gaunt00.png) | [![](/figures/gaunt10.png)](/figures/gaunt10.png) | [![](/figures/gaunt11.png)](/figures/gaunt11.png) | [![](/figures/gaunt20.png)](/figures/gaunt20.png)



### Building the Legendre polynomials
In practice, the $P_l^m/\sin\theta$ and $\dP_{\theta}P_l^m$ are found recursively using :
\begin{eqnarray}
  m\frac{P_l^m(\cos\theta)}{\sin\theta} &=& 0
    ~\mbox{,}~ m=0\\
  \frac{P_{l+1}^m(\cos\theta)}{\sin\theta} &=&\frac{1}{l-m+1}\Big[
    (2l+1)\cos\theta\frac{P_l^m(\cos\theta)}{\sin\theta}\\
    &&-(l+m)\frac{P_{l-1}^m(\cos\theta)}{\sin\theta}
    \Big]
    ~\mbox{,}~ 0<m<l\\
  \frac{P_{l+1}^{m+1}(\cos\theta)}{\sin\theta} &=&
    -(2l+1)\sin\theta\frac{P_l^m(\cos\theta)}{\sin\theta}
    ~\mbox{,}~ m=l\\
  \dP_{\theta}P_l^m(\cos\theta) &=&
      m\cos\theta\frac{P_l^m(\cos\theta)}{\sin\theta} \\
      &&+\sin\theta\frac{P_l^{m+1}(\cos\theta)}{\sin\theta}
      ~\mbox{,}~ m\geq l
\end{eqnarray}
where $P_l^{m+1}=0$ if $m\geq l$.
