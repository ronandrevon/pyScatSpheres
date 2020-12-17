import time
import numpy as np
import scipy.special as spe
import sage.all
import py3nj
import importlib as imp
import utils.displayStandards as dsp ;imp.reload(dsp)

gaunt = sage.functions.wigner.gaunt
w3j = lambda n1,n2,n3,m1,m2,m3 : py3nj.wigner3j(2*n1,2*n2,2*n3,2*m1,2*m2,2*m3)
yn   = lambda n,z : np.sqrt(np.pi/(2*z))*spe.yv(n+0.5,z)
jn   = lambda n,z : np.sqrt(np.pi/(2*z))*spe.jv(n+0.5,z)
hn1  = lambda n,z : np.sqrt(np.pi/(2*z))*spe.hankel1(n+0.5,z)
jnp  = lambda n,z : np.sqrt(np.pi/(2*z))*(spe.jvp(n+0.5,z)  - spe.jv(n+0.5,z)/(2*z) )
hn1p = lambda n,z : np.sqrt(np.pi/(2*z))*(spe.h1vp(n+0.5,z) - spe.hankel1(n+0.5,z)/(2*z) )
Pl = lambda l,theta: np.sqrt(4*np.pi/(2*l+1))*spe.sph_harm(0,l,0,theta)
# psi_lm = lambda z_l,l,m,r,theta,phi:z_l(l,r)*spe.sph_harm(m,l,phi,theta)

def get_MNnm(kr,theta,phi, Nmax, zn=hn1,znp=hn1p,pol='y'):
    '''Computes Mnm and Nmn in spherical coordinates up to Nmax'''
    Plms,dtPlm = get_jmdtPlm(Nmax,theta,phi,phi_opt=False,Ymn_opt=False)
    st = np.sin(theta)
    Ns = dtPlm.shape[0]
    ss = (Ns,3)+kr.shape
    Mnm,Nnm = np.zeros(ss,dtype=complex),np.zeros(ss,dtype=complex)
    kr[kr==0]=1e-3
    if pol=='xy':
        for n in range(1,Nmax+1):
            zn1kr  = zn(n,kr)           #;print(hn1kr)
            dzn1kr = znp(n,kr)+zn1kr/kr #;print(dhn1kr)
            for m in range(n+1):
                idm = int(n*(n+1)/2)+abs(m)
                ep = np.exp(1J*m*phi)   #; print(ep)
                Mnm[idm,1] =  zn1kr*1J*m*Plms[idm]*ep
                Mnm[idm,2] = -zn1kr*dtPlm[idm]*ep
                Nnm[idm,0] =  n*(n+1)*zn1kr/kr*Plms[idm]*st*ep #does not work for m=0 !!!!!!
                Nnm[idm,1] =  dzn1kr*dtPlm[idm]*ep
                Nnm[idm,2] =  dzn1kr*1J*m*Plms[idm]*ep
    #### wiki signs Mo,Ne - x polarization
    if pol=='x':
        for n in range(1,Nmax+1):
            zn1kr  = zn(n,kr)           #;print(hn1kr)
            dzn1kr = znp(n,kr)+zn1kr/kr #;print(dhn1kr)
            for m in range(n+1):
                idm = int(n*(n+1)/2)+abs(m)
                cp,sp = np.cos(m*phi), np.sin(m*phi)  #;print(cp)
                Mnm[idm,1] =  m*cp*Plms[idm]*zn1kr
                Mnm[idm,2] = -sp*dtPlm[idm]*zn1kr
                Nnm[idm,0] =  zn1kr/kr*cp*n*(n+1)*Plms[idm]*st #does not work for m=0 !!!!!!
                Nnm[idm,1] =  cp*dtPlm[idm]*dzn1kr
                Nnm[idm,2] = -m*sp*Plms[idm]*dzn1kr
    ### wiki signs Me,No -y polarization
    if pol=='y':
        for n in range(1,Nmax+1):
            zn1kr  = zn(n,kr)           #;print(hn1kr)
            dzn1kr = znp(n,kr)+zn1kr/kr #;print(dhn1kr)
            for m in range(n+1):
                idm = int(n*(n+1)/2)+abs(m)
                cp,sp = np.cos(m*phi), np.sin(m*phi)  #;print(cp)
                Mnm[idm,1] = -m*sp*Plms[idm]*zn1kr
                Mnm[idm,2] = -cp*dtPlm[idm]*zn1kr
                Nnm[idm,0] =  zn1kr/kr*sp*n*(n+1)*Plms[idm]*st #does not work for m=0 !!!!!!
                Nnm[idm,1] =  sp*dtPlm[idm]*dzn1kr
                Nnm[idm,2] =  m*cp*Plms[idm]*dzn1kr

    return Mnm,Nnm

def get_jmdtPlm(l,theta,phi=None,Ymn_opt=False,phi_opt=False,v=0):
    r'''Computes Plm/sin(cos(theta)) and d_theta Plm(cos(theta)) for :
    - l     : order so that Plm for Plmsin = [Pl00,Pl10,Pl11,Pl20,Pl21,Pl22,Pl30,...]/sin(theta)
    - theta : 2d-array
    - phi,phi_opt,Ymn_opt : obsolete
    returns :
        - Plmsin : $P_l^m(\cos\theta)/\sin\theta$
        - dtPlm  : $\Partial_{\theta}P_l^m(\cos\theta)$
    '''
    nls = int((l+1)*(l+2)/2)    #;print(nls)
    ct,st = np.cos(theta),np.sin(theta)
    if isinstance(theta, int) or isinstance(theta, float):
        ones,zeros,ss = 1,0,(nls)
    else:
        ones,zeros,ss = np.ones(ct.shape,dtype=float),np.zeros(ct.shape,dtype=float),((nls,)+theta.shape)

    if v:print('...Plm/sin...')
    Plmsin = np.zeros(ss,dtype=complex)
    Plmsin[0] = zeros
    Plmsin[1] = zeros
    Plmsin[2] = -ones
    for n in range(2,l+1):
        idn0 = int((n-2)*(n-1)/2)
        idn1 = int((n-1)*(n+0)/2)
        idn2 = int((n-0)*(n+1)/2)
        Plmsin[idn2] = zeros
        for p in range(1,n-1):
            Plmsin[idn2+p] = ((2*(n-1)+1)*ct*Plmsin[idn1+p] - (n-1+p)*Plmsin[idn0+p])/(n-1-p+1); # print('idnp:',idn2+p)
        Plmsin[idn2+n-1] = (2*(n-1)+1)*ct*Plmsin[idn1+n-1]      #; print('idnp:',idn2+n-1, 'using ', idn1+n-1)
        Plmsin[idn2+n] = -(2*(n-1)+1)*st*Plmsin[idn1+n-1]       #;print('idnp:',idn2+n, 'using ', idn1+n)

    if v:print('...derivative...')
    dtPlm  = np.zeros(ss,dtype=complex)
    for n in range(1,l+1):
        idn = int(n*(n+1)/2)
        for p in range(n):
            dtPlm[idn+p] = p*ct*Plmsin[idn+p] + st*Plmsin[idn+p+1]  #;print('idnp:',idn+p)
        dtPlm[idn+n] = n*ct*Plmsin[idn+n]                           #;print('idnp:',idn+n)

    if v:print('...u,v...')
    if phi_opt:
        print('warning : Trying to use Nlm with obsolete phi_opt')
        # for n in range(1,l+1):
        #     idn = int(n*(n+1)/2)
        #     for p in range(n+1):
        #         snm = 1#[1,np.sqrt((2*l+1)/(4*np.pi)*np.product(np.arange(1,l-p+1))/np.product(np.arange(1,l+p+1)))][Ymn_opt]
        #         # print('idnp:',idn+p)
        #         ep = np.exp(1J*p*phi)
        #         Plmsin[idn+p]   *= snm*1J*p*ep
        #         dtPlm[idn+p]    *= snm*ep
    # else:
    #     # print('here')
    #     for n in range(1,l+1):
    #         idn = int(n*(n+1)/2)
    #         for p in range(n+1):
    #             Plmsin[idn+p] *= p

    return Plmsin,dtPlm


def get_MNn1_t(kdp,n,kr,theta,phi,Nmax,zn=jn,znp=jnp,pol='y',fz=np.real):
    Annu,Bnnu = get_vec_trans_coeff_lin0(kdp,n,Nmax)
    Mn1,Nn1   = get_MNn1(kr,theta,phi, Nmax, zn,znp,pol)
    ss = (3,)+kr.shape
    Mr,Nr = np.zeros(ss,dtype=complex),np.zeros(ss,dtype=complex)
    em = 1J#{}
    for n in range(Nmax):
        # Ql1*Mlm[idlm] + im*Pl1*Nlm[idlm]
        Mr +=    Annu[n]*Mn1[n]+em*Bnnu[n]*Nn1[n]
        Nr += em*Annu[n]*Nn1[n]+   Bnnu[n]*Mn1[n]
    if fz:
        return fz(Mr),fz(Nr)
    else:
        return Mr,Nr

def get_MNn1(kr,theta,phi, Nmax, zn=hn1,znp=hn1p,pol='y'):
    '''Computes Mn1 and Nm1 in spherical coordinates up to Nmax'''

    Pl1s,dtPl1 = get_Pl1sdtPl1(Nmax,theta)
    st    = np.sin(theta)
    cp,sp = np.cos(phi), np.sin(phi)  #;print(cp)

    Ns = dtPl1.shape[0]
    ss = (Ns,3)+kr.shape
    kr[kr==0]=1e-3
    Mn1,Nn1 = np.zeros(ss,dtype=complex),np.zeros(ss,dtype=complex)
    m=1
    #### wiki signs Mo,Ne - x polarization
    if pol=='x':
        for n in range(1,Nmax+1):
            zn1kr  = zn(n,kr)
            dzn1kr = znp(n,kr)+zn1kr/kr
            Mnm[n-1,1] =  m*cp*Pl1s[n]*zn1kr
            Mnm[n-1,2] = -sp*dtPl1[n]*zn1kr
            Nnm[n-1,0] =  zn1kr/kr*cp*n*(n+1)*Pl1s[n]*st
            Nnm[n-1,1] =  cp*dtPl1[n]*dzn1kr
            Nnm[n-1,2] = -m*sp*Pl1s[n]*dzn1kr
    ### wiki signs Me,No -y polarization
    if pol=='y':
        for n in range(1,Nmax+1):
            zn1kr  = zn(n,kr)
            dzn1kr = znp(n,kr)+zn1kr/kr
            Mn1[n-1,1] = -m*sp*Pl1s[n]*zn1kr
            Mn1[n-1,2] = -cp*dtPl1[n]*zn1kr
            Nn1[n-1,0] =  zn1kr/kr*sp*n*(n+1)*Pl1s[n]*st
            Nn1[n-1,1] =  sp*dtPl1[n]*dzn1kr
            Nn1[n-1,2] =  m*cp*Pl1s[n]*dzn1kr

    return Mn1,Nn1

def get_Pl1sdtPl1(l,theta,v=0):
    r'''Computes Pl1/sin(cos(theta)) and dP_theta Pl1(cos(theta)) for :
    - l     : order so that Pl1 for Pl1sin = [Pl00,Pl10,Pl20,Pl30,...]/sin(theta)
    - theta : 2d-array
    - phi,phi_opt,Ymn_opt : obsolete
    returns :
        - Plmsin : $P_l^1(\cos\theta)/\sin\theta$
        - dtPlm  : $\Partial_{\theta}P_l^1(\cos\theta)$
    '''
    nls = l+1
    ct,st = np.cos(theta),np.sin(theta)
    if isinstance(theta, int) or isinstance(theta, float):
        ones,zeros,ss = 1,0,(nls)
    else:
        ones,zeros,ss = np.ones(ct.shape,dtype=float),np.zeros(ct.shape,dtype=float),((nls,)+theta.shape)

    if v:print('...Plm/sin...')
    Pl1s = np.zeros(ss,dtype=complex)
    Pl1s[1] = -ones
    for n in range(1,l):
        Pl1s[n+1] = ((2*n+1)*ct*Pl1s[n] - (n+1)*Pl1s[n-1])/n; # print('idnp:',idn2+p)

    if v:print('...derivative...')
    dtPl1  = np.zeros(ss,dtype=complex)
    dtPl1[1] = ct*Pl1s[1]
    for n in range(2,l+1):
        dtPl1[n] = n*ct*Pl1s[n] - (n+1)*Pl1s[n-1]  #;print('idnp:',idn+p)

    return Pl1s,dtPl1

################################################################################
#                           Plane wave
################################################################################
def get_plane_wave_scalar(y,z,lam=1,d_p=2.5,alpha=0,Nmax=10,fz=np.real,timeOpt=False):
    r'''Plane Wave Ex = exp(jkz) expansion upon Spherical waves using :
            $$e^{jkz} = e^{jkr\cos(\theta)}=sum_{n=0}^{\infty} a_nJ_n(kr)P_n(\cos(\theta))$$
    - lam : wavelength
    - d_p : distance to origin
    - alpha : incident angle degrees
    - Nmax  : max order expansion
    -fz     : np.real,np.imag,np.abs or None
    '''
    k,alpha,npts = 2*np.pi/lam, np.deg2rad(alpha),y.shape[0]
    r,theta = np.sqrt(y**2+z**2),np.arctan(y/z)
    Ex = np.zeros(z.shape,dtype=complex)
    n  = np.arange(Nmax+1)
    if alpha:
        t = time.time()
        for n in range(Nmax+1):
            ymn = np.zeros(Ex.shape,dtype=complex)
            for m in range(-n,n+1):
                ymn += spe.sph_harm(m,n,0,theta)*np.conj(spe.sph_harm(m,n,0,alpha))
            Ex += 1J**n*spe.spherical_jn(n,k*r)*ymn
        if timeOpt:print('loop:',time.time()-t)
        Ex *= 4*np.pi*np.exp(1J*k*d_p*np.cos(alpha))
    else:
        t = time.time()
        Pn = np.zeros((npts,npts,Nmax+1))
        for i in range(npts): #enumerate(np.cos(theta.flatten()):
            for j in range(npts):
                Pn[i,j,:] = spe.lpn(Nmax,np.cos(theta[i,j]))[0]
        if timeOpt:print('poly:',time.time()-t)

        t=time.time()
        Jn = spe.spherical_jn
        an = lambda n:1J**(-n)*(2*n+1)
        for n in range(Nmax+1):
            Ex += an(n)*Jn(n,k*r)*Pn[:,:,n]
        if timeOpt:print('sum Jn:',time.time()-t)
    if fz:
        return fz(Ex)
    else:
         return Ex


def get_plane_wave_vector(r,theta,phi=np.pi/2,lam=1,alpha=0,d_p=0,Nmax=10,fz=np.real,pol='x',v=False):
    r'''Plane Wave expansion upon Vector Spherical waves
        - r,theta,phi : grid coordinates
        - lam : wavelength
        - d_p : distance to origin
        - alpha : incident angle (degrees)
        - Nmax  : max order expansion (inclusive)
        -fz     : np.real,np.imag,np.abs or None
    returns :
        - E = Er,Et,Ep for plane wave
    '''
    k,alpha = 2*np.pi/lam, np.deg2rad(alpha)

    if v:print('\t\t plane wave vector')
    E = np.zeros((3,)+r.shape,dtype=complex)

    if alpha>0:
        if v:print('...Mlm,Nlm,Plm(alpha)... : ',end='');t = time.time()
        Mlm,Nlm = get_MNnm(k*r,theta,phi, Nmax, zn=jn, znp=jnp,pol=pol)
        jmPlm,dtPlm = get_jmdtPlm(Nmax,alpha,phi=phi,Ymn_opt=False,v=0)
        if v:print('%.2f'%(time.time()-t))

        if v:print('..summation.. : ',end='');t = time.time()
        for l in range(1,Nmax+1):
            for m in range(-l,l+1):
                slm = np.product(np.arange(1,l-abs(m)+1))/np.product(np.arange(1,l+abs(m)+1))
                idlm = int(l*(l+1)/2)+abs(m)
                Plm = -1J**l*(2*l+1)/(l*(l+1))*slm*jmPlm[idlm]
                Qlm = -1J**l*(2*l+1)/(l*(l+1))*slm*dtPlm[idlm]
                # print('l=%-3d,m=%-3d, slm=%.1E,' %(l,m,slm))
                # print('\t P=%.1E, Q=%.1E' %(Plm,Qlm))#, jmPlm[idlm.max(),dtPlm[idlm].max()))
                E += Plm*Nlm[idlm]+Qlm*Mlm[idlm]
    else:
        # Ml1,Nl1 = get_MNn1(k*r,theta,phi, Nmax, zn=jn, znp=jnp)
        Mlm,Nlm = get_MNnm(k*r,theta,phi, Nmax, zn=jn, znp=jnp,pol=pol)
        if v:print('..summation.. : ',end='');t = time.time()
        im = {'x':-1J,'y':1J,'xy':1}[pol]
        for l in range(1,Nmax+1):
            idlm = int(l*(l+1)/2)+1
            Pl1 = -1J**l*(2*l+1)/(2*l*(l+1))
            Ql1 = -1J**l*(2*l+1)/(2*l*(l+1))
            # print('l=%-3d,m=%-3d, slm=%.1E,' %(l,m,slm))
            # print('\t P=%.1E, Q=%.1E' %(Plm,Qlm))#, jmPlm[idlm.max(),dtPlm[idlm].max()))
            E +=  Ql1*Mlm[idlm] + im*Pl1*Nlm[idlm]
        if d_p:E *= np.exp(1J*k*d_p)

    if v:print('%.2f'%(time.time()-t))

    if fz:
        return fz(E)
    else:
         return E


########################################################################################
#                               Translation addition theorem
########################################################################################
# def get_vec_trans_coeff(l,m,n,p,**kwargs):
#     '''vectorial translational addition theorem coefficients'''
#     Almnp = 1/(n*(n+1))*(
#         0.5*np.sqrt((l-m)*(l+m+1)*(n-p)*(n+p+1))*a_lmnp(l,m+1,n,p+1,**kwargs)+
#         0.5*np.sqrt((l+m)*(l-m+1)*(n+p)*(n-p+1))*a_lmnp(l,m-1,n,p-1,**kwargs)+
#         m*p*a_lmnp(l,m,n,p,**kwargs)
#         )
#     Blmnp =(
#     -0.5*np.sqrt((l-m)*(l+m+1))*(
#         1/(n+1)*np.sqrt((n+p+1)*(n+p+2)/((2*n+1)*(2*n+3)))*a_lmnp(l,m+1,n+1,p+1,**kwargs)+
#         1/n    *np.sqrt((n-p-1)*(n-p  )/((2*n-1)*(2*n+1)))*a_lmnp(l,m+1,n-1,p+1,**kwargs)
#         )
#     +0.5*np.sqrt((l+m)*(l-m+1))*(
#         1/(n+1)*np.sqrt((n-p+1)*(n-p+2)/((2*n+1)*(2*n+3)))*a_lmnp(l,m-1,n+1,p-1,**kwargs)+
#         1/n    *np.sqrt((n+p-1)*(n+p  )/((2*n-1)*(2*n+1)))*a_lmnp(l,m-1,n-1,p-1,**kwargs)
#         )
#     +m*(
#         1/(n+1)*np.sqrt((n+p+1)*(n-p+1)/((2*n+1)*(2*n+3)))*a_lmnp(l,m,n+1,p,**kwargs)+
#         -1/n   *np.sqrt((n+p  )*(n-p  )/((2*n-1)*(2*n+1)))*a_lmnp(l,m,n-1,p,**kwargs)
#         )
#     )
#     return Almnp,Blmnp
#

def a_ln(lmax,nmax, fz_q,r_d,theta_d):
    '''Scalar translational addition theorem coefficients
    - Nmax : order of expansion
    - fz_q : type of expansion (bessel or hankel)
    - r_d,theta_d : translation vector
    '''
    l = np.arange(lmax+nmax+2)
    Yq  = np.sqrt((2*l+1)/(4*np.pi))*spe.lpn(lmax+nmax+1,np.cos(theta_d))[0]
    # Yq = np.array([spe.sph_harm(0,q,0,theta_d) for q in range(lmax+nmax)])
    zq = np.array([fz_q(q,r_d) for q in range(lmax+nmax+2)])

    aln = np.zeros((lmax+1,nmax+1),dtype=complex)
    for l in range(lmax+1):
        # print('l=%d' %l)
        for n in range(nmax+1):
            q = np.arange(abs(l-n),n+l+1)   #;cput = time.time()
            Glnq = np.array([np.float(gaunt(l,n,q, 0,0,0)) for q in q])
            # Glnq = spu.w3j(l,n,q, 0,0,0)
            aln[l,n] = 4*np.pi*np.sum((-1J)**(l-n-q)*zq[q]*Yq[q]*Glnq) #;print('n=%d' %n, q,zq[q],Yq[q],Glnq)
    return aln

def a_lmnp(l,m,n,p, fz_q,r_d,theta_d,phi_d):
    '''Scalar translational addition theorem coefficients
    - l,m,n,p : index to pass from psi_np to psi_lm
    - fz_q : type of expansion (bessel or hankel)
    '''
    q = np.arange(abs(l-n),n+l+1)
    cput = time.time()
    # Glmnpq = np.array([np.float(gaunt(l,n,q, m,p,m-p)) for q in q])                           #;printf('Glmnpq:',t-cput) = time.time()
    Glmnpq = np.array([np.float(gaunt(l,n,q, m,-p,-m+p)) for q in q])                           #;printf('Glmnpq:',t-cput) = time.time()
    z_q    = np.array([fz_q(q,r_d) for q in q])                                                 #;printf('zq    :',t-cput) = time.time()
    Ympq   = np.array([ [0,spe.sph_harm(m-p,q,phi_d,theta_d)][bool(abs(m-p)<=q)] for q in q ])  #;printf('Ympq  :',t-cput) = time.time()
    # a_lmnp = 4*np.pi*np.sum(
    #     (-1J)**(n-l-q)*z_q*Ympq*Glmnpq*(-1)**m)
    a_lmnp = 4*np.pi*np.sum(
        (-1J)**(l-n-q)*z_q*Ympq*Glmnpq*(-1)**m)
    return a_lmnp

def get_trans_coeff(n,p,r_d,theta_d,phi_d, N=5):
    a_lmnp,i = np.zeros((N**2),dtype=complex),0
    for l in range(N):
        q = np.arange(abs(l-n),n+l+1)
        Jq = np.array([spe.spherical_jn(q,r_d) for q in q])
        # print('l=%-3d' %l,q,Jq )
        for m in range(-l,l+1):
            Glmnpq = np.array([np.float(gaunt(l,n,q, m,-p,-m+p)) for q in q])
            Ympq = np.array([ [0,spe.sph_harm(m-p,q,phi_d,theta_d)][bool(abs(m-p)<=q)] for q in q ])
            # print('\t m=%-3d : ' %m,Glmnpq,Ympq)
            a_lmnp[i] = 4*np.pi*np.sum(
                (-1J)**(l-n-q)*(-1)**m*Glmnpq*Jq*Ympq)
            i+=1
    return a_lmnp

def translation_addth_scalar(r,theta,phi, l,m,r_d,theta_d,phi_d,Nmax=10,zl=jn,zq=jn):
    '''Scalar translation addition theorem for psi_lm
    - r,theta,phi       : coordinates at which output is calculated
    - r_d,theta_d,phi_d : centre of spherical coordinates reference
    - l,m   : index of output function
    - Nmax  : max order of expansion
    - zl : function of the expansion jn('in') or hn1('out')
    - zq : function jn(out-out or in-in) or hn1(out-in or in-out)
    returns:
    - psi_lm : psi_lm(r,theta,phi) = almnp(r_d,theta_d,phi_d)*psi_np(r,theta,phi)
    '''
    # print(r_d,theta_d,phi_d,Nmax)
    ns = np.hstack([ [n]*(2*n+1) for n in np.arange(Nmax)])
    ps = np.hstack([np.arange(-n,n+1) for n in np.arange(Nmax)])
    psi_lm = np.zeros(r.shape,dtype=complex)
    a=[]
    for n,p in zip(ns,ps):
        almnp = a_lmnp(l,m,n,p, zq,r_d,theta_d,phi_d)
        psi_lm += almnp*zl(n,r)*spe.sph_harm(p,n,phi,theta) #;print(n,p,almnp)
        a+=[almnp]
    a=np.abs(np.array(a));ida=a>0;a=a[ida]; print(a,ns[ida],ps[ida])
    # print(['%.2E' %i.real for i in a])

    # almnp = get_trans_coeff(l,m,r_d,theta_d,phi_d, N=Nmax)
    # for i in range(ns.size):
    #     psi_lm += almnp[i]*zl(ns[i],r)*spe.sph_harm(ps[i],ns[i],phi,theta)

    # print(['%.2E' %i for i in almnp.real])
    return psi_lm



def get_vec_trans_coeff_lin0(kdp,n,Nmax=5):
    '''vectorial translational addition coefficients for endfired (normal incidence alpha=0) linear configuration
        - kdp  : normalized translation distance
        - n    : int - order(n>=1) of function to translate
        - Nmax : max included order
    returns :
        - Annu,Bnnu : coefficients nu=1..Nmax with m=mu=1
    '''
    nus = np.arange(Nmax)+1
    hp = np.array([hn1(p, kdp) for p in range(n+Nmax)])
    Annu,Bnnu = np.zeros((Nmax),dtype=complex),np.zeros((Nmax),dtype=complex)
    for nu in nus:
        p = np.arange(abs(n-nu),n+nu)
        Annup,Annupq = (2*p+1)*np.sqrt(n*(n+1)/(nu*(nu+1)))*w3j(n,nu,p,1,-1,0)*np.array([w3j(n,nu,p,0,0,0),w3j(n,nu,np.maximum(p-1,0),0,0,0)])
        annup = (1J)**(nu-n+p)/(2*nu*(nu+1))*(
            2*nu*(nu+1)*(2*nu+1) + (nu+1)*(n+nu-p)*(n+p-nu+1) - nu*(n+nu+p+2)*(nu+p-n+1))
        bnnup = -(1J)**(nu-n+p)*(2*nu+1)/(2*nu*(nu+1))*np.sqrt(
            (n+nu+p+1)*(nu+p-n)*(n+p-nu)*(n+nu-p+1))

        Annu[nu-1] = np.sum( Annup *annup*hp[p])
        Bnnu[nu-1] = np.sum(-Annupq*bnnup*hp[p])

        #f=lambda x:np.abs(x).max()
        #print('n=%d, nu=%d, p =' %(n,nu), p);print('hp=%.1f,Anp=%.1f,anp=%.1f,Anq=%.1f,bnp=%.1f' %(f(hp),f(Annup),f(annup),f(Annupq),f(bnnup)))

    return Annu,Bnnu


##################################################################################################
#                   MESH, COORDS
##################################################################################################
def sphere2cart(r,theta,phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x,y,z

def cart2sphere(x,y,z):
    r   = np.sqrt(x*x+y*y+z*z)
    rho = np.sqrt(x**2+y**2)
    theta = np.arccos(z/r)
    phi   = np.arctan2(y,x)
    return r,theta,phi

def Esphere2cart(Esph,theta,phi):
    ''' Convert field map from spherical coordinates to cartesian
    - Esp : spherical field map
    - theta,phi : angular spherical coordinates
    returns :
    - E : cartesian field map
    '''
    E = np.zeros(Esph.shape)
    ct,st = np.cos(theta),np.sin(theta)
    cp,sp = np.cos(phi),np.sin(phi)

    E[0] =  st*cp*Esph[0]  + ct*cp*Esph[1] - sp*Esph[2] #E_x
    E[1] =  st*sp*Esph[0]  + ct*sp*Esph[1] + cp*Esph[2] #E_y
    E[2] =  ct   *Esph[0]  - st   *Esph[1] + 0 *Esph[2] #E_z
    return E

def polar_mesh2D(cart_opt=False,npts=100,r=(0,1)):
    ''' 2D polar mesh :
        - cart_opt : True(cartesia mesh) ot False(polar mesh)
        - r     :
            - cart_opt=True  (4-tuple) - y[0],z[0]
            - cart_opt=False (2-tuple) - range of radius
        - npts  : int - discretization
    returns :
        - r,theta,y,z : polar mesh and corresponding cartesian
    '''
    if cart_opt:
        if len(r)==2:r=r*2
        y,z = np.meshgrid( np.linspace(r[0],r[1],npts), np.linspace(r[2],r[3],npts))
        r,theta = np.sqrt(y**2+z**2),np.arctan2(y,z)
    else:
        r,theta = np.meshgrid(np.linspace(r[0],r[1],npts),np.linspace(0,np.pi,npts))
        z,y = r*np.cos(theta), r*np.sin(theta)
    return r,theta,y,z

def mesh_3D(npts=50,r=(0,1)):
    if len(r)==2:r=r*3
    u = np.linspace(r[0],r[1],npts)
    v = np.linspace(r[2],r[3],npts)
    w = np.linspace(r[4],r[5],npts)
    x,y,z = np.meshgrid(u,v,w)
    return x,y,z

##################################################################################################
#                   Display_utils
##################################################################################################
def plot_E3d_plane(fE,cmpts=[0,1,2],cart_opt=True,r=5,npts=50,phi=0,name='E',fieldname='$E$',caxis=None,**kwargs):
    r'''plots 3d field map in phi-planes
    - fE : function to compute E as fE(r,theta,phi)
    - cmpts     : components to plot
    - cart_opt  : 'M'(cartesian mesh), 'E'(Fieldmap)
    - r,npts    : mesh options
    - phi       : phi-plane
    - field : str - field name for title
    - name  : str - fig name prefix - name+'%s.png' %cmpt
    '''
    if isinstance(cart_opt,bool):['','ME'][cart_opt]

    r,theta,y,z = polar_mesh2D('M' in cart_opt,npts,r)
    E = fE(r,theta,phi)

    if 'E' in cart_opt:
        E = Esphere2cart(E,theta,phi)
        scmpt,ncmpt=['x','y','z'],['x','y','z']
    else:
        scmpt,ncmpt = ['r',r'\theta',r'\phi'],['r','t','p']

    if not caxis:
        mm = np.abs([E.min(),E.max()]).max()
        caxis=[-mm,mm]
    #display
    for cmpt in cmpts:
        tle  = 'Field %s \n $e_{%s}$ component' %(fieldname,scmpt[cmpt])
        nameE = '_%s.png' %(ncmpt[cmpt])

        dsp.stddisp(im=[y,z,E[cmpt]],labs=['$y$','$z$'],title=tle,
            imOpt='cv',caxis=caxis,axPos='V',
            name=name+nameE,**kwargs)

class sphere_array:
    def biRCS_show(self,npts=100,is_3d=True,phi=[0,90],bargs={},im='',**kwargs):
        ''' plots bistatic RCS
        - npts  : int - nb points in theta and phi(if is_3d=1)
        - is_3d : bool - True=3d plot
        - phi   : list or np.ndarray - discrete values of phi(degrees) for 2d case(is_3d=0)
        '''
        theta = np.linspace(0,np.pi,2*npts)
        if is_3d:
            phi = np.linspace(0,2*np.pi,npts)
        else:
            phi = np.deg2rad(phi)
        theta3d,phi3d = np.meshgrid(theta,phi)
        rcs = self.compute_biRCS(theta3d,phi3d,**bargs)

        tle = r'bistatic RCS $\lambda$=%.1f, $a_p$=%.1f, $ka=$%.1f  $\alpha$=%d$^{\circ}$, $N=%d$' %(self.lam,self.a_p[0],2*np.pi*self.a_p[0]/self.lam, np.rad2deg(self.alpha),self.N)
        if is_3d:
            x,y,z =sphu.sphere2cart(np.log10(rcs),theta3d,phi3d)
            dsp.stddisp(scat=[x,y,z,np.abs(z)],rc='3d',labs=['$x$','$y$','$z$'],**kwargs)
        else:
            mono = self.compute_monoRCS(**bargs)
            cs   = dsp.getCs('Spectral',phi.size); # print(theta.shape,rcs.shape)
            plts = [[np.rad2deg(theta),rcs[i,:],cs[i],r'$\phi=%d^{\circ}$' %phi] for i,phi in enumerate(np.rad2deg(phi))]
            plts += [[180,mono,'bo' ,'$mono$']]
            labs = [r'$\theta(^{\circ})$',r'$\sigma(\theta,\phi)/\pi a_0^2$']
            if im:
                dsp.image_bg(im,plots=plts,lw=2,labs=labs,title=tle,
                    **kwargs)
            else:
                dsp.stddisp(plts,lw=2,labs=labs,title=tle,name=name,
                    **kwargs)

    def fields_show(self,ftype='t',cmpts=[0,1,2],idp=None,phi=np.pi/2,name='',
        cart_opt=False,npts=200,r=(1,10),**kwargs):
        ''' Plot the near field electric field
        - ftype : str - 't'(total), 's'(scattered), 'i'(incident)
        - cmpts : list- indices of components to plot 0(E_r) 1(E_theta) 2(E_phi)
        - idp   : int - sphere index (None => all)
        - phi   : float - angle(rad) for the phi plane
        - name  : prefix for name of figure where naming convention will be :
            name+'E%s%s_%s.png' %(ftype,cmpt,idp)
        - cart_opt,npts,r : mesh options, see spherical_utils.polar_mesh2D
        '''
        print('...Computing fields...')
        r,theta,y,z = sphu.polar_mesh2D(cart_opt,npts,r)
        E = self.compute_fields(r,theta,phi,ftype,idp)

        print('...Plotting fields...')
        ftypes = {'t':'Total','s':'Scattered', 'i':'Incident'}
        str_p,str_pn='',''
        if isinstance(idp,int):str_p, str_pn=' sphere %d' %(idp), '_%d' %(idp)

        #spheres
        t = np.linspace(-np.pi/2,np.pi/2,100)
        ct,st = np.cos(t),np.sin(t)
        plts = [ [ap*ct, dp+ap*st,'k-',''] for ap,dp in zip(self.a_p,self.d_p)]
        for cmpt in cmpts:
            mm = np.abs([E[cmpt].min(),E[cmpt].max()]).max()
            tle = r'%s Field $E^{%s}_{%s}$ '%(ftypes[ftype],ftype,['r',r'\theta', r'\phi'][cmpt])
            tle+=str_p
            dsp.stddisp(plots=plts,im=[y,z,E[cmpt]],labs=['$y$','$z$'],title=tle,
                imOpt='cv',caxis=[-mm,mm],axPos='V',lw=2,
                name=name+'E%s%s%s.png' %(ftype,['r','t','p'][cmpt],str_pn),**kwargs)
        return y,z,E

    def save(self,file=''):
        if not file:file='PECspheres%d%d.pkl' %(self.N,self.lam)
        with open(file,'wb') as out :
            pickle.dump(self, out, pickle.HIGHEST_PROTOCOL)
        if v:print(colors.green+"object saved\n"+colors.yellow+file+colors.black)


# spe.lmnp
# mnify = lambda x:np.hstack([np.hstack([np.flipud(x[1:n+1,n]),x[:n+1,n]]) for n in range(Nmax+1)])
# m = np.hstack([np.arange(-n,n+1) for n in np.arange(Nmax+1)])
# n = np.hstack([[n]*(2*n+1) for n in np.arange(Nmax+1)])


# if __name__== "__main__":
#     plt.close('all')
    # print(get_jmdtPlm(l=3,theta=0,phi=np.pi,v=0))
    # print(get_jmdtPlm(l=1,theta=np.pi/6))
