from utils import *
import matplotlib.patches as patches
import scipy,copy,time
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import spherical_utils as spu   ;imp.reload(spu)
from scipy.integrate import trapz,quad
import scipy.special as spe
plt.close('all')
ka,kp=6,1.01
N=10
nmax = 10
kd = 10
kds = np.arange(25,30,5)
opts = 'c'

qdot1 = qsa.QdotSphereArray(N=1,ka=ka,kp=kp,nmax=nmax,solve=0)
qdot1.solve(opt2=1)
# qdot1.test_convergence([10,11,12,13,15])
t0 = time.time()
qdot1 = qsa.QdotSphereArray(N=10,ka=ka,kp=kp,nmax=nmax,solve=1)
print(time.time()-t0)

if 'f' in opts:
    qdot2   = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=1)
    qdot2_0 = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=0)
    qdot2.test_convergence([5,7,9])#,11])

    theta_d = np.linspace(0,60,2000)
    theta = np.deg2rad(theta_d)
    ff =  qdot2.get_ff(theta)
    ffu = qdot2_0.get_ff(theta)
    ffmax = abs(ff).max()
    plts = [[theta_d,abs(ff)/ffmax,'b-','exact'],[theta_d,abs(ffu)/ffmax,'r--','approx']]
    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$f(\theta)$"],lw=2)
    idx = abs(ff)/ffmax>1e-3
    err = abs((ff[idx]-ffu[idx])/ff[idx]).mean()
    print(err)
    # if err<0.05:qdot2.show_f()

if 'c' in opts:
    qdot2u = [qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd,nmax=nmax,solve=1,copt=1) for kd in kds]

    cs = dsp.getCs('Spectral',nmax+1)

    ll = np.arange(nmax+1)
    jll = 1#1J**ll
    bp0 = qdot1.bp/jll
    bp1 = np.array([qdot.bp[:nmax+1] for qdot in qdot2u])/jll
    bp2 = np.array([qdot.bp[1*(nmax+1):2*(nmax+1)] for qdot in qdot2u])*np.exp(-1J*kds)[:,None]/jll
    bp3 = np.array([qdot.bp[2*(nmax+1):3*(nmax+1)] for qdot in qdot2u])*np.exp(-2J*kds)[:,None]/jll

    # for qdot in qdot2u:qdot.solve(optT=1)
    # bpT = np.array([qdot.bp[nmax+1:] for qdot in qdot2u])*np.exp(-1J*kds)[:,None]/jll

    # alnu = spu.a_ln(nmax,nmax,spu.hn1,kd,0,0)
    # bps = np.sum(qdot1.bp)
    vl = qdot1.vl
    bp2a = np.array([
        qdot1.bp*np.exp(1j*kd)+vl*np.sum(qdot1.bp*spu.a_ln(nmax,nmax,spu.hn1,kd,0,0).T,axis=1) for kd in kds])*np.exp(-1J*kds)[:,None]/jll
    bp3a = np.array([
        qdot1.bp*np.exp(2j*kd)+vl*np.sum(qdot1.bp*(spu.a_ln(nmax,nmax,spu.hn1,2*kd,0,0).T+spu.a_ln(nmax,nmax,spu.hn1,kd,0,0).T),axis=1) for kd in kds])*np.exp(-2J*kds)[:,None]/jll

    plts,legElt=[],{}
    plts+= [[bp0.real[i],bp0.imag[i],[cs[i],'o'],''] for i in range(nmax+1)]     ;legElt['$b_p^{(0)}$']=['k','o' ]
    # plts+= [[bp1.real[:,i],bp1.imag[:,i],[cs[i],'v-'],''] for i in range(nmax+1)];legElt['$b_1}$'     ]=['k','v-']
    # plts+= [[bp2.real[:,i],bp2.imag[:,i],[cs[i],'s-' ],''] for i in range(nmax+1)];legElt['$b_2}$'      ]=['k','s-']
    # plts+= [[bp2a.real[:,i],bpa.imag[:,i],[cs[i],'x--'],''] for i in range(nmax+1)];legElt['$b_2^{(a)}$' ]=['k','x--']
    plts+= [[bp3.real[:,i],bp3.imag[:,i],[cs[i],'s-' ],''] for i in range(nmax+1)];legElt['$b_3}$'      ]=['k','s-']
    plts+= [[bp3a.real[:,i],bp3a.imag[:,i],[cs[i],'x--'],''] for i in range(nmax+1)];legElt['$b_3^{(a)}$' ]=['k','x--']
    # plts+= [[bpT.real[:,i],bpT.imag[:,i],[cs[i],'d--'],''] for i in range(nmax+1)];legElt['$b_2^{(T)}$' ]=['k','d--']
    # legElt = {'$b_p^{(0)}$' :['k','o'],'$b_1$' :['k','d-'],
    #     '$b_2$' :['k','s-'],'$b_2^{(a)}$' :['k','x--']}
    dsp.stddisp(plts,labs=[r"$\Re(b_p)$",r"$\Im(b_p)$"],lw=2,legElt=legElt)

if 'C' in opts:
    nu = 1
    theta0  = np.linspace(0,np.pi,10)
    theta_p = np.hstack([theta0,np.flipud(theta0)])
    phi_p   = np.hstack([np.pi/2*np.ones(theta0.shape),-np.pi/2*np.ones(theta0.shape)])
    x,y,z   = spu.sphere2cart(ka*np.ones(theta_p.shape),theta_p,phi_p)
    kr,theta,phi =spu.cart2sphere(x,y,z+kd)

    psi_nu0 = spu.hn1(nu,kr)*spe.sph_harm(0,nu,phi,theta)

    ul,vl = qdot1.ul,qdot1.vl
    n_p  = qdot1.kp
    alnu = spu.a_ln(nu,nmax,spu.hn1,kd,0,0)[-1]
    # Yl0 = np.array([ spe.sph_harm(0,l,phi_p,theta_p) for l in range(nmax)])
    psi_out,psi_in=np.zeros(psi_nu0.shape,dtype=complex),np.zeros(psi_nu0.shape,dtype=complex)
    psi_nu0p=np.zeros(psi_nu0.shape,dtype=complex)
    for l in range(nmax):
        Yl0 = spe.sph_harm(0,l,phi_p,theta_p)
        psi_in   += alnu[l]*ul[l]*spu.jn(l,n_p*ka)*Yl0
        psi_out  += alnu[l]*vl[l]*spu.hn1(l,ka)*Yl0
        psi_nu0p += alnu[l]      *spu.jn(l,ka) *Yl0
    # scat = ([y,z+kd,psi_nu0.real,'o'],
    #         [y,z,psi_out.real,'s'])
            # (y,z,psi_in.real,psi_in.real,'d')]
    # dsp.stddisp(scat=scat)
    print(abs(psi_nu0+psi_out-psi_in).sum())
    # print(abs(psi_nu0-psi_nu0p).sum())
